#!/bin/bash -l
#SBATCH --account=m3748_g
#SBATCH --constraint=gpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --job-name=warpx-coupled-audit
#SBATCH --output=build/reaudit-2026-09-21/run-%j.log

set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
root=$PWD
audit="$root/build/reaudit-2026-09-21"
source "$audit/perlmutter_gpu_warpx.profile"
export OMP_NUM_THREADS=16
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export AMREX_DEFAULT_INIT='amrex.use_gpu_aware_mpi=1'
export MPLCONFIGDIR="$audit/mpl-cache"
phase=${1:-regression}
build_name=${WARPX_AUDIT_BUILD_NAME:-build_pm_gpu_py}
validation_build="$root/$build_name"
if [[ "${2:-branch}" == stock ]]; then
    validation_build="$root/build/stock-warpx/$build_name"
fi
export PYTHONPATH="$validation_build/lib/site-packages"
study="$audit/$phase${WARPX_AUDIT_STUDY_TAG:+-$WARPX_AUDIT_STUDY_TAG}-$SLURM_JOB_ID"
mkdir -p "$study"
git rev-parse HEAD > "$study/revision.txt"
git diff --binary > "$study/source.patch"
python3 Tools/Algorithms/PrescribedFluids/record_sources.py "$study/source-sha256.json"
cp "$validation_build/CMakeCache.txt" "$study/CMakeCache.txt"
sha256sum "$validation_build"/lib/site-packages/pywarpx/*.so > "$study/python-library-sha256.txt"
printf '%s\n' "$validation_build" > "$study/validation-build.txt"
if [[ "${2:-branch}" == stock ]]; then
    git -C "$root/build/stock-warpx" rev-parse HEAD > "$study/stock-revision.txt"
fi
module list > "$study/modules.txt" 2>&1
python3 -m pip freeze > "$study/python-packages.txt"
nvidia-smi > "$study/device.txt"
for name in AMREX_DEFAULT_INIT CUDA_MODULE_LOADING MPICH_GPU_SUPPORT_ENABLED \
    MPICH_OFI_NIC_POLICY OMP_NUM_THREADS SLURM_JOB_ID SLURM_JOB_NUM_NODES SLURM_NTASKS; do
    printf '%s=%s\n' "$name" "${!name:-unset}"
done > "$study/runtime.txt"
rank="$root/Tools/Algorithms/PrescribedFluids/perlmutter_gpu_rank.sh"
launcher="srun -N 1 -n 1 -c 32 --cpu-bind=cores bash $rank"
srun -N "$SLURM_JOB_NUM_NODES" -n "$SLURM_NTASKS" -c 32 --cpu-bind=cores \
    bash "$rank" python3 Tools/Algorithms/PrescribedFluids/gpu_topology.py \
    > "$study/topology.json"

if [[ "$phase" == topology ]]; then
    exit 0
elif [[ "$phase" == physics ]]; then
    status=0
    physics_build=$build_name
    if [[ "${2:-double}" == float ]]; then physics_build="${build_name}_float"; fi
    printf '%s\n' "$root/$physics_build" > "$study/physics-build.txt"
    sha256sum "$physics_build"/{fluid-tests,pjg,mcc-physics}/test_* \
        > "$study/physics-executable-sha256.txt"
    if [[ "${2:-double}" == float ]]; then
        cp "$audit/amrex-float-$build_name/CMakeCache.txt" \
            "$study/physics-AMReX-CMakeCache.txt"
    fi
    for directory in fluid-tests pjg mcc-physics; do
        $launcher ctest --test-dir "$physics_build/$directory" --output-on-failure \
            --output-junit "$study/$directory.xml" > "$study/$directory.log" 2>&1 || status=1
    done
    python3 -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v \
        > "$study/python-physics.log" 2>&1 || status=1
    python3 Tools/Algorithms/PrescribedFluids/gaussian_reference.py \
        > "$study/gaussian.log" 2>&1 || status=1
    exit "$status"
elif [[ "$phase" == preflight ]]; then
    status=0
    ctest --test-dir "$validation_build" -j 1 --output-on-failure --timeout 180 \
        --output-junit "$study/ctest.xml" \
        -R 'test_rz_fluid_interfaces|test_rz_rigid_source_Yee_picmi|test_1d_background_mcc_bad_inputs|test_3d_proton_impact_ionization_bad_inputs' \
        -E checksum > "$study/ctest.log" 2>&1 || status=1
    (cd "$study" && $launcher python3 \
        "$root/Examples/Tests/collision/inputs_test_3d_proton_impact_ionization_bad_inputs_picmi.py" \
        --rz > rigid-inputs.log 2>&1) || status=1
    exit "$status"
elif [[ "$phase" == runtime ]]; then
    fixture="$root/Examples/Tests/collision"
    status=0
    (cd "$study" && $launcher python3 \
        "$fixture/inputs_test_1d_background_mcc_bad_inputs_picmi.py" > mcc.log 2>&1) || status=1
    for geometry in particle rigid; do
        options=()
        if [[ "$geometry" == rigid ]]; then options=(--rz); fi
        (cd "$study" && $launcher python3 \
            "$fixture/inputs_test_3d_proton_impact_ionization_bad_inputs_picmi.py" \
            "${options[@]}" > "$geometry.log" 2>&1) || status=1
    done
    exit "$status"
elif [[ "$phase" == mixed_source ]]; then
    fixture="$root/Examples/Tests/collision"
    status=0
    for weight in 200 50 12.5; do
        for seed in 41 42 43 44 45 46; do
            output="$study/weight$weight-seed$seed"
            mkdir -p "$output"
            (cd "$output" && $launcher python3 \
                "$fixture/inputs_test_3d_proton_impact_ionization_picmi.py" \
                --mixed --seed "$seed" --product-weight "$weight" --max-products 500000 \
                > run.log 2>&1) || status=1
            (cd "$output" && python3 "$fixture/analysis_proton_impact_ionization.py" \
                > analysis.log 2>&1) || status=1
        done
    done
    exit "$status"
elif [[ "$phase" == poynting ]]; then
    cmake --build "$validation_build/fluid-tests" --parallel 4 > "$study/build.log" 2>&1
    $launcher ctest --test-dir "$validation_build/fluid-tests" -R poynting \
        --output-on-failure --output-junit "$study/ctest.xml" > "$study/ctest.log" 2>&1
elif [[ "$phase" == regression ]]; then
    options=()
    if [[ -n "${3:-}" ]]; then options+=(--include "$3"); fi
    python3 Tools/Algorithms/PrescribedFluids/run_regression.py "$validation_build" "$study" "${options[@]}" \
        > "$study/ctest.log" 2>&1
elif [[ "$phase" == convergence ]]; then
    comparison=${2:-source}
    python3 Tools/Algorithms/PrescribedFluids/convergence.py "$comparison" \
        --output "$study/$comparison" --mcc-all --launcher "$launcher"
    python3 Tools/Algorithms/PrescribedFluids/analyze_convergence.py "$comparison" \
        "$study/$comparison"
elif [[ "$phase" == dcs ]]; then
    fixture="$root/Examples/Tests/collision"
    (cd "$study" && $launcher python3 \
        "$fixture/inputs_test_1d_background_mcc_iaa_elastic_picmi.py" \
        --n2-dcs "$audit/reference-data/DCS.e-N2" \
        --o2-dcs "$audit/reference-data/DCS.e-O2" > run.log 2>&1)
    (cd "$study" && python3 "$fixture/analysis_background_mcc_iaa_elastic.py" \
        > analysis.log 2>&1)
elif [[ "$phase" == probes || "$phase" == flux ]]; then
    fixture="$root/Examples/Tests/reduced_diags"
    mkdir -p "$study/full" "$study/restart"
    options=(--record-only)
    analysis_options=()
    if [[ "$phase" == flux ]]; then
        options+=(--nonzero-flux --solver "${3:-Yee}")
        analysis_options+=(--nonzero-flux)
    fi
    if [[ "${2:-branch}" == stock ]]; then
        export PYTHONPATH="$root/build/stock-warpx/$build_name/lib/site-packages"
        options+=(--stock)
    fi
    (cd "$study/full" && srun -N 1 -n 4 -c 32 --cpu-bind=cores bash "$rank" \
        python3 "$fixture/inputs_test_rz_fluid_diagnostics_picmi.py" "${options[@]}" \
        > run.log 2>&1)
    (cd "$study/restart" && srun -N 1 -n 2 -c 32 --cpu-bind=cores bash "$rank" \
        python3 "$fixture/inputs_test_rz_fluid_diagnostics_picmi.py" "${options[@]}" \
        --continue-output --restart "$study/full/diags/chk000002" > run.log 2>&1)
    (cd "$study" && python3 "$fixture/analysis_rz_fluid_diagnostics.py" full restart "${analysis_options[@]}" \
        > analysis.log 2>&1)
elif [[ "$phase" == profile ]]; then
    for ranks in 1 4; do
        for suite in source coupled; do
            output="$study/$suite-ranks$ranks"
            options=()
            if [[ "$suite" == coupled ]]; then options=(--mcc-all); fi
            python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite "$suite" \
                --output "$output" --cells 64 256 --beams fluid quiet --ppc 64 \
                --seeds 1 --steps 20 --subcycles 2 --profile \
                --launcher "srun -N 1 -n $ranks -c 32 --cpu-bind=cores bash $rank" \
                "${options[@]}"
        done
    done
elif [[ "$phase" == noise ]]; then
    suite=${2:-coupled}
    options=()
    if [[ "$suite" == coupled ]]; then options=(--mcc-all); fi
    python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite "$suite" \
        --output "$study/$suite" --cells 64 256 --ppc 16 64 256 \
        --seeds 6 --steps 20 --subcycles 2 --launcher "$launcher" "${options[@]}"
    python3 Tools/Algorithms/PrescribedFluids/analyze.py "$study/$suite"
    if [[ "$suite" == coupled ]]; then
        python3 Tools/Algorithms/PrescribedFluids/assess_ensemble.py "$study/$suite"
    fi
elif [[ "$phase" == coupled ]]; then
    solver=${2:-Yee}
    for subcycles in 1 2 4; do
        output="$study/$solver-subcycles$subcycles"
        python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite coupled \
            --output "$output" --cells 64 256 --max-grid-size 32 \
            --beams fluid quiet --ppc 64 --seeds 6 --steps 20 --solver "$solver" \
            --mcc-all --subcycles "$subcycles" --launcher "$launcher"
        python3 Tools/Algorithms/PrescribedFluids/analyze.py "$output"
        python3 Tools/Algorithms/PrescribedFluids/assess_ensemble.py "$output"
    done
elif [[ "$phase" == equal_error ]]; then
    python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite fields \
        --output "$study/fields" --cells 256 1024 \
        --radial-sigmas 32 --longitudinal-sigmas 48 --ppc 4 16 64 256 \
        --seeds 6 --steps 20 --dt 1e-13 --launcher "$launcher"
    python3 Tools/Algorithms/PrescribedFluids/analyze_equal_error.py "$study/fields"
elif [[ "$phase" == coupled_restart ]]; then
    solver=${2:-Yee}
    ranks=$((4 * SLURM_JOB_NUM_NODES))
    for seed in 41 42 43 44 45 46; do
        output="$study/$solver-seed$seed"
        mkdir -p "$output/full" "$output/restart"
        options=(--mode coupled --mcc-all --beam fluid --ions fluid --cells 32 256
                 --max-grid-size 32 --steps 20 --subcycles 2 --solver "$solver"
                 --seed "$seed" --checkpoint --checkpoint-period 10 --snapshot-particles)
        (cd "$output/full" && srun -N "$SLURM_JOB_NUM_NODES" -n "$ranks" -c 32 \
            --cpu-bind=cores bash "$rank" python3 \
            "$root/Tools/Algorithms/PrescribedFluids/benchmark.py" "${options[@]}" \
            > run.log 2>&1)
        (cd "$output/restart" && srun -N 1 -n $((ranks / 2)) -c 32 \
            --cpu-bind=cores bash "$rank" python3 \
            "$root/Tools/Algorithms/PrescribedFluids/benchmark.py" "${options[@]}" \
            --restart "$output/full/diags/chk000010" --check-restored "$output/full/result.npz" \
            > run.log 2>&1)
    done
    python3 Tools/Algorithms/PrescribedFluids/analyze_restart.py "$study"
elif [[ "$phase" == ownership ]]; then
    status=0
    ranks=$((4 * SLURM_JOB_NUM_NODES))
    resumed=$((ranks / 2))
    launch_full="srun -N $SLURM_JOB_NUM_NODES -n $ranks -c 32 --cpu-bind=cores bash $rank"
    launch_restart="srun -N 1 -n $resumed -c 32 --cpu-bind=cores bash $rank"
    fixture="$root/Examples/Tests/collision"
    for solver in Yee PSATD semi_implicit_em semi_implicit_mm; do
        for process in source attachment; do
            output="$study/$solver-$process"
            mkdir -p "$output/full" "$output/restart"
            # Eight axial boxes occupy all ranks on two nodes and exceed
            # PSATD's eight-cell guards. RZ PSATD does not partition r.
            options=(--solver "$solver" --max-grid-size 16 --cells-z 128)
            if [[ "$process" == source ]]; then
                input="$fixture/inputs_test_rz_rigid_source_picmi.py"
                options+=(--subcycles 2)
                if [[ "$solver" == Yee ]]; then options+=(--openpmd); fi
            else
                input="$fixture/inputs_test_rz_immobile_products_picmi.py"
                options+=(--kind attachment)
            fi
            (cd "$output/full" && $launch_full python3 "$input" "${options[@]}" > run.log 2>&1) || status=1
            (cd "$output/restart" && $launch_restart python3 "$input" "${options[@]}" \
                --restart "$output/full/diags/chk000002" > run.log 2>&1) || status=1
            if [[ "$process" == source ]]; then
                analysis_options=(--cells-z 128 --solver "$solver")
                if [[ "$solver" == Yee ]]; then analysis_options+=(--openpmd); fi
                for directory in full restart; do
                    (cd "$output/$directory" && python3 \
                        "$fixture/analysis_rz_prescribed_outputs.py" "${analysis_options[@]}" > analysis.log 2>&1) || status=1
                done
            fi
        done
    done
    exit "$status"
elif [[ "$phase" == storage ]]; then
    for steps in 20 80; do
        output="$study/steps$steps"
        python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite coupled \
            --output "$output" --cells 64 256 --beams fluid --seeds 3 \
            --steps "$steps" --mcc-all --subcycles 2 --checkpoint \
            --launcher "$launcher"
        python3 Tools/Algorithms/PrescribedFluids/analyze.py "$output"
    done
elif [[ "$phase" == weak ]]; then
    counts='1 2 4'
    if [[ "$SLURM_JOB_NUM_NODES" -gt 1 ]]; then counts='8'; fi
    for ranks in $counts; do
        nodes=$(((ranks + 3) / 4))
        # Fixed dx and electron work per rank; the prescribed beam is held
        # fixed while the uniform kinetic-electron domain grows in z.
        output="$study/push-ranks$ranks"
        python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite push \
            --output "$output" --cells 64 $((256 * ranks)) \
            --longitudinal-sigmas $((12 * ranks)) --beams fluid quiet --ppc 64 \
            --seeds 3 --steps 20 --electron-ppc 16 \
            --launcher "srun -N $nodes -n $ranks -c 32 --cpu-bind=cores bash $rank"
        python3 Tools/Algorithms/PrescribedFluids/analyze.py "$output"
    done
elif [[ "$phase" == scaling ]]; then
    counts='1 2 4'
    if [[ "$SLURM_JOB_NUM_NODES" -gt 1 ]]; then counts='8'; fi
    for ranks in $counts; do
        nodes=$(((ranks + 3) / 4))
        for suite in fields push coupled; do
            output="$study/$suite-ranks$ranks"
            options=()
            if [[ "$suite" == coupled ]]; then options=(--mcc-all); fi
            python3 Tools/Algorithms/PrescribedFluids/ensemble.py --suite "$suite" \
                --output "$output" --cells 256 1024 --beams fluid quiet --ppc 64 \
                --seeds 3 --steps 20 --dt 1e-13 --weight 10 --subcycles 2 \
                --launcher "srun -N $nodes -n $ranks -c 32 --cpu-bind=cores bash $rank" \
                "${options[@]}"
            python3 Tools/Algorithms/PrescribedFluids/analyze.py "$output"
        done
    done
else
    echo "Unknown re-audit phase: $phase" >&2
    exit 2
fi
