#!/bin/bash -l
#SBATCH --account=m3748_g
#SBATCH --constraint=gpu
#SBATCH --qos=debug
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=none
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --job-name=warpx-psatd-failures
#SBATCH --output=build/failure-audit-2026-09-21/psatd-%j.log

set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
root=$PWD
source build/reaudit-2026-09-21/perlmutter_gpu_warpx.profile
export PYTHONPATH="$root/build_pm_gpu_latest/lib/site-packages"
export OMP_NUM_THREADS=16
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export AMREX_DEFAULT_INIT='amrex.use_gpu_aware_mpi=1'
scripts="$root/Tools/Algorithms/PrescribedFluids"
study="$root/build/failure-audit-2026-09-21/psatd-$SLURM_JOB_ID"
reference="$study/test_rz_rigid_source_PSATD_picmi"
baseline_ranks=${1:-8}
cells_z=${2:-256}
baseline_nodes=$(((baseline_ranks + 3) / 4))
restart_ranks=$((baseline_ranks / 2))
mkdir -p "$study"
printf 'baseline_ranks=%s\ncells_z=%s\nrestart_ranks=%s\n' \
    "$baseline_ranks" "$cells_z" "$restart_ranks" > "$study/parameters.txt"
git rev-parse HEAD > "$study/revision.txt"
git diff --binary > "$study/source.patch"
python3 "$scripts/record_sources.py" "$study/source-sha256.json"
cp build_pm_gpu_latest/CMakeCache.txt "$study/CMakeCache.txt"
sha256sum "$PYTHONPATH"/pywarpx/*.so > "$study/python-library-sha256.txt"
srun -N "$baseline_nodes" -n "$baseline_ranks" -c 32 --cpu-bind=cores \
    bash "$scripts/perlmutter_gpu_rank.sh" \
    python3 "$scripts/gpu_topology.py" > "$study/topology.json"
status=0
for case in test_rz_rigid_source_PSATD_picmi repeat same_rank_restart fewer_rank_restart load_balance fewer_rank_full; do
    ranks=$baseline_ranks
    nodes=$baseline_nodes
    options=()
    if [[ "$case" == fewer_rank* ]]; then
        ranks=$restart_ranks
        nodes=$(((ranks + 3) / 4))
    fi
    if [[ "$case" == *restart ]]; then options+=(--restart "$reference/diags/chk000002"); fi
    if [[ "$case" == load_balance ]]; then options+=(--load-balance); fi
    mkdir -p "$study/$case"
    result=0
    (cd "$study/$case" && srun -N "$nodes" -n "$ranks" -c 32 --cpu-bind=cores \
        bash "$scripts/perlmutter_gpu_rank.sh" python3 \
        "$root/Examples/Tests/collision/inputs_test_rz_rigid_source_picmi.py" \
        --solver PSATD --cells-z "$cells_z" --max-grid-size 16 --subcycles 2 \
        "${options[@]}" > run.log 2>&1) || result=$?
    printf '%s %s\n' "$case" "$result" >> "$study/status.txt"
    if [[ "$result" != 0 ]]; then status=1; fi
    if [[ "$case" != test_rz_rigid_source_PSATD_picmi ]]; then
        for step in 2 6; do
            if [[ -f "$study/$case/state_$step.npz" ]]; then
                python3 "$scripts/report_state_difference.py" \
                    "$reference/state_$step.npz" "$study/$case/state_$step.npz" \
                    --output "$study/$case-step$step-state-difference.json" \
                    --rigid-source-criteria \
                    > "$study/$case-step$step-analysis.log" 2>&1 || status=1
            fi
        done
    fi
done
exit "$status"
