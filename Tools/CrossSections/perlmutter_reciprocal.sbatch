#!/bin/bash
#SBATCH --constraint=gpu&hbm40g
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00

# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# Supply --account=<GPU allocation>. Production and cumulative-reference
# bundles must already exist; this job performs no scientific data generation.
set -euo pipefail
ulimit -c 0
cd "$(git rev-parse --show-toplevel)"
: "${WARPX_ROTATION_PYTHON:?Set the prepared Python environment}"
: "${WARPX_ROTATION_DATA:?Set the prepared production data root}"
: "${WARPX_ROTATION_CUMULATIVE:?Set the offline cumulative-reference data root}"
: "${WARPX_ROTATION_SOURCE_DATA:?Set the ordinary cross-section source root}"
rotation_build="${WARPX_ROTATION_BUILD:-${PWD}/build}"
rotation_output="${PWD}/build/reciprocal-gpu-${SLURM_JOB_ID}"
mkdir -p "${rotation_output}"
export CCACHE_DIR="${PWD}/build/ccache"
export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export CRAY_ACCEL_TARGET=nvidia80
export PYTHONPATH="${rotation_build}/lib/site-packages${PYTHONPATH:+:${PYTHONPATH}}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
nvidia-smi --query-gpu=name --format=csv,noheader |
    awk '!/A100/ {bad=1} END {exit bad}'
nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits |
    awk '$1 < 38000 || $1 > 45000 {bad=1} END {exit bad}'
rotation_mpi=OFF
if [[ "${SLURM_NTASKS:-1}" -gt 1 ]]; then rotation_mpi=ON; fi
rotation_dependencies=()
if [[ "${rotation_mpi}" == ON ]]; then
    rotation_dependencies+=("-DCMAKE_CUDA_HOST_COMPILER=$(command -v CC)")
    export AMREX_DEFAULT_INIT="${AMREX_DEFAULT_INIT:-} amrex.use_gpu_aware_mpi=1"
fi
if [[ ! -f "${rotation_build}/CMakeCache.txt" && -d "${PWD}/build/_deps/fetchedamrex-src" ]]; then
    rotation_dependencies+=("-DWarpX_amrex_src=${PWD}/build/_deps/fetchedamrex-src")
fi
if [[ ! -f "${rotation_build}/CMakeCache.txt" && -d "${PWD}/build/_deps/fetchedpyamrex-src" ]]; then
    rotation_dependencies+=("-DWarpX_pyamrex_src=${PWD}/build/_deps/fetchedpyamrex-src")
fi
if [[ ! -f "${rotation_build}/CMakeCache.txt" && -d "${PWD}/build/_deps/fetchedpybind11-src" ]]; then
    rotation_dependencies+=("-DFETCHCONTENT_SOURCE_DIR_FETCHEDPYBIND11=${PWD}/build/_deps/fetchedpybind11-src")
fi
cmake -S . -B "${rotation_build}" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=CC -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DWarpX_DIMS=1 -DWarpX_COMPUTE=CUDA -DWarpX_MPI="${rotation_mpi}" \
    -DWarpX_PRECISION="${WARPX_ROTATION_REAL_PRECISION:-DOUBLE}" \
    -DWarpX_PARTICLE_PRECISION="${WARPX_ROTATION_PARTICLE_PRECISION:-DOUBLE}" \
    -DWarpX_PYTHON=ON -DWarpX_PYTHON_IPO=OFF -DpyAMReX_IPO=OFF \
    -DWarpX_FFT=OFF -DWarpX_OPENPMD=OFF -DWarpX_QED=OFF \
    -DPython_EXECUTABLE="${WARPX_ROTATION_PYTHON}" "${rotation_dependencies[@]}"
cmake --build "${rotation_build}" -j 8 --target test_reciprocal_rotation \
    benchmark_reciprocal_mcc test_rotation_thresholds pyWarpX_1d \
    pyWarpX_python_sources pyAMReX_python_sources
srun -n 1 --gpus-per-task=1 "${rotation_build}/bin/test_rotation_thresholds"
if [[ "${rotation_mpi}" == OFF ]]; then
    srun -n 1 --gpus-per-task=1 "${WARPX_ROTATION_PYTHON}" \
        Examples/Tests/collision/analysis_reciprocal_inputs.py \
        --program "${rotation_build}/bin/benchmark_reciprocal_mcc" \
        --sampler "${rotation_build}/bin/test_reciprocal_rotation" \
        --data-dir "${WARPX_ROTATION_DATA}" --source-dir "${WARPX_ROTATION_SOURCE_DATA}" \
        --output "${rotation_output}/bounds"
fi
for rotation_gas in N2 O2; do
    for rotation_mode in alias cumulative; do
        rotation_root="${WARPX_ROTATION_DATA}"
        rotation_cumulative=0
        if [[ "${rotation_mode}" == cumulative ]]; then
            rotation_root="${WARPX_ROTATION_CUMULATIVE}"
            rotation_cumulative=1
        fi
        srun -n 1 --gpus-per-task=1 "${rotation_build}/bin/test_reciprocal_rotation" \
            "file=${rotation_root}/${rotation_gas}/IAA/reciprocal_hybrid_300K/thermal_rotation.rot" \
            "cross_section=${WARPX_ROTATION_SOURCE_DATA}/${rotation_gas}/IAA/elastic.txt" \
            "output=${rotation_output}/${rotation_gas}-${rotation_mode}.txt" \
            "cumulative=${rotation_cumulative}" samples=1048576 timing_repetitions=8 \
            amrex.the_arena_init_size=8388608
        "${WARPX_ROTATION_PYTHON}" Tools/CrossSections/reciprocal_rotation/sampling_reference.py \
            --bundle "${rotation_root}/${rotation_gas}/IAA/reciprocal_hybrid_300K/thermal_rotation.rot" \
            --sampler "${rotation_output}/${rotation_gas}-${rotation_mode}.txt" --samples 1048576
    done
    # This isolated sampler has no device MPI traffic. Cray MPICH probes host
    # pointers through cuPointerGetAttribute, which memcheck reports as API
    # errors. Disable those probes here; timed MPI runs remain GPU-aware.
    srun -n 1 --gpus-per-task=1 env MPICH_GPU_SUPPORT_ENABLED=0 DARSHAN_DISABLE=1 \
        AMREX_DEFAULT_INIT=amrex.use_gpu_aware_mpi=0 \
        compute-sanitizer --tool memcheck --error-exitcode=99 \
        "${rotation_build}/bin/test_reciprocal_rotation" \
        "file=${WARPX_ROTATION_DATA}/${rotation_gas}/IAA/reciprocal_hybrid_300K/thermal_rotation.rot" \
        "cross_section=${WARPX_ROTATION_SOURCE_DATA}/${rotation_gas}/IAA/elastic.txt" \
        "output=${rotation_output}/${rotation_gas}-memcheck.txt" samples=4096 \
        amrex.the_arena_init_size=8388608
done
if [[ "${rotation_mpi}" == OFF ]]; then
    srun -n 1 --gpus-per-task=1 compute-sanitizer --tool memcheck --error-exitcode=99 \
        "${WARPX_ROTATION_PYTHON}" Examples/Tests/collision/inputs_test_1d_reciprocal_rotation_picmi.py \
        --data-dir "${WARPX_ROTATION_DATA}" --source-dir "${WARPX_ROTATION_SOURCE_DATA}" \
        --target air --pic --particles 8192 --steps 8 --density 1e28 --dt 1e-16 \
        --output "${rotation_output}/mcc-memcheck.json"
fi
if [[ -n "${WARPX_ROTATION_IONIZATION_DATA:-}" ]]; then
    for rotation_gas in N2 O2; do
        srun -n "${SLURM_NTASKS:-1}" --gpus-per-task=1 "${WARPX_ROTATION_PYTHON}" \
            Examples/Tests/collision/inputs_test_1d_reciprocal_rotation_picmi.py \
            --data-dir "${WARPX_ROTATION_DATA}" --source-dir "${WARPX_ROTATION_SOURCE_DATA}" \
            --target "${rotation_gas}" --mode mono --energy 100 --translational-temperature 0 \
            --ionization-table "${WARPX_ROTATION_IONIZATION_DATA}/${rotation_gas}.txt" \
            --particles 131072 --steps 1 --density 2e27 --dt 1e-14 \
            --output "${rotation_output}/combined-ionization-${rotation_gas}.json"
    done
fi
for rotation_gas in N2 O2 air; do
    for rotation_electron_temperature in 75 300 1200; do
        srun -n "${SLURM_NTASKS:-1}" --gpus-per-task=1 "${WARPX_ROTATION_PYTHON}" \
            Examples/Tests/collision/inputs_test_1d_reciprocal_rotation_picmi.py \
            --data-dir "${WARPX_ROTATION_DATA}" --source-dir "${WARPX_ROTATION_SOURCE_DATA}" \
            --target "${rotation_gas}" --electron-temperature "${rotation_electron_temperature}" \
            --particles 131072 --steps 512 --dt 1e-17 --density 1e30 --check-thermal \
            --output "${rotation_output}/thermal-${rotation_gas}-${rotation_electron_temperature}.json"
    done
done
if [[ "${WARPX_ROTATION_BENCHMARK:-1}" == 0 ]]; then exit 0; fi
rotation_comparison=(--disabled)
if [[ -n "${WARPX_ROTATION_REFERENCE:-}" ]]; then
    rotation_comparison=(--reference-dir "${WARPX_ROTATION_REFERENCE}" \
        --samplings alias cumulative reference)
fi
"${WARPX_ROTATION_PYTHON}" Tools/CrossSections/benchmark_reciprocal_rotation.py \
    --data-dir "${WARPX_ROTATION_DATA}" --cumulative-dir "${WARPX_ROTATION_CUMULATIVE}" \
    --source-dir "${WARPX_ROTATION_SOURCE_DATA}" --output "${rotation_output}/timesteps" \
    --mpi-ranks "${SLURM_NTASKS:-1}" --particles 262144 1048576 \
    --targets air --steps 64 --repeats 5 "${rotation_comparison[@]}" \
    --mcc-program "${rotation_build}/bin/benchmark_reciprocal_mcc"
