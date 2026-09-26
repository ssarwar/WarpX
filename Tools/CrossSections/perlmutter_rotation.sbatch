#!/bin/bash
#SBATCH --constraint=gpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00

# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# Submit from the repository root with --account=<GPU allocation> and a Python
# environment containing NumPy, SciPy, CuPy and the WarpX build dependencies:
# WARPX_ROTATION_PYTHON=/path/to/python sbatch ... this_file [before|after]
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
: "${WARPX_ROTATION_PYTHON:?Set the Python interpreter for this validation}"
rotation_label=${1:-after}
rotation_output="${PWD}/build/rotation-${rotation_label}"
rotation_data_args=()
if [[ -n "${WARPX_ROTATION_DATA:-}" ]]; then
    rotation_data_args=(--data-dir "${WARPX_ROTATION_DATA}")
fi
mkdir -p "${rotation_output}"
export CCACHE_DIR="${PWD}/build/ccache"
export PYTHONPATH="${PWD}/build/lib/site-packages${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1

git rev-parse HEAD
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC \
    -DCMAKE_CUDA_ARCHITECTURES=80 \
    -DWarpX_DIMS=1 -DWarpX_COMPUTE=CUDA -DWarpX_MPI=OFF \
    -DWarpX_PYTHON=ON -DWarpX_PYTHON_IPO=OFF -DpyAMReX_IPO=OFF \
    -DWarpX_FFT=OFF -DWarpX_OPENPMD=OFF -DWarpX_QED=OFF \
    -DPython_EXECUTABLE="${WARPX_ROTATION_PYTHON}"
cmake --build build -j "${SLURM_CPUS_PER_TASK:-16}"
ctest --test-dir build --output-on-failure \
    -R 'rotation|background_mcc_iaa_elastic|background_mcc_bad_inputs|background_mcc_relativistic|background_mcc_nonrelativistic' \
    | tee "${rotation_output}/physics.log"
srun -n 1 "${WARPX_ROTATION_PYTHON}" Tools/CrossSections/benchmark_rotation.py \
    "${rotation_data_args[@]}" \
    --output "${rotation_output}/timesteps" --particles 524288 --steps 32 --repeats 5
for rotation_mode in 0 1; do
    srun -n 1 build/bin/test_mcc_rotation \
        reference="${PWD}/build/Tools/CrossSections/rotation_reference/reference.txt" \
        output="${rotation_output}/sampler-${rotation_mode}.txt" \
        cumulative="${rotation_mode}" samples=1048576
done
