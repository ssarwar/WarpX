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
#SBATCH --time=00:20:00
#SBATCH --job-name=warpx-restart-failures
#SBATCH --output=build/failure-audit-2026-09-21/restarts-%j.log

set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source build/reaudit-2026-09-21/perlmutter_gpu_warpx.profile
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=GPU
export AMREX_DEFAULT_INIT='amrex.use_gpu_aware_mpi=1'
build="$PWD/build_pm_gpu_latest"
if [[ "${1:-branch}" == stock ]]; then build="$PWD/build/stock-warpx/build_pm_gpu_latest"; fi
export PYTHONPATH="$build/lib/site-packages"
output="build/failure-audit-2026-09-21/acceleration-${1:-branch}-$SLURM_JOB_ID"
mkdir -p "$output"
git rev-parse HEAD > "$output/revision.txt"
git diff --binary > "$output/source.patch"
python3 Tools/Algorithms/PrescribedFluids/record_sources.py "$output/source-sha256.json"
cp "$build/CMakeCache.txt" "$output/CMakeCache.txt"
sha256sum "$build"/lib/site-packages/pywarpx/*.so > "$output/python-library-sha256.txt"
if [[ "${1:-branch}" == stock ]]; then
    git -C build/stock-warpx rev-parse HEAD > "$output/stock-revision.txt"
fi
python3 Tools/Algorithms/PrescribedFluids/investigate_restarts.py "$build" "$output/cases"
