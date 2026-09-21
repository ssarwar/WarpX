#!/bin/bash
#SBATCH --account=m3748_g
#SBATCH --constraint=gpu
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --job-name=warpx-fluid-validation
#SBATCH --output=fluid-validation-%j.log

set -euo pipefail
module load gcc-native/13.2 cudatoolkit/13.2 cmake/3.30.2
module load cray-python/3.11.7 cray-hdf5-parallel/1.14.3.9
cd "${SLURM_SUBMIT_DIR}"
source .venv/bin/activate
export PYTHONPATH="$PWD/build/lib/site-packages${PYTHONPATH:+:$PYTHONPATH}"
export MPLCONFIGDIR="$PWD/build/.mpl-cache"
export OMP_NUM_THREADS=1
# These tests exercise ordinary MPI transport, including two-rank restarts.
# Cray's Python MPI binding is not linked to the CUDA GTL library.
export MPICH_GPU_SUPPORT_ENABLED=0

study="build/fluid-a100-${SLURM_JOB_ID}"
mkdir -p "$study"
git rev-parse HEAD > "$study/revision.txt"
nvidia-smi > "$study/device.txt"
module list > "$study/modules.txt" 2>&1

ctest --test-dir build --output-on-failure --timeout 90 -j 1 \
    -R 'test_rz_(rigid_source|rigid_beam|rigid_train|immobile|fluid_charge)' \
    > "$study/ctest.log" 2>&1

# Timing runs are serial and synchronize each step. The field suite measures
# noise versus cost; push uses identical deterministic electron populations.
for suite in fields push source coupled; do
    python Tools/Algorithms/PrescribedFluids/ensemble.py \
        --suite "$suite" --output "$study/$suite" --cells 64 256 \
        --ppc 16 64 256 --seeds 6 --steps 20 \
        --launcher 'srun -n 1 -c 32 --gpus-per-task=1 --cpu-bind=cores'
    python Tools/Algorithms/PrescribedFluids/analyze.py "$study/$suite"
done
