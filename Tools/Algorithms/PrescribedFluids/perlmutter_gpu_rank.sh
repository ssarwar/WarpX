#!/usr/bin/env bash
# Match the GPU/NIC placement in Tools/machines/perlmutter-nersc/perlmutter_gpu.sbatch.
set -euo pipefail
export CUDA_VISIBLE_DEVICES=$((3 - SLURM_LOCALID))
exec "$@"
