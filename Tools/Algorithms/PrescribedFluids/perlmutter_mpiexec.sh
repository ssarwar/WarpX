#!/usr/bin/env bash
# CTest's MPI launcher: select the device before MPI initializes the GPU transport.
set -euo pipefail
if [[ "$1" != -n ]]; then
    echo 'Usage: perlmutter_mpiexec.sh -n RANKS COMMAND [ARGUMENTS...]' >&2
    exit 2
fi
ranks=$2
shift 2
directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec srun -N $(((ranks + 3) / 4)) -n "$ranks" -c 32 --cpu-bind=cores \
    bash "$directory/perlmutter_gpu_rank.sh" "$@"
