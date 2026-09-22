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
#SBATCH --job-name=warpx-audit-suite
#SBATCH --output=build/reaudit-2026-09-21/suite-%j.log

# Run independent studies in one allocation while retaining individual status.
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
status=0
run_study () {
    local tag
    tag=$(IFS=-; echo "$*")
    local result=0
    WARPX_AUDIT_STUDY_TAG="$tag" \
        bash Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch "$@" || result=$?
    printf '%s %s\n' "$tag" "$result" \
        >> "build/reaudit-2026-09-21/suite-$SLURM_JOB_ID.status"
    if [[ "$result" != 0 ]]; then status=1; fi
}

case "${1:?Specify convergence, noise, restart, scaling, physics or coupled_mm}" in
    convergence)
        run_study convergence joint
        run_study convergence solvers
        run_study convergence deposition
        ;;
    noise)
        for suite in fields push source coupled; do run_study noise "$suite"; done
        ;;
    restart)
        for solver in Yee PSATD semi_implicit_em semi_implicit_mm; do
            run_study coupled_restart "$solver"
        done
        ;;
    scaling)
        run_study scaling
        run_study weak
        if [[ "$SLURM_JOB_NUM_NODES" == 1 ]]; then run_study storage; fi
        ;;
    physics)
        run_study physics float
        run_study dcs
        ;;
    coupled_mm)
        run_study coupled semi_implicit_mm
        ;;
    *) echo "Unknown audit suite: $1" >&2; exit 2 ;;
esac
exit "$status"
