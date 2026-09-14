# Saved validation evidence

This directory records past numerical evidence as well as the checks made
when archiving it. A passing old log is not a claim that the current branch's
new dependency revisions were rebuilt, or that a CPU result validates a GPU.
The production model and collision code are unchanged by this archival update.

## Historical payloads

| Archive | Contents and interpretation |
| --- | --- |
| `validation/legacy-numerical-data.tar.gz` | September 6 total-only-refit lookup table, initial physics/performance NPZ files, three smoke-run OpenPMD datasets and used inputs; obsolete model evidence |
| `validation/production-pic-results.tar.gz` | Saved N2/O2 production-source physics, low-energy, cold, mixed-parent, budget, checkpoint/restart and performance NPZ diagnostics and used inputs |
| `validation/historical-logs.tar.gz` | CPU double/particle-float/Apple-sanitizer CTest logs, unsuccessful conda-sanitizer attempt, integration collision-test summaries and Sphinx warning logs |

Exact members, byte sizes and hashes are in [manifest.json](manifest.json).
The original paths and log timestamps are retained. Production PIC output
can be associated with its used inputs and the test definitions in Git;
the latest integration test additions are in `27cb037d8`. Not every original
run embedded a full Git/dependency fingerprint, so the archive does not
retroactively invent one. Original per-run overwrites are not recoverable.

`development-sync-final-tests.log` reports 79 passing run/analysis tests.
Earlier `development-sync-collisions.log` and `development-sync-2d.log`
include incomplete-build or stale-build-tree-PICMI failures that were resolved
before that final run; they are deliberately not discarded. Sphinx warnings
include documentation-framework issues and are not failed physics checks.
The conda Clang sanitizer attempt stalled before `main`; only the subsequent
AppleClang address/undefined-sanitizer run is counted. Its linked AMReX
library was a release build and was not fully instrumented.

Saved standalone tests cover relativistic/free and molecular kinematics,
the independent Dirac trace, actual table interpolation in both precisions,
moments, conditional energy/angle sampling and phase precision. PIC tests
cover represented yield, unchanged beam particles, equal pair weights and
positions, kinetic spectra including above-free tails, neutral-thermal ions,
capped counts, cold parser/constant inputs, mixed incident energies,
fractional-weight carryover, bare-alpha scaling, and checkpoint budgets.
The [current reference guide](../README.md) gives the commands and exact
acceptance budgets. Numerical tolerances are not data uncertainty estimates.

Performance measurements in that guide are Apple M3 Pro CPU results. The
five-repeat PIC timing medians/ranges survive as a summary; the archive's
NPZ timing files are the saved latest runs, not all five raw repeats. Model
startup is separate from steady-state sampling. Capped collision-only smoke
throughput is not full electromagnetic-simulation throughput. There was no
CUDA/HIP/SYCL compiler/device result to preserve or claim.

## Archival checks

`verify_research_archive.py` checks every listed payload and tar member,
rejects unsafe paths and non-regular archive entries, and reports unlisted
payloads. The Python archive tests also check the final coefficient identity,
the final/source-refresh input-hash chain, finite JSON values and preserved
table dimensions. Reproduction uses the archived inputs without network
access. The check results for this archival update are recorded in
`validation/archive-validation.json`; they are separate from historical logs.

These checks establish preservation and numerical reproducibility, not new
experimental agreement. Limitations such as the N2 low-energy mean discrepancy,
approximate optical closure, unvalidated backward angular spectrum and
relativistic molecular extrapolation remain explicit in the final theory.

The HTML documentation rebuild exited successfully without warnings for the
proton-impact chapter. It was not a warning-free full API build: the local
documentation environment lacks Doxygen/pywarpx and could not fetch optional
Doxygen tag files. The warning output is preserved in
`validation/archive-sphinx-warnings.txt`. Ruff and the archive verifier were
run directly; the full pre-commit runner was not installed. A local pre-commit
hook now verifies the archive whenever its payloads or verifier change.
