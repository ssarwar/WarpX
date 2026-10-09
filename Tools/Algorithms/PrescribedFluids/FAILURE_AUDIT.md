# Perlmutter failure investigation, 22 September 2026

The latest GitHub branch fetched for this investigation was
`codex/rigid-beam-immobile-ions` at
`f5177036b8fbaf8e75d12f734165719a1a0fe3f4`. The Perlmutter checkout was
fast-forwarded against that remote branch; it was already at the same revision.
Build job `58731325` successfully reconfigured and built the CUDA application
and Python libraries for 1D, RZ and 3D. Its elapsed time was 1 minute 58 seconds.
The build directory is
`/pscratch/sd/s/ssarwar/warpx-rigid-fluid-validation/build_pm_gpu_latest`.

The build follows the repository's [Perlmutter instructions](../../../Docs/source/install/hpc/perlmutter.rst)
through [the recorded build recipe](perlmutter_stock_build.sh): stock modules,
dependency versions, Cray compiler wrappers, CUDA architecture 80 and CMake.
The account and software prefix are private to this validation checkout.
MPI, Python, FFT and openPMD are enabled; QED and embedded boundaries are
disabled. The built Python package is selected through `PYTHONPATH` to keep
branch and stock libraries separate. These scope and installation choices
are described in [the build record](REAUDIT.md#build-environment).

All simulations, numerical analysis and lint checks in this investigation run
on Perlmutter. No production equation, collision implementation or regression
tolerance was changed in this follow-up. The changes add failure controls,
complete discrepancy reports and cross-section documentation.

Follow-up GPU job `58732047` passes all ten standalone C++ physics checks,
all 64 Python model/archive tests and the independent Gaussian-field reference.
The raw results and executable hashes are in
`build/reaudit-2026-09-21/physics-failure-followup-58732047`.
CPU job `58732373` recomputes the plotfile reports, verifies immediate restored
fields and passes Ruff and shell-syntax checks. The Sphinx build succeeds with
the previously recorded macro-generated `FieldType` warning; the updated
parameter and PICMI pages render successfully.

## Acceleration restart failures occur in branch and stock

The failing tests are:

- `test_3d_acceleration_restart.analysis`
- `test_3d_acceleration_psatd_restart.analysis`
- `test_3d_acceleration_psatd_time_avg_restart.analysis`

These are upstream **3D boosted-frame accelerator** fixtures: a moving window,
kinetic beams and plasma, and a laser, on a `32 x 32 x 256` mesh. They run
10 steps, checkpoint at step 5, and restart to step 10. They use CKC with PML,
Galilean PSATD, and Galilean PSATD with time averaging, respectively. They do
not use the new rigid fluid beam, immobile ions, PJG or electron MCC.

Jobs `58731776` (branch) and `58731777` (stock `471191e3f`) each run 15
simulations: the three fixtures, each with C++ full/restart runs and Python
native-field full/repeat/restart controls. All 30 simulations finish normally.
Both audit jobs return failure because the unchanged numerical comparisons
fail. The Python observer uses one `evolve(-1)` call, as the application does;
splitting it into `evolve(1)` calls could add momentum synchronization steps.

The comparison matches particle records by their persistent creation-CPU/ID
pair and verifies that the identities are unique and unchanged. It reports
every field instead of stopping at the first assertion. Output times, domain
bounds and mesh dimensions agree. The earlier ordering-only test defect is
already repaired, without changing the original `1e-12` tolerance.

| Control | CKC | PSATD | Time-averaged PSATD |
| --- | --- | --- | --- |
| Branch full versus restart | Fails | Fails | Fails |
| Stock full versus restart | Fails | Fails | Fails |
| Branch independent full versus full | Fails | Fails | Fails |
| Stock independent full versus full | Fails | Fails | Fails |
| Native E/B immediately after loading step 5 | Exact in both builds | Exact in both builds | Exact, including all averaged E/B components, in both builds |

For CKC and ordinary PSATD, the failing quantities are principally the very
small transverse proton momentum and sometimes transverse current. In the
earlier complete-field reanalysis (`58731626`), the branch maximum proton
`p_y` discrepancies are `8.32e-36` and `2.86e-35 kg m/s`, respectively.
Their component-normalized errors exceed `1e-12`, while their errors relative
to the peak momentum vector are `1.67e-18` and `5.74e-18`.

Time-averaged PSATD also has larger field discrepancies and must be considered
separately. The new branch C++ continuation reaches an `E_z` discrepancy of
`4.88e-9` of its component peak; a branch independent uninterrupted repeat
reaches `4.67e-9`. Stock repeats show the same class of failures. Native
averaged fields already differ between independent runs at step 5, before
either run has restarted. The full records include instantaneous and averaged
fields at steps 5, 6 and 10, as well as every plotfile quantity.

These controls exclude missing saved E/B fields as an explanation of the
observed valid-cell discrepancies and show that restart is not necessary to
trigger the failures. GPU deposition/reduction ordering and sensitivity of
the averaged update are consistent with the evidence. The exact kernel-level
amplification mechanism has not been isolated; these remain numerical
reproducibility failures, not a claim of error-free stock numerics. No particle
sorting change or broader field tolerance is used to hide them.

## RZ prescribed-source PSATD ownership failures

The source fixture combines N2/O2 production, capped creation, pending yield,
quiet sampling, kinetic and immobile ion destinations, and subcycling by two.
It checks particle charge footprints, physical totals, source energy budgets,
native fields and immediate restart state. Kinetic products are frozen and
their simulation deposition is disabled in this isolated fixture; immobile
ions still contribute charge. Its solver-corrected current therefore includes
the PSATD continuity correction and is distinct from the stored analytic beam
current.

The original failure was a current comparison after load balancing or MPI
redistribution. The new controls add a second uninterrupted run, a restart
with unchanged rank count, a restart with fewer ranks, load balancing, and a
fresh run with fewer ranks. Every comparison uses the fixture's existing
criteria: `2e-14` pointwise relative error for current and persistent source
state; its existing vector-scale roundoff criterion for E/B. Both the
simulation assertion and supplemental all-field analysis are retained.

| Job and mesh | Repeat full run | Same-rank restart | Fewer-rank restart | Load balancing | Fresh run with fewer ranks |
| --- | --- | --- | --- | --- | --- |
| `58731928`, 16 x 256, 8 GPUs / 2 nodes; fewer = 4 | Current assertion fails | Pass | Pass | Ownership changes; current assertion fails | Current assertion fails |
| `58731946`, 16 x 64, 2 GPUs; fewer = 1 | Pass | Pass | Pass | Ownership changes; pass | Current assertion fails |

For the eight-GPU repeat, 33 small `J_z` values exceed the pointwise criterion.
The largest absolute difference is `1.13687e-13 A/m^2`, or `4.38625e-18` of
peak `J_z`. The load-balanced and fresh four-GPU comparisons have maximum
absolute differences of `5.68434e-14 A/m^2`, or `2.19313e-18` of peak.
The original two-GPU layout's fresh one-GPU comparison fails at three small
`J_r` values; its maximum difference is `7.10543e-15 A/m^2`, or
`6.84425e-18` of peak. Its load-balanced continuation passes in this repeat,
although the earlier regression run failed.

The immediately restored source state passes before any resumed step in all
these restart controls. Particle coordinates/momenta/weights, source budgets,
sampling counters and analytic beam density/current agree exactly in the
native comparisons; accumulated ion densities differ at floating-point
roundoff during evolution. MPI-summed electron populations also show reduction
roundoff, while the individual particle weights agree exactly. Similar
differences appear without checkpointing.
This supports an accumulation/transform ordering explanation rather than lost
or duplicated chemistry. The strict current checks remain recorded failures;
the new passes do not erase the earlier failed runs.

An initial eight-GPU control (`58731699`) used only eight mesh boxes. Its
load-balancing case failed the assertion requiring an actual ownership change,
so it does not establish load-balancing acceptance. The revised two-node case
above uses 16 boxes and reaches the intended ownership comparison. The
initial eight-to-four and same-rank restart runs themselves passed.

## Other failures retained in the audit record

[REAUDIT.md](REAUDIT.md) and its two earlier JSON archives retain the original
failures, investigations and successful follow-up jobs. They are not replaced
by the narrower controls above.

| Failure class | Disposition |
| --- | --- |
| Missing external-field and helium reference data | Public datasets installed at the stock CI paths; all 20 repeated stages pass in `58718189` |
| QED-dependent inputs in a QED-disabled build | Configuration mismatch; explicitly excluded from this simulation's acceptance |
| PSATD output analysis inferred Yee from directory name | Explicit solver argument; ownership/plotfile/openPMD checks repeated |
| Particle restart comparison assumed unchanged storage order | Match and validate persistent identities; five stock/branch failures resolved, three acceleration failures retained above |
| Mixed-energy PJG percentile sampling | Independent weight/seed convergence; finer product weight and adequate creation cap pass the unchanged spectrum bounds |
| Native-float monoenergetic PJG table mismatch | Same device row selection/fraction as the particle executor; strict float/double boundary and tail tests pass |
| Release GPU parser/majorant validation | Host-visible device error flags replace disabled assertions; negative-input tests pass |
| Per-fluid macroparticle-count diagnostic indexed kinetic species | Fluid macroparticle counts return zero; kinetic histogram, plotfile, openPMD and restart checks pass |
| FieldProbe/Poynting initial integration and cumulative restart defects | Reproduced on stock and repaired; exact constant-field/surface-flux checks pass |
| RZ/1D Poynting component-index and cell-centered interpolation defects | Independent manufactured-field tests expose and verify repairs |
| Stale incremental object files and MPI launcher flags | Rebuild with changed-source mtimes and clear cached flags; invalid runs are not physics acceptance |
| CUDA captures in earlier fixes | Corrected and rebuilt before runtime acceptance |
| Initial SYCL dependency/link configuration | Corrected; compilation passes, runtime coverage remains unavailable |
| Sphinx macro-generated enum warning | Existing Doxygen `FieldType` warning remains; relevant parameter/PICMI pages build |

The first new acceleration harness submissions (`58731621`, `58731622`)
duplicated AMReX's input prefix and stopped before initialization. The corrected
jobs above retain the CTest-relative input path. Submission `58731929` was
cancelled while pending after identifying a topology-probe node-count mismatch;
`58731946` is its corrected replacement. None is interpreted as a physics
failure or successful acceptance run.

## Automatic proton-impact cross sections

`collision_type = proton_impact_ionization` uses the built-in calibrated
PJG-type N2/O2 SDCS. `PJGModel` numerically integrates that SDCS to obtain
the total inclusive electron-production cross section and interpolates it on
its incident-energy grid. No external total proton-impact cross-section table
is required for either projectile representation.

- Particle projectiles evaluate the total at each projectile's kinetic energy.
- A fixed-energy rigid fluid projectile builds only the two required energy
  rows, using the same interpolation as the full particle table, and caches
  the total and sampling state.
- Separate N2/O2 collision instances supply their molecular number densities.
  For protons, the combined primary source is
  `n_b * abs(v_b) * (n_N2*sigma_N2(E_b) + n_O2*sigma_O2(E_b))`.
- Bare heavy-ion projectiles apply the implemented mass mapping and `Z_b^2`
  factor. The IAA secondary-angle closure changes sampled directions, not
  the total-rate normalization.
- Electron MCC total-rate tables are separate user inputs. Its RBEQ energy
  sharing and IAA angular models do not replace those total-rate tables.

The SDCS represents inclusive electron yield through effective electron/ion
pairs. The 800 MeV use is an extrapolation of the calibrated molecular model,
as documented in [the theory](../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst).
Automatic numerical integration does not supply new experimental validation.

## Evidence and reproduction

[The machine-readable archive](results/perlmutter-failure-audit-2026-09-22.json)
contains all reported field differences, both build revisions, commands,
library hashes, source-manifest hashes, topology and scheduler status. Raw
logs, fields and checkpoints remain under the approved Perlmutter root in
`build/failure-audit-2026-09-21`; previous physics/performance outputs remain in
`build/reaudit-2026-09-21`.

Run from the validation checkout after following the recorded build recipe:

```bash
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_restart_failures.sbatch branch
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_restart_failures.sbatch stock
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_psatd_failures.sbatch 8 256
sbatch --nodes=1 Tools/Algorithms/PrescribedFluids/perlmutter_psatd_failures.sbatch 2 64
```

`report_plotfile_difference.py` checks all quantities using the original
acceleration tolerance. `report_state_difference.py --rigid-source-criteria`
evaluates the existing source-fixture assertions. The latter is also applied
to fresh-run controls that the original fixture does not compare.
`collect_failure_audit.py --help` describes archival arguments for new job
directories. Its `--verify-plotfiles` option recomputes the complete reports
and checks time/domain metadata before archiving. A successful simulation
process alone is not acceptance when the associated comparison fails.
