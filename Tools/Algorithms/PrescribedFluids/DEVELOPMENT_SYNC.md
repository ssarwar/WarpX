# Fork development integration, 2026-09-22

## Source provenance and scope

The integration branch is
`codex/rigid-beam-immobile-ions-development-sync`. Its first parent is the
user's fork `ssarwar/WarpX:development`, at
`cee2ca2fd5d514dca1c8096571e4b213a7e5a125`; its merged parent is
`ssarwar/WarpX:codex/rigid-beam-immobile-ions`, at
`a681a47c7f76be5f356c2d0e269c45b3933f8b76`. Their merge base is
`8d71e576744b961cc709a79d32b7f912a97235ee`. No upstream repository branch
was substituted for the fork's development branch. Neither source branch
is changed by this integration. No pull request is opened.

The resolved source tree, before adding this report and follow-up validation
changes, is `a671deb9e1ee4d0bd59a5edaf4eaa2422160309a`. The same indexed tree
was transferred to the authorized Perlmutter checkout at
`/pscratch/sd/s/ssarwar/warpx-rigid-fluid-validation` and verified with
`git write-tree` before CUDA compilation. Raw logs are retained under
`build/development-sync-cee2ca2fd` on each machine.

## Merge conflicts

The four conflicts come from the overlap between the feature branch's
collision work and development's Legendre fusion scattering commit
`a9811a54a`, rather than the mass-matrix commit:

- `BackgroundMCCCollision.cpp`: development extends the shared two-product
  scattering call with a Legendre coefficient table and range-status pointer;
  the feature branch replaced the surrounding collision selection,
  kinematics, attachment, and product-destination logic. Retain that complete
  feature implementation and pass an empty table and null status for its
  non-Legendre scattering call. MCC explicitly rejects Legendre input.
- `ImpactIonization.H`: development adds an enum include to a header the
  feature branch deleted after replacing its physics/product implementation.
  Retain the deletion and the replacement helpers; restoring this header
  would not integrate the new collision architecture.
- `WarpXAlgorithmSelection.H`: both branches append to the scattering enum
  and add neighboring enums. Retain IAA's existing ordinal, append Legendre,
  and keep both ionization-energy-sharing and fusion-table-format enums.
- `parameters.rst`: both branches document additions at the same location.
  Retain both the Legendre table/product-order requirements and the
  proton-impact/PJG particle/fluid documentation.

Review of the automatically merged `ImplicitSolver.cpp` confirms that the
new mass-matrix completion helper is retained unchanged. Prescribed-fluid
current remains outside the field-dependent matrices, is added after
kinetic current scaling, and charge is not resummed during matrix-only
residual evaluations.

## What the development mass-matrix change addresses

Commit `cee2ca2fd` deposits exactly one half of each diagonal 2D mass-matrix
stencil, and completes the other half by the symmetry
`S(i,d) = S(i+d,-d)`. The deposited half satisfies
`di+dj < 0 || (di+dj == 0 && dj <= 0)`; the complementary half reads those
separate components. The new common folding routine also serves 1D.

The independent stock tests compare `S dE` with particle push/current
response, for three particle shapes and both synchronization schemes, with
nonzero magnetic field and periodic/multibox boundaries. All six cases in
each of 1D and 2D pass locally at the existing tolerance. This fixture does
not cover RZ; coupled RZ mass-matrix simulations require separate validation.

The previous acceleration/restart reproducibility failures occurred with
explicit CKC/PSATD solvers. This commit does not modify those paths. It is
therefore not evidence that those failures are fixed. Fresh GPU controls
(`58773853` branch, `58775987` unchanged development) reproduce failures in
all three C++ full/restart comparisons and in independent uninterrupted
repeats. All observed native E/B values at checkpoint step 5 are restored
exactly, before the first resumed step, including averaged fields. For
time-averaged PSATD, the relative Ez differences are 4.668e-9/3.974e-9 for
branch restart/repeat and 4.882e-9/2.450e-9 for stock restart/repeat, against
the original 1e-12 criterion. This is the same class of reproducibility issue
as before the mass-matrix change. The kernel-level amplification remains
unisolated; no assertion tolerance is changed. Every component, particle-ID
comparison and native-state check is retained in
[the fresh acceleration archive](results/development-sync-acceleration.json).

## Initial validation and build environment

The merged source compiles locally in 1D, 2D, RZ, and 3D with MPI, Python,
FFT, and openPMD enabled; QED and embedded boundaries are disabled. It uses
the repository-pinned AMReX `d6d1aa11f38cab90fc592ee193ae574cdc001aaf`
and pyAMReX 26.09. The local build directory is `build-sync`.

Local setup issues were distinguished from physics failures:

- The conda Clang linker cannot read the current macOS SDK's `arm64e.x1`
  stubs. The build uses Apple's `/usr/bin/clang{,++}`.
- Vendored openPMD 0.17.1 fails to compile against the current libc++ in
  `ADIOS2IOHandler`. The compatible installed conda openPMD 0.17.1 is used
  via `WarpX_openpmd_internal=OFF`.
- MPICH 5.0.1's default sockets provider hangs in `MPI_Finalize`, including
  in a minimal mpi4py program with no WarpX import. A sampled stack locates
  the hang in the OFI/libfabric sockets finalization path. `FI_PROVIDER=tcp`
  makes that control initialize and finalize successfully and is used for
  local MPI validation. No simulation code is changed for this workaround.

Before the merge commit, validation passes comprise 64 Python model/archive
tests, independent Gaussian quadrature checks, both Cartesian mass-matrix
unit suites, and eight MCC/fluid/interface/source/restart/load-balance CTest
stages. No physics assertion tolerances or checksum reference files were
changed.

Perlmutter uses the repository's documented module profile and CMake build
recipe, through `perlmutter_stock_build.sh`, in separate branch/stock
`build_pm_gpu_sync` directories. Both include the added 2D geometry needed
for mass-matrix validation. CUDA builds are in progress at this checkpoint;
this report must not be interpreted as completed GPU acceptance.

## Angular-model compatibility follow-up

The shared enum allowed IAA to be selected for fusion or DSMC even though
those kernels do not implement the electron-MCC closure. A zero-step fusion
initialization with IAA returned success before the fix. Explicit guards
now reject IAA in both constructors, before any scattering can fall through
to a different model. A new test checks forward/backward/isotropic success
and IAA/Legendre-without-a-table rejection for fusion and DSMC. The MCC
invalid-input test also checks the development branch's Legendre rejection
for both elastic scattering and ionization. Both CTest stages pass in the
local Release build; their tolerances are unchanged.

## Reproducible complete-suite driver

`run_full_regression.py BUILD OUTPUT` inventories and runs every configured
CTest stage, preserving the dependency ordering through checksum stages.
Checksum failures are reported separately, not hidden by removing those
stages. The build-configuration exclusions include the nonlinear QED tests, collider
diagnostics, beam-beam, and nodal-electrostatic groups requiring QED diagnostics or photon emission
when `WarpX_QED=OFF`. Embedded-boundary injection groups require `WarpX_EB=ON`.
The manifest, raw log, JUnit XML, and classified summary are retained.
An optional MPI launcher runs the Python unit suites inside an allocated
GPU step, with their original CTest environment and working directory.

The stock Perlmutter build driver accepts `WARPX_AUDIT_DIMS` (defaulting to
its previous `1;RZ;3` configuration), and uses CMake's actual `BUILD_TESTING`
option. This campaign sets `WARPX_AUDIT_DIMS='1;2;RZ;3'` and
`WARPX_AUDIT_BUILD_NAME=build_pm_gpu_sync`. The restart-failure driver also
accepts that build-name setting, so the new comparison cannot accidentally
load an older library. Shell syntax and Python lint/compilation checks pass.

The same six mass-matrix identities in each Cartesian geometry also pass
with two local MPI ranks. The GPU driver includes one-, two-, and four-rank
runs, with separate JUnit files per rank.

The complete local run exposed missing `lasy`, `axiprop`, and `dill` test
packages. After installing these test-only packages without replacing the
active MPI/openPMD stack, all six laser-file preparation stages pass.
Perlmutter's missing generator and pytest dependencies are also installed
before GPU execution. Three stock flux-from-embedded-boundary groups are
registered even in an EB-disabled build; their required capability is absent.
The initial local failures are retained (the processes additionally hung in
MPI_Abort). The driver now identifies all nine stages as configuration
exclusions for `WarpX_EB=OFF`, alongside the QED-dependent process and diagnostic stages. The nodal-electrostatic test
requests QED emission and checks the QED chi diagnostic, so removing its
zero-valued Schwinger option alone does not make it a non-QED test; its
original input deck is retained.


The three embedded-boundary injection registrations now use the same
`if(WarpX_EB)` guard as the neighboring EB test directories. The supplemental
SYCL driver accepts `WARPX_AUDIT_SYCL_BUILD_NAME` so this campaign can use
`build_pm_sycl_sync` while retaining the earlier build. SYCL compilation
is distinct from GPU runtime acceptance; Perlmutter has no Intel or AMD GPU.

The QED configuration filter also excludes the nonlinear Breit-Wheeler,
quantum-synchrotron, and Schwinger groups when QED support is absent. The
hybrid Maxwell vacuum-polarization test remains enabled: its field solver
is available independently of the particle-QED build option. The initial
local run retained the unavailable-feature failures; the corrected driver
records the explicit excluded-test list before launching the GPU suite.
The beamsize/virtual-photon groups also request compiled QED functionality
and are listed with those exclusions. The manifests record exact names;
none of the rigid-fluid, PJG, MCC, implicit, or fluid diagnostic/restart
checks is excluded by these capability filters.


## Follow-up diagnostic findings

The 116-stage OpenMP selection completes with 111 passes and five strict
`BeamRelevant` restart-comparison failures. All five are the same near-zero
centroid issue: the uninterrupted run reports exactly zero mean z, while
MPI redistribution gives approximately -1.10436e-18 m. The comparison divides
by its 1e-30 floor. A separate run with every fluid feature disabled also
fails the unchanged centroid comparison. The original diagnostic implementation
was identical to stock development. These failures remain recorded; no assertion
or tolerance has been relaxed. Other field, charge, energy, and momentum
comparisons in these cases agree at roundoff.

Review of that diagnostic identifies a separate dimensional error predating
this branch: its dimensionless Twiss alpha contains an extra speed-of-light
factor. For x=x0+a*s and u=u0+b*(h*s+k*t), with independent equiprobable signs
s,t, the exact alpha is -h/abs(k). The new 64-particle covariance test expects
-2/3 and obtains -199861638.666667 from the old implementation, exactly the
extra factor c. Alpha must be -Cov(x,u)/sqrt(Var(x)Var(u)-Cov(x,u)**2), or
-Cov(x,p)/(m*c*normalized_emittance). This follows the covariance definition
of Twiss parameters in the [USPAS ellipse notes](https://uspas.fnal.gov/materials/18ODU/Fund/some-notes-on-ellipses.html).

A correction removes that factor and places the empty-population guard before
normalizing by total weight. The associated documentation retains normalized
emittance and beta definitions. Exact populated- and empty-beam tests pass
in the 2D/RZ OpenMP build and all four NOACC unit suites. The same covariance
reference fails on untouched development by the factor c, confirming the
error predates this merge. GPU acceptance of the correction is pending.
This correction is independent of the centroid reduction-order sensitivity.


The documented CUDA build now passes for both the merged source (58773767)
and unchanged fork development (58773769). Both revisions pass the independent
1D/2D mass-matrix identities on one, two, and four A100s (58774151/58774153).
Native float-particle physics tests pass on A100 (58774663); the supplemental
SYCL RZ build and standalone kernels compile (58774662). These portability
results precede the diagnostic-only correction above; the corrected diagnostic
is being rebuilt separately before the remaining runtime studies are released.

Two local analysis timeouts were caused by interactive `plt.show()` calls.
The complete-suite driver now sets `MPLBACKEND=Agg` for its subprocesses.
The Compton analysis then passes; the effective-potential analysis reaches
its existing assertion and reports 0.07277108% against a 0.07% bound. The
unchanged stock comparator reproduces those same effective-potential values.
Neither analysis assertion was changed. The EB-dependent ion-extraction
case is also excluded when EB support is not compiled.

For local stock comparisons, a detached checkout at `cee2ca2fd` uses the
same AppleClang, AMReX, pyAMReX, and openPMD as the branch. Installing the
full original build into a local prefix exposed an absolute `/usr/local`
WarpX alias in its install script; the install stopped there and no alias
was created. The already installed dependency packages are sufficient.
Their Python modules needed an explicit local library RPATH before they
could load; those isolated copied modules were adjusted and re-signed.
Setup failures are retained separately from the subsequent physics controls.

GitHub pushes from the activated conda shell initially used conda Git's
credential-cache helper rather than macOS Git's Keychain helper. Using
`/usr/bin/git` fixes the push path. This is a local tooling difference,
not a GitHub authorization or source-state change.

The corrected diagnostic compiles in CUDA (58775901) and SYCL (58775908).
The Perlmutter checkout had a restricted Git fetch refspec: fetching a branch
updated `FETCH_HEAD` without advancing its tracking ref. An explicit branch
refspec and a checked `08124f33a` HEAD resolved this before either rebuild ran.
The subsequent runtime jobs depend on that verified CUDA build.

The first fresh stock acceleration comparison (58773855) timed out on
`nid008341`; its logs report NVML "GPU requires reset" before the first step.
This is hardware execution failure, not checkpoint or physics evidence.
A separate rerun excludes that node and retains the original failed attempt.

The original local helium discharge checks give 9.62% (MCC) and 7.40%
(DSMC) density RMS error against the unchanged 6.50% Turner-profile bound.
The fresh stock comparator gives 3.03% and 5.30%. These are not classified as
stock failures. `discharge_convergence.py` retains the original analysis and
its assertion, records every seed and input/library hash, and varies particle
count, mesh, and timestep independently before joint refinement. Physical
simulation time and averaging time remain fixed. Its results must be assessed
before attributing the difference to sampling or finite-step collision effects.
It does not replace failed individual assertions with an ensemble pass.

`export_reaudit.py` now archives nested JUnit files using relative paths,
including complete-suite and per-rank mass-matrix/diagnostic evidence, without
silently overwriting identically named XML files from different configurations.

The initial joint-refinement driver enlarged the mesh after allocating its
ion-density averaging array, producing a 33-versus-65-element callback error.
The driver now resizes that array consistently and also preserves the DSMC
neutral rethermalization interval in physical time when reducing the timestep.
The affected attempts are retained separately and repeated; this was an audit
driver error, not a change to either simulation's collision implementation.
MPI launches are bounded and their process groups are terminated on timeout.

## Complete local regression disposition

The original complete run attempted 1,141 CTest stages. After accounting for
the disabled QED/EB capabilities, 1,090 stages are eligible. Six missing laser
generator dependencies and the interactive plotting backend were repaired
and their affected stages repeated. The corrected four-geometry Python unit
suites pass. The reconciled result retains 16 physics-analysis failures and
162 platform-dependent checksum failures; the latter are ignored according
to this repository's instructions, without modifying benchmark files.

| Remaining local analysis failures | Investigation |
| --- | --- |
| Two differential-luminosity tests | Untouched stock reproduces the same failures: total luminosity errors 0.2173% and 0.2636%, versus 0.2%. |
| Three focusing/rotated Gaussian tests | Untouched stock also fails their existing slice-width assertions. |
| RZ Ohm modes, 1D Ohm ion-beam instability | Untouched stock reproduces the same amplitude/growth-rate reference failures. |
| Effective-potential electrostatic | Both versions give a maximum density RMS error of 0.07277108%, versus 0.07%. |
| Time-averaged PSATD acceleration restart | Fresh stock/GPU repeat controls above also fail; saved native fields restore exactly. |
| Five fluid diagnostic continuations | Only the near-zero BeamRelevant z centroid fails; the corrected 26-stage diagnostic rerun retains those five failures. A no-fluid OpenMP control changes by 1.767e-18 m. |
| Helium MCC and DSMC discharge | Seed and resolution studies are separate from the failed original assertions; they are not classified as stock failures. |

[The local evidence archive](results/development-sync-local.json) retains all
original complete-suite outcomes, explicit capability exclusions, follow-up
provenance, failure text, stock controls, absolute centroid differences and
the corrected diagnostic/unit runs. The initial stock Python loading failures
are setup attempts; the subsequent `stock-python-controls.xml` provides the
valid physical comparison. The legacy collision assertions are not superseded
by a favorable seed or a differently resolved run.

The first GPU reference batch (58773847) passes all ten standalone C++ tests,
64 Python model/archive checks, independent Gaussian quadrature, external IAA
DCS comparisons, invalid-input guards, stationary probes and Yee/PSATD flux
diagnostic continuations from four ranks to two. Initial integrated power is
exact in both flux cases. Their worst continuation discrepancy is 4.98e-15;
these GPU runs pass the unchanged diagnostic criterion, unlike the local
near-zero centroid cases.

## Helium discharge resolution controls

The local ensemble uses two MPI ranks for both MCC and DSMC, four explicit
WarpX seeds, and the unchanged input and 6.50% analysis bound. Three of four
branch seeds pass for each method; all four stock seeds pass. RMS errors of
the mean density profiles are 5.02%/5.00% for branch MCC/DSMC and 3.50%/3.92%
for stock. These means are descriptive and do not replace individual failures.
The original DSMC CTest uses one rank and remains a separately recorded case.

| Seed 1 control | Branch MCC | Stock MCC | Branch DSMC | Stock DSMC |
| --- | ---: | ---: | ---: | ---: |
| 32 cells, 256 particles/cell | 9.622% | 3.032% | 5.853% | 5.287% |
| Half timestep | 6.647% | 4.210% | 5.879% | 5.747% |
| Four times particles/cell | 4.856% | 3.487% | 5.476% | 4.043% |
| 64 cells, 512 particles/cell, half timestep | 3.741% | 3.876% | 4.070% | 3.648% |

Both versions pass the original density criterion after particle refinement
and joint refinement. Their jointly refined density profiles differ by
1.36% (MCC) and 1.82% (DSMC) in 1D RMS relative to stock
at the interior nodes. This supports sampling/discretization sensitivity in
the coarse regression, rather than a demonstrated converged density bias.
It does not prove exact equality of their stochastic trajectories or turn the
original assertion into a pass. Every attempted and corrected run is retained
in [the discharge archive](results/development-sync-discharge-local.json).

The equation review also makes an existing profile limitation explicit in
the input documentation: a finite rigid train is clipped in laboratory z,
not periodically wrapped. Entry/exit tests use nonperiodic longitudinal
boundaries. The periodic particle-reference benchmarks require the profile
to remain inside the domain; their driver rejects beam transit through a
periodic boundary. No new periodic-beam model is introduced by this merge.

## Implicit collision substep timing

The second temporal audit finds another pre-existing stock error. Implicit
collisions run after the field/particle push, with `cur_time` already at the
end of the PIC step. Adding `i_sub*dt_sub` to that time evaluates later
substeps in the future. An exact control defines zero gas throughout the
simulated interval [0, 1 ps], with gas appearing only after 1.125 ps.
Four implicit substeps incorrectly scatter all 512 test particles on the
branch, and three on untouched development. The different counts reflect
the collision algorithms; the analytical answer is exactly zero in both.

The handler now samples left endpoints for explicit evolution and right
endpoints for implicit evolution within the same elapsed physical interval.
For an implicit step ending at t, the placement is
`t - (N-1-i)*dt/N`, for i=0,...,N-1; N=1 preserves the original time exactly.
The rigid-source integration interval remains independently defined, so this
also corrects time-dependent neutral temperature sampling without changing
its Gaussian space-time integration. Constant gas/temperature benchmarks
are unaffected by this change.

Nine short MPI controls cover future gas at one/four substeps and gas present
only early in the step, for explicit Yee, semi-implicit EM and mass matrices.
The early-gas control prevents accidentally fixing the future evaluation by
sampling every substep at the final endpoint. All nine pass in the rebuilt
local MPI library. The affected regression selection retains only the known
coarse helium-MCC assertion failure. GPU acceptance is in progress; the
pre-fix controls and their logs remain retained.

The repeated four-GPU RZ covariance unit test exposed a separate test setup
error: its manually constructed RZ simulation omitted the common helper's
on-demand arena setting. With repeated initialization, it aborted in
`Arena::Initialize` with an out-of-memory error before the empty-beam case
(uncaptured reproduction 58779101). It now uses the same zero initial arena
reservation and exception settings as the Cartesian helper. Both populated
and empty cases pass on four local MPI ranks. No covariance expectation or
tolerance is changed. The earlier NVML reset failures remain separate
device-execution incidents, with their logs preserved.

The post-fix OpenMP/MPI selection completes 121 stages with 116 passes and
the same five near-zero-centroid comparisons described above. The three
scatter-deposition fixtures explicitly use four OpenMP threads on two MPI
ranks; ordinary CTests retain the repository's one-thread default. The
[local follow-up archive](results/development-sync-followups-local.json)
retains all nine exact timing cases, the 77-stage collision selection, the
121-stage OpenMP selection, and pre-fix branch/stock timing controls.

CUDA acceptance job `58779585` completes successfully: populated/empty exact
covariance tests in 1D, 2D, RZ, and 3D on one and four GPUs, all nine timing
cases, and all ten standalone fluid/PJG/MCC physics tests. CUDA and SYCL
rebuilds `58779552`/`58779553` also succeed. Subsequent performance studies
keep checkout `242559e71` and its compiled libraries fixed for their duration.

An interrupted allocation can leave an empty topology JSON or partial XML.
The evidence exporter now records the filename, parse error, content hash,
length and final bytes explicitly, while retaining other valid records.
Such output is never counted as a successful check. Its regression uses a
mixture of passing, failing and incomplete nested records. Future source
manifests also include `tests/unit`; earlier manifests identify those tests
through their Git revision and JUnit results instead of per-file hashes.

## Distributed studies and retained unsuccessful attempts

The [GPU study archive](results/development-sync-gpu-studies.json) preserves
completed studies and incomplete attempts separately, with each study's
own source/library fingerprint. It includes the following earlier controls:

- `58773847`: standalone physics, 64 Python model/archive tests, Gaussian
  quadrature, measured N2/O2 elastic DCS, input guards, stationary probes and
  nonzero Yee/PSATD Poynting flux. Four-to-two-GPU diagnostic continuations
  pass, with maximum reported normalized difference 4.98e-15.
- `58773852`: 12 particle-deposition comparisons, three joint refinements
  and eight solver/timestep settings pass their unchanged analyses.
- `58773849`: two-node/eight-GPU ownership cases cover both source and
  attachment in all four solver configurations, including four-GPU restarts.
  Coupled six-pair restarts pass for Yee, PSATD and semi-implicit EM. The
  mass-matrix continuation fails only its O-minus population criterion:
  mean resumed-minus-uninterrupted difference -1583.33, standard error 258.74,
  versus the predeclared five-standard-error bound 1293.68. Immediate saved
  fields, particles, counters, primary budgets and charge checks pass. A
  separately scheduled 48-pair follow-up retains this original failure.
- `58773848`: all ten completed coupled solver/subcycling ensembles pass
  population/energy and simultaneous spectral-band checks. Yee at two
  substeps stops because a documentation/tool commit changed the recorded
  checkout revision during the ensemble; the four-substep Yee case was not
  run. The same metadata issue stops `58773850` during four-GPU scaling.
  A fixed-checkout 1/2/4-GPU repeat `58777473` completes successfully.
- `58773851`: fields and equal-electron-work noise ensembles complete, but
  the source ensemble is interrupted by the allocation limit; the coupled
  phase never acquires a usable topology record. This is incomplete
  evidence, not an accepted noise campaign. A longer allocation repeats it.
- `58776972`: independent continuum-error qualification and 16 synchronized
  source/coupled profiles complete. Instrumented timings include profiler
  synchronization and source startup; they are separate from ordinary
  timestep timing ensembles.

The four-seed, two-GPU helium control `58776386` retains one coarse MCC
failure: branch seed 1 has 8.719% RMS error against the original 6.5%
criterion. Branch seeds 2--4 give 3.913%, 3.684%, 4.957%; untouched stock
gives 3.712%, 3.954%, 5.007%, 5.386%. All eight DSMC controls pass
(branch 4.093--5.556%, stock 4.461--6.129%). The planned independent
timestep, particle and joint refinements are separate from these original
assertions, as in the local study.

[The submission archive](results/development-sync-submissions.json) stores
the actual Slurm-retrieved scripts and job descriptions for the final
campaign. Its descriptions were captured before those jobs completed and
are allocation provenance, not test outcomes. Command-line and scheduler
overrides take precedence over the saved script's default directives.

The same-rank counterpart also passes: exact attachment survival and native
ion/electron charge-footprint checks run uninterrupted and through checkpoint
restoration on two local MPI ranks for all four solver configurations. GPU
job `58779962` repeats the eight runs on four GPUs with unchanged rank count
and passes. These supplement the rank-redistribution controls; neither
expects bitwise replay of the stochastic MCC history.

## Focused source review

The equation/time review rechecks Gaussian normalization with the annular
measure, signed axial velocity/current, the order-p charge versus order-(p-1)
cumulative current primitive, cylindrical axis normalization, and time
placement versus source integration. It verifies that the fixed-energy PJG
source uses the same two neighboring incident-energy rows and interpolation
as the particle implementation. The total cross section is the code's
integral of the calibrated N2/O2 PJG spectrum; a user-supplied proton total
cross-section table is not required. Electron MCC rate tables remain separate
inputs; `RBEQ` names its energy-sharing model and `IAA` its angular model.

The execution/lifecycle review rechecks that scatter loops use `amrex::For`
and host/device atomics where necessary, each creation scan owns disjoint
particle output, and fluid-ion footprints use the electron's stored position
and weight. Only fresh increments undergo volume scaling and `SumBoundary`.
Persistent populations, remainders, budgets and quiet counters are mandatory
restart state; analytic caches rebuild without chemistry or self-field replay.
The implicit residual adds prescribed current once, after kinetic scaling,
and never commits chemistry inside a nonlinear/Jacobian evaluation. The
development stencil-folding helper remains unchanged by conflict resolution.
The timing and Twiss fixes above are the additional production corrections
found during these passes; all other findings retain their measured evidence
and limitations instead of changing assertion thresholds.

## Final-library performance controls

Jobs `58779668` and `58779669` complete strong scaling on 1/2/4/8 A100s with
the same compiled libraries from checkout `242559e71`. Four GPUs use a full
node and eight use two nodes; topology records check distinct devices and
a GPU-buffer MPI reduction. Each number below is the mean of three per-run
median ordinary timestep times, after two warm-up steps. Timers synchronize
the device and report the maximum rank duration. The mesh is 256 by 1024;
quiet particle beams use 64 particles/cell.

| Timestep, milliseconds | 1 GPU | 2 GPUs | 4 GPUs | 8 GPUs |
| --- | ---: | ---: | ---: | ---: |
| Coupled, fluid beam / fluid ions | 28.718 | 18.411 | 13.219 | 10.712 |
| Coupled, quiet beam / fluid ions | 91.802 | 49.955 | 36.468 | 22.650 |
| Coupled, fluid beam / frozen ions | 23.351 | 14.287 | 10.855 | 8.955 |
| Coupled, quiet beam / frozen ions | 85.995 | 46.034 | 34.086 | 20.705 |
| Identical electrons, fluid beam | 10.440 | 6.798 | 4.894 | 3.583 |
| Identical electrons, quiet beam | 31.430 | 17.782 | 12.679 | 7.636 |

For the coupled cases with fluid ions, changing the beam representation gives
speedups of 3.20, 2.71, 2.76 and 2.11 respectively. The equal-electron-work
cases use 1,048,576 identical kinetic electrons with chemistry disabled.
These are fixed-resolution comparisons, not universal speedup estimates.
Sparse frozen ions remain faster here because fluid updates communicate mesh
data; their particle storage instead grows with accumulated ion events.

Independent continuum-error qualification in `58779671` retains the
predeclared Er/B-theta 1% and Ez 2% targets, checking six seeds and all saved
outputs on the expanded domain. The fluid case qualifies at 4.554 ms/step.
Quiet particle beams qualify already at 4 particles/cell and 3.808 ms/step;
random particle beams first qualify among the tested resolutions at 64
particles/cell and 4.969 ms/step. Thus the fluid beam removes beam sampling
noise, but is not the fastest representation for every accuracy requirement.
The separate all-channel plasma noise ensemble is still being completed;
beam determinism does not imply zero kinetic-electron sampling noise.

The same-build storage comparison retains 682,112 bytes for four persistent
ion-density fields at both 20 and 80 steps, including guards. Matched frozen
ion particle payload grows from 3.504 to 14.975 MB. Complete fluid/frozen-ion
checkpoints grow from 8.211/11.033 MB to 19.063/33.364 MB because electrons
remain kinetic in both representations. Inclusive checkpoint-step times are
46.59/60.46 ms and 66.18/96.35 ms; these include the physical step and variable
filesystem cost. Density increments and aggregate-charge scratch are excluded
from the persistent-density figure and also scale with mesh size.
The numerical results and confidence intervals are preserved in
[the follow-up GPU archive](results/development-sync-gpu-followups.json).

The final Sphinx build `58780066` succeeds after the timing/boundary input
documentation changes, retaining the previously observed single Doxygen
`warpx::fields::FieldType` lookup warning. Ruff lint and format checks pass
for all ten Python files changed by this integration. No generated stubs or
checksum reference files are updated.

## CUDA memory-check investigation

The first sanitizer job `58777655` completes its three physics programs but
returns the requested error status 99 for each (776, 180 and 124 reported
errors). Printed reports are CUDA API diagnostics: Cray MPI calls
`cuPointerGetAttribute` on ordinary host addresses, and CUDA's extended
logging reports internal symbol probes in `libcublasLt`. The initial print
limit truncates later records, so this run is not memory-check acceptance.

The MPI-only program in `58780012` reproduces the pointer-query reports
without importing WarpX or CuPy. Its GPU-disabled counterpart completes the
host collective but the sanitizer returns 255 because no CUDA API was used;
the follow-up control explicitly allows this intentional host-only case.
The CUDA driver documents `CUDA_ERROR_INVALID_VALUE` for this query on
unregistered pointers, and the sanitizer distinguishes returned API errors
from device memory-access failures. See the
[CUDA pointer-query reference](https://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__UNIFIED.html)
and [sanitizer reporting/suppression reference](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html#cuda-api-error-checking).

The next attempt, `58780113`, fails before running physics: its stack-qualified
suppression does not match, and a test executable named `positive` collides
with the intended output directory. Both setup issues are retained in the
archive. The corrected attempt uses a separate executable name and an exact
API/result-code suppression (`cuPointerGetAttribute`, error 1), with explicit
API reporting and unrestricted device-memory checks. GPU-aware MPI remains
enabled. A deliberate out-of-bounds CUDA write must still produce status 99,
while the suppressed MPI-only control and all three physics programs must
return zero. This acceptance is pending; no physical tolerance is changed.
