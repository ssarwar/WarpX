# Perlmutter coupled-physics re-audit

This record covers the re-audit requested after the upstream merge at
`97078e27b`. The stock comparison is the merged upstream revision
`cb5672fae0b9099d2c7631e93ffff3404a553821`, not the older feature-branch base.
All new compilation, physics tests, regression tests and timing measurements
run on Perlmutter. Results remain pending until their logs are recorded below.

## Previously committed records

- [Fluid validation](VALIDATION.md), its compact JSON measurements and figures
  record the implementation assumptions, independent references, convergence,
  CPU/A100 measurements and known coverage limits through `3ac32229a`.
- [Proton model theory](../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst)
  and the [PJG research archive](../ProtonImpactIonization/Research/README.md)
  preserve the equations, calibration, input provenance, superseded approaches
  and uncertainty of the 800 MeV extrapolation.
- [Collision theory](../../../Docs/source/theory/multiphysics/collisions.rst)
  and [MCC checks](../BackgroundMCC/README.md) document the event selector,
  relativistic recoil, RBEQ energy sharing, IAA angular models and attachment.
- [Application examples](../../../Examples/Physics_applications/proton_beam_air/README.md)
  and the parameter reference document text/PICMI configuration and limits.

Raw build trees, checkpoints and most runtime logs were not committed. Compact
results and reproducible drivers were committed. This re-audit will retain
machine-readable acceptance summaries and provenance in Git; large raw outputs
remain in the authorized Perlmutter validation directory.

## Acceptance matrix

| Area | Independent check and regression coverage | Status |
| --- | --- | --- |
| Stock build | Current stock Perlmutter profile, documented CMake and dependency recipes, matching branch/stock configurations | Pending |
| Rigid beam | Gaussian normalization, exact space-time yield, continuity, self-fields, axis/walls, pulse overlap, shape orders | Pending |
| Immobile ions | Identical-event frozen-particle footprints, signed charge transfer, nonnegativity, no repeated synchronization/filtering | Pending |
| Proton impact | PJG N2/O2 spectra/totals, IAA angular closure, bare-ion scaling, tail precision, caps/remainders | Pending |
| Electron MCC | Elastic/excitation IAA DCS, RBEQ SDCS and IAA ionization angles, two-/three-body attachment, recoil and null events | Pending |
| Coupled operation | All channels together, both particle/fluid representations, subcycling and all four solver configurations | Pending |
| Diagnostics/restart | Native/plotfile/openPMD/reduced outputs, new upstream per-species particle counts, immediate restored state and continuation | Pending |
| GPU ownership | One/two/four GPUs on one node, eight GPUs on two nodes, CUDA-aware MPI and changed restart decomposition | Pending |
| Performance/noise | Synchronized strong/weak scaling, matched electron work, particle references and repeat-seed statistics | Pending |
| Stock regressions | Collision, fluid, implicit, diagnostics and restart; reproduce failures on the merged upstream revision | Pending |

The previous coupled timing driver included ionization and attachment but did
not include elastic or excitation. The new combined fixture must include both.
In the input API, the electron SDCS model is named `RBEQ`; `IAA` names the
ionization angular closure and the tabulated elastic/excitation DCS option.
Synthetic tables establish numerical correctness, not experimental air-chemistry
accuracy. No assertion tolerance will be relaxed to accommodate a failure.

## Equation and lifecycle review

The source review checks the following invariants against the independent
tests; runtime results are recorded separately below.

- The Gaussian normalization uses the cylindrical measure `2*pi*r*dr*dz`,
  `sigma_z = abs(v_z)*sigma_t`, and the finite-support Gaussian integral.
  Current and bunch charge are positive input magnitudes. The deposited axial
  current has the sign of `charge*v_z`; clipping by the domain does not change
  the prescribed normalization. `sigma_r` is the Cartesian transverse RMS
  width, while untruncated radial RMS is `sqrt(2)*sigma_r`.
- Yee charge is projected with the particle B-spline shape; current uses its
  cumulative primitive over the physical step. The discrete spatial
  difference of that primitive cancels the discrete temporal charge change.
  RZ axis folding and the Verboncoeur volume are included. PSATD uses its
  cell-centered radial/axial layout and native spectral current correction.
- Primary production integrates `n_b*n_gas*Z_b^2*sigma(E_b)*abs(v_b)` over
  the annular cell and physical collision interval. Constant gas has a
  separable exact-measure quadrature; parser-valued gas uses positive composite
  quadrature. Placement time and integration interval are distinct, including
  the end-of-step collision placement in semi-implicit integration.
- A fractional product yield carries no charge until an electron is emitted.
  A creation cap raises the product weight; it does not discard the yield.
  Emitted electron and ion density use the same stored position, weight and
  deposition shape. Attachment deposits the removed electron's footprint
  before compaction. Rejected recoil events do not update ion density.
- Fresh ion increments are volume-scaled and synchronized once before being
  added to persistent density. Axis mirrors are consumed once, and physical
  wall shape support is retained. The persistent density is never treated as
  a fresh deposition. Only owned charge copies enter subsequent synchronization.
- Semi-implicit chemistry is outside residual and Jacobian evaluations.
  The field-independent prescribed current is added after kinetic current
  scaling, so it enters neither stored kinetic current nor mass matrices.
  Charge is not re-summed or re-filtered in a mass-matrix-only evaluation.
- Electron elastic/excitation kinematics include molecular recoil;
  ionization solves kinetic-plus-binding energy conservation for the chosen
  energy share and angles. Immobile destinations discard ion motion only
  after those kinematics. The rigid beam is externally prescribed, so its
  source energy is an external input rather than a conserved particle loss.
- MCC subcycling repeats each collision operator with the smaller collision
  timestep. Distinct gas/source operators retain WarpX's sequential splitting;
  subcycling alone does not remove the error from their ordering. Global
  timestep, emission weight and sampling refinement are separate studies.
- Checkpoints validate prescribed model/configuration metadata and require
  the density, remainder, budget and quiet-sampling counter files. Analytic
  caches are reconstructed without replaying collisions or initial fields.
  The restart driver checks restored state before stepping, then checks output
  timestamps, cadence and continuation.

The IAA angular closure is a physical approximation, not an exact molecular
double-differential cross section. The heavy-projectile version adapts its
free-electron cosine and conditions the permitted interval. PJG's extrapolation
to 800 MeV and the underlying molecular data retain the uncertainties documented
in the theory/research archive. Numerical acceptance cannot remove those model
uncertainties.

## Build environment

The saved user profile refers to older Python/HDF5 modules. A fresh copy of the
current repository's `perlmutter_gpu_warpx.profile.example` is used with account
`m3748_g`. Existing BLAS++ links CUDA 12 and existing ADIOS2 links an unavailable
HDF5 library, while the stock profile now loads CUDA 13.2/HDF5 1.14.3.9. These
dependencies must be rebuilt with the stock recipes into an isolated validation
prefix. The user's shared installation and home profile are preserved.

The build recipe follows the current
[stock WarpX Perlmutter instructions](https://warpx.readthedocs.io/en/latest/install/hpc/perlmutter.html).
Its installation/build paths are isolated within this checkout; the destructive
cleanup of the shared user prefix in the stock dependency installer is omitted.

The relevant optional features are MPI, Python, FFT and openPMD/HDF5/ADIOS2.
The validation build covers 1D, RZ and 3D, with double fields/particles. QED and
embedded boundaries are outside this simulation's scope. GPU-aware MPI and
rank-to-GPU placement follow the stock Perlmutter batch example. The earlier
single-GPU timings do not establish full-node or multi-node performance.

## Review findings and results

The merged upstream per-species macroparticle-count diagnostic accepted fluid
indices from the extended diagnostic registry and then indexed the kinetic
container with them. The fix returns an exact zero field for fluid species.
Regression coverage includes kinetic cell counts independently histogrammed
from saved particle positions, fluid-only PICMI/text output, mixed-species
plotfile/openPMD output and restart output. Runtime acceptance is pending.

Stock dependency preparation was submitted as CPU build job `58708323`. The
script's shell syntax was checked on Perlmutter before submission.
It completed successfully in 5 minutes 8 seconds. The fresh environment uses
Python 3.12.12 and rebuilt `mpi4py` against Cray MPICH 9.1.0.794. The clean
three-geometry CUDA build was submitted as job `58708497`.

The source-document review confirms Schmalzried Eq. (11.132)'s incident-energy
denominator, and Table 11.4's screening radii 0.6052/0.5677 for N2/O2. The beam
example now distinguishes Roy et al.'s illustrative 1.4 pC CST calculation,
two-sigma support, and rounded experimental beam parameters from this input's
explicit 0.6 A, 25 ps, eight-sigma normalization (37.5994 pC). A Gaussian
transverse RMS width remains a user-supplied assumption.

Release GPU runtime guards are also under review: ordinary device assertions
are compiled out by `NDEBUG`, including checks on parser-valued gas properties
and a user-supplied collision majorant. Invalid-input rejection must be verified
in the same Release configuration used for timings.

The expanded run script uses exclusive GPU nodes, the stock inverse local-rank
GPU mapping and GPU-aware MPI. `gpu_topology.py` records each rank's node/PCI
device and reduces device buffers across the communicator. The ownership phase
uses enough mesh boxes for every GPU and changes the restart decomposition.
Strong scaling holds the full problem fixed. The weak study holds cell spacing
and uniform electron work per rank fixed while extending the axial domain;
the beam is held fixed, so this is an equal-electron-work weak study, not a
weak-scaling claim for the primary source or its Gaussian support.

The baseline profile (`58709463`) completed 16 particle/fluid source and combined
MCC cases on one and four GPUs. In the four-GPU fluid/fluid combined case,
separate full-density validity scans consumed about 0.040 s of 0.218 s for
20 steps. A fused update/validity reduction is implemented; its speedup and
conservation acceptance remain pending. Timings include all MCC channels in
the new fixture, but synthetic DCS are still numerical fixtures.

### Fresh-build and test provenance

- The clean branch build finished in `58709119` after the 30-minute build limit
  interrupted `58708497`; the separate stock build finished in `58709446` after
  `58708830`. This is a scheduler limit, not a compiler failure.
- `58709466`: all seven original double-particle C++ physics tests passed.
  A strict Python archive comparison exposed a platform-dependent one-ulp
  coefficient recomputation. Commit `e2529276d` preserves the archived literals
  and separately verifies the defining expression; all 64 Python tests passed
  without changing the archive tolerance.
- `58710270`: six of seven native float-particle C++ tests passed. The cached
  monoenergetic PJG cross section differed from device evaluation. Row selection,
  interpolation fraction and cached total must be evaluated on the same backend
  as the particle reference, including energies adjacent to row boundaries.
  The expanded strict comparison is pending.
- `58710512`: the FieldProbe repair passed a four-GPU run restarted on two GPUs.
  The exact solution is constant axial B: energy is `B^2*pi*R^2*L/(2*mu0)` and
  the integrated probe is `B*t`. Maximum energy relative error was `1.12e-15`;
  cumulative-probe continuation was exact. Other reduced continuations differed
  by less than `4e-15` relative to column scales.
- `58711837`: the separate stock build reproduced both probe defects (an extra
  initial `dt`, and lost cumulative history on restart). `58712061` also
  reproduced an extra initial/restart Poynting integration interval using a
  nonzero initial radial flux with exact surface power `2*pi*R*L*E_theta*B_z/mu0`.
- `58712062`: the first stock PSATD flux fixture was invalid because its axial
  boxes were no larger than the eight-cell guards. The revised fixture uses
  16-cell boxes; this failure is not evidence of a production PSATD defect.
- `58712279`: the measured N2/O2 DCS fixture passed its angular-moment,
  backward-probability and recoil checks for elastic/excitation scattering,
  spanning 0.1 eV to 1 GeV. The source data are the user's `elmolcs` DCS files;
  their license is retained beside the raw results. A CuPy shutdown warning
  followed successful completion; it did not change the process exit status.

Two validation infrastructure errors were found and corrected. The first
incremental upload preserved source mtimes older than freshly built baseline
objects, so several fixes were not recompiled. All modified source was explicitly
touched on Perlmutter before rebuild; subsequent uploads use checksums without
preserving mtimes. Second, cached CMake MPI flags were passed as executable names
to the new CTest launcher. The recipe now clears both MPI pre- and postflags.
The failed broad run `58710513` therefore establishes no physics acceptance.
The verified rebuild `58711820` exposed real CUDA portability issues in the
new fixes: constructor-local device lambdas and host/device capture-list
differences. Those are being corrected before repeating runtime acceptance.

The clean two-node topology test `58709445` verified eight distinct A100 devices
and device-buffer MPI communication. This is communication coverage, not yet
acceptance of an eight-GPU coupled WarpX simulation. Early ownership runs
`58710514`/`58710515` used stale diagnostic objects and failed; they are not
reported as successful multi-GPU physics tests.

Record subsequent fixed defects, rejected hypotheses, unresolved limitations,
test revisions, compiler/dependency versions, job IDs and measured results here.
