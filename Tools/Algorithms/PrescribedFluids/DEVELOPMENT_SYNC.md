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
therefore not evidence that those failures are fixed. Fresh stock and branch
comparisons are pending; the earlier results remain in `FAILURE_AUDIT.md`.

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
stages. The only build-configuration exclusions are the collider diagnostics
and beam-beam groups that require photon emission when `WarpX_QED=OFF`.
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
exclusions for `WarpX_EB=OFF`, alongside the five QED-dependent stages.
