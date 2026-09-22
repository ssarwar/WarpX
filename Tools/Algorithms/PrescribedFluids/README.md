# Prescribed-fluid validation and performance studies

The short CTests under `Examples/Tests/collision` and
`Examples/Tests/langmuir_fluids` check continuity, shapes 1–4, exact collision
charge footprints, explicit/implicit solvers, diagnostic totals and checkpoint
restoration before the first resumed step. This directory contains larger
particle-reference studies, separate from routine CI.

Recorded CPU/A100 results, convergence limits and backend coverage are in
[VALIDATION.md](VALIDATION.md).
The subsequent upstream-merge review, all-channel MCC comparisons and
distributed GPU acceptance are recorded separately in [REAUDIT.md](REAUDIT.md).

`test_profiles.cpp` compares the analytic projection with an independent million
particle quadrature using WarpX's actual particle shape functions. It also
checks normalization and the Yee continuity equation on CPU or GPU.
Closed Gaussian antiderivatives independently check the annular and space-time
source integrals for both beam directions, overlapping pulses and finite or
unbounded support. The attachment CTests also check exact exponential survival,
with the weighted Bernoulli variance accounting for cylindrical particle weights.
`gaussian_reference.py` integrates the Coulomb Green function for a Gaussian in its rest
frame and Lorentz transforms its fields. Its own check uses the closed spherical
Gaussian solution, its exact Lorentz-contracted moving solution, and doubled
quadrature order. No production projection code
is used by this reference.

## Build and run

Use a WarpX RZ build with Python, MPI and FFT support. Set `PYTHONPATH` to that
build's `lib/site-packages`. NumPy, SciPy, mpi4py and Matplotlib are required;
CUDA/HIP measurements also need CuPy. The scripts use the existing build and do
not modify installed packages or source files.

```bash
python Tools/Algorithms/PrescribedFluids/gaussian_reference.py
python Tools/Algorithms/PrescribedFluids/ensemble.py \
    --suite fields --output build/beam-fields --cells 32 128 \
    --ppc 4 16 64 256 --seeds 8 --steps 20
python Tools/Algorithms/PrescribedFluids/analyze.py build/beam-fields
```

Repeat with `--suite source`, `--suite coupled` and `--suite push`.
`source` freezes emitted electrons to isolate the integrated PJG yield and
spectrum. `coupled` evolves electrons with PJG production, RBEQ energy sharing,
IAA ionization angles and attachment. Add `--mcc-all` to include tabulated IAA
elastic/excitation scattering and three-body attachment in the same run.
**The MCC rate tables are synthetic regression
fixtures; these runs compare numerical representations, not measured air
chemistry.** Both suites compare fluid ions with identical frozen kinetic ions.
`benchmark.py --ions thermal` provides a moving-ion displacement reference.
`push` uses the same deterministic electron population in every run, without
chemistry, to isolate the cost of the beam representation at equal electron
work. Particle beams are ballistic and do not gather fields.

`--solver` accepts `Yee`, `PSATD`, `semi_implicit_em` and `semi_implicit_mm`.
Semi-implicit particle references use charge-conserving Villasenor deposition,
including mass matrices, to match the prescribed beam's continuity equation.
Use `--implicit-deposition direct` to study the native direct-deposition option;
its current shape differs at finite mesh spacing and requires mesh convergence.
Use `--dt`, `--cells`, `--weight`, `--cap`, `--source-resolution` and
`--subcycles` independently for convergence. Additional arguments passed to
`ensemble.py` reach each individual benchmark. Quiet particle counts per cell
must be perfect squares. No beam wraps around the periodic longitudinal domain
in these comparison runs; the driver rejects runs long enough to do so.

`convergence.py` writes and runs five fixed study manifests. Analyze each with
`analyze_convergence.py` using the same study name and output directory:

```bash
python Tools/Algorithms/PrescribedFluids/convergence.py deposition --output build/deposition-study
python Tools/Algorithms/PrescribedFluids/analyze_convergence.py deposition build/deposition-study
python Tools/Algorithms/PrescribedFluids/convergence.py source --output build/source-study
python Tools/Algorithms/PrescribedFluids/analyze_convergence.py source build/source-study
python Tools/Algorithms/PrescribedFluids/convergence.py joint --output build/joint-study
python Tools/Algorithms/PrescribedFluids/analyze_convergence.py joint build/joint-study
python Tools/Algorithms/PrescribedFluids/convergence.py solvers --output build/solver-study
python Tools/Algorithms/PrescribedFluids/analyze_convergence.py solvers build/solver-study
python Tools/Algorithms/PrescribedFluids/convergence.py continuum --output build/continuum-study
python Tools/Algorithms/PrescribedFluids/analyze_convergence.py continuum build/continuum-study
```

The source study varies one parameter at a time over six seeds and includes
moving thermal ions. At fixed product weight, mesh refinement increases the
number of cells retaining fractional yields. Include pending production in the
primary-yield budget, then reduce the weight to converge the emission delay and
secondary chemistry. A finer mesh alone does not remove that sampling error.
The joint study reduces product weight by eight when halving each mesh spacing,
so the global un-emitted fraction tends to zero during mesh refinement.
The solver study compares moving electrons with ionization and attachment across
Yee, PSATD and both semi-implicit configurations, at two resolved timesteps.

For a mesh study, hold domain extents fixed and double both cell counts. For a
domain study, double `--radial-sigmas`, `--longitudinal-sigmas` and cell counts
together. The finite-domain initial self-field solve has zero potential at the outer
radius and periodic axial boundaries in **all** representations. Yee retains a
conducting outer wall during evolution; PSATD uses its native radial boundary. Its axial
electric field can differ substantially from the unbounded continuum result;
domain convergence must accompany comparison with `gaussian_reference.py`.

`--dry-run` writes the complete command manifest before execution. `--resume`
skips completed runs only if that manifest is unchanged. Each case saves its
configuration, revision, backend, MPI count, timings, populations, native field
snapshots, and electron spectrum. `analyze.py` writes `summary.json` and standalone
noise/cost plots. Preserve the individual JSON/NPZ files to permit reanalysis.

## Interpreting results

- Density and field RMS errors compare with a fluid beam on the same mesh.
  They include particle quadrature error; they are not continuum error estimates.
  Seed noise is computed around each representation's ensemble mean.
- Independent unbounded self-field errors are reported separately in the core
  region `r < 3 sigma_r`, `|z| < 3 sigma_z`, at initialization.
- Primary yield includes the fractional pending population. Compare it against
  the independent PJG quadrature using the existing table accuracy bound of
  `1e-3`, plus independently established spatial/temporal quadrature error.
- Frozen fluid/particle ion destinations should agree to accumulated floating
  point roundoff when supplied identical events. Charge footprints and rejected
  events are checked directly by the short CTests.
- Means include two-sided 99% Student confidence intervals across seeds. Spectrum
  and rare-tail comparisons require more samples than total-yield checks.
  Compare seed distributions and convergence; do not apply an IID count-error
  formula to correlated quiet samples or unequal particle weights.
- Timings exclude Python measurements, synchronize the device and MPI ranks,
  and omit two warm-up steps. Report the actual electron work as well as wall
  time. `--profile` enables synchronized TinyProfiler instrumentation for kernel
  attribution; use separate uninstrumented runs for overall timings.
- `--checkpoint` records the final checkpoint's bytes and the inclusive
  checkpoint-step time. This is not an isolated I/O timing. Initialization can
  include an initial checkpoint. Particle payload and density bytes are logical
  storage, distinct from reserved device memory and allocator capacity.
- An equal-error speedup requires a stated error target and a converged
  reference. An equal-particle-count timing alone does not establish that result.

## Current Perlmutter re-audit

`perlmutter_stock_build.sh` uses the repository's stock Perlmutter GPU profile
and documented dependency/CMake recipes. It isolates rebuilt dependencies and
the Python environment under `build/reaudit-2026-09-21/software`, with only the
account and installation path changed in the profile. Run `prepare`, then
`configure`, `build`, `physics`, `float` and `test-data` as separate phases on
Perlmutter. Compilation does not require a GPU. `stock` builds the selected
upstream comparator in a separate checkout. The environment variables
`WARPX_AUDIT_BUILD_NAME` and `WARPX_AUDIT_STOCK_REVISION` select an additional
build directory and comparator without replacing an active build.

`perlmutter_reaudit.sbatch` requests exclusive GPU nodes, enables GPU-aware MPI
and uses the stock local-rank GPU placement through `perlmutter_gpu_rank.sh`.
Its CTest launcher is `perlmutter_mpiexec.sh`; configure with empty
`MPIEXEC_PREFLAGS` and `MPIEXEC_POSTFLAGS`. Every study records its command
manifest, source-file hashes, loaded Python-library hashes, build configuration,
software versions and rank-to-device topology. A one-node allocation has four
A100s, but the `noise` phase deliberately runs one rank/GPU at a time.

```bash
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch physics
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch regression
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch coupled_restart semi_implicit_em
sbatch --nodes=2 Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch ownership
sbatch --nodes=2 Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch coupled_restart PSATD
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch noise coupled
sbatch Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch scaling
sbatch --nodes=2 Tools/Algorithms/PrescribedFluids/perlmutter_reaudit.sbatch scaling
```

The `coupled` phase accepts each of the four solver names and compares all four
fluid/particle beam-ion combinations at subcycle counts 1, 2 and 4. `dcs` uses
the measured N2/O2 reference DCS files; combined representation comparisons use
the synthetic rates described above. `ownership` checks source/attachment on
four or eight GPUs and restores on half as many. `coupled_restart` verifies
restored fields, particles and source state before stepping, then applies paired
ensemble bounds to stochastic continuation. `scaling` holds the problem fixed
at 1/2/4 GPUs or 8 GPUs on two nodes. `weak` fixes electron work per rank while
extending the axial domain; it does not scale the prescribed Gaussian beam.
Use `convergence` followed by one of the five study names above for independent
mesh, timestep, weight, cap, sampling and subcycling studies.

`perlmutter_reaudit_suite.sbatch` groups finite studies into one allocation and
records each exit status, continuing independent studies after a failure. It
accepts `convergence` (joint, solver and deposition sweeps), `noise` (all four
modes), `restart` (all four solvers), `scaling` (strong/weak scaling and one-node
storage), `physics` (native float and measured DCS), or `coupled_mm`. For example:

```bash
sbatch --export=ALL,WARPX_AUDIT_BUILD_NAME=build_pm_gpu_latest \
    --nodes=2 Tools/Algorithms/PrescribedFluids/perlmutter_reaudit_suite.sbatch restart
```

The native-float physics build uses `<build-name>_float`, with a separate AMReX
build, so dependency updates cannot replace an active validation library.
`perlmutter_docs_build.sh` renders Sphinx/PICMI documentation in an isolated
environment. `perlmutter_sycl_build.sh` is a supplementary oneAPI compilation
check based on the repository's Intel CI and Aurora BLAS++/LAPACK++ recipes;
it is separate from the stock Perlmutter CUDA build. A successful compilation
does not establish SYCL GPU runtime acceptance on Perlmutter.

`export_reaudit.py` archives compact job, test, configuration and measurement
records while retaining failures. `report_state_difference.py` quantifies native
snapshot discrepancies without altering an acceptance test or its tolerance.
Large checkpoints and arrays remain in the validation directory.

The `equal_error` phase compares fluid, quiet and random beams on the expanded
256 by 1024 mesh (`r_max = 32 sigma_r`, `|z_max| = 48 sigma_z`). It uses six seeds
and 4, 16, 64 and 256 particles per cell. `analyze_equal_error.py` compares the
initial and final fields with the independently integrated translating Gaussian
solution in `r < 3 sigma_r`, `|z-v*t| < 3 sigma_z`. Its predefined targets are
1% cylindrical RMS error for Er/Btheta and 2% for Ez, chosen from the separate
domain-convergence study. A representation qualifies when the upper end of its
99% interval for the worse of the two outputs meets every target. The report
retains unsuccessful particle resolutions and compares the fastest qualifying
particle setting with the fluid case. This measures matched vacuum-field
accuracy, not matched nonlinear plasma accuracy. Existing physics assertions
are unchanged.

## Historical one-GPU measurements

`perlmutter.sbatch` records hardware, modules, Python packages, the CMake cache,
revision and any tracked source changes, runs GPU regressions,
then executes four ensembles on one A100 at a time. It requests one Perlmutter
GPU node for up to 30 minutes; edit the account when using another project.
Run it from the isolated validation checkout after compiling the CUDA build.
Build the standalone tests in `build/fluid-tests` and `build/pjg` against that
build's AMReX package. To fit the debug queue's time limit, submit individual
phases, for example `sbatch Tools/Algorithms/PrescribedFluids/perlmutter.sbatch tests`
and then `fields`, `push`, `source`, `coupled`, `storage` or `profile`. Omitting
the phase runs everything. The storage phase compares checkpoint cost at two
durations; the profile phase measures source cost separately from timing ensembles.
The test phase collects failures from each independent group before returning
an unsuccessful status; a failed group still prevents subsequent timing phases
in an `all` run.
The build's CTest MPI launcher should be `srun` with
`MPIEXEC_PREFLAGS='--cpu-bind=cores;--gpus-per-task=1'`.

TinyProfiler's inclusive source tables include first-call setup even though
the separate timestep statistic omits two warm-up steps. To distinguish CUDA
module loading from source execution, repeat the profiling phase with
`sbatch --export=ALL,CUDA_MODULE_LOADING=EAGER Tools/Algorithms/PrescribedFluids/perlmutter.sbatch profile`.
The script records this setting; keep both profiles and the uninstrumented
ensembles. Eager loading is a profiling control, not a required simulation setting.

`export_report.py` produces the compact JSON record and combined noise/cost plot.
Gather the phase directories into the layout described in its help, including
the acceptance CMake cache and optional `profile-eager/profile/` directory:

```bash
python Tools/Algorithms/PrescribedFluids/export_report.py build/fluid-comparisons/a100 \
    --output Tools/Algorithms/PrescribedFluids/results/a100-2026-09-21 \
    --date 2026-09-21 --hardware 'Perlmutter A100-SXM4-40GB' \
    --compiler 'CUDA 13.2, GCC 13.2'
```

Use the same MPI library for WarpX and `mpi4py`. In the recorded CUDA 13.2
environment, Cray Python 3.11.7's bundled `mpi4py` loads MPICH 9.0.1, while the
compiler wrapper links WarpX with MPICH 9.1.0. Mixing them stalls AMReX
initialization. The validation virtual environment uses `mpi4py` 4.1.2 built
from source with `MPICC=cc` and `MPI4PY_BUILD_MPICC=cc` after loading the build's
compiler and MPI modules. Verify `MPI.Get_library_version()` and the linked
libraries before running the suite. NumPy 2.4.6 and SciPy 1.17.1 are installed
together in that environment; the older system SciPy has an incompatible ABI.

The RZ FFT dependencies must also match the accelerator toolchain. These runs
use BLAS++ and LAPACK++ v2024.05.31 built for CUDA 13.2 in
`build/backend-deps/cuda13`; the batch script selects that library directory
when it exists. The Python package list and CMake cache saved with each study
record the precise environment used for its measurements.

The standalone SYCL complex-transform check and its build instructions are in
[SpectralRZ](../SpectralRZ/README.md). It compares the production strided FFT with
a direct long-double DFT, including all complex modes and changes of queue.
Compilation alone is not runtime acceptance: the test needs a SYCL GPU that
AMReX and oneMKL both support.
