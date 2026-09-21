# Prescribed-fluid validation and performance studies

The short CTests under `Examples/Tests/collision` and
`Examples/Tests/langmuir_fluids` check continuity, shapes 1–4, exact collision
charge footprints, explicit/implicit solvers, diagnostic totals and checkpoint
restoration before the first resumed step. This directory contains larger
particle-reference studies, separate from routine CI.

Recorded CPU results, convergence limits and outstanding GPU acceptance are in
[VALIDATION.md](VALIDATION.md).

`test_profiles.cpp` compares the analytic projection with an independent million
particle quadrature using WarpX's actual particle shape functions. It also
checks normalization and the Yee continuity equation on CPU or GPU.
`gaussian_reference.py` integrates the Coulomb Green function for a Gaussian in its rest
frame and Lorentz transforms its fields. Its own check uses the closed spherical
Gaussian solution and doubled quadrature order. No production projection code
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
IAA scattering and attachment. **The MCC rate tables are synthetic regression
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
The solver study compares moving electrons and all collision channels across
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

`perlmutter.sbatch` records hardware, modules and revision, runs GPU regressions,
then executes four ensembles on one A100 at a time. It requests one Perlmutter
GPU node for up to 30 minutes; edit the account when using another project.
Run it from the isolated validation checkout after compiling the CUDA build.
The build's CTest MPI launcher should be `srun` with
`MPIEXEC_PREFLAGS='--cpu-bind=cores;--gpus-per-task=1'`.
