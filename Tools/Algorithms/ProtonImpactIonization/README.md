# Proton-impact ionization reference and validation

The production model is the calibrated PJG-type inclusive electron source for
N2 and O2. Its total is the integral of its SDCS. The complete equations,
original PJG comparison (including the 1977 erratum), fitted coefficients,
data provenance, residuals, stopping powers and limitations are in the
[theory chapter](../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst).
The original printed model is retained only as an offline comparison; there
is no runtime selector for a superseded calibration.

The source emits effective electron/singly charged molecular-ion pairs.
It does not resolve exclusive single ionization, fragmentation or capture.
The beam is rigid and ions are neutral-thermal. The proton incident range is
5 keV–10 GeV, with no lower cut on secondary-electron energy. The angular
closure is outside this calibration.

## Independent Python reference

Install NumPy and SciPy, then run from the repository root:

```sh
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
```

- `calibrated_pjg.py` exposes the frozen SDCS and integrated total, in
  cm2/eV and cm2 per molecule.
- `pjg_model.py` contains the parameterized positive Lorentzian model and
  molecular gate. `target_parameters.py` preserves the printed Table III.
- `original_pjg.py` evaluates printed proton PJG with the erratum, preserving
  signed values rather than clipping negative cross sections.
- `reference.py` independently implements Rudd recommendations and a
  relativistic longitudinal/transverse PWBA kernel. A molecular response
  is needed to turn the latter into a cross section.
- `experimental.py` and `source_datasets.py` retain tabulated and digitized
  data, PDF digests, axes, units and normalization conventions.
- `pjg_moments.py`, `pjg_properties.py` and `optical_reference.py` provide
  moment, stopping, optical-response and interpolation checks.
  `pstar_reference.json` records the NIST electronic-stopping rows.

The tests cover independent adaptive versus segmented quadrature, free and
molecular limits, positivity, moment inequalities, hard-secondary matching,
the relativistic dipole logarithm, data conversions and source semantics.
Numerical tolerances are not experimental uncertainty estimates.

## Reproduce the calibration and figures

The numerical proton constraints are included in the Python sources.
The two external optical inputs must be obtained separately; PDFs are not
redistributed with WarpX. Install matplotlib and pypdf for this workflow.

```sh
mkdir -p build/pjg-data
curl -L --fail -o build/pjg-data/NIFS-DATA-109.pdf \
  https://nifs-repository.repo.nii.ac.jp/record/11706/files/NIFS-DATA-109.pdf
curl -L --fail -o build/pjg-data/leiden-o2.txt \
  https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/O2/O2.txt
python Tools/Algorithms/ProtonImpactIonization/extract_nifs_oscillators.py \
  build/pjg-data/NIFS-DATA-109.pdf --output build/pjg-data/nifs-optical.json
python Tools/Algorithms/ProtonImpactIonization/fit_pjg.py \
  --nifs build/pjg-data/nifs-optical.json \
  --o2-leiden build/pjg-data/leiden-o2.txt \
  --output build/pjg-data/fit-results.json \
  --figures build/pjg-data/figures
```

The fit compares six-parameter and reduced five-parameter models for both
gases, including an alternate starting point. Its output distinguishes
measurements, authors' fitted totals, recommended SDCS curves and evaluated
optical data. It reports absolute and shape residuals, mean/second moments,
quantiles, above-free tails, optical comparisons and NIST stopping ratios.
No production coefficients are rewritten automatically.

Add `--frozen` to evaluate and plot the committed coefficients instead of
reoptimizing them. An alternate-start optimization remains a validation
check, without changing those coefficients. Add `--read-results` to redraw
an existing result without recomputing it.

The extraction used in the calibration has these SHA-256 digests:

| Input | SHA-256 |
| --- | --- |
| NIFS-DATA-109 PDF | `432332e8835a00c1e0fcee32d1cc8a85d4a983f8841c799abd832f9e19af76d7` |
| NIFS extracted JSON | `46da92de977acbd605a0ba6faf3296e685b1140686c67b2898f6323a0571b87a` |
| Leiden O2 text | `3787a1dd6f474f6c0bc691dfb7dfbb3a76901afdde892826c36b4b098a15d72c` |

The JSON digest includes extraction metadata. For a refreshed download,
compare numerical rows and source provenance as well as the digest.
Newer repository dates do not imply new measurements. The theory chapter
explains the NIFS absorption versus Leiden ionization distinction and the
channel-resolved checks against Gallagher (1988) and Mahla–Mehnen (2025).

## Production C++ checks

Use the AMReX package from this worktree's build:

```sh
cmake -S Tools/Algorithms/ProtonImpactIonization -B build/proton-impact-physics \
  -DAMReX_DIR="$PWD/build/_deps/fetchedamrex-build/lib/cmake/AMReX" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/proton-impact-physics -j 8
ctest --test-dir build/proton-impact-physics --output-on-failure
```

The kinematics test compares free Tmax with the invariant expression,
the complete Bhabha factor with an independent Dirac trace, and molecular
endpoints with an independent recoil-mass quadratic. It covers float/double,
several projectile masses and threshold/ultrarelativistic limits.

The sampling test uses analytic spectra with different parent energies to
check conditional quantile uniformity, moments, stratification and
parent-order invariance. An ordered energy sequence paired with ordered
parents is not a valid mixed-energy sampler.
It also checks fixed-point angular/thermal phases at large sequence indices,
including float precision and wraparound near 2^32. Floating multiplication
followed by a fractional part loses these phases at large indices.

The angular test runs the production closure on the device for float/double,
5 keV–10 GeV protons, zero/soft/hard/above-free-endpoint secondaries and
several bindings. It checks the conditional CDF, first/second angular moments,
unit directions, rotation of the incident axis and the independent free
projectile mass-shell equation. Conditioning the uniform IAA-inspired interval
on physical cosines removes the former clipping-induced point mass at the
forward direction (25% at and above free Tmax for nonzero binding). This
change does not alter the calibrated SDCS or its parameters.

The model test initializes AMReX, builds the actual production tables and
runs the executor through `amrex::ParallelFor`. It compares host SDCS and
moments with independent Python/Simpson fixtures and checks float/double
table interpolation, monotone quantiles, molecular support and invalid input.
The fixture can be regenerated with:

```sh
python Tools/Algorithms/ProtonImpactIonization/generate_reference.py
```

The double-build test converts its actual tables to both precisions.
A particle-single-precision WarpX build is still needed for end-to-end
single-precision source coverage. Link an accelerator-enabled AMReX package
to run the device test on CUDA/HIP/SYCL; a CPU run does not establish GPU
correctness or performance.

Full PIC tests include both gases at 5 and 50 keV, mixed 50/500-keV parents,
fractional-weight carryover, and capped many-cell production at 800 MeV:

```sh
cmake --build build --target pyWarpX_python_sources
ctest --test-dir build -R proton_impact_ionization --output-on-failure
```

The first command refreshes the build-tree Python package; building only a
selected C++ or Python-extension target can leave an older PICMI wrapper there.

The tests check represented yield, paired weights/positions, unchanged beams,
electron-energy distributions, above-free tails, thermal-ion velocities
and bounded product counts. Added cases check 1e-8 K neutrals (constant and
parser inputs), rejection of projectile/product species aliasing, and the
emitted-weight-plus-remainder budget over seven steps. The latter compares
coarse/fine product weights, a density pulse that switches off, and bare-alpha
Z^2 scaling at fixed speed, allowing for the finite-projectile-mass correction.
The angular analysis checks conditional CDFs and second moments as well as
the mean direction. The checkpoint/restart pair reuses the eight budget cases
and restarts after step two, before either coarse-weight source has emitted a
pair. It checks that both gases restore their positive fractional remainders
before advancing, then compares final emitted weights, remainders and counts
with the uninterrupted run. It does not require subsequent random angles to
match across restarts or compute backends.

## Performance and numerical results

Separate startup from steady-state sampling:

```sh
build/proton-impact-physics/benchmark_proton_ionization_model
build/proton-impact-physics/benchmark_proton_ionization_sampling
```

The model benchmark synchronizes the execution stream, runs one warmup and
seven timed passes over 2^20 events with interleaved incident energies,
and reports the median/range plus a checksum. It includes total lookup,
energy/binding sampling, input-energy generation and output stores.
The quantile-only benchmark is a separate host diagnostic; its ordered
baseline is not a valid source algorithm. It also compares floating and
fixed-point angular/thermal phases at large indices and prints their second
moments. The floating float sequence can collapse to zero; the fixed-point
sequence retains the uniform second moment of 1/3 at comparable host cost.

The September 2026 Apple M3 Pro CPU run used Clang 19, release optimization,
one OpenMP thread, AMReX particle-double precision, NumPy 2.5.1 and SciPy 1.18.0.
The measured model-kernel medians were 22.85/29.09 ns per event for N2
(float/double) and 22.68/28.10 ns for O2. These are not GPU timings or a
controlled speedup comparison with the superseded implementation.

For the 512-cell PIC smoke cases (32768 beam particles, five steps and
20480 capped pairs), five repeated runs gave median startup times of 0.364 s
for N2 and 0.378 s for O2. The five timesteps took 0.00920 s and 0.00902 s,
respectively, or 17.8/18.2 million parent-particle collision calls per second.
The timestep ranges were 0.00880–0.00961 s and 0.00878–0.00937 s. These cases disable
particle pushing, gathering and deposition to isolate the source workflow;
they are not full electromagnetic-simulation throughput measurements.

The source audit hoists the neutral thermal speed out of the product loop,
reuses the selected parent's frame and free endpoint, bypasses parsers for
constant backgrounds, and skips product work in empty-source tiles. Counts
are scanned in 64 bits before checking particle-tile indexing limits.
Fixed-point angular/thermal phases remove a float accuracy failure at roughly
the same host cost (about 0.75–0.78 ns per phase in the microbenchmark).
These measurements do not isolate a causal speedup for each optimization.

Across 49 incident energies and both table precisions, maximum relative
errors were below 0.030% (total), 0.013% (mean), 0.077% (second moment)
and 0.0013% (binding). Numerical acceptance budgets are 0.1%, 0.1%, 0.2%
and 0.1%, respectively. Each target stores 525056 scalar entries:
4.20 MB in double or 2.10 MB in float. There is no event-local quadrature,
inverse-CDF search, rejection loop or allocation.

The many-cell performance regression is a catastrophic-regression guard,
not a portable throughput requirement: it bounds initialization/timestep
time and requires the configured capped product count. Measure occupancy,
table locality and cell/product-count variation on the intended GPU.

For local MPICH builds whose libfabric sockets provider stalls at shutdown,
`FI_PROVIDER=tcp` was required in this environment. This is a test-host
setting, not a WarpX input or a model change.

All four standalone tests also passed with native single-precision particles,
and with AppleClang 21 and
`-fsanitize=address,undefined -fno-omit-frame-pointer`. The model and test
sources were instrumented; the linked release AMReX library was not.
The conda Clang 19 AddressSanitizer runtime stalled in shadow-memory
initialization before `main` on this macOS host, so that run is not counted
as coverage. No CUDA/HIP/SYCL compiler or GPU was available for this run.
