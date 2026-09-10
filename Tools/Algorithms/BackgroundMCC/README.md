# Background MCC physics checks

These standalone tests exercise the production host/device kernels against
independent conservation equations and manufactured interpolation tables.
They complement the full PIC collision tests; they do not establish the
experimental accuracy of the RBEQ spectrum or IAA angular closure.

## Build and run

Use an AMReX package built inside this worktree:

```sh
cmake -S Tools/Algorithms/BackgroundMCC -B build/mcc-physics \
  -DAMReX_DIR="$PWD/build/_deps/fetchedamrex-build/lib/cmake/AMReX" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/mcc-physics -j 8
ctest --test-dir build/mcc-physics --output-on-failure
```

Link a CUDA, HIP or SYCL AMReX package to exercise the same device code on
that backend. Also repeat with `AMReX_PARTICLES_PRECISION=SINGLE`: merely
testing float tables in a particle-double build does not cover the native
particle-float kinematics and RBEQ probability conversion.

- `test_kinematics.cpp` checks electron elastic/excitation energy and momentum
  conservation from the recoil-shifted threshold to 1 GeV, for N2/O2 masses
  and forward, transverse and backward draws. It checks ionization through
  100 GeV, five energy shares and all four angular models, and verifies that
  inadmissible near-threshold ionization leaves every output unchanged.
- `test_interpolation.cpp` uses a decreasing cross-section table spanning
  17 orders of magnitude. Both selectors must retain its positive endpoint
  and constant high-energy extrapolation. A separate manufactured RBEQ
  inverse table checks probabilities from `1-1e-8` through `1-1e-11`; forming
  `1-q` after narrowing `q` to float incorrectly loses this tail.

The CTest timeout is 30 seconds per standalone executable. Neither test
uses wall time as a physics acceptance criterion.

## Full PIC regressions

Build WarpX with dimensions `1;3` and Python bindings, then run:

```sh
ctest --test-dir build -R background_mcc \
  -E 'background_mcc_picmi\.' --output-on-failure
```

The focused cases cover process frequencies and ordering, thermal attachment,
relativistic velocities and recoil, RBEQ spectra, bad inputs, and many-channel
product creation. The loose-majorant case increases the bound by 10000 while
keeping the same physical rates, timestep and analysis as the tight-bound
case. Its one-event probability must remain `-expm1(-nu*dt)`, not
`-expm1(-nu_max*dt)*nu/nu_max`. The derivation, thermal-target limitations and
subcycling requirements are in the
[collision theory chapter](../../../Docs/source/theory/multiphysics/collisions.rst).

The existing helium discharge is an additional independent application check:

```sh
ctest --test-dir build -R 'background_mcc_picmi\.(run|analysis)$'
```

It requires the external `warpx-data` repository at the path expected by its
input. This longer test is distinct from the quick manufactured regressions.
Do not adjust its density-error tolerance to accommodate a collision change.

## September 2026 audit

The CPU verification used an Apple M3 Pro, Clang 19, release optimization and
one OpenMP thread. Both standalone executables passed with native double and
float particles. They also passed with AppleClang 21 address/undefined-behavior
sanitizers; the test/kernel sources were instrumented, but linked AMReX was
not. `FI_PROVIDER=tcp` was needed for this host's MPICH shutdown; it is not a
WarpX setting. No CUDA/HIP/SYCL compiler or GPU was available.

The audit corrected a spurious backward excitation branch near threshold,
rest-energy cancellation, unchecked ionization recoil residuals, rare-tail
probability rounding and decreasing-table endpoint cancellation. Failed recoil
draws are now null events and their reserved product slots are compacted before
another collision operator sees them. This can reduce yield in the very narrow
region where a prescribed angle/share is kinematically inadmissible; it does
not recalibrate a tabulated cross section.

Removing the loose-majorant probability bias restored the independent helium
discharge density comparison: 5.35% RMS error versus its unchanged 6% limit.
The pre-correction run had 99.91% RMS error. Reference data and assertion
tolerances were not changed.

## Performance considerations

The selector caches its incident-energy interval through rejection and channel
selection. Equal energy sharing skips unnecessary random draws. The common
ionization recoil case retains three Newton updates, with a checked, bounded
fallback only when needed; accepting an unconverged solve is not an optimization.
Shared direction construction uses AMReX's combined sine/cosine operation.

The corrected event probability needs one `expm1` per attempted physical state,
in particle precision. Invalid ionization products require a reduction and,
only when present, compaction. These correctness costs are explicit; no
unmeasured GPU speedup is claimed. The many-channel PIC tests retain their
catastrophic-regression timing guards. Profile launch/synchronization cost,
register pressure, table locality and fallback frequency on the intended GPU,
using identical physical rates and product counts for comparisons.
