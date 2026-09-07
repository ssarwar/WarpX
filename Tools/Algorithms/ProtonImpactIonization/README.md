# Proton-impact ionization reference checks

These tests exercise the actual host/device kinematics and sampling headers independently of
the PIC time step and the empirical cross-section fit. They do not initialize
MPI. Use the AMReX package from the current worktree's build, for example:

```sh
cmake -S Tools/Algorithms/ProtonImpactIonization -B build/proton-impact-physics \
  -DAMReX_DIR="$PWD/build/_deps/fetchedamrex-build/lib/cmake/AMReX" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/proton-impact-physics -j 4
ctest --test-dir build/proton-impact-physics --output-on-failure
```

The checks cover both `float` and `double`, independently of the configured
particle precision:

- The free stationary-electron endpoint is compared with its invariant-energy
  expression for four projectile masses and 361 incident energies each.
- The complete pointlike spin-1/2 Bhabha factor is compared with the invariant
  tree-level Dirac trace, including the endpoint cancellation.
- The three-body molecular endpoint is compared with an independently solved
  recoil-mass quadratic for two neutral masses, three projectile masses, three
  binding energies, and 100 incident energies each. Subthreshold exclusion and
  the available-energy bound are also checked.
- A two-parent analytic spectrum checks the mixed-energy source mean, second
  moment, invariance under parent-order reversal, conditional quantile
  uniformity, and stratification within an aligned parent block.

The free and molecular endpoints have different meanings. At 5 keV a proton
can transfer only about 10.88 eV to a free stationary electron. Molecular
ionization can emit an electron above that value because the residual ion can
recoil. The molecular endpoint only establishes what is allowed; it does not
give the probability of that emission.

These are host execution tests. Linking an accelerator-enabled AMReX package
does not establish GPU execution coverage or CUDA/HIP/SYCL performance.
Cross-section validation, sampler statistics, and accelerator benchmarks must
be reported separately.

An independent host microbenchmark isolates the quantile-generation cost:

```sh
build/proton-impact-physics/benchmark_proton_ionization_sampling
```

It reports the median of seven runs for each precision and a checksum to keep
the calculations observable. The ordered sequence is a timing baseline, not
a correct mixed-parent sampler. The benchmark does not time a full collision
call and has no platform-dependent pass/fail timing threshold.

## Research cross-section references

With NumPy and SciPy installed, run the independent Python checks with:

```sh
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
```

`reference.py` implements the Rudd recommended total and differential curves
and the separate longitudinal/transverse relativistic PWBA integration kernel.
Tests compare the latter against an analytic longitudinal integral, independent
adaptive transverse quadrature, and the dipole limit including both
`ln(gamma^2)` and the `-beta^2` term. The kernel requires a target response; it
is not by itself a molecular SDCS.

`experimental.py` separately preserves the original Table I totals for N2 and
O2 from Crooks and Rudd (1971), and the N2 totals and mean electron energies
from Rudd (1979). These measurements are not silently normalized to the
recommended curves. The 1979 N2 absolute scale was tied to the 1971 data, so
the two data sets must not be treated as independent absolute calibrations.
At 5 keV the 1979 N2 total is 23.8% above the 1985 recommended total; the
original paper also warns of much larger uncertainties for some lowest-energy
measurements. The mean electron energies provide an additional spectral
check, not another measurement of the total yield.

See [RECONSTRUCTION.md](RECONSTRUCTION.md) for the derivations, source
conventions, rejected model trial, and remaining validation requirements.

The current PJG-preserving direction is described in [PJG_REPAIR.md](PJG_REPAIR.md).
It retains the two-Lorentzian structure while separating its soft and hard
limits. These are research candidates, not production replacements. Reproduce
the recorded diagnostics, or repeat the staged fits, with:

```sh
python Tools/Algorithms/ProtonImpactIonization/fit_pjg_repair.py
python Tools/Algorithms/ProtonImpactIonization/fit_pjg_repair.py --refit
```

The original experimental totals and mean electron energies are reported
separately from the recommended curves used in the diagnostic fit. The
nonrelativistic trial rejects incident energies outside 5--4000 keV; it must
not be used as a relativistic bound-tail implementation.
