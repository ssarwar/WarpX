# Proton-impact ionization reference checks

These tests exercise the actual host/device kinematics header independently of
the PIC time step and the empirical cross-section fit. They do not initialize
MPI. Use the AMReX package from the current worktree's build, for example:

```sh
cmake -S Tools/Algorithms/ProtonImpactIonization -B build/proton-impact-physics \
  -DAMReX_DIR="$PWD/build/_deps/fetchedamrex-build/lib/cmake/AMReX"
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

The free and molecular endpoints have different meanings. At 5 keV a proton
can transfer only about 10.88 eV to a free stationary electron. Molecular
ionization can emit an electron above that value because the residual ion can
recoil. The molecular endpoint only establishes what is allowed; it does not
give the probability of that emission.

These are host execution tests. Linking an accelerator-enabled AMReX package
does not establish GPU execution coverage or CUDA/HIP/SYCL performance.
Cross-section validation, sampler statistics, and accelerator benchmarks must
be reported separately.
