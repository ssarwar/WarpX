# Electron collision tables

The RBEQ exporter shares its host source definition with the WarpX sampler.
Build it inside the WarpX worktree:

```sh
cmake -S Tools/CrossSections -B build/cross-sections
cmake --build build/cross-sections
build/cross-sections/export_rbeq N2 iaa_thesis_2023 N2_ionization.txt
python Examples/Tests/collision/analysis_rbeq_sources.py build/cross-sections/export_rbeq
```

Use `O2` for oxygen and `elmolcs_b8643810` for the archived fitted model.
An optional final `raw` argument exports signed source totals for comparison;
these are not production MCC tables. Production tables sum positive partials,
span the outer-shell threshold to 1 GeV and target 0.02% linear-interpolation
error (with an absolute floor of 1e-8 of the peak). Their metadata are checked
by WarpX. The independent check requires NumPy and SciPy.

See `Docs/source/theory/mcc_iaa_sources.rst` for the source distinctions,
near-threshold continuation, and rotational production gates.

## Thermal rotation

`rotation_reference.py` constructs state-resolved forward kernels, finite-mass
COM detailed-balance reverse kernels, and an explicit unchanged contribution.
`audit_iaa_rotation.py` tests the supplied IAA constructions. Invoke it with
the archived package source on `PYTHONPATH`, the archive via `--archive`, and
an audit directory via `--output`. It records failed gates without producing
a production bundle. Reference temperature zero is an explicit assumption
for this audit; the inclusive measurements do not establish that temperature.
The spectator/Born angular shapes are approximations. Interpreting laboratory
elementary shapes as COM shapes introduces a separate heavy-target angular
approximation; the threshold, detailed balance and recoil factors use finite
target masses. No source uncertainty is hidden by the numerical error targets.

Bundles have three ASCII header lines followed by a little-endian payload:

```text
WARPX_THERMAL_ROTATION_V1
N2 iaa_sudden_spectator J_MAX T_REFERENCE_K
N_ENERGY N_ANGLE N_TRANSITION
```

The arrays are energies (binary64), cosine-bin edges (binary64), transition
pairs `(J_initial,J_final)` (int32), then angular-bin rate coefficients
(binary64, C order `[energy,angle,component]`). Component zero is the unchanged
contribution; the others correspond to the transition pairs. Rates have units
m³/s and precede Boltzmann weighting. Rigid-rotor constants and spin weights
are canonical, rather than taken from rounded table labels. The source's
maximum initial J must converge both the population and transfer moments.
Indistinguishable float32 threshold knots are merged upward; loss labels are
never averaged. The runtime checks finite-mass kinematic accessibility too.

The host rejects negative components and inconsistent inclusive rates. It
constructs the angle marginal and sparse conditional tables; the recycled
alias variate supplies a uniform position within the selected cosine bin.
This represents the tabulated joint histogram exactly, without smoothing a
zero-loss atom or adding inverse-CDF interpolation error. Angular binning and
incident-energy interpolation still require independent convergence checks.
Data arrays and final sampler storage each have a 512 MiB limit per copy.
The source arrays are released after initialization.

`refine_grid` checks rates and absolute first/second transfer moments at
validation temperatures before exporting a production grid. The default
target is 0.1%, with a 1e-4-of-peak floor for insignificant moment tails;
independent midpoint checks require 0.2%. Four-float32-ulp threshold bands
are treated separately because their energy quantization limits relative
accuracy. The first interval uses a `sqrt(E)` mixture to preserve the
finite-elastic-cross-section limit next to zero. The test kernels also check
state-sum moment convergence and thermal power balance at 100, 300 and 1000 K.

In a 1D WarpX test build, CTest registers independent quadrature, state-sum,
detailed-balance, joint-sampler, signed-recoil and complete MCC tests. Use:

```sh
ctest --test-dir build -R 'rotation|rbeq_source' --output-on-failure
```

The portable `test_mcc_rotation` executable runs AMReX device kernels and
reports initialization time, device table bytes, sampling time, and eight
distribution moments. Add `cumulative=1` for the explicit cumulative reference.
The same code builds for CUDA, HIP and SYCL; successful CPU tests do not
constitute GPU validation. Full-timestep comparisons use
`inputs_test_1d_background_mcc_rotation_picmi.py --steps N --particles N`,
with `--cumulative` for the reference. Use several independent seeds and
report variance times computational cost, alongside throughput and memory.
