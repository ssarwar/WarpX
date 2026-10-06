# Electron collision tables

For production N₂/O₂ rotational scattering, start with the
[complete rotational DCS specification](../../Docs/source/theory/multiphysics/rotational_scattering/rotational_dcs.rst)
and the [continuing development guide](../../Docs/source/theory/multiphysics/rotational_scattering/development.rst).
They record the physics, data provenance, implementation map, reproducible
commands and validation requirements needed in a fresh task or checkout.
The tools below also cover ionization and alternative rotational models.


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

See `Docs/source/theory/multiphysics/mcc_iaa_sources.rst` for the source distinctions,
near-threshold continuation, and rotational production gates.

## Thermal rotation

Physical cross sections and bundles are prepared offline and stored in
`warpx-data`. WarpX reads existing files; it never invokes elmolcs or a
cross-section generator. Initialization applies the chosen Boltzmann populations
and packs shared in-memory sampling tables. This supports arbitrary fixed
rotational temperatures without generating physical data during a simulation.

`rotation_reference.py` constructs integral forward/reverse rates in the
heavy-target limit, with relativistic electron phase space and canonical
rigid-rotor level spacings as thresholds. `export_elmolcs.py` converts real
source tables; `test_elmolcs.py` independently checks the conversion and
prepares sampler references in the build directory.

The angle comes from the existing elastic DCS and is independent of the
unchanged/excitation/de-excitation draw. An energy-row mixture followed by one
alias lookup selects the internal change. There is no angular-accessibility
search, rejection, rotational-state loop, or special function on the device.
A subthreshold excitation caused by grid rounding becomes an unchanged event;
the other channels are not renormalized. Loss labels retain double precision.
The optional cumulative reference uses the same probabilities. Only that
reference stores a cumulative table; the alias path omits it.

Bundles have three ASCII header lines and a little-endian payload:

```text
WARPX_THERMAL_ROTATION_RATES
N2 elastic_dcs J_MAX T_REFERENCE_K
N_ENERGY N_TRANSITION
```

Arrays contain energies (binary64), transition pairs `(J_initial,J_final)`
(int32), and integral rate coefficients (binary64, C order
`[energy,component]`). Component zero is unchanged. Other components correspond
to transition pairs, before Boltzmann weighting. The rate format permits `J_initial =
J_final`: these state-dependent unchanged rates are summed to one unchanged
outcome on the host. Older scalar-unchanged input remains readable. Rates have units m³/s.
There are no angular bins. The thresholds and source-rate detailed balance omit the molecular recoil shift.
The state sum and rate/transfer moments must converge. Input and sampler
storage each have a 512 MiB limit per copy.

The signed recoil solve preserves the sampled angle. For excitation in the
small energy-only band `0 <= E-loss <= 2*m_e/M*E`, it neglects recoil energy
in the electron update, setting `E_out = E-loss`. The virtual neutral receives
the momentum difference. The energy defect is bounded by `2*m_e/M*E` for
N2/O2 and is explicitly tested. Outside that band, two-body recoil conserves
four-momentum. Unchanged outcomes reuse ordinary elastic recoil. Thus nominal
thresholds and independent angles do not require angle projection or rejection.

The source rates satisfy heavy-target integral detailed balance. The selected
energy-dependent elastic angular shape does not enforce differential detailed
balance. This approximation and the threshold continuation are separate from
numerical interpolation accuracy.

`refine_grid` checks rates and absolute first/second transfer moments. Its
0.1% target uses a 1e-4-of-peak floor for insignificant tails; independent
midpoints must stay below 0.2%. Float32 threshold bands are tested separately.
Analytic bundles check population and moment convergence, thermal power
balance at 100/300/1000 K, thresholds and angular independence.

Export real source cross sections offline with elmolcs on `PYTHONPATH`:

```sh
python Tools/CrossSections/export_elmolcs.py --output /path/to/warpx-data/MCC_cross_sections
```

The N2 and O2 `IAA` directories contain source cross sections and short usage
notes. Keep synthetic fixtures, validation output and candidate bundles in
the build directory. N2's residual elastic and elementary rotational tables
are incompatible near resonance; a combined production family needs an
explicit source reconciliation. `--test-bundles` prepares source-driven
software tests, including all unchanged `J->J` channels. It does not resolve
that physical normalization problem. Tests and simulations only load the
prepared files; they never evaluate elmolcs or generate source data at runtime.

```sh
python Tools/CrossSections/export_elmolcs.py --output build/elmolcs-data --test-bundles
python Tools/CrossSections/test_elmolcs.py --data-dir build/elmolcs-data --output build/elmolcs-checks
cmake -S . -B build -DWarpX_ELMOLCS_TEST_DATA="$PWD/build/elmolcs-data"
ctest --test-dir build -R 'elmolcs|rotation|rbeq_source|secondary_angles' --output-on-failure
python Tools/CrossSections/benchmark_rotation.py --elmolcs --data-dir build/elmolcs-data --output build/rotation-benchmark --repeats 5
```

`test_mcc_rotation` measures table initialization, memory and isolated lookup
throughput. `benchmark_rotation.py` measures complete timesteps and variance
per computational cost using identical physics for alias and cumulative paths.
It performs one separate warmup run per sampler before measuring fresh seeded
ensembles, and records the timing mean, median and standard deviation.
`perlmutter_rotation.sbatch` builds and runs the portable checks and benchmarks
on one Perlmutter GPU; provide the GPU allocation and Python environment when
submitting. CPU results alone do not establish CUDA/HIP/SYCL performance.

## Kinetic IAA spectator model

`rotation_model="iaa_spectator"` is a separate N2 model. It keeps the ordinary
inclusive elastic cross section, samples its elastic DCS, and then selects a
discrete rotational outcome conditional on that angle. The transition
probabilities follow Eqs. (11.24), (11.35) and (2.48); their angular averages
are not forced to reproduce the input elementary integral cross sections.
No matrix balancing or additive rotational rate is used.

Prepare its source data offline in an environment with elmolcs:

```sh
python Tools/CrossSections/thermiaa_spectator.py \
    --output /path/to/warpx-data/MCC_cross_sections/N2/IAA/thermal_spectator.rot
```

The spectator data file has three ASCII header lines followed by a little-endian payload:

```text
WARPX_IAA_SPECTATOR
N2 iaa_spectator J_MAX 0
N_ENERGY N_ANGULAR_NODES
```

Arrays are energies (binary64), angular-row offsets (int32, `N_ENERGY+1`),
`sin(theta/2)` knots (binary32), angular basis values (binary32,
`[angular_node,4]`), and transition weights (binary32,
`[energy,J_initial,7,4]`). The seven choices are `J_final = J_initial +
2*(-3,-2,-1,0,1,2,3)`; negative final states have zero weight. The four ranks
are 0, 2, 4 and 6. A common scale at each energy cancels in normalization.
The fourth field of the second header line is reserved and must be zero.
The data are temperature independent and shared on the host and device.

Initialization prepares only the small Boltzmann alias. The device selects
a virtual initial state in constant time and evaluates seven outcomes from
the tabulated weights. This factorization avoids storing a large table of
normalized probabilities and performs no Bessel functions or state-population
loop on the device. The virtual state is not stored or evolved. Canonical
energy changes remain binary64. `rotation_sampling="cumulative"` uses a
cumulative reference for the initial-state draw with identical physics.
Energy knots also remain binary64 when particles use binary32. Row mixtures
use a square-root coordinate immediately above zero and the first excitation
threshold; this resolves the channel onset inside the smallest energy interval.

The data use the real elmolcs elementary tables and cover their 0–1000 eV
range. The initial state sum resolves a 1000 K bath. Rank truncation, the
spectator model's low-energy inaccuracies, and its lack of detailed balance
are explicit physical limitations. Independent checks compare the first moment
with Eq. (2.48),
and also test the second moment, signed outcomes and angle–energy correlation.
The angle-independent model preserves integral detailed balance in its reference rates;
that alone does not establish differential detailed balance either.

The independent reference requires NumPy/SciPy and the published physical
tables, without an elmolcs installation:

```sh
cmake -S . -B build \
    -DWarpX_SPECTATOR_TEST_DATA=/path/to/warpx-data/MCC_cross_sections
cmake --build build -j 8
ctest --test-dir build -R spectator --output-on-failure
python Tools/CrossSections/benchmark_rotation.py --spectator \
    --data-dir /path/to/warpx-data/MCC_cross_sections \
    --spectator-reference build/Tools/CrossSections/spectator_reference/reference.txt \
    --output build/spectator-benchmark --repeats 5
```

Add `--thermal-balance` to `test_thermiaa_spectator.py` to measure equilibrium
heating and cooling independently of the sampler. This reports the physical
model's imbalance rather than treating a nonzero result as a tabulation error.
For the Perlmutter script, set `WARPX_SPECTATOR_DATA` to the same data directory.
Add `--broad-spectrum` to the benchmark driver to initialize a log-uniform
0.0021--900 eV electron ensemble at each of the four bath temperatures. This
exercises divergent energy lookups and table-cache access. Its recorded mean
electron energy change is an ensemble difference, independent of particle
reordering; the fixed-energy tests provide the separate physics checks.
The benchmark distributes the total particle count across 128 cells by default
(`--cells`). WarpX's within-cell injection loop is serial, so putting the entire
ensemble in one cell exaggerates GPU startup time. The fixed-energy physics
tests retain their minimal one-cell setup.
