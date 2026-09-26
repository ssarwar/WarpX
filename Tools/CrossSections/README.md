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

Physical cross sections and bundles are prepared offline and stored in
`warpx-data/IAA`. WarpX reads existing files; it never invokes elmolcs or a
cross-section generator. Initialization applies the chosen Boltzmann populations
and packs shared in-memory sampling tables. This supports arbitrary fixed
rotational temperatures without generating physical data during a simulation.

`rotation_reference.py` constructs integral forward/reverse rates in the
heavy-target limit, with relativistic electron phase space and canonical
rigid-rotor level spacings as thresholds. `audit_iaa_rotation.py` audits the
pinned source data and accepts `--reference-temperature` (default zero, an
explicit assumption). It records unresolved physics checks without exporting
an unvalidated production bundle.

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
WARPX_THERMAL_ROTATION_V3
N2 elastic_dcs J_MAX T_REFERENCE_K
N_ENERGY N_TRANSITION
```

Arrays contain energies (binary64), transition pairs `(J_initial,J_final)`
(int32), and integral rate coefficients (binary64, C order
`[energy,component]`). Component zero is unchanged. Other components correspond
to transition pairs, before Boltzmann weighting. Rates have units m³/s.
There are no angular bins. V1/V2 bundles must be regenerated; V3 omits the
molecular recoil shift from thresholds and source-rate detailed balance.
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

Export the synthetic verification inputs offline:

```sh
python Examples/Tests/collision/analysis_rotation_reference.py --export-only /path/to/warpx-data/MCC_cross_sections/IAA/rotation/verification
```

These files are also stored on `warpx-data/IAA`; they are not production
molecular cross sections. CTest regenerates verification fixtures as a separate
setup step. The PICMI input and performance driver only read prepared files.

```sh
ctest --test-dir build -R 'rotation|rbeq_source|secondary_angles' --output-on-failure
python Tools/CrossSections/benchmark_rotation.py --data-dir /path/to/warpx-data/MCC_cross_sections/IAA/rotation/verification --output build/rotation-benchmark --repeats 5
```

`test_mcc_rotation` measures table initialization, memory and isolated lookup
throughput. `benchmark_rotation.py` measures complete timesteps and variance
per computational cost using identical physics for alias and cumulative paths.
It performs one separate warmup run per sampler before measuring fresh seeded
ensembles, and records the timing mean, median and standard deviation.
`perlmutter_rotation.sbatch` builds and runs the portable checks and benchmarks
on one Perlmutter GPU; provide the GPU allocation and Python environment when
submitting. CPU results alone do not establish CUDA/HIP/SYCL performance.
