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

`rotation_reference.py` constructs integral forward rates, finite-mass
COM detailed-balance reverse rates, and an unchanged integral contribution.
`audit_iaa_rotation.py` audits the pinned source rates. Invoke it with the
archived package source on `PYTHONPATH`, the archive via `--archive`, and an
output directory via `--output`. It records unresolved physics gates without
exporting an unvalidated production bundle. Reference temperature zero is an
explicit audit assumption; the inclusive measurements do not establish it.

The collision draws its angle from the existing elastic DCS, then draws a
signed rotational energy change. No separate rotational DCS is used. The
thesis's spectator rainbow parameter is not a hard quantum cutoff. Only
exact energy/recoil accessibility conditions the rotational outcome.

Integral bundles have three ASCII header lines and a little-endian payload:

```text
WARPX_THERMAL_ROTATION_V2
N2 elastic_dcs J_MAX T_REFERENCE_K
N_ENERGY N_TRANSITION
```

The arrays are energies (binary64), transition pairs `(J_initial,J_final)`
(int32), then integral rate coefficients (binary64, C order
`[energy,component]`). Component zero is unchanged; the others correspond to
transition pairs. Rates have units m³/s and precede Boltzmann weighting.
There are no angular bins. V1 bundles must be regenerated. Rigid-rotor constants,
spin weights and thresholds are canonical, and the state sum must converge.

At initialization, WarpX weights the rotational states and constructs one
sparse alias distribution and one cumulative distribution per energy row.
Outcomes are ordered by increasing energy loss, with descending initial
internal energy for equal losses. A conservative mass-dependence bound checks
that this ordering also orders the accessibility thresholds. Every nonempty
row must retain an unchanged or de-excitation outcome. Reference input and
sampler storage each have a 512 MiB limit.

The usual device path samples one alias. When a rotational excitation is
inaccessible at the elastic angle, a bounded bisection finds the accessible
prefix. Both the energy-row mixture and the outcome distribution are then
conditioned on that prefix. There is no rejection loop, state loop, or
special-function evaluation. The cumulative reference uses the same physics.
The angular draw is unchanged, and loss labels are never interpolated.

The neutral-rest-frame energy and the signed laboratory-angle recoil solver
use relativistic transformations. The solver chooses the higher outgoing
energy if two forward solutions exist immediately above threshold. Unit tests
compare accessibility with independent numerical minimization of the final
energy budget, preserve the sampled angle, and check conditioned mixtures
against explicit probabilities.

Integral detailed balance is enforced in the reference rates. Reusing an
energy-dependent elastic DCS does not impose exact differential detailed
balance. Conditioning can change integral rotational rates inside the narrow
recoil-accessibility bands. These are explicit approximations; they must be
checked separately from interpolation and thermal power balance.

`refine_grid` checks aggregate rates and absolute first/second transfer moments.
The default target is 0.1%, with a 1e-4-of-peak floor for insignificant tails;
independent midpoint checks require 0.2%. Four-float32-ulp threshold bands are
tested separately. The first interval uses a `sqrt(E)` mixture to preserve
finite elastic cross sections next to zero. Analytic verification bundles
check population/moment convergence and thermal power balance at 100, 300
and 1000 K. They do not establish production N2/O2 source accuracy.

```sh
ctest --test-dir build -R 'rotation|rbeq_source|secondary_angles' --output-on-failure
```

`test_mcc_rotation` runs portable AMReX kernels and reports startup time, table
bytes, sampling cost and distribution moments. Full timestep measurements use
`benchmark_rotation.py --output build/rotation-benchmark --repeats 5` in a
WarpX Python environment. Compare `alias` and `cumulative` with identical
physics and report variance times cost. CPU measurements do not establish
CUDA, HIP or SYCL throughput.
