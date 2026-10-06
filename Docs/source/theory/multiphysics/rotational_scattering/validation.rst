.. _rotation-validation:

Validation and performance
==========================

The 300 K production data use real elmolcs inputs and the published constraints
listed in :ref:`rotation-sources`. Numerical accuracy refers to the prescribed model;
it does not establish experimental accuracy of the low-energy O2 bridge or
the high-J continuation. Detailed logs and temporary data belong in build
folders, outside warpx-data.

Physics and precision
---------------------

Independent source evaluations and decoded distributions gave maximum
low-energy interpolation errors of 0.120% for N2 and 0.100% for O2, including
angle-weighted first and second energy-transfer moments. Equilibrium power
imbalance was below 0.0065% of heating plus cooling. The high-energy
rate/angular interpolation error was below 0.103%; the separate
momentum-transfer approximation was below 0.038%.

CPU particle tests cover heating, cooling and equilibrium at 0, 100, 250, 300,
350 and 1000 K for N2, O2 and their air mixture. The additional-temperature
fixtures have a 20 eV domain and are prepared offline for thermal tests.
Equilibrium checks include paired changes in the first three energy moments
and several energy-distribution cutoffs. Unequal bath temperatures,
timestep refinement, subcycling and restart are tested separately.

Tests of the real combined family with real RBEQ ionization tables compare
ionization counts to the independently computed one-event channel probability.
They check that the created ion and electron counts agree. This exercises
competition between the combined family and ordinary MCC processes at a
large optical depth, where double counting the elastic contribution would
be readily detectable.

Endpoint checks cover direct lookup, cached selection and the fallback
selector, finite majorants, malformed metadata and bundle temperature.
Double, single-particle and all-single CPU builds pass the sampler and
threshold tests. Double and all-single CUDA runs also pass the 28-energy
sampler checks, heating/cooling/equilibrium cases and RBEQ competition checks.
The all-single CUDA run additionally passes the direct/cached/fallback endpoint
checks. A seeded device test verifies that collision and sampling uniforms
resolve the upper half interval more finely than binary32.

Double particle precision is needed to retain meV transfers in MeV electrons.
For example, a forward 10 meV loss at 2.5 MeV was stored as
0.01000000071 eV with double particles and as zero with single particles.
This is distinct from the accuracy of the double-energy sampling tables.
Moving-target recoil tests check laboratory four-momentum conservation;
stationary-target tests include tiny deflections, backward scattering, signed
losses, and the documented narrow nominal-threshold continuation.

A100 measurements
-----------------

The Perlmutter campaign uses A100-SXM4-40GB devices, CUDA compilation for
``sm_80``, and double field/particle precision. Five independent seeds are
measured after warmup, alternating alias, cumulative and elastic-only controls.
Each GPU has one MPI process. Both gases are resident. Workloads include
thermal, 2.47 eV resonance, 50 eV, 2.5 MeV, and log-uniform mixed energies from
0.002 eV to 3 MeV within warps. Two global particle counts, 262144 and 1048576,
are measured. The grid has 128 cells per MPI process.

The full PIC timing includes gather, charge/current deposition, field advance,
particle advance and collisions in a homogeneous periodic test. Its small
physical electron density keeps self-fields negligible while retaining those
operations. The complete MCC operator is timed separately, including its
normal selection, recoil and synchronization costs. Particle observations
and file loading are outside the timed step region.

For 1048576 particles and 64 measured steps, the original prepared binary alias tables
without search indices gave the following median times:

=====================  ==============  =============  ==============  =============
Distribution           MCC, one GPU    PIC, one GPU   MCC, four GPUs  PIC, four GPUs
=====================  ==============  =============  ==============  =============
Thermal, 300 K          35.14 ms        1300.13 ms     18.51 ms        480.52 ms
Resonance, 2.47 eV      47.00 ms        1312.27 ms     22.72 ms        486.39 ms
Intermediate, 50 eV     50.29 ms        1315.40 ms     24.91 ms        487.41 ms
Relativistic, 2.5 MeV   37.17 ms        1302.81 ms     19.26 ms        480.93 ms
Mixed energies         53.88 ms        1318.99 ms     24.43 ms        486.17 ms
=====================  ==============  =============  ==============  =============

Alias sampling was faster than the identical-physics cumulative reference in
all these cases. At one million particles it reduced complete MCC time by
about 1--9%, and complete PIC time by about 1--2% on one GPU. Relative to the
elastic-only control, the added rotational physics increased complete PIC
time by at most 1.4% across the one-GPU cases. Its additional complete-MCC
cost was 11--41%, depending on workload; joint-distribution searches, discrete
outcomes and signed recoil account for additional work. Comparisons of
isolated sampler timings alone do not determine this acceptance.

Compute Sanitizer checks the sampler without suppressing memory errors.
For the isolated one-rank sampler in an MPI build, GPU-aware MPI pointer
classification is disabled: Cray MPICH otherwise probes host addresses with
``cuPointerGetAttribute`` and reports handled API errors to the sanitizer.
Timed MPI simulations retain GPU-aware communication. Sampler and full-MCC
memcheck runs reported zero errors, including the all-single CUDA build.
HIP and SYCL toolchains/hardware were unavailable for this campaign; those
backends have not been compiled or executed here.

``benchmark_reciprocal_rotation.py`` records the variance of energy and
longitudinal momentum transfer multiplied by elapsed step time. Transfer
observables subtract the initial state to remove initial-sample noise.
Five seeds give only a coarse estimate of the variance ratio; an apparent
noise advantage should not be inferred from that small sample alone.

Exact search indices
~~~~~~~~~~~~~~~~~~~~

The production bundles include optional quantile lookup bounds. They add
7,849,816 bytes, including alignment, bringing both gases to 751,195,072 bytes
(716.4 MiB). All existing physical arrays remain bitwise unchanged. The CPU
and GPU tests compare indexed and full searches with identical random draws
and require identical angles and discrete outcomes for each event.

At one million particles, the indices reduced complete MCC time by about
4--9% for the thermal, resonance, intermediate and mixed-energy cases.
The 2.5 MeV case changed by less than 1%, because it uses the analytic
Rutherford angle and the momentum-transfer bank. These measurements use prepared binary inputs. Text input constructs the same
search indices once during initialization.

Full PIC timing depends on particle layout. In the 128-cell test, indexed
alias runs took about 1--2% longer overall, with the difference in current
deposition; sorting time was unchanged. Sorting uses atomic linked lists,
so separate GPU runs do not promise identical particle ordering or seeded
particle trajectories. A second test kept four particle tiles but used 4096
cells (256 particles per cell). With five alternating seeds, indexed/full
search PIC-time ratios were 0.993 (thermal), 0.989 (resonance), 0.992
(intermediate), 1.003 (2.5 MeV), and 0.978 (mixed). Thus the collision search
improvement is measured, while an end-to-end speedup is workload dependent.
No unexplained regression above 5% remained against the same-physics
reference. The original full-grid lookup remains supported for prepared binary bundles
without the optional indices.

Readable-data validation
------------------------

The readable probability encoding has its own accuracy budget and is checked
against the decoded, previously validated probability cells. It preserves all
scalar rates, energy grids, angular quantiles, deflections, changing
probabilities, discrete labels and threshold supports. Only the storage of
conditional spectra changes. Reconstruction is checked after decimal rounding
and after float32 alias packing. This comparison is additional to, and does
not replace, the independent source and thermal checks above.

A deterministic native check of every cell after alias packing found maximum
relative errors of 1.996e-6 for both gases in unchanged/up/down probabilities
and the separate positive/negative first, second and fourth moments. Both
float64 and all-single CPU builds exercise the text reader and full MCC tests.

The initial local CPU measurements of text reading, probability reconstruction,
alias preparation and table validation were 4.31 seconds for N2 and 2.74 seconds
for O2; they are initialization times, not GPU collision throughput. Measure
cold filesystem/cache conditions separately when deploying on a new system.
These trials are within the allowed minute-scale startup budget. Text inputs
remain substantial: approximately 140 MiB for N2 and 95 MiB for O2, compared
with 716.4 MiB combined binary data. Human readability does not imply that
every probability is a two-column cross section; their representation is
specified in :ref:`rotation-data-format`.

On a 40 GB Perlmutter A100, both readable alias and cumulative inputs passed
independent 28-energy sampler checks, exact-search equivalence, threshold/range
and malformed-input checks, and air heating/equilibrium/cooling cases. Sampler
Compute Sanitizer runs for both gases reported zero errors. Deterministic
GPU decoding of every probability cell reproduced the CPU encoding-error
bounds above. This text-input campaign used double precision; the all-single CUDA and
four-GPU results in the earlier section are from the binary-input campaign.

First-read initialization measurements were 6.89 seconds for N2 and 18.40
seconds for O2. Filesystem state matters: subsequent reads were much faster.
In the repeated complete-PIC benchmark with both gases resident, median
initialization was 6.27--6.44 seconds for text input, versus 0.84--0.93 seconds for
prepared binary input. No initialization work is repeated per collision.

With 1048576 electrons, 4096 cells, 64 timed steps and three alternating seeds,
the text-input/binary-input median full-MCC time ratios were 0.9992 (thermal), 1.0024 (resonance),
0.9997 (50 eV), 1.0012 (2.5 MeV) and 0.9993 (mixed). The corresponding complete
PIC ratios were 0.9906, 0.9911, 1.0013, 0.9899 and 0.9996. These measurements
show no material collision-throughput regression from preparing aliases at
startup. They do not establish a noise advantage: three seeds give a weak
variance estimate, and a changed alias layout changes seeded trajectories.
The packed second and fourth moments provide the stronger deterministic check
that the encoding preserves physical energy diffusion and its rare tails.

The isolated device allocation increments were 471859200 bytes for N2 aliases
and 281018368 for O2, including allocator granularity: about 718 MiB together,
below 1 GiB. The complete million-electron benchmark's measured peak device
allocation was 1415643136 bytes for both text and binary cases; that includes
particles, fields and other allocations as well as rotational tables.

Reproducing numerical checks
----------------------------

Use the commands in ``Tools/CrossSections/reciprocal_rotation/README.rst``.
The independent checks have different responsibilities:

* ``prepare_jung.py`` reproduces the inferred angular fractions from the
  readable manual figure readings, with their stated reading uncertainty.
* ``source_anchor_checks.py`` checks Jung branch fractions, the O2 published
  integral and momentum-transfer constraints, and the zero-energy limit.
* ``test_transition_normalization.py`` checks the implicit normalization root
  against an independent scalar solve.
* ``verify.py`` evaluates additional source energies and angularly weighted
  rate/energy-transfer moments, including equilibrium heating and cooling.
* ``verify_high.py``, ``q_transition_check.py`` and ``phase_average_check.py``
  separately check angular continuation, momentum-transfer binning and the
  oscillatory tail approximation.
* ``check_text_cells.py`` compares every native decoded probability cell to the original
  validated probability distribution, including rare tails and float32 packing.
* ``test_reciprocal_rotation`` samples the actual C++/GPU implementation;
  ``sampling_reference.py`` compares its first and higher moments with an
  independently decoded distribution. This comparison alone is not a source
  validation.
* ``test_rotation_thresholds`` exercises signed recoil, moving targets,
  tiny/large angles, threshold continuation and particle precision.
* The collision example tests exercise complete MCC rates, thermal evolution,
  ionization competition, direct/cached/fallback selectors, invalid inputs,
  subcycling and restart. Benchmarks measure complete MCC and PIC costs.

Grid selection and error measures
---------------------------------

The source reference uses mixed logarithmic/linear energy grids with physical
thresholds and source joins inserted explicitly, plus angular quadrature in
:math:`y=\sin(\theta/2)`. Reference refinement solves new energies before
comparing with interpolated old ones. Angular quadrature refinement and
initial-state/rank refinement are separate checks. The subsequent exporter
refines sampling energies and conditional quantiles against rates and the
positive/negative first and second transfer moments, including angular
weighting. It does not choose one fixed point count for all gases.

For linear differential density :math:`D(y)=D_0+s(y-y_0)` on a segment of
width :math:`h`, its exact angular mass is

.. math::

   8\pi\left[y_0D_0h+\tfrac12(D_0+y_0s)h^2+\tfrac13sh^3\right].

This is the integral used before inversion; summing unweighted angular samples
would be wrong. Above the cold join the exporter evaluates the actual IAA
angular CDF, avoiding a second approximate quadrature of it on the rotational
grid. The conditional-spectrum compression checks both total variation and
energy moments. Probabilities below :math:`10^{-15}` are pruned only with a
checked bound on the summed omitted probability; energetic rare tails are
also checked by their transfer moments.

The primary acceptance targets are less than 0.2% error in combined rates
and first/second transfer moments, and less than 0.1% equilibrium power
imbalance relative to heating plus cooling. Near a zero moment, absolute
error and the separate positive/negative powers are used; a relative error
in a nearly canceled signed mean is not meaningful. Statistical particle
tests use sampling uncertainty rather than weakening deterministic accuracy
requirements. Float32 alias rounding, float64 physical labels, interpolation,
state truncation, q-bin approximation and finite timestep errors are distinct
contributions. The printed experimental uncertainty is not part of a
numerical discretization tolerance.

Open physical limits
--------------------

An accurate numerical sampler cannot supply missing experimental physics.
The O2 1--20 eV bridge, the unmeasured N2 angular tails, unresolved O2 spin,
fixed vibrational state, bound rigid-rotor high-J closure, cold scattering
lengths and high-energy sudden limit remain assumptions. No exact
finite-energy reciprocity claim is made for the high-energy bank, or exact
recoil conservation inside the documented nominal-threshold continuation.
Room-temperature source agreement and tests support the intended use; they
do not establish universal accuracy at arbitrary temperature or electron
energy. The complete 0--1 GeV domain is a supported modeled range, not a
measurement range.
