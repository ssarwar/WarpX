.. _beam-air-validation:

Physics verification and performance
====================================

This page separates independent physics checks, numerical precision, physical
model uncertainty and backend measurements. Test commands and recorded limits
apply to the feature versions identified in each campaign.

Independent analytical references
---------------------------------

The analytical fixtures use synthetic data under WarpX build directories, never
under warpx-data. They isolate a known solution from source-model uncertainty.
State the fixed-energy, stationary/heavy-neutral assumptions when using angular
or fixed-loss references; finite recoil changes those idealized problems.

* **Angular sampling.** Isotropy has :math:`p(\mu)=1/2`,
  :math:`\langle\mu\rangle=0` and :math:`\langle\mu^2\rangle=1/3`.
  For the N2 :math:`d\pi` law, the exact CDF is

  .. math::

     F(\mu)=\frac12+\frac{15}{28}
         \left(\mu-\mu^3+\frac{14}{15}\mu^5\right),

  with :math:`\langle\mu^2\rangle=3/7`,
  :math:`\langle\mu^4\rangle=137/441`, and Legendre moments
  :math:`a_2=1/7`, :math:`a_4=8/63`. Mixture CDFs and moments are the
  probability-weighted sums of component CDFs and moments.

* **Angular relaxation.** In a fixed-energy Lorentz gas with collision rate
  :math:`\nu` and single-event moment :math:`a_\ell`, continuous relaxation is
  :math:`A_\ell(t)=A_\ell(0)e^{-\nu(1-a_\ell)t}`. With one event permitted per
  substep, the exact discrete expectation after :math:`N` steps is instead
  :math:`A_\ell(0)[1-p(1-a_\ell)]^N`, where
  :math:`p=1-e^{-\nu\Delta t}`. Check this expression first, then demonstrate
  convergence to the continuous solution. A repeated sampler test does not
  exercise the complete MCC dispatch, gas sampling and majorant path.

* **Attachment.** For constant competing rates :math:`\nu_j`, electron survival
  is :math:`e^{-\nu t}`, with :math:`\nu=\sum_j\nu_j`, and product fraction
  :math:`(\nu_j/\nu)(1-e^{-\nu t})`. Raw three-body data require
  :math:`\nu_3=n_{O_2}n_Mv\sigma_3`. Change each density factor independently.
  Unequal particle weights require weighted Bernoulli variance, not a count-only
  Poisson error estimate.

* **Discrete cooling.** For loss :math:`Q`, constant rate above threshold and
  terminal level :math:`K=\lfloor E_0/Q\rfloor`, levels :math:`k<K` follow
  :math:`{N\choose k}p^k(1-p)^{N-k}` after :math:`N` one-event steps. The
  absorbing level has the remaining probability. The continuous limit is the
  corresponding Poisson ladder with mean :math:`\nu t`. The constant-rate,
  heavy-target fixture must not be confused with real energy-dependent cooling.

* **Kinematics and sources.** Independent invariant endpoints and energy,
  momentum and charge accounting check accepted collision products. Gaussian
  antiderivatives check beam normalization and space-time yield; a separate
  Coulomb Green-function calculation checks unbounded Gaussian fields and
  their Lorentz transform. Discrete continuity and exact electron/ion charge
  footprints are checked independently of that continuum field reference.

Permanent test entry points
---------------------------

With a configured 1D Python build, enable real production data explicitly:

.. code-block:: bash

   cmake -S . -B build -DWarpX_DIMS=1 -DWarpX_PYTHON=ON \
       -DWarpX_IAA_AIR_TEST_DATA=/path/to/warpx-data/MCC_cross_sections \
       -DWarpX_RECIPROCAL_TEST_DATA=/path/to/warpx-data/MCC_cross_sections
   cmake --build build -j 8
   ctest --test-dir build -N -R 'iaa|reciprocal|rotation_thresholds'
   ctest --test-dir build --output-on-failure \
       -R 'test_1d_(iaa|reciprocal|rotation_thresholds)'

Select the compiler, MPI and accelerator options appropriate to the machine.
The test listing is important: an omitted external-data option can omit the
production checks entirely. Repeat in a separate particle-single build using
``WarpX_PARTICLE_PRECISION=SINGLE``; record field precision independently.

``Tools/CrossSections/check_iaa_angular.py`` checks the actual sampler against
independent angular quadrature, analytic distributions and mixtures.
``test_iaa_air.py`` checks source anchors, channel accounting, repairs, RBEQ
pairing, units, continuation and manifest selections. The complete-MCC
``inputs_test_1d_iaa_analytical_picmi.py`` checks attachment and cooling.
The committed relaxation fixture is a repeated angular-kernel test; the wider
full-MCC relaxation audit has not yet become a permanent CTest.

Separate standalone CMake projects under ``Tools/Algorithms/BackgroundMCC``,
``ProtonImpactIonization`` and ``PrescribedFluids`` test the production C++
helpers with independent references. Their READMEs give build commands and
dependencies. In an RZ Python build, tests under ``langmuir_fluids``,
``collision``, ``implicit`` and ``reduced_diags`` cover beam continuity,
immobile products, solver coupling, source budgets and restart guards.
Distributed and long performance studies are separate from short CPU CI.

Current production-air coverage and limits
------------------------------------------

The October 13, 2026 audit exercised WarpX ``30a2685e67267289320ddbeb6bb01aaf7fa32dc0``
and warpx-data ``f39ddb790bb9cb5ce6ff7d4695b2d6a6f117f04d`` on Perlmutter
**A100-SXM4-40GB** GPUs in double and single particle precision. The subsequent
gas-directory reorganization preserved numerical data and process ordering.
These production-air results are distinct from the earlier elastic/rotation
campaigns below. HIP and SYCL GPU physics and performance remain untested.
Earlier SYCL compilation of prescribed-fluid support does not close that gap.

The campaign covered analytical attachment/cooling, full-MCC angular relaxation,
production cooling/swarm/pulse checks, and compact versus dense representations
of the same physical interpolants. The audit's full-MCC isotropic, d-pi and
mixture relaxation cases checked P1--P4 with analytical estimator variances.
Private audit fixtures are additional evidence, not maintained regression tests.

Independent host quadrature of all 2,992 actual angular-table rows found maximum
CDF error :math:`5.245\times10^{-4}` and P1--P4 moment errors below
:math:`2.8\times10^{-4}` in both particle precisions. These are numerical row
errors. Endpoint completion, angular transition choices and uncertain
high-energy continuations have separate, larger model uncertainties.

At 1,048,576 particles and five paired ensembles, compact-table full-PIC median
costs were about 1.1--3.4% lower in double precision and 1.4--5.5% lower in single
precision across the measured workloads. The double-precision selector storage
was 13,663,752 versus 66,317,592 bytes across the gases; angular storage was
25,691,432 versus 98,797,400 bytes. These figures exclude unchanged rotational
families and are not peak process memory. Perlmutter jobs ``59745672`` and the
single-precision continuation campaign supplied the full-step timings; isolated
lookup measurements were run separately. No new timing measurement is implied
by this documentation reorganization.

The following gaps prevent describing the full beam--air model as validated:

* Production cooling/swarm tests check charge, cooling and drift direction.
  They do not establish mobility, diffusion, ionization/attachment coefficients
  or an energy distribution against independent experiments or a Boltzmann solver.
* In the A100 dry-air cooling refinement (job ``59754199``), 20 eV initial
  electrons evolved for 72 steps of 1e-13 s. Changing 16 to 64 coupled subcycles
  changed final mean energy by 0.489%, with paired 95% interval 0.372--0.607%.
  The 32-to-64 difference was still resolved. Convergence at 64 is not established.
* The reduced 800 MeV pulse regression is homogeneous and uses different pulse
  widths/spacing from the beam paper. A geometry- and time-matched experimental
  yield comparison remains incomplete. Older prescribed-fluid studies with
  synthetic electron rates validate representation choices, not IAA chemistry.
* Single-particle native lookup exceeds the offline 0.05% pointwise rate budget
  near some rounded zeros/thresholds. O2 RBEQ reaches about 0.0726%, above its
  0.02% export target. Observable sensitivity and a separate runtime precision
  budget are needed. A malformed negative raw value can also vanish before the
  resampled reader validates it; production files do not contain negative rates.
* A 7 eV loss at 1 GeV can disappear when the final momentum is stored in float.
  Tiny oblique deflections can disappear as well. Double intermediate
  ``1-cos(theta)`` does not overcome final-storage precision. Use double
  particles when these transfers must be resolved.
* Ordinary MCC rejects a positive-loss event below its recoil-shifted threshold
  before calling the excitation helper. Direct helper calls below that threshold
  can clamp inaccessible phase space and must not be interpreted as conservative
  events. Permanent checks around the boundary, moving neutrals and oblique
  stored-momentum changes remain incomplete.
* Completion/transition/tail sensitivity has not been propagated through a
  converged beam or swarm calculation. Source-table loss moments alone do not
  validate the complete recoil and ionization transfer moments.
* The new complete-data GPU campaign used one GPU and 1D without MPI, FFT or
  openPMD. Earlier distributed rotational/fluid tests do not substitute for
  multi-GPU, RZ or application-scale validation of the complete IAA inventory.
  The default GPU CI smoke-test selection does not run this production-data suite.

Keep these limits separate from implementation pass/fail results. Do not loosen
physics assertions to turn an unresolved discrepancy into a passing check.

Prescribed-fluid and integration evidence
-----------------------------------------

``Tools/Algorithms/PrescribedFluids/README.md`` gives independent Gaussian
references, quiet/ordinary particle controls, equal-error studies, source-weight
and mesh refinement, and Perlmutter submission scripts. Its ``VALIDATION.md``,
``REAUDIT.md`` and ``DEVELOPMENT_SYNC.md`` retain dated hardware, solver and
restart evidence. ``Tools/Algorithms/BackgroundMCC/development_sync_2026_10.rst``
and ``dsmc_branch_audit_2026_10.rst`` record unchanged-stock controls and
integration failures. Historical results must be read with their recorded
revision, physics inputs, GPU memory variant and unresolved limitations.

GPU performance validation
--------------------------

CPU timing and accelerator compilation are necessary checks, but neither proves
a GPU speedup. Validate a performance change with reference and candidate
builds on the target cluster. Use the same compiler, optimization flags, GPU
architecture, AMReX options, MPI layout, particle decomposition, input data and
random seed. Discard a warm-up run, collect at least 20 timed repetitions and
compare medians together with a confidence interval rather than a single wall
time.

The ``inputs_test_1d_background_mcc_many_attachment_picmi.py`` test accepts
``--particle-count``, ``--process-count``, ``--steps``, ``--subcycles``,
``--cell-count``, ``--max-grid-size`` and ``--target-acceptance`` arguments for
stand-alone performance runs. ``--particle-count`` is the global source count
and must be divisible by ``--cell-count``. Multiple cells and a smaller maximum
grid size create enough boxes to occupy multiple MPI ranks and GPUs. Pass
``--mpi-timing`` under ``mpiexec`` with a ``WarpX_MPI=ON`` build so only rank
zero writes the generated table and the reported time is the maximum over
ranks. The benchmark aborts if the MPI world sizes do not agree. When specified,
``--target-acceptance`` chooses an explicit conservative majorant and, with one
step and one subcycle, an optical depth whose expected product-event fraction
is the requested value. Large values are useful lifecycle stress tests, but are
not time-converged physical simulations.

Exercise at least the following matrix:

* 1, 8, 32, 64 and 128 processes, using both shared and different energy grids;
* attachment-only, ionization-only and representative elastic/excitation/product
  mixtures;
* approximately 0, 1, 10, 50 and 90 percent product-event fractions;
* both small and large particle tiles and single- and double-precision particles;
* one GPU, every GPU on one node and a multi-node weak- and strong-scaling case.

Use the ABLASTR profiler regions
``BackgroundMCCCollision::selectAndScatter()``,
``BackgroundMCCCollision::createProducts()``,
``BackgroundMCCCollision::compactAttachedElectrons()`` and
``BackgroundMCCCollision::doCollisions()`` to separate selection, creation,
removal and total MCC time. Also inspect a GPU timeline and kernel metrics with
an appropriate vendor profiler. In particular, check kernel-launch count,
device-to-host synchronization, allocation time, atomic contention, achieved
occupancy, register spills, cache and memory throughput, and MPI idle time.

Accept an optimization only if it preserves particle counts, charge, weights,
discrete losses and the statistical angle and energy distributions, reduces
the total MCC median for the intended workload, does not regress the
low-acceptance production regime and keeps peak device memory bounded. This
measurement is also what determines the next optimization: for example,
profiles dominated by per-tile allocation or count synchronization motivate
reusable scratch storage or batched count transfers, while atomic-dominated
high-acceptance runs motivate a scan-based event queue. Such changes should not
be selected from CPU timings alone.


.. _rotation-validation:

Elastic and rotational validation campaigns
-------------------------------------------

The 300 K production data use real elmolcs inputs and the published constraints
listed in :ref:`rotation-sources`. Numerical accuracy refers to the prescribed model;
it does not establish experimental accuracy of the low-energy O2 bridge or
the high-J continuation. Detailed logs and temporary data belong in build
folders, outside warpx-data.

Physics and precision
~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~

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
^^^^^^^^^^^^^^^^^^^^

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

Subcycled MCC runtime optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The October 2026 integration adds exact energy-search bounds, resolves alias
alternates to outcome identifiers at initialization, and specializes kernels
for scalar/parser backgrounds and mutually exclusive rotational models.
Pure scattering retains its device error flag across successful calls and
omits redundant iterator exit synchronization. The host still checks the
flag after every call; the RNG kernel synchronization and iterator entry
synchronization remain. Empty ranks skip only pure scattering. Product
operators still participate in fluid-increment commits on every rank.

On an A100-SXM4-80GB, with both gases resident, 4096 cells, four particle
boxes and 64 subcycles per gas, the following times are medians of three
seed medians. Each seed has 20 synchronized blocks of eight PIC steps,
after four warmup steps. Variant order alternates between seeds. Positions
and momenta are restored outside each block. The baseline is integration
commit ``305e203ff``; the acceptance law and physical tables are identical.

.. list-table:: Time per PIC step, milliseconds
   :header-rows: 1
   :widths: 46 18 18 18

   * - Workload
     - Baseline
     - Optimized
     - Reduction
   * - 1048576 electrons, thermal MCC
     - 35.005
     - 32.038
     - 8.5%
   * - 32768 electrons, thermal MCC
     - 14.879
     - 13.080
     - 12.1%
   * - 1048576 electrons, mixed energies, full PIC
     - 51.860
     - 47.843
     - 7.8%
   * - Thermal MCC with inactive ionization product channels
     - 43.751
     - 42.158
     - 3.6%
   * - Mixed-energy MCC with parser-defined backgrounds
     - 55.264
     - 43.126
     - 22.0%

The product-channel control uses zero ionization cross sections to retain a
fixed population while exercising product buffers. The parser expressions
have constant values. The MCC workloads disable particle push, gather and
deposition in the normal PIC-step driver; field-driver work and bookkeeping
remain in the timed region. These controls isolate implementation costs; they are
not evolving air-plasma performance predictions. The full-PIC case includes
push, gather, deposition and field advance, with negligible self-fields.
Mixed energies are log-uniform from 0.002 eV to 3 MeV within warps. The PIC
timestep is 6.4e-13 s. A separate ordinary-elastic control improved by 2.4%.
The isolated ablation campaign also covers 2.47 eV, 50 eV and 2.5 MeV.

Alias remains the fastest tested sampler. The cumulative option additionally
uses exact cell-search bounds; it remains a useful independent search path,
but has greater memory use. Both-gas table storage is 716.425 MiB for alias
and 909.993 MiB for indexed cumulative sampling, below the combined 1 GiB
budget. The energy index adds 30728 bytes across the two gases; the optional
cumulative-cell indices add 27.738 MiB. In the repeated sampler comparison,
median initialization changed from 6.95 to 7.07 s for aliases and 7.67 to
7.87 s for cumulative tables. These are warm-filesystem measurements.

CPU and CUDA comparisons reproduce every decoded cell and the seeded
28-energy sample output byte for byte, comparing each sampler to its own
unoptimized implementation. Indexed energy interpolation agrees exactly
with the original search at 50439 N2 and 49539 O2 queries, including every
physical knot and adjacent floating-point values. Cumulative lookup checks
edge draws in all 223733 cells and retains the original strict upper-bound
search. Signed recoil, moving-target four-momentum conservation, analytic
Rutherford transport moments, thermal evolution, subcycling and product
creation remain separately tested. The implementation introduces no new
physical approximation or change in random-draw count.

The reproducible timing driver and the complete audit are in
``Tools/Algorithms/BackgroundMCC/benchmark_rotational_subcycling.py`` and
``Tools/Algorithms/BackgroundMCC/rotation_optimization_2026_10.rst``.
HIP and SYCL execution remain unverified. Nsight Systems captured CUDA
kernel/API activity; Nsight Compute could not acquire the driver's profiling
resource, so these results do not claim measured cache-hit rates or occupancy.

Readable-data validation
~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~

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
