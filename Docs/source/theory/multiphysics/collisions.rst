.. _multiphysics-collisions:

Collisions
==========

WarpX includes several different models to capture collisional processes
including collisions between kinetic particles (Coulomb collisions, DSMC,
nuclear fusion) as well as collisions between kinetic particles and a fixed
(i.e. non-evolving) background species (MCC, background stopping).

.. _multiphysics-collisions-mcc:

Background Monte Carlo Collisions (MCC)
---------------------------------------

Several types of collisions between simulation particles and a neutral
background gas are supported including elastic scattering, back scattering,
charge exchange, excitation collisions, impact ionization and electron
attachment.

The so-called null collision strategy is used in order to minimize the
computational burden of the MCC module. This strategy is standard in PIC-MCC and
a detailed description can be found elsewhere, for example in :cite:t:`b-Birdsall1991`.
In short the maximum collision probability is found over a sensible range of
energies and is used to pre-select the appropriate number of macroparticles for
collision consideration. Only these pre-selected particles are then individually
considered for a collision based on their energy and the cross-sections of all
the different collisional processes included.

The MCC implementation assumes that the background neutral particles are
**thermal** and move at non-relativistic velocities in the laboratory frame.
For each simulation particle considered for a collision, a neutral ordinary
velocity :math:`\boldsymbol{V}_n` is sampled from the user-specified classical
Maxwellian distribution. WarpX particle momentum components are normalized
momenta, or proper velocities,
:math:`\boldsymbol{u}=\gamma\boldsymbol{v}`.

For electron projectiles, Background MCC uses the fast approximate relative
proper velocity

    .. math::

       \widetilde{\boldsymbol{u}}
       = \boldsymbol{u}_e-\boldsymbol{V}_n,
       \qquad
       \widetilde{\gamma}
       = \sqrt{1+\frac{\widetilde{u}^2}{c^2}}.

This subtraction is exact to leading order when the electron is
non-relativistic, because :math:`\boldsymbol{u}_e\simeq\boldsymbol{v}_e`, and
is exact for a stationary neutral at any electron energy. When the electron is
relativistic, the neglected correction due to neutral motion has relative size
of order :math:`V_n/c`, which is negligible for a classical gas. This
approximation is specific to relativistic electrons colliding with
non-relativistic neutrals; it is not a general relativistic relative-velocity
formula.

The cross-section lookup energy and physical collision-rate speed are

    .. math::

       E_{\mathrm{lookup}}
       = \frac{m_e\widetilde{u}^2}
              {e(\widetilde{\gamma}+1)},
       \qquad
       g_{\mathrm{coll}}
       = \frac{\lvert\widetilde{\boldsymbol{u}}\rvert}
              {\widetilde{\gamma}},

where :math:`E_{\mathrm{lookup}}` is in electronvolts. In the
stationary-neutral limit, :math:`g_{\mathrm{coll}}` is the ordinary electron
speed, not the stored proper speed :math:`\lvert\boldsymbol{u}_e\rvert`.

The lookup energy approximates the electron kinetic energy in the neutral rest
frame. For an electron incident on an atomic or molecular neutral, its
difference from the total center-of-momentum kinetic energy is of relative
order :math:`m_e/M` until extreme relativistic energies. WarpX therefore does
not distinguish these two energies for cross-section lookup in this model.
Using ``ParticleUtils::getCollisionEnergy()`` would evaluate the exact two-body
center-of-momentum energy and add another square root without a useful increase
in accuracy. A full three-vector Lorentz transformation is likewise not needed.

Thus, the frequency for process :math:`i` is
:math:`\nu_i=n_n\sigma_i(E_{\mathrm{lookup}})g_{\mathrm{coll}}`.

All configured processes, including impact ionization and attachment, compete
in one draw. A particle therefore undergoes at most one accepted process in
each Background MCC collision substep. In-place elastic and excitation outcomes
are applied immediately. Product-changing outcomes are recorded and then
created from one compact event queue. Channels are grouped by process type and
destination species for allocation, but all groups are populated by one
particle-creation kernel.

By default, the null-collision majorant is constructed from the union of all
cross-section table knots. WarpX reuses the total-cross-section row of the
cumulative process table described below. Between consecutive union knots the
summed cross section is linear. WarpX evaluates both endpoints and, on a
decreasing segment, the one possible stationary point of
:math:`\sigma_{\mathrm{tot}}(E)g_{\mathrm{coll}}(E)`. The electron stationary
point is analytic; the initialization-only non-electron path uses a fixed
bisection. This is a tighter bound than multiplying the larger endpoint cross
section by the upper-endpoint speed. If the cumulative table is disabled by its
memory limit, a k-way merge constructs the same union intervals without
materializing the table. Neither path steps through the total energy span using
the smallest input spacing. For electrons, the constant high-energy table
extrapolation is also bounded using :math:`g_{\mathrm{coll}}<c`.

A user can bypass automatic construction with a collision-level majorant in
:math:`\mathrm{s}^{-1}`::

    mcc.nu_max = 1.0e12

or with ``picmi.MCCCollisions(..., nu_max=1.0e12)``. The supplied value must
bound the sum of all configured process frequencies for every particle state
and background density encountered by that MCC object; an underestimated value
biases the collision probabilities. This includes the relative speeds produced
by thermal-neutral sampling. WarpX checks the bound for sampled collision
candidates and aborts if it is violated. Each MCC object retains its own
majorant.
The candidate probability is recalculated from the current collision timestep
as :math:`P_{\max}=1-\exp(-\nu_{\max}\Delta t_{\mathrm{coll}})`, so
collision subcycling and variable timesteps use the appropriate probability.

For non-electron projectiles, the existing center-of-momentum collision-energy
convention from ``ParticleUtils::getCollisionEnergy()`` is retained. Its inverse
uses both projectile and neutral masses when constructing the automatic
majorant. The collision rate uses the corresponding ordinary relative speed,
not proper speed. This distinction is negligible in the intended
non-relativistic ion regime, but keeps the rate physical and bounded by
:math:`c`. The automatic non-electron majorant covers both the tabulated range
and the endpoint-clamped high-energy continuation. It bounds the latter by the
last total cross section multiplied by :math:`c`; a user ``nu_max`` is not
needed solely because a non-electron projectile can exceed the last table
energy.

The configured ``background_mass`` always denotes the neutral target mass and
is kept separate from the mass of any ionization product species. If it is
omitted for an ionizing electron-neutral collision, the neutral mass is inferred
from the positive-ion product mass plus one electron mass. Attachment-only MCC
objects must specify ``background_mass`` because a dissociative negative-ion
product does not uniquely determine the original neutral mass.

Ionization and attachment process names may have unique suffixes, for example
``ionization_N2`` and ``attachment_O2_dissociative``. Each such process names
its destination species with ``<process>_species``. An ionization event keeps
the incident electron, creates one additional electron and one singly charged
positive ion, while an attachment event consumes the incident electron and
creates one singly charged negative ion. Created macroparticles inherit the
incident electron weight and receive new particle IDs. Positive and negative
product species must have charge ``+q_e`` and ``-q_e``, respectively.

Attachment cross-section units are explicit. Set
``<process>_cross_section_units = m2`` when the table already contains an
effective two-body cross section. For a raw three-body table, set
``<process>_cross_section_units = m5`` and provide a positive
``<process>_third_body_density`` in :math:`\mathrm{m}^{-3}`. WarpX multiplies
the raw table by this density exactly once before constructing the majorant.
The unscaled values are retained in double precision so that small
:math:`\mathrm{m}^{5}` values do not underflow in single-precision builds.

Ionization energy sharing and kinematics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user-supplied integral cross-section table always determines the
ionization event rate. By default, WarpX subtracts the configured threshold
``<process>_energy`` and shares the remaining energy equally between the two
outgoing electrons. The per-process ``RBEQ`` energy-sharing model instead uses
the relativistic binary-encounter Bethe, or RBEQ, singly differential cross
section described by :cite:t:`b-Schmalzried2023`::

    mcc.ionization_N2_energy_sharing_model = RBEQ
    mcc.ionization_N2_rbeq_target = N2
    mcc.ionization_N2_energy = 15.58

``RBEQ`` is available for ``N2`` and ``O2``. WarpX first selects one of five
:math:`\mathrm{N}_2` or six :math:`\mathrm{O}_2` subshells from its positive
RBEQ partial cross section. The selected binding energy :math:`B_i`, rather
than only the outer-shell threshold, is then removed from the incident energy.
The lower-energy electron is sampled from the conditional RBEQ distribution on

    .. math::

       0 \leq T_s \leq \frac{T-B_i}{2},

and the other electron receives the remainder, apart from molecular-ion recoil.
The configured process threshold must match the target outer-shell binding
energy within 0.05 eV: 15.58 eV for :math:`\mathrm{N}_2` and 12.07 eV for
:math:`\mathrm{O}_2`.

Some unit-oscillator-strength RBEQ partials become slightly negative immediately
above their subshell thresholds. A negative partial is unphysical and cannot
define a probability, so WarpX clamps it to zero until the analytic expression
becomes positive; it is never reflected with an absolute value. WarpX stores
the analytic zero-crossing coordinate and interpolates the unnormalized,
non-negative partial cross sections. Consequently, interpolation cannot activate
a shell below its crossing, and normalization occurs only after interpolation
at the collision energy.

For several of these shells, the integrated partial becomes positive slightly
before the published SDCS is non-negative over its complete kinematic interval.
A non-monotone cumulative function cannot be sampled as a probability. In that
narrow interval, WarpX uses a uniform conditional energy distribution and
switches to RBEQ once the complete SDCS is non-negative. This continuation is
also used within 0.1 percent of every binding threshold, where direct evaluation
is ill-conditioned. It preserves a finite, symmetric threshold limit without
turning negative values into artificial positive probability.

The partial cross sections and conditional inverse CDFs are precomputed on
logarithmic energy grids during initialization. Subshell selection uses 2,049
energy points, while the much larger inverse-CDF data use 257. The two grids
cover the same range, so the shell coordinate is obtained from the inverse-CDF
coordinate by one fixed scale factor rather than another logarithm or search.
This independently resolves sharp shell onsets without multiplying the
inverse-CDF memory footprint. The inverse CDF uses 513 samples of the symmetric
probability map

    .. math::

       q(x) = \frac{x^4}{x^4+(1-x)^4}, \qquad 0 \leq x \leq 1,

which resolves both high-energy probability tails without linearly
interpolating from an interior quantile to the physical endpoint. The tables
cover incident energies through at least 1 GeV. Each accepted event then
requires one logarithm and fixed-size interpolation, with no device allocation
or iterative root solve.

The angular model is selected independently from energy sharing. With
``<process>_scattering_angle_model = IAA``, the collision is evaluated in the
sampled neutral rest frame. Let :math:`T_a=T-B_i`, and let :math:`T_p` and
:math:`T_s` be the unrecoiled primary and secondary energies. The primary
electron follows the relativistic binary-encounter relation

    .. math::

       \cos\theta_p =
       \sqrt{\frac{T_p(T_a+2m_ec^2)}{T_a(T_p+2m_ec^2)}}.

The bound/free interpolation for the secondary electron follows
Schmalzried Eq. (11.132) directly:

    .. math::

       \cos\theta_s =
       \frac{T_s}{T_s+B_i}
       \frac{T_s+B_i/2}{\sqrt{T_s T}}
       + \frac{B_i}{T_s+B_i}\xi,
       \qquad \xi\sim\mathcal{U}[-1,1].

Here :math:`T` is the incident electron energy before the binding loss. The
sampled result is restricted to the physical cosine interval
:math:`[-1,1]`; this affects only the upper edge of the crude near-threshold
square-window model.

The electron azimuths differ by :math:`\pi`; the product ion receives the
remaining momentum. WarpX solves for the total electron energy after ion recoil
with three fixed Newton updates. The analytic derivative includes both outgoing
electron momenta, while a stable :math:`pc` representation avoids subtracting
the ion rest energy. Thus even deliberately adverse backward and isotropic
events retain energy and momentum conservation through 1 GeV without a
device-side convergence loop. The three products are then Lorentz transformed
back to the simulation frame.

For every electron process with a positive discrete loss :math:`Q`, finite
target recoil raises the stationary-target laboratory threshold above
:math:`Q` to

    .. math::

       T_{\mathrm{thr}} = Q\left(1+\frac{m_e}{M}\right)
       + \frac{Q^2}{2Mc^2}.

An integral cross-section table can rise immediately above :math:`Q`, leaving
a narrow interval in which interpolation gives a nonzero value but the final
state is kinematically forbidden. WarpX treats a selection in
:math:`Q<T<T_{\mathrm{thr}}` as a null event. This guard applies to excitation,
ionization and any other positive-loss electron channel, independently of its
angular or energy-sharing model.

IAA/elmolcs differential scattering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For elastic and excitation electron scattering, ``IAA`` selects an angular
differential cross-section (DCS) table such as the :math:`\mathrm{N}_2` and
:math:`\mathrm{O}_2` tables from the
`elmolcs project <https://codeberg.org/jesenek/elmolcs>`_::

    mcc.elastic_N2_scattering_angle_model = IAA
    mcc.elastic_N2_differential_cross_section = /path/to/DCS.e-N2
    mcc.excitation_N2_scattering_angle_model = IAA
    mcc.excitation_N2_differential_cross_section = /path/to/DCS.e-N2
    mcc.excitation_N2_energy = 6.17

This DCS table is separate from ``<process>_cross_section``: the latter controls
the event rate, while the DCS controls only the conditional scattering angle.
WarpX reads the elmolcs ``DCS.e-N2`` and ``DCS.e-O2`` files directly. It ignores
header, metadata and separator rows whose first token is not numeric. Each
numeric row contains a positive, strictly increasing energy in eV followed by
at least three non-negative angular values. All numeric rows must have the same
number of values, uniformly spaced from :math:`0` to :math:`\pi`. The elmolcs
tables contain 361 values at 0.5-degree spacing and extend from their low-energy
endpoint through 1 GeV. Their DCS units are
:math:`10^{-20}\,\mathrm{m}^2\,\mathrm{sr}^{-1}`, although a common positive
scale cancels from angular sampling.

For a DCS :math:`D(E,\theta)`, WarpX constructs the polar-angle density

    .. math::

       p(\theta\mid E) =
       \frac{D(E,\theta)\sin\theta}
            {\int_0^\pi D(E,\vartheta)\sin\vartheta\,d\vartheta}.

Below 10 keV, WarpX follows the IAA interpolation variables
:math:`x=\log E` and :math:`y=\sin(\theta/2)`. It treats each DCS row as
piecewise linear in :math:`y` and integrates the solid-angle density
:math:`yD(E,y)` exactly on every angular interval. A tail-resolving inverse CDF
is precomputed for every table energy on the host and copied to the device. The
inverse table stores :math:`1-\cos\theta=2y^2` rather than
:math:`\cos\theta`, retaining small forward deflections in single-precision
particle builds.

At runtime, the DCS linear in :math:`x` is sampled as an exact mixture of its
two bracketing row distributions. If :math:`f` is the logarithmic energy
fraction and :math:`I_0,I_1` are the precomputed row integrals, the row weights
are :math:`(1-f)I_0` and :math:`fI_1`. One bisection, one logarithm and one
uniform variate therefore select a row and its conditional quantile. Only two
adjacent values from that row's inverse CDF are loaded; interpolation does not
blend inverse angles from different energies.

The 0.5-degree elmolcs grid cannot resolve the increasingly narrow forward
lobe at relativistic energies. Consequently, increasing only the inverse-CDF
resolution cannot recover the missing sub-grid probability. For files whose
``SPECIES:`` metadata identifies :math:`\mathrm{N}_2` or :math:`\mathrm{O}_2`,
WarpX follows the IAA high-energy prescription at and above 10 keV and samples
the screened-Rutherford continuation analytically:

    .. math::

       \eta = \frac{1}{(2ak)^2}, \qquad
       1-\cos\theta = \frac{2\eta\xi}{1-\xi+\eta},
       \qquad \xi\sim\mathcal{U}[0,1),

where :math:`k=\sqrt{\tau(\tau+2)}/\alpha` in inverse Bohr radii,
:math:`\tau=T/(m_ec^2)`, and the fitted screening radii are
:math:`a=0.6052\,a_0` for :math:`\mathrm{N}_2` and
:math:`a=0.5677\,a_0` for :math:`\mathrm{O}_2`. This path is valid through the
1 GeV IAA endpoint and avoids an energy bisection and inverse-CDF table reads
for accepted high-energy events. WarpX retains only the tabulated rows through
the first row at or above 10 keV, so unused relativistic rows consume neither
initialization time nor device memory.

The sampled DCS angle is the outgoing-electron angle for a stationary target,
not a center-of-momentum angle. For elastic scattering, WarpX applies exact
relativistic two-body recoil in the sampled neutral rest frame. For excitation,
the target's final rest energy is increased by the configured discrete loss
:math:`Q`. The outgoing momentum magnitude is the analytic two-body solution of

    .. math::

       T + m_ec^2 + Mc^2 = T' + m_ec^2
       + \sqrt{(Mc^2+Q)^2 + c^2\lvert\boldsymbol{p}-\boldsymbol{p}'\rvert^2},

for the sampled angle. This conserves energy and momentum, including molecular
recoil, before transforming the electron back to the simulation frame. The
model therefore requires an electron projectile and a neutral target heavier
than the electron. Above the recoil-shifted physical threshold, an angular draw
outside the allowed two-body laboratory range is projected to its nearest
kinematically allowed value.

Many-channel selection and DCS reuse
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cumulative selector covers **all** processes listed in one Background MCC
object: elastic, excitation, ionization, attachment, charge exchange and
two-product reactions. It is not an excitation-only table. The original input
order defines process indices :math:`i=0,\ldots,P-1`. During initialization,
WarpX constructs the sorted union

    .. math::

       \mathcal{E} = \operatorname{unique}\!\left(
       \bigcup_{i=0}^{P-1}\mathcal{E}_i\right)
       = \{E_0,\ldots,E_{U-1}\}

of the knots from every integral cross-section grid, then stores the prefix
rows

    .. math::

       C_j(E_k) = \sum_{i=0}^{j}\sigma_i(E_k),
       \qquad j=0,\ldots,P-1.

Each :math:`\sigma_i` is evaluated with the same piecewise-linear interpolation
and endpoint clamping as an ordinary ``ScatteringProcess`` lookup. Every input
breakpoint is present in :math:`\mathcal{E}`, so all :math:`\sigma_i` and all
:math:`C_j` are linear between adjacent union knots. Linear interpolation of a
prefix row is therefore exact relative to the stored cross-section
representation, apart from floating-point roundoff; resampling onto the union
does not smooth a threshold or otherwise approximate a channel.

For each particle that passes the global null-collision preselection, WarpX
samples one neutral velocity and calculates one collision energy. Process
selection then consists of:

#. one bisection of :math:`\mathcal{E}` to obtain the bracketing index
   :math:`k` and interpolation fraction;
#. one interpolation of the last prefix row
   :math:`C_{P-1}(E)=\sigma_{\mathrm{tot}}(E)`;
#. one uniform draw converted from collision-probability units to
   cross-section units; and
#. one binary search over :math:`j` for the first interpolated prefix with
   :math:`C_j(E)` above that draw.

Thus, yes: when every process was supplied on a common grid, that grid is also
the union and the cross-section energy interval is found only once per
candidate particle. The same remains true when the input grids differ, because
they were combined at initialization. The selector cost is
:math:`O(\log U+\log P)` rather than :math:`P` separate energy bisections plus a
linear process scan. For example, 64 processes require at most six prefix
comparisons after the energy search.

"One energy search" refers specifically to the integral cross sections used
for event acceptance and process choice. A selected IAA angular DCS has its own
independent energy grid and performs one additional DCS bisection, but only for
the chosen process. RBEQ energy sharing likewise uses its own fixed logarithmic
coordinate. WarpX never searches the angular or RBEQ tables for processes that
were not selected.

The prefix table contains cumulative cross-section weights only; it does not
contain or average energy losses. The binary search returns the original
process index :math:`j`.
WarpX then reads that process's unchanged discrete ``<process>_energy``, type,
angular model, ionization model and product-species group. An excitation event
therefore subtracts exactly the loss configured for its selected channel. An
ionization or attachment event records the same process index for the grouped
particle-creation pass, so cumulative selection cannot disconnect a product
from its channel.

The storage in one host or device copy is

    .. math::

       (P+1)U\,\mathrm{sizeof}(\mathtt{ParticleReal}),

where the extra row is the union energy grid and the :math:`P` other rows are
the prefixes. WarpX enables the selector only when this total is at most
64 MiB. Consequently,

    .. math::

       U_{\max} = \left\lfloor
       \frac{64\,\mathrm{MiB}}
            {(P+1)\,\mathrm{sizeof}(\mathtt{ParticleReal})}
       \right\rfloor.

For 64 processes this permits 129,055 union points with double-precision
particles or 258,111 with single-precision particles. For 128 processes the
limits are 65,027 and 130,055 points, respectively. Typical cross-section sets
are far smaller: if 64 channels each use the same 100-point grid, then
:math:`U=100`, not 6,400. Only completely distinct knots can make :math:`U`
approach the sum of the individual grid sizes. A GPU build retains one host
copy and one device copy, each subject to the 64 MiB limit. Original process
grids, angular DCS tables and RBEQ tables are separate allocations.

If the limit would be exceeded, WarpX does not truncate or coarsen any grid. It
uses the exact legacy fallback, which evaluates each process cross section and
scans the cumulative probabilities at runtime. The first-step information
message reports the active selection path, union point count and bytes per
copy, making an accidental fallback visible in production logs.

Prefix rows are stored process-major. Particles in a GPU warp begin each
process binary search at the same middle prefix, while nearby particle energies
usually address nearby entries within that row. All selector arrays are fixed,
contiguous device data; selection performs no allocation, virtual dispatch or
iterative solve in the particle kernel.

Elastic and excitation channels that name the same DCS file share one
precomputed inverse-CDF table and one device allocation. This avoids duplicating
initialization work and device memory when many excitation channels use the
same target DCS. Use the same path string for channels that should share a
table.

The elmolcs data are distributed separately from WarpX and are not copied into
the BSD-licensed source tree. Users must obtain the tables separately and comply
with their license. The construction and intended range of the IAA database are
documented in :cite:t:`b-Schmalzried2023`.

Product creation and attachment removal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ionization and attachment preserve the original selected process index in a
compact product-event record. The record also contains the source-particle
index, sampled neutral velocity, collision energy and destination-group offset.
Consequently, the creation kernel still reads the selected channel's discrete
energy loss, ionization model and angular model; grouping never averages channel
physics. Carrying the already calculated collision energy also avoids repeating
its square root for every accepted ionization event.

On a GPU, the selection kernel reserves one slot in the compact event queue and
one offset in the destination group for each product-changing event. These are
device-local atomic increments. At the recommended
:math:`\nu_{\max}\Delta t_{\mathrm{coll}}\lesssim 0.1`, fewer than about ten
percent of the source particles can pass even the null-collision preselection,
and product-changing events are a subset of those candidates. The atomics are
therefore sparse in the converged operating regime. On a CPU, selection first
writes independent per-source records; one host pass both compacts those
records and computes the group offsets, without unsafe updates inside an
``amrex::ParallelFor``.

The host receives only the small vector of per-group event counts. Each source
or destination particle tile is then resized once. All destination groups and
all secondary electrons are populated by one kernel over the number of compact
events, not by a kernel over every source particle and not by one source scan
per process or product group. Contiguous particle-ID ranges are reserved before
that kernel and assigned within it, avoiding separate ID kernels. Immutable
smart-copy and product-group metadata are constructed once and reused across
collision substeps.

An attachment marks its incident electron invalid. WarpX compacts only the
source tile that contained an attachment and only when that tile actually had
an attachment event. It does not rescan every tile merely because an attachment
channel is configured. Each MPI rank and GPU performs this creation and removal
on its local particle tiles; the algorithm introduces no inter-rank collective
or cross-GPU synchronization. A device-to-host count transfer remains necessary
before resizing a particle tile, because the new allocation size is a host-side
container operation.

GPU performance validation
^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Collision timestep
^^^^^^^^^^^^^^^^^^

Background MCC permits at most one event per particle per collision substep.
Use ``ndt_subcycle`` and verify convergence when the collision optical depth is
not small. This remains important for quantitative air simulations.

Once a particle is selected for a specific collision process, that process determines how the particle is scattered as outlined below.

.. _multiphysics-collisions-proton-impact-ionization:

Proton and bare-ion impact ionization
-------------------------------------

The ``proton_impact_ionization`` collision type creates an electron and a
singly charged molecular ion when a kinetic proton or bare-ion beam traverses
a prescribed :math:`\mathrm{N}_2` or :math:`\mathrm{O}_2` background. The
initial production model is a source model: the projectile is copied only to
initialize product positions and attributes, and its momentum and weight are
left unchanged. The neutral background is not depleted. For example::

    collisions.collision_names = n2_ionization
    n2_ionization.type = proton_impact_ionization
    n2_ionization.species = beam
    n2_ionization.product_species = electrons n2_ions
    n2_ionization.ionization_target = N2
    n2_ionization.background_density = 2.5e25
    n2_ionization.background_temperature = 300.0
    n2_ionization.fixed_product_weight = 1.0e7
    n2_ionization.max_products_per_cell = 64

The corresponding PICMI class is
``picmi.ProtonImpactIonizationCollisions``.

Corrected Porter--Jackman--Green model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The integral ionization cross section and the singly differential cross
section (SDCS) in ejected-electron kinetic energy :math:`T` use the
Porter--Jackman--Green (PJG) model :cite:t:`b-Porter1976`. With
:math:`E` and :math:`M` the projectile kinetic energy and mass,
:math:`E_e=m_e\beta^2c^2/2`, and continuum threshold :math:`I_j`, the
implemented form is

.. math::

   \begin{aligned}
   S(E,T) = \sum_j \frac{F_j(E)}{E_e}\Bigg\{&
   K\Gamma^2 L_j
   \left[
   \frac{1}{(T-T_0)^2+\Gamma^2}
   -\frac{B(E_e)}{(T-T_1)^2+\Gamma_1^2}
   \right] \\
   &+N_e\pi e_{\mathrm c}^4\left[
   \frac{1}{2(E+Mc^2)^2}
   -\frac{\beta^2}
   {(T_{\max}+I_j)(T+I_j)}
   \right]\Bigg\},
   \\
   L_j ={}& \ln\left(
   \frac{4E_eC_j}{I_j(1-\beta^2)}+\mathrm e
   \right)-\beta^2,
   \\
   F_j(E) ={}& f_j\frac{E_e^{\nu+1}}
   {J^{\nu+1}+E_e^{\nu+1}}.
   \end{aligned}

Here :math:`e_{\mathrm c}` is the elementary charge in the Gaussian units of
the original fit and :math:`\mathrm e` is Euler's number. Except for the
refitted parameters described below, the energy-dependent functions and the
:math:`\mathrm{N}_2` and :math:`\mathrm{O}_2` parameters are those in PJG
Table III. WarpX converts the result to SI units. It analytically integrates
every continuum over the full range :math:`0\leq T\leq T_{\max}` to construct
the total cross section; no low-energy ejected-electron cutoff is imposed.

The maximum transferable kinetic energy printed below PJG Eq. (16) is not the
relativistic two-body result. WarpX instead uses

.. math::

   T_{\max} =
   \frac{2m_ec^2\beta^2\gamma^2}
   {1+2\gamma m_e/M+(m_e/M)^2}.

This has the required non-relativistic limit
:math:`T_{\max}\rightarrow 4m_eME/(M+m_e)^2\simeq 4m_eE/M`; the printed PJG
expression is smaller by approximately a factor of two in that limit. For an
800 MeV proton, the exact value is 2.4807396 MeV, whereas the printed expression
gives 0.670701 MeV.

The last bracket in the SDCS is also corrected against Bhabha's relativistic
heavy-particle result :cite:t:`b-Bhabha1938`. After extracting the common PJG
factor :math:`1/E_e`, the spin term is
:math:`1/[2(E+Mc^2)^2]` and the inverse-energy term contains :math:`\beta^2`.
The published PJG expression has neither dependence correctly. The empirical
:math:`\delta` offset in the upper denominator is not part of the free-electron
Bhabha result, and removing it changes the total cross section by less than
:math:`3.1\times10^{-5}` in relative value over 2 keV--1 MeV. It is therefore
set to zero rather than retained as an unidentifiable parameter. The threshold
shift :math:`I_j` remains the PJG bound-electron continuation.

The 1977 correction notice is contained in
:cite:t:`b-Garvey1977`. The complete correction list was checked: in Eq. (14),
:math:`\arctan\alpha_2` becomes :math:`\arctan(\alpha_2/2)` and the denominator
:math:`\beta_1` becomes :math:`\beta_2`; Fig. 1(a) has :math:`w=58.2`; Table II
adds :math:`b_5` to the :math:`\beta_2` equation and :math:`g_5` to the
:math:`\gamma` equation, with :math:`b_5=0.654` for :math:`\mathrm{N}_2`,
:math:`b_5=0.619` for :math:`\mathrm{O}_2`, and :math:`g_5=1.0` for both;
the :math:`\mathrm{O}_2` value is :math:`a_2=4.107\times10^{-1}`; the mark in
Eq. (16) is an ordinary plus sign, not a superscript; and the
:math:`\mathrm{O}_2` Table III value is :math:`E_\Gamma=-129.1` eV. The
Eq. (14), Table II, and :math:`E_\Gamma` corrections belong to the
electron-impact fits and do not enter the proton Eq. (16) implementation. The
Eq. (16) plus-sign correction does enter its parsing.

As an additional source audit, a direct evaluation of Eq. (16) with the Table
III proton parameters does not reproduce the low-energy curve in PJG Fig. 6.
The discrepancy is not repaired by restoring the erroneous printed Bhabha
factors or maximum-transfer expression. The most likely source is therefore
an undocumented difference between the calculation used for the figure and
the final equation or parameter table, rather than the 1977 correction notice.

Refit to recommended proton data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The exact :math:`T_{\max}` and Bhabha factors above are fixed, not fitted. With
:math:`\delta=0`, the target-wide distortion scale :math:`J` was first fitted
by itself to the recommended total cross sections of Rudd et al.
:cite:t:`b-Rudd1983,b-Rudd1985`. The fit minimizes the mean squared logarithmic
residual at 301 uniformly log-spaced energies over 5--4000 keV, the range of
the underlying Rudd et al. measurements. The reference curve is

.. math::

   \sigma_{\mathrm R}(E) = 4\pi a_0^2
   \left[
   \frac{U}{A\ln(1+U)+B}+\frac{1}{CU^D}
   \right]^{-1},
   \qquad
   U=\frac{E}{(m_p/m_e)R_{\mathrm y}},
   \quad R_{\mathrm y}=13.6057\ \mathrm{eV}.

Rudd et al. give :math:`(A,B,C,D)=(3.82,2.78,1.80,0.70)` for
:math:`\mathrm{N}_2` and :math:`(4.77,0,1.76,0.93)` for
:math:`\mathrm{O}_2`. Their recommended curves result from a complete
literature survey whose fit weights account for quoted and estimated
uncertainties, independent normalization, number of points, and covered
energy span.

For :math:`\mathrm{N}_2`, changing only :math:`J` is sufficient and gives
:math:`J=59.65` eV; :math:`K=7.58\times10^{-16}\ \mathrm{cm}^2` is retained.
For :math:`\mathrm{O}_2`, a :math:`J`-only fit leaves a 47.5 percent maximum
deviation because :math:`F_j\rightarrow f_j` at high energy. The minimal
additional identifiable change is the soft-continuum normalization, giving
:math:`J=19.28` eV and
:math:`K=4.581\times10^{-16}\ \mathrm{cm}^2`. All :math:`\nu`, continuum
fractions, thresholds, Bethe constants, and line-shape parameters remain
unchanged. The resulting root-mean-square logarithmic residuals are 0.1245 for
:math:`\mathrm{N}_2` and 0.0628 for :math:`\mathrm{O}_2`; the largest
symmetric multiplicative errors are 23.0 and 25.8 percent, respectively.
These residuals are comparable to the reliability assigned to the low- and
peak-energy molecular data by Rudd et al.; no extra shape parameter is
justified by those data.

Reference unit-charge proton totals used by the independent physics test are:

.. list-table:: Corrected PJG total ionization cross sections (:math:`\mathrm{m}^2`)
   :header-rows: 1
   :widths: 25 25 25

   * - Proton energy
     - :math:`\mathrm{N}_2`
     - :math:`\mathrm{O}_2`
   * - 50 keV
     - :math:`4.9571136\times10^{-20}`
     - :math:`5.3533138\times10^{-20}`
   * - 200 keV
     - :math:`3.4782052\times10^{-20}`
     - :math:`3.9221672\times10^{-20}`
   * - 500 keV
     - :math:`2.2309605\times10^{-20}`
     - :math:`2.4668388\times10^{-20}`
   * - 1 MeV
     - :math:`1.4438120\times10^{-20}`
     - :math:`1.5732801\times10^{-20}`
   * - 10 MeV
     - :math:`2.3912617\times10^{-21}`
     - :math:`2.4865961\times10^{-21}`
   * - 100 MeV
     - :math:`3.5000091\times10^{-22}`
     - :math:`3.2936513\times10^{-22}`
   * - 800 MeV
     - :math:`1.0497733\times10^{-22}`
     - :math:`8.9033261\times10^{-23}`

Energy and angle sampling
^^^^^^^^^^^^^^^^^^^^^^^^^

Initialization tabulates the total cross section and conditional SDCS inverse
CDF on 256 logarithmically spaced projectile energies and 513 transformed
quantiles. The quantile coordinate uses

.. math::

   q(x)=\frac{x^4}{x^4+(1-x)^4},\qquad 0\leq x\leq1,

so both endpoints, particularly the long relativistic Bhabha tail, are
resolved. Ejected energy is interpolated in :math:`\log(1+T)`. A device event
then uses one projectile-energy logarithm, one exponential, and fixed-size
interpolation, with no search or root solve. The double-precision tables use
about 2.1 MB per collision object, and single-precision particle builds use
about half that amount.

The ejected-electron angle uses a heavy-projectile extension of the practical
IAA bound/free closure described by :cite:t:`b-Schmalzried2023`. The exact
free-electron binary-encounter polar cosine is

.. math::

   \cos\theta_{\mathrm{free}} =
   \sqrt{\frac{T(T_{\max}+2m_ec^2)}
   {T_{\max}(T+2m_ec^2)}}.

For effective continuum binding energy :math:`I`, WarpX samples

.. math::

   \cos\theta = \operatorname{clamp}_{[-1,1]}\left[
   \cos\theta_{\mathrm{free}}\frac{T+I/2}{T+I}
   +\frac{I}{T+I}\xi\right],
   \qquad \xi\sim\mathcal U[-1,1].

The effective :math:`I` is the SDCS-weighted mean of the positive continuum
contributions at the sampled energy. This closure recovers the free
binary-encounter direction as :math:`I/T\rightarrow0` and an isotropic bound
component as :math:`T/I\rightarrow0`; it is not a molecular angular DCS. The
azimuth is uniform about the projectile direction.

Low-noise macroparticle source and GPU execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In cell :math:`c`, one collision call produces the physical pair weight

.. math::

   W_c = n_{n,c}\Delta t\sum_{p\in c}
         w_p Z_p^2\sigma_{\mathrm{PJG}}(E_p)v_p,

where :math:`Z_p` is the positive integer projectile charge state. A
checkpointed cell residual carries fractional product weight between calls.
WarpX creates
:math:`\lfloor(W_c+R_c)/w_{\mathrm{fixed}}\rfloor` electron--ion pairs and
samples their parent projectiles systematically with probability proportional
to :math:`w_p\sigma(E_p)v_p`. Ejected-energy quantiles are stratified within
each cell. This removes Bernoulli event-count noise and is especially useful
when the expected yield per beam macroparticle is small.

``max_products_per_cell`` bounds work and memory growth in one call. When the
bound is reached, WarpX creates exactly that many pairs with equal increased
weights, preserving the complete accumulated physical weight rather than
discarding events. The electron and molecular ion in a pair have the same
position and weight. The ion velocity is sampled from the zero-drift neutral
Maxwellian at ``background_temperature``; it is deliberately not assigned the
event recoil. Every sampled electron, including arbitrarily low-energy ones,
is represented kinetically.

The implementation uses dense particle bins, one independent thread per cell,
a device exclusive scan, and one product-creation kernel. It uses no scatter
updates, atomics, device-side allocation per event, or host callback per event.
Per-cell scratch uses two device allocations per tile and collision call. A
host-visible product count remains necessary before particle-tile resizing.
The performance regression input
``inputs_test_3d_proton_impact_ionization_performance_picmi.py`` accepts
``--cells-per-direction``, ``--particles-per-cell-direction``, ``--steps`` and
``--max-products-per-cell`` for accelerator profiling. CPU timings are only a
regression check; GPU speedups must be measured on the intended CUDA, HIP, or
SYCL system.

Scope and limitations
^^^^^^^^^^^^^^^^^^^^^

Only single total ionization of :math:`\mathrm{N}_2` and
:math:`\mathrm{O}_2` is implemented. The PJG table is a unit-charge proton
model. Bare ions use the usual :math:`Z_p^2` scaling with the actual projectile
mass and velocity; electron screening, charge exchange, and corrections for
large :math:`Z_p` are outside this model. The default lookup interval is 1 keV
to 1 GeV, and the cross section is zero outside the configured interval.

Because projectile momentum and energy, neutral depletion, and ion recoil are
all omitted, total energy and momentum of the represented particles are not
closed: the beam and background act as external reservoirs. The background
density and temperature expressions are evaluated at cell centers. Resolve
their spatial variation and verify convergence with collision timestep, cell
size, ``fixed_product_weight``, and ``max_products_per_cell``.

.. _multiphysics-collisions-pulseddecay:

Pulsed Decay
------------

This collision module can be used to have a parent species decay into two product
species with a user-defined decay rate. Mathematically, it solves the following
rate equations on a cell-by-cell basis:

    .. math::

       \begin{aligned}
        \frac{dn_1}{dt} &= -\nu(t)n_1, \\
        \frac{dn_A}{dt} &= +\nu(t)n_1 = \frac{dn_B}{dt}, \\
       \end{aligned}

where :math:`n_1` is the parent species density, :math:`n_A` and :math:`n_B` are the product species densities,
and :math:`\nu(x,y,z,t)` is the user-specified decay rate.

This can be used, for example, to represent ionization of a parent species by an externally applied laser pulse.

.. _multiphysics-collisions-dsmc:

Direct Simulation Monte Carlo (DSMC)
------------------------------------

The algorithm by which binary collisions are treated is outlined below. The
description assumes collisions between different species.

1. Particles from both species are sorted by grid-cells.
2. The order of the particles in each cell is shuffled.
3. Within each cell, particles are paired to form collision partners. Particles
   of the species with fewer members in a given cell is split in half so that
   each particle has exactly one partner of the other species.
4. Each collision pair is considered for a collision using the same logic as in
   the MCC description above.
5. Particles that are chosen for collision are scattered according to the
   selected collision process.

Scattering processes
--------------------

Charge exchange
^^^^^^^^^^^^^^^

This process can occur when an ion and neutral (of the same species) collide
and results in the exchange of an electron. The ion velocity is simply replaced
with the neutral velocity and vice-versa.

Elastic scattering
^^^^^^^^^^^^^^^^^^

The ``elastic`` option uses isotropic scattering, i.e., with a differential
cross section that is independent of angle.
This scattering process as well as the ones below that relate to it, are all
performed in the center-of-momentum (COM) frame. Designating the COM velocity of
the particle as :math:`\boldsymbol{u}_c` and its labframe velocity as :math:`\boldsymbol{u}_l`,
the transformation from lab frame to COM frame is done with a general Lorentz
boost (see function ``ParticleUtils::doLorentzTransform()``):

    .. math::
            \begin{bmatrix}
                \gamma_c c \\
                u_{cx} \\
                u_{cy} \\
                u_{cz}
            \end{bmatrix}
         = \begin{bmatrix}
                \gamma & -\gamma\beta_x & -\gamma\beta_y & -\gamma\beta_z \\
                -\gamma\beta_x & 1+(\gamma-1)\frac{\beta_x^2}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} \\
                -\gamma\beta_y & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & 1 +(\gamma-1)\frac{\beta_y^2}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} \\
                -\gamma\beta_z & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} & 1+(\gamma-1)\frac{\beta_z^2}{\beta^2} \\
            \end{bmatrix} \begin{bmatrix}
                \gamma_l c \\
                u_{lx} \\
                u_{ly} \\
                u_{lz}
            \end{bmatrix}

where :math:`\gamma` is the Lorentz factor of the relative speed between the lab frame and the COM frame, :math:`\beta_i = v^{COM}_i/c` is the i'th component of the relative velocity between the lab frame and the COM frame with

    .. math::

        \boldsymbol{v}^{COM} = \frac{m \boldsymbol{u}_c}{\gamma_u m + M}

The particle velocity in the COM frame is then isotropically scattered using the function ``ParticleUtils::RandomizeVelocity()``. After the direction of the velocity vector has been appropriately changed, it is transformed back to the lab frame with the reversed Lorentz transform as was done above followed by the reverse Galilean transformation using the starting neutral velocity.

Back scattering
^^^^^^^^^^^^^^^

The process is the same as for elastic scattering above except the scattering angle is fixed at :math:`\pi`, meaning the particle velocity in the COM frame is updated to :math:`-\boldsymbol{u}_c`.

Excitation
^^^^^^^^^^

The process is also the same as for elastic scattering except the excitation energy cost is subtracted from the particle energy. This is done by reducing the velocity before a scattering angle is chosen.

Forward scattering
^^^^^^^^^^^^^^^^^^

This process operates in two ways:

1. If an excitation energy cost is provided, the energy cost is subtracted from the particle energy and no scattering is performed.
2. If an excitation energy cost is not provided, the particle is not scattered and the velocity is unchanged (corresponding to a scattering angle of :math:`0` in the elastic scattering process above).

See :cite:t:`b-Janssen2016` for a recommended use of this process.

Benchmarks
----------

See the :ref:`MCC example <examples-capacitive-discharge>` for a benchmark of the MCC
implementation against literature results.

Particle cooling due to elastic collisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is straight forward to determine the energy a projectile loses during an elastic collision with another body, as a function of scattering angle, through energy and momentum conservation.
See for example :cite:t:`b-Lim2007` for a derivation. The result is that given a projectile with mass :math:`m`, a target with mass :math:`M`, a scattering angle :math:`\theta`, and collision energy :math:`E`, the post collision energy of the projectile is given by

    .. math::

       \begin{aligned}
       E_{final} = E - &[(E + mc^2)\sin^2\theta + Mc^2 - \cos(\theta)\sqrt{M^2c^4 - m^2c^4\sin^2\theta}] \\
       &\times \frac{E(E+2mc^2)}{(E+mc^2+Mc^2)^2 - E(E+2mc^2)\cos^2\theta}
       \end{aligned}

The impact of incorporating relativistic effects in the MCC routine can be seen in the plots below where high energy collisions are considered with both a classical and relativistic implementation of MCC. It is observed that the classical version of MCC reproduces the classical limit of the above equation but especially for ions, this result differs substantially from the fully relativistic result.

.. figure:: https://user-images.githubusercontent.com/40245517/170900079-74e505a5-2790-44f5-ac84-5847eda954e6.png
   :alt: Comparison of classical and relativistic MCC collision results
   :width: 96%

   Classical v. relativistic MCC.

.. bibliography::
    :keyprefix: b-
