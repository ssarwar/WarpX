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
:math:`c`. The automatic non-electron majorant covers the tabulated energy
range. Supply ``nu_max`` when non-electron particles can exceed that range and
the endpoint-clamped total cross section is nonzero.

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
becomes positive; it is never reflected with an absolute value. Shell
probabilities and conditional inverse CDFs are precomputed on a logarithmic
energy grid during initialization. The inverse CDF uses 513 samples of the
symmetric probability map

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
remaining momentum. Two fixed recoil corrections conserve energy through
:math:`O((m_e/M)^3)` without a device-side convergence loop, and the three
products are Lorentz transformed back to the simulation frame.

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
than the electron. In the extremely narrow interval between the configured
loss and the recoil-shifted physical threshold, WarpX projects the sampled
angle to the nearest kinematically allowed value.

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
index, sampled neutral velocity and destination-group offset. Consequently,
the creation kernel still reads the selected channel's discrete energy loss,
ionization model and angular model; grouping never averages channel physics.

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

Collision timestep
^^^^^^^^^^^^^^^^^^

Background MCC permits at most one event per particle per collision substep.
Use ``ndt_subcycle`` and verify convergence when the collision optical depth is
not small. This remains important for quantitative air simulations.

Once a particle is selected for a specific collision process, that process determines how the particle is scattered as outlined below.

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
