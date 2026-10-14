.. _multiphysics-collisions-mcc:

Background Monte Carlo collisions
=================================

For molecular nitrogen and oxygen, see :ref:`theory-rotational-scattering`
for the combined vibrationally elastic and rotational model, and
:ref:`theory-iaa-air` for production excitation and attachment data. RBEQ energy
sharing is described in :ref:`electron-impact-ionization`; implementation costs and product
handling are described in :ref:`beam-air-runtime`.

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
cumulative process table in :ref:`beam-air-runtime`. Between consecutive union knots the
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

For a sampled collision state, let
:math:`\nu=n_n g_{\mathrm{coll}}\sigma_{\mathrm{tot}}(E)`. Conditional on a
candidate, the acceptance probability is

.. math::

   A(E)=\frac{1-\exp(-\nu\Delta t_{\mathrm{coll}})}{P_{\max}}.

Thus the unconditional one-event probability is
:math:`P_{\max}A=1-\exp(-\nu\Delta t_{\mathrm{coll}})`, independent of how
loose the majorant is. After acceptance, process :math:`j` has probability
:math:`\sigma_j/\sigma_{\mathrm{tot}}`. One uniform draw performs both decisions:
reject if :math:`u\geq A`, otherwise select the cross-section prefix containing
:math:`(u/A)\sigma_{\mathrm{tot}}`. ``expm1`` evaluates the small-optical-depth
limit without cancellation. Rate accumulation, conditional acceptance and
selection uniforms use double precision on the device, including with
single-precision particle storage.

The former acceptance :math:`\nu/\nu_{\max}` instead gave
:math:`(1-\exp(-\nu_{\max}\Delta t))\nu/\nu_{\max}`. Its ratio to the
linearized physical probability :math:`\nu\Delta t` was
:math:`(1-\exp(-\tau_{\max}))/\tau_{\max}`, which can be arbitrarily small
for a loose majorant, even when the physical collision rate is well resolved.
The corrected rule removes that additional finite-step bias. It does **not**
represent multiple collisions within one call, changes of particle energy
after an earlier collision in that call, or continuous resampling of a thermal
target within the step. Subcycling/timestep convergence remains necessary
when the physical collision frequency times the step is not small. The
warning based on :math:`\nu_{\max}\Delta t` is a conservative sufficient
criterion, not evidence that the actual physical rate is that large.

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
Different positive-ion product masses require an explicit neutral mass;
specifying it bypasses mass inference, not product-charge validation. This
does not make an effective dissociative channel a resolved fragmentation model.

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


.. _mcc-finite-recoil:

Finite-recoil thresholds
------------------------

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
-----------------------------------

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
Near threshold that allowed region is a **forward** cone. The squared energy
equation also admits a backward cone with negative outgoing momentum; it must
not be used. WarpX solves for the outgoing momentum directly and evaluates
kinetic energy as :math:`(p'c)^2/(\sqrt{(m_ec^2)^2+(p'c)^2}+m_ec^2)`.
The backward-angle momentum root is rationalized separately. These forms
avoid subtracting the electron rest energy when the outgoing kinetic energy
approaches zero.


.. _mcc-collision-timestep:

Collision scheduling and timestep
---------------------------------

Background MCC permits at most one event per particle per collision substep.
Use :pp:param:`<collision_name>.ndt_subcycle` and verify convergence when the
physical collision optical depth is not small. The majorant warning is a
sufficient resolution check, not an estimate of the physical rate.

The handler orders due collision substeps by their evaluation endpoint in the
elapsed PIC interval: left endpoints for explicit evolution and right endpoints
for implicit evolution after the push. Equal endpoints retain the input order.
For equal N2 and O2
subcycle counts this alternates gases at each substep; it does not advance one
gas through the entire PIC interval before the other gas acts. Product creation,
attachment removal and fluid increments complete before the next operator runs.
The operators remain split, so interleaving does not remove timestep error.
``ndt_supercycle`` retains the elapsed interval and delayed-call convention.
Explicit and implicit calls use the same physical interval, including parser
evaluation times. See :ref:`beam-air-validation` for scheduling tests.

Input options are listed in :ref:`running-cpp-parameters-collision`.
The data reader contract is in :ref:`usage-collision-data`, and the complete
beam setup is in :ref:`usage-beam-air`.
