.. _theory-rotational-scattering:

Electron--molecule elastic and rotational scattering
====================================================

``rotation_model = reciprocal_hybrid`` treats a collision of an electron with
N2 or O2 as one vibrationally elastic event. The event can leave the molecule's
rotation unchanged, excite it, or de-excite it. Its scattering angle and its
discrete rotational energy change are sampled jointly. Ordinary elastic recoil
is also applied. Vibrational and electronic excitation, ionization, attachment,
proton-impact collisions and the beam source remain separate physical models.

This is a source-constrained hybrid model, not a complete ab initio molecular
scattering calculation. Its main purpose is transport and thermal energy
exchange in air near room temperature, including a high-energy continuation
for electrons produced by a relativistic beam. The published production data
have a fixed rotational temperature of 300 K and a supported collision-energy
domain of 0--1 GeV. Experimental support does not span this entire domain.
The distinction between measurements, theoretical inputs, fitted closures and
numerical accuracy is part of the specification below.

.. toctree::
   :maxdepth: 1

   angular_models
   sources
   sampling
   alternatives

Reading and using the model
---------------------------

Start with :ref:`rotation-physics` for the physical quantities and
:ref:`rotation-dcs` for the complete piecewise N2/O2 angular construction.
:ref:`rotation-sources` records the evidence supporting each energy range.
:ref:`rotation-reciprocity` derives the reconciliation of the IAA inclusive
elastic data with state-resolved rotational data. :ref:`rotation-continuations`
states what is assumed outside the measured ranges. :ref:`rotation-algorithm`
connects those equations to MCC selection, sampling and recoil.
:ref:`rotation-data-format` specifies the production files and initialization.
:ref:`rotation-validation` gives reproducible checks and measured limitations.

Configure one elastic process per gas with its ordinary ``elastic/cross_section.txt`` as
``cross_section``, ``scattering_angle_model = IAA``, and
``rotation_model = reciprocal_hybrid``. Set ``rotation_file`` to
``MCC_cross_sections/<gas>/IAA/rotational/elastic_rotation_300K/thermal_rotation.rot``
in warpx-data, ``rotational_temperature = 300``, and normally
``rotation_sampling = alias``. The combined bundle contains its angular law;
there is no separate elastic DCS input for this option. The ordinary elastic
file checks source consistency. Its rate is replaced, rather than added, by
the bundle's rate. Do not also configure the elementary rotational files as
ordinary MCC processes. Both proton-beam examples resolve the bundle relative
to their cross-section manifest.

The neutral translational temperature may differ from the rotational
temperature. Omission of ``rotational_temperature`` uses a constant
translational temperature, which must match the bundle temperature. Spatially
or temporally varying rotational baths, evolving neutral state populations,
spin-resolved O2, isotope mixtures and vibrationally excited rotors are outside
this model. Additional fixed temperatures require separately prepared data.

The alternative ``elastic_dcs`` model samples a rotational outcome independently
of the elastic angle from state-resolved integral rates;
its integral balance does not guarantee differential balance.
``iaa_spectator`` is the N2 kinetic counterpart of the thesis mean-loss
prescription and does not impose the coupled reciprocal balance constraint.
See :ref:`rotation-alternatives` before using those compatibility models.

.. _rotation-physics:

Physical variables and the collision operator
---------------------------------------------

Inclusive elastic scattering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this model, *vibrationally elastic* means that the target's electronic and
vibrational states remain unchanged. Rotation need not remain unchanged.
The IAA residual integral cross section is constructed by subtracting the
other inelastic channels from total scattering (thesis Eq. 11.10). Its
interpretation, the unresolved DCS in Eq. 8.129, and the database comparisons
in Chapters 11--12 are consistent with an inclusive rotational sum.
Consequently, adding rotational excitation cross sections to that residual
would count some collisions twice.

Write :math:`D_{if}(E,\Omega)=d\sigma_{i\to f}/d\Omega` for a
state-resolved DCS and :math:`b_i` for the equilibrium population of initial
rotational state :math:`i`. The desired inclusive angular marginal is

.. math::
   :label: rotation-inclusive

   D_{\rm inc}(E,\Omega)=\sum_i b_i\sum_f D_{if}(E,\Omega),\qquad
   \sigma_{\rm inc}(E)=\int D_{\rm inc}\,d\Omega.

The supplied IAA elastic DCS and the residual integral have independently
chosen normalizations. The production construction therefore uses

.. math::

   F(E,\Omega)=\frac{D_{\rm IAA}(E,\Omega)}
                       {\int D_{\rm IAA}(E,\Omega')d\Omega'},\qquad
   D_{\rm inc}=\sigma_{\rm residual}(E)F(E,\Omega).

This preserves the source's chosen total rate and angular shape separately.
It does not assume that integrating the unmodified DCS reproduces the
residual integral. The cold continuation is the explicit exception described
in :ref:`rotation-continuations`.

Notation, units and molecular states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Energies :math:`E,\epsilon_J,\Delta` are in eV, cross sections in m2,
DCS in m2/sr, rates :math:`K` in m3/s and temperatures in kelvin. In equations
below, :math:`P(E)=p(E)c=\sqrt{E(E+2m_ec^2)}` is momentum times light speed,
in eV; using :math:`p` instead leaves all momentum ratios unchanged.

.. math::

   v(E)=c\frac{P(E)}{E+m_ec^2},\qquad
   y=\sin(\theta/2),\quad d=1-\cos\theta=2y^2,\quad
   d\Omega=8\pi y\,dy.

The azimuth is uniform. Polar angles are those of an electron scattered from
a stationary molecule in its rest frame, corresponding to the laboratory
DCS measured with an initially stationary target.

A rigid rotor has :math:`H_{\rm rot}=\boldsymbol J^2/(2I)` and hence
:math:`B=\hbar^2/(2I)`. We use ground-relative levels

.. math::
   :label: rotation-levels

   \epsilon_J=B[J(J+1)-J_g(J_g+1)],\qquad
   \Delta_{if}=\epsilon_f-\epsilon_i.

Positive :math:`\Delta` is an electron energy loss. Negative :math:`\Delta`
is superelastic energy gain. The internal excitation threshold is
:math:`E=\Delta`, without the small molecular recoil shift. The target
isotope/statistical model and constants are:

.. list-table::
   :header-rows: 1
   :class: rotational-scattering-table
   :widths: 35 32 33

   * - Quantity
     - N2
     - O2
   * - Rotational constant, cm-1
     - 1.998
     - 1.438
   * - Rotational constant, eV
     - 0.0002477204284695341
     - 0.00017828927734694198
   * - Ground rotational level
     - :math:`J_g=0`
     - :math:`N_g=1`
   * - Statistical weight
     - :math:`g_J=(2J+1)(6\text{ if even, }3\text{ if odd})`
     - :math:`g_N=2N+1` for odd :math:`N`; zero for even :math:`N`
   * - Molecular mass
     - 28.0134 atomic mass units
     - 31.9988 atomic mass units
   * - Scattering separation :math:`R/a_0`
     - 2.068
     - 2.281
   * - Screening radius :math:`a/a_0`
     - 0.6052
     - 0.5677
   * - Bound-rotor energy ceiling
     - 9.759 eV
     - 5.116 eV

O2's :math:`N` is nuclear rotation, with electron-spin fine structure
unresolved. The implementation calls this index ``J`` in common routines;
it does not mean that the resolved total electronic-plus-rotational angular
momentum has been included. N2's 6:3 weights describe the predominant
14N2 nuclear-spin statistics. These weights are not a model of isotopic
abundances. The listed masses are the gas masses used by the implementation.
Scattering separation and :math:`B` are independent adopted source
parameters; recomputing :math:`B` from :math:`R` would change the model.

The equilibrium bath is

.. math::
   :label: rotation-populations

   b_J(T)=\frac{g_J\exp[-\epsilon_J/(k_BT)]}{Z(T)},\qquad
   Z(T)=\sum_K g_K\exp[-\epsilon_K/(k_BT)].

At zero temperature, all population is in :math:`J_g`. No neutral states
are advanced in time. The precomputation sums initial states in increments
of eight until a conservative upper bound on omitted population is below
:math:`10^{-10}`. With :math:`a=B/(k_BT)` and first omitted index :math:`j`,
the bound used for the unnormalized tail is

.. math::

   6(2j+1+a^{-1})\exp[-a(j(j+1)-J_g(J_g+1))].

It combines a decreasing first-term bound with the integral of
:math:`(2J+1)e^{-aJ(J+1)}`; it is used only on the decreasing tail.
Initial-state limits at 300 K are 56 for N2 and 64 for O2 (only allowed parity
states contribute). The energy-transfer moments are checked separately;
a small missing population alone does not bound a large energy-loss tail.
The final-state space is discussed in :ref:`rotation-continuations`.

Why a rate coefficient is tabulated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The frequency for a target density :math:`n_n` is :math:`\nu=n_nK(E)`, where

.. math::
   :label: rotation-rate

   K(E)=v(E)\int\sum_{if}b_iD_{if}(E,\Omega)\,d\Omega.

Detailed balance can give :math:`\sigma_{fi}(E)\sim E^{-1/2}` as
:math:`E\to0` while :math:`v\sim E^{1/2}`. Thus the *rate* is finite at
rest even when the cross section diverges. Storing a finite zero-energy
cross section would set that heating rate to zero. Storing an infinity
would make ordinary interpolation and :math:`0\times\infty` ambiguous.
The combined family stores the physically meaningful finite :math:`K`.
Ordinary channels continue to use their existing cross-section interface.

Kinetic extension of the thesis mean loss
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The thesis Eq. 2.48 averages :math:`\Delta_{if}` over the bath and the
angle-dependent rotational probabilities. Our discrete conditional law is

.. math::
   :label: rotation-conditional

   \Pr(i,f\mid E,\Omega)=\frac{b_iD_{if}(E,\Omega)}{D_{\rm inc}(E,\Omega)}.

Draw :math:`\Omega` from :math:`F`; then draw :math:`(i,f)` using this law.
The average of the drawn energy loss is the analogous first moment. Unlike
applying that mean deterministically, a discrete draw also retains
:math:`\langle\Delta^2\rangle`, the probability of heating, and correlation
with angle. For a fixed energy the rotational drift and diffusion moments are

.. math::

   A_E=-n_n v\sum_{if}b_i\int\Delta_{if}D_{if}\,d\Omega,\qquad
   B_E=n_n v\sum_{if}b_i\int\Delta_{if}^{\,2}D_{if}\,d\Omega.

These are moments of a jump process, not a replacement Fokker--Planck model.
Using the mean alone eliminates its physical energy diffusion. Collapsing
all :math:`i=f` outcomes into one zero-change outcome preserves this
rotational law. The existing elastic recoil uses the ground molecular mass
for that collapsed outcome; the negligible rotational mass increment is
retained for resolved changing outcomes only.

.. _rotation-reciprocity:

Reciprocity and reconciliation of the inclusive rate
----------------------------------------------------

The source mismatch
~~~~~~~~~~~~~~~~~~~

The IAA residual is a thermally unresolved vibrationally elastic rate; the
Kutz--Meyer elementary curves are calculations for a specified initial
rotational state. They also use different scattering approximations. At
2.3 eV the supplied N2 excitation curves alone sum to approximately
:math:`2.40\times10^{-19}` m2, exceeding the IAA residual of
:math:`1.73\times10^{-19}` m2. Boltzmann weighting is required, but it does
not resolve the mismatch: the earlier full sudden reconstruction at 300 K
was approximately :math:`5.09\times10^{-19}` m2. These numerical comparisons
refer to the unreconciled source reconstruction, not the production kernel.

A negative *angular residual* means that at some :math:`E,\theta`

.. math::

   D_{\rm inc}(E,\theta)-\sum_{i\ne f}b_iD_{if}(E,\theta)<0.

It is impossible to interpret this as a nonnegative unchanged probability.
Clipping that residual to zero changes the inclusive rate; renormalizing
incoming rows independently generally spoils detailed balance. The source
construction instead specifies which constraints are retained, and solves
for a positive common forward/reverse kernel under the inclusive constraint.

Detailed balance from time reversal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an unpolarized target in the heavy-target scattering approximation,
time-reversal invariance relates a transition and its reverse at equal total
energy. The incoming/outgoing state densities and incident flux give

.. math::
   :label: rotation-differential-balance

   g_iP(E)^2 D_{if}(E,\theta)
       =g_fP(E-\Delta)^2D_{fi}(E-\Delta,\theta),
   \qquad \Delta=\epsilon_f-\epsilon_i>0.

The reverse experiment reverses both asymptotic electron momenta and hence
has the same angle between incoming and outgoing directions. Integrating this
identity gives integral balance. The converse is false: integral balance
places no constraint on angular shapes. A model with
:math:`D_{if}=\sigma_{if}F(E,\theta)` and independently sampled rotational
outcomes needs :math:`F(E,\theta)=F(E-\Delta,\theta)` for differential balance,
which is not generally true of the elastic DCS. An angular dependence is not
by itself evidence for or against balance; the paired identity is the test.

For a Boltzmann electron gas, the isotropic density per energy is proportional
to :math:`P(E)(E+m_ec^2)e^{-E/(k_BT)}`. Multiplication by the speed leaves
:math:`P(E)^2 e^{-E/(k_BT)}`. Meanwhile
:math:`b_i/b_f=(g_i/g_f)e^{\Delta/(k_BT)}`. These factors and
:eq:`rotation-differential-balance` make the upward and downward thermal
fluxes equal after the change of variable :math:`E'=E-\Delta`.
The nonrelativistic limit gives the usual Maxwellian proof with
:math:`P^2\propto E`. Integral balance therefore suffices for zero net
rotational power for an isotropic thermal distribution. Differential balance
also constrains angular relaxation and angle--energy correlations away from
that equilibrium.

A shared forward primitive
~~~~~~~~~~~~~~~~~~~~~~~~~~

For each excitation pair define a positive quantity :math:`X_{if}` by

.. math::
   :label: rotation-primitive

   D_{if}(E,\theta)&=\frac{P(E-\Delta)}{P(E)^2}X_{if}(E,\theta),\\
   D_{fi}(E-\Delta,\theta)&=\frac{g_i}{g_f}
                                \frac{X_{if}(E,\theta)}{P(E-\Delta)}.

Both sides of :eq:`rotation-differential-balance` are then
:math:`g_iP(E-\Delta)X_{if}`. In the reference, :math:`X_{if}` is a positive
sum of rank primitives with the coefficients in :eq:`rotation-coupling`.
O2 also has explicit transition-specific Born primitives. The same stored
primitive and interpolation are used for both directions.

Multiplying by :math:`v(E)` makes the finite limits explicit. With
:math:`\kappa(E)=c/(E+m_ec^2)`, the upward and downward rate densities at
incident energy :math:`E` are

.. math::
   :label: rotation-paired-rates

   R^+_{if}(E,\theta)&=b_i\kappa(E)
            \frac{P(E-\Delta)}{P(E)}X_{if}(E,\theta),\\
   R^-_{fi}(E,\theta)&=b_f\kappa(E)\frac{g_i}{g_f}
                                      X_{if}(E+\Delta,\theta).

Excitation is zero below threshold. The reverse term remains finite at
:math:`E=0`. Changing a forward primitive at :math:`E` must also change its
reverse partner at :math:`E-\Delta`; multiplying all incoming outcomes at
one energy by an unrelated factor does not do this.

Inclusive constraint and the common correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At fixed temperature the sum of unchanged, forward and reverse rate densities
is constrained to :math:`v(E)\sigma_{\rm IAA}(E)F(E,\theta)` above the cold
join. A common multiplier :math:`c(E,\theta)` changes the forward rank
primitives. The reverse term automatically uses :math:`c(E+\Delta,\theta)`
through the shared primitive. This is an angle-dependent, shifted-energy
normalization problem, not a separate probability normalization for every
incoming rotational state.

The equation is solved from high to low energy because reverse contributions
refer to higher-energy forward kernels. When :math:`E+\Delta` falls between
the current and next grid row, its interpolation includes a contribution
from the unknown current row. That implicit contribution must be retained;
omitting it produces an energy-grid-dependent result.

The reference makes two different choices where data support them:

* N2 below 1 eV and O2 through 200 eV protect the changing primitives and
  solve for the rank-zero unchanged remainder. A materially negative
  remainder fails the construction; it is not clipped.
* Above those ranges a common multiplier adjusts all ranks while preserving
  their forward relative strengths at a given energy and angle. Transitions
  over 1--1.25 eV for N2 and 200--220 eV for O2 release the protected
  background smoothly.

The multiplier does not preserve all original *integral* channel rates or
all thermally averaged channel ratios: the reverse factors are evaluated at
shifted energies, and angular integration weights the corrected kernels.
It preserves the specified forward relative rank content before those sums.
This is the explicitly chosen reconciliation of incompatible sources, not a
consequence uniquely required by time reversal. Other positive reciprocal
kernels can share the same inclusive marginal.

Derivation of the implicit transition solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At one energy and angle let :math:`S` be the prescribed inclusive rate minus
known future-row contributions, :math:`D` the current unscaled forward rate
excluding rank zero (including nonzero-rank diagonal terms), :math:`Q` its
original rank-zero rate and :math:`H` the coefficient of the implicit reverse
contribution. All four have units of rate per solid angle. Let :math:`w` be
the smooth transition fraction. The protected rank-zero contribution before
the common multiplication is :math:`S-cH-D`, so the final equation is

.. math::

   c\{D+(1-w)(S-cH-D)+wQ+H\}=S.

Equivalently,

.. math::
   :label: rotation-quadratic

   (1-w)Hc^2-[(1-w)S+w(D+Q)+H]c+S=0.

The stable smaller root is

.. math::

   c=\frac{2S}{A+\sqrt{A^2-4(1-w)HS}},\qquad
   A=(1-w)S+w(D+Q)+H.

The discriminant is evaluated as a sum of nonnegative terms to avoid
cancellation. The reconstructed background uses :math:`S-cH-D`, not
:math:`S-H-D`. At the protected endpoint this recovers :math:`c=1`; at the
fully released endpoint it gives :math:`c=S/(D+Q+H)`. Independent bisection
checks the physical root. All resulting primitives and backgrounds must be
finite and nonnegative.

What is and is not exactly reciprocal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The paired primitive reference obeys differential reciprocity algebraically
in its heavy-target model. Numerical quadrature, finite grids, interpolation
of sampled distributions and the high-energy sudden approximation have
separate errors. Below 1 keV the runtime approximates this reference through
mixtures of prepared rows; it does not evaluate the paired identity exactly
at every sampled off-grid energy.

Above 1 keV the compact momentum-transfer bank uses the bounded sudden
spectrum. It neglects finite :math:`\Delta/E` corrections in its channel
weights and does not explicitly evaluate a reverse spectrum at
:math:`E+\Delta`. It is an asymptotic continuation, not a proof of exact
finite-energy differential balance. The molecular recoil update also uses
finite masses while the source reciprocity relation uses a heavy target.
Numerical equilibrium and transfer-moment checks quantify selected
consequences of these approximations; they do not turn these assumptions
into exact microscopic identities. Room-temperature equilibrium lies far
below the high-energy switch.
