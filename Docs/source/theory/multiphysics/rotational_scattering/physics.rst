.. _rotation-physics:

Physical variables and the collision operator
=============================================

Inclusive elastic scattering
----------------------------

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
------------------------------------

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
-----------------------------------

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
-----------------------------------------

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
