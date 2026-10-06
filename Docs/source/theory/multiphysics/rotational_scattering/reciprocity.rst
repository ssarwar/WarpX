.. _rotation-reciprocity:

Reciprocity and reconciliation of the inclusive rate
====================================================

The source mismatch
-------------------

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
-----------------------------------

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
--------------------------

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
----------------------------------------------

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
-------------------------------------------

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
-------------------------------------

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
