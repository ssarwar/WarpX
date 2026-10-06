.. _rotation-continuations:

Cold limits, high-energy physics and supported ranges
=====================================================

Source ranges and the runtime contract
--------------------------------------

The ordinary residual elastic data end at 6 keV. Elementary N2 rotational
curves end at 1 keV; the O2 ground-state Born integral ends at 20 eV. These
are source endpoints, not justification for a constant tail. The production
combined family supports incident relative electron energies from zero to
1 GeV using the continuations specified here. An evaluation above 1 GeV
fails with a range error. Additional shifted-energy support needed to
construct reverse channels is generated during preparation and does not
extend the runtime domain.

For new ordinary tables, comment metadata declares ``energy_min_eV``,
``energy_max_eV`` and ``outside_energy_range = error``. Known zeros below
excitation thresholds are included in that supported domain. Background MCC
intersects all declared ranges during initialization and checks the result
before either its cached or fallback selection path. A relative tolerance of
eight ``ParticleReal`` machine epsilons at the upper endpoint accommodates
roundoff only; accepted roundoff is clamped to the endpoint. Negative energies
and genuinely out-of-range energies are not silently extended.
Unmarked legacy tables retain their established endpoint clamping.

No arbitrary distant zero is appended. Such a point would define a linear
taper where there was no source or model, followed by a permanent cutoff.
An asymptotic plateau that follows from an explicitly specified model is
also distinct from accidentally holding the last tabulated value.

Cold electrons and the unchanged background
-------------------------------------------

At and above 1 meV the reference uses the IAA inclusive rate and angular
marginal. Below that energy it retains the reciprocal changing kernels and
constructs a rank-zero s-wave background. Set :math:`E_c=0.001` eV and let
:math:`A_c(\theta)` be the positive unchanged background solved at that
energy. The continuation is

.. math::

   X_0(E,\theta)=P(E)\left[
      (1-s(\sqrt{E/E_c}))a_s^2a_0^2
      +s(\sqrt{E/E_c})A_c(\theta)\right],
   \qquad s(u)=u^2(3-2u),\quad 0\le u\le1.

The adopted dimensionless scattering-length magnitudes are 0.44 for N2 and
0.30 for O2. This is a specified cold background, not an additional fitted
rotational transition. Its DCS tends to :math:`a_s^2a_0^2`, giving a total
:math:`4\pi a_s^2a_0^2`, while its rate vanishes with speed. Superelastic
rates retain their finite limits from :eq:`rotation-paired-rates`. The
reference's first positive energy is continued with the analytic
:math:`P(E)/P(E_{\min})` factor rather than a linear extrapolation of a
nonzero primitive to zero. At exactly zero relative momentum there is no
incident axis, so the emitted electron's direction is isotropic. Angular
conditioning at that degenerate point is a convention, not a measured DCS.

Inclusive Born continuation
---------------------------

Above :math:`E_b=6000` eV the inclusive rate uses the elmolcs fitted Born
formula with a continuous scale match:

.. math::

   \sigma(E)=\sigma_{\rm residual}(E_b)
                 \frac{\sigma_B(E)}{\sigma_B(E_b)},\qquad
   \sigma_B=\pi a_0^2\frac{\alpha^2}{\beta^2}
          (A+B_1r+C_1r^2),
   \quad r=\frac{\alpha^2}{\beta^2\gamma^2},\quad
   \gamma=1+\frac{E+D}{m_ec^2},\quad\beta^2=1-\gamma^{-2}.

The elmolcs fit coefficients :math:`(A,B_1,C_1,D)` are
:math:`(129.243,-270.926,74854.382,342.492)` for N2 and
:math:`(148.738,-380.676,166556.802,412.32)` for O2, with :math:`D` in eV.
The fit is an elastic-scattering Born expansion, not an ionization Bethe
formula. Its relativistic limit is finite; that model limit must not be
confused with an uncontrolled endpoint hold. Continuity at 6 keV is imposed;
equality of derivatives to the residual interpolation is not assumed.

Angular continuation
--------------------

The native IAA DCS is interpolated linearly in :math:`y=\sin(\theta/2)` and
in :math:`\log E`, then normalized over solid angle. Source angular knots
and the integral rate are treated separately. With
:math:`k=P(E)/(\alpha m_ec^2)` in inverse Bohr radii, screening radius
:math:`a` in Bohr radii and :math:`\eta=(4a^2k^2)^{-1}`, the normalized
screened-Rutherford shape and its CDF are

.. math::
   :label: rotation-rutherford

   F_{\rm SR}(E,y)=\frac{\eta(1+\eta)}{4\pi(\eta+y^2)^2},\qquad
   U(y)=\frac{(1+\eta)y^2}{\eta+y^2},\qquad
   d(U)=\frac{2\eta U}{1-U+\eta}.

Between 8 and 10 keV, the angular density and CDF are mixed with fraction
:math:`s(\log(E/8000)/\log(10000/8000))`. Above 10 keV the analytic inverse
uses the actual incident energy. The continuation therefore keeps becoming
more forward-directed with increasing energy; it does not reuse the last
source angular row.

Bounded spectator spectrum
--------------------------

The two-centre spectator approximation transfers momentum to one atom during
a rapid collision. Expanding the orientation-dependent phase in spherical
harmonics gives elementary strengths proportional to

.. math::
   :label: rotation-spectator

   s_L(z)=(2L+1)j_L(z)^2,\qquad
   z=kR\sin(\theta/2),\qquad L=0,2,4,\ldots.

Here :math:`j_L` is a spherical Bessel function. The scale :math:`z`, called
the rotational rainbow scale in the thesis, is not a hard quantum upper
bound on :math:`\Delta J`. The production continuation recouples these
strengths to thermally occupied initial states. Define

.. math::

   w_{if}(z)=b_i\sum_L C_{ifL}s_L(z),\qquad
   P_{if}(z)=\frac{w_{if}(z)}{\sum_{i,f\in\mathcal B}w_{if}(z)}.

The set :math:`\mathcal B` includes only final levels with
:math:`\epsilon_f` below the declared dissociation ceiling. The largest final
levels are 197 for N2 and odd 167 for O2. This is conditioning on the bound,
ground-vibrational rotor manifold; it is not a high-J spectroscopic model
or a dissociation cross section. A rigid rotor without this restriction
would permit arbitrarily large rotational excitation. Missing dissociation
and non-rigid-rotor effects remain physical uncertainties.

Between 200 and 1000 eV the reference maps the 200 eV data to equal
:math:`z`: :math:`y_{200}=yP(E)/P(200)`. Outside :math:`y_{200}\le1`, a smooth
blend over :math:`1<y_{200}<1.05` returns to the spectator prior. The entire
angular-rank distribution then blends to the spectator spectrum with
:math:`s(\log(E/200)/\log(1000/200))`. The low-energy reference uses even
ranks through 48; its convergence is checked independently. The high-energy
bank uses all ranks required by the bounded final-state basis. These are
different finite representations, with a tested transition at 1 keV.

Momentum-transfer bins and oscillatory tails
--------------------------------------------

The high-energy bank tabulates probabilities in :math:`z`, not at one fixed
incident energy. The incident energy and drawn angle determine :math:`z`
at runtime. For the Rutherford law, with :math:`b=(R/(2a))^2`, the density
in :math:`z` is proportional to

.. math::

   w_z(z)=\frac{2bz}{(z^2+b)^2}.

Each bin stores :math:`\int_{z_l}^{z_h}w_zP_{if}\,dz/\int w_z\,dz`.
This averages probabilities, preserving the discrete energy spectrum and
its variance. It does not replace the outcome by a bin's mean loss.
The same bank is used between 1 and 10 keV under the blended elastic angular
law; independent angularly weighted checks bound that bin approximation.

The grid begins at zero, uses approximately 1% relative spacing from
:math:`10^{-5}` to :math:`10^4`, and 100 logarithmic tail intervals through
:math:`10^6`. The extreme tail uses a probability phase average of the
oscillatory spherical Bessel functions. For
:math:`V_{L,:}=\sqrt{2L+1}(j_L,-y_L)` and
:math:`M=V^T\operatorname{diag}(c_L)V`, where
:math:`c_L=\sum_{if}b_i C_{ifL}`, the identity

.. math::

   \left\langle\frac{uu^T}{u^TMu}\right\rangle_\phi
       =\frac{M^{-1/2}}{\operatorname{tr}(M^{1/2})},\qquad
   u=(\cos\phi,\sin\phi)

gives a positive phase-averaged rank spectrum. The slowly varying Hankel
envelope is held fixed during that phase average. Resolved quadrature checks
both ordinary and :math:`z^2`-weighted moments. Extending the last bin, if
needed by another supported energy range, requires an explicit tail bound;
the current runtime maximum is covered by the supplied bank.
