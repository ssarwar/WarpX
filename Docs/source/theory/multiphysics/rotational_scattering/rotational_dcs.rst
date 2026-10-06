.. _rotation-dcs:

Rotational differential cross sections for N2 and O2
====================================================

This page specifies the angular kernels used by ``rotation_model =
reciprocal_hybrid``. It distinguishes the imported elastic angular marginal,
the rotational source kernels, the reciprocal normalization, and the
conditional channel probabilities sampled by MCC. The source bibliography
and data qualifications are in :ref:`rotation-sources`; the normalization
solve is derived in :ref:`rotation-reciprocity`.

A rank :math:`L` below denotes an elementary angular-momentum-transfer tensor.
It is not an initial molecular level, and is not generally the final molecular
level. Only for initial :math:`J=0` does its label coincide with the elementary
final level. Recoupling each tensor gives transitions between arbitrary initial
and final rotor levels, including unchanged :math:`i=f` contributions from
nonzero :math:`L`. A sampled channel is a pair :math:`i\to f`, with its
own discrete energy gap, rather than a mean loss assigned to a rank.

Definitions and the elastic marginal
------------------------------------

Use :math:`E` in eV, :math:`y=\sin(\theta/2)`, :math:`\mu=1-2y^2`, and
:math:`P(E)=\sqrt{E(E+2m_ec^2)}` in eV. Define

.. math::
   :label: rotation-dcs-definitions

   I(E,\theta)&=\sigma_{\rm IAA}(E)F(E,\theta),
       \qquad \int F(E,\theta)d\Omega=1,\\
   C_{ifL}&=|\langle i0,L0\mid f0\rangle|^2,
       \qquad \Delta_{if}=B[f(f+1)-i(i+1)].

The source reference uses :math:`m_ec^2=510998.95069` eV,
:math:`E_h=27.211386245988` eV (one Hartree),
:math:`a_0=5.29177210544\times10^{-11}` m (the Bohr radius), and
:math:`\alpha=7.2973525693\times10^{-3}` (the fine-structure constant).
The Boltzmann factor uses :math:`k_B=8.617333262145\times10^{-5}` eV/K.

Here :math:`I` is the prescribed thermally inclusive, vibrationally elastic
DCS, in m2/sr. Its total comes from the residual integral and its normalized
shape from the elastic DCS; their original source normalizations need not
agree. Neither :math:`I` nor :math:`F` is the rotationally unchanged DCS.
The canonical constants, statistical weights and populations :math:`b_i`
are specified in :ref:`rotation-physics`.

An angular function written in :math:`\theta` is evaluated at the corresponding
:math:`y=\sin(\theta/2)` when passed to a code routine that uses :math:`y`.
Degrees appear only in the measured-angle interpolation and its explicit
boundary fractions.

All joins below use the clipped cubic smoothstep

.. math::
   :label: rotation-dcs-smoothstep

   s(x)=u^2(3-2u),\qquad u=\min(1,\max(0,x)).

Its derivative vanishes at either end. When the argument contains
:math:`\log E`, the blend is smooth in logarithmic energy. Source interpolation
within a regime can still have derivative changes at individual tabular knots;
a smooth join does not promise global differentiability of the source data.

The normalized spectator *rank prior* used while preparing the low/intermediate
energy kernels is

.. math::
   :label: rotation-dcs-rank-prior

   S_L(E,y)=\frac{(2L+1)j_L(z)^2}
                 {\sum_{K\in\mathcal L}(2K+1)j_K(z)^2},\qquad
   z=\frac{R P(E)}{\alpha m_ec^2}\,y,
   \qquad \mathcal L=\{0,2,\ldots,48\}.

In the sudden, two-centre picture the orientation-dependent amplitude of
identical scatterers contains the phases
:math:`e^{i\boldsymbol q\cdot\boldsymbol R/2}+e^{-i\boldsymbol q\cdot\boldsymbol R/2}`.
The Rayleigh expansion gives

.. math::

   2\cos\!\left(\frac{qR}{2}\cos\gamma\right)
     =2\sum_{L\ {\rm even}}(2L+1)i^L j_L(qR/2)P_L(\cos\gamma),

where :math:`\gamma` is the molecular-axis angle relative to momentum transfer.
Projection onto unpolarized rotor states yields elementary weights proportional
to :math:`(2L+1)j_L(qR/2)^2`. With the sudden elastic approximation
:math:`q=2k\sin(\theta/2)`, the Bessel argument is
:math:`qR/2=kR\sin(\theta/2)`, not :math:`qR`. A common atomic scattering
amplitude cancels in the conditional probabilities; the IAA data supply the
absolute inclusive rate and angular marginal. This derivation explains the
angular dependence but does not incorporate finite rotational-gap corrections
or establish its accuracy in a molecular resonance.

:math:`R` is 2.068 Bohr for N2 and 2.281 Bohr for O2. :math:`j_L` is a
spherical Bessel function. The high-energy collision spectrum uses a larger
bound-rotor basis and a different, thermally recoupled normalization, specified
below. Confusing these two normalizations changes the relative channels.

N2: complete construction by energy range
-----------------------------------------

The following :math:`A_L(E,\theta)` are reduced angular strengths before
reciprocal reconciliation. Multiplication by :math:`P(E)` gives source
primitives. They are not yet the final normalized state-to-state DCS.

.. list-table:: N2 source-kernel regimes
   :header-rows: 1
   :class: rotational-scattering-table
   :widths: 18 42 40

   * - Electron energy
     - Kernel used
     - Evidence and applicability
   * - Up to 1 eV
     - Isotropic elementary rotational strengths, canonical threshold opening
     - Real elmolcs/Itikawa and Kutz--Meyer integrals; leading low-energy angular approximation
   * - 1--1.25 eV
     - Cubic blend into the resonance construction
     - Chosen join; not a separately measured DCS
   * - 1.25--4 eV
     - Kutz--Meyer energy strengths with Jung-constrained, Read-completed angular multipliers
     - Anchors at 2.22 and 2.47 eV; direct-background and missing-angle assumptions are explicit
   * - 4--10 eV
     - Log-energy blend into the Gote construction
     - Chosen connection to the 10 eV angular fractions
   * - 10--200 eV
     - Gote fractions under the IAA elastic marginal
     - Measured angular range 10--160 degrees; spectator completion outside it
   * - 200--1000 eV
     - Equal-momentum-transfer continuation of the 200 eV shape, blended to the spectator prior
     - High-energy model continuation beyond the measured fractions
   * - At least 1000 eV
     - Bounded, thermally recoupled spectator outcome spectrum
     - Sudden high-energy approximation; finite-gap reciprocity is approximate

Elementary strengths and the sub-eV kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For :math:`L=2,4,6`, reduce the elementary integral by its threshold factor:

.. math::

   U_L(E)=\sigma_{0\to L}(E)
             \frac{P(E)}{P(E-\Delta_{0L})},\qquad U_0(E)=\sigma_{0\to0}(E).

Positive knots of :math:`U_L` use log--log interpolation. Between its canonical
threshold and first positive datum the first reduced value is held. The
primitive can also be evaluated below that elementary threshold when
recoupling a different initial state; the *actual* :math:`i\to f` channel
is opened only by its own :math:`\Delta_{if}`. This continuation is a
threshold-corrected sudden model, not a new independent measurement.

The low-energy prior is

.. math::

   A^{\rm low}_0=I(E,\theta),\qquad
   A^{\rm low}_L=U_L(E)/(4\pi)\quad(L=2,4,6),

with other ranks zero. The isotropic nonzero-rank strengths do not make the
total elastic DCS isotropic: rank zero is subsequently replaced by the
positive background required by the inclusive constraint.

Resonance shapes, measured fractions and missing angles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The normalized Read--Andrick tensors used to complete the angular coverage are

.. math::

   h_0(\mu)&=\frac{5(3\mu^2-1)^2}{16\pi},\\
   h_2(\mu)&=\frac{5(9\mu^4-9\mu^2+4)}{56\pi},\\
   h_4(\mu)&=\frac{5(\mu^2+3)^2}{224\pi},
   \qquad \int h_Ld\Omega=1.

The N2 temporary-anion shape resonance is treated as a predominantly
:math:`{}^2\Pi_g` electronic configuration, with capture/emission dominated
by :math:`\ell=2`, :math:`|\Lambda|=1` partial waves. Coupling two d waves
allows elementary even tensor ranks 0, 2 and 4. Projection onto rotor states,
averaging initial magnetic substates and summing unobserved final substates
gives the polynomials above. The single-resonance construction assumes a
lifetime short compared with the rotational period and separable electronic
and vibrational motion. Rank six retained below is an additional source
contribution, not a prediction of this pure d-wave resonance algebra.

They describe the resonance contribution near the N2 shape resonance; they
are not an accurate universal rotational law from threshold to high energy.
The actual construction also includes a direct/background component and uses
Jung's vibrationally elastic data, rather than identifying the whole elastic
DCS with a pure Read tensor.

At each anchor :math:`E_a\in\{2.22,2.47\}` eV, the processed source file
``n2_jung.json`` supplies rank fractions at 15, 30, ..., 105 degrees. They are
inferred from thermal branches at 500 K using the recoupling coefficients,
not raw state-to-state measurements. ``prepare_jung.py`` reproduces that
inference from ``n2_jung_digitization.json``. Interpolate each fraction
linearly in :math:`\theta`. Outside the measured interval continue its edge
value by :math:`h_L(\theta)/h_L(\theta_{\rm edge})`, and normalize the three
fractions to sum to one.

Call those continued fractions :math:`r_L^a`. Define

.. math::

   m(\theta)=s\left(\frac{15^\circ-\theta}{15^\circ}\right)
             +s\left(\frac{\theta-105^\circ}{15^\circ}\right),
   \qquad (a_2,a_4)=(0.52823198,0.84013224),\\
   \widetilde r_L^a=r_L^a[1-(1-a_L)m]\quad(L=2,4),\qquad
   \widetilde r_0^a=r_0^a+\sum_{L=2,4}(r_L^a-\widetilde r_L^a).

These fitted missing-angle factors reduce changing intensity outside the
measured range and put it into rank zero. The same factors are used at both
anchors. They are adopted fit parameters, not universal constants from
Read's theory. Their calibration to Jung's integral branch sums and the
remaining reading/model uncertainty are discussed in :ref:`rotation-sources`.

The angular multipliers and energy-dependent resonance prior are

.. math::

   H_L^a(\theta)&=\frac{I(E_a,\theta)\widetilde r_L^a(\theta)}{U_L(E_a)},\\
   t_R(E)&=\min\left(1,\max\left(0,\frac{E-2.22}{0.25}\right)\right),\\
   H_L(E,\theta)&=(1-t_R)H_L^{2.22}+t_R H_L^{2.47},\\
   A_L^{\rm res}(E,\theta)&=U_L(E)H_L(E,\theta)\quad(L=0,2,4),\\
   A_6^{\rm res}&=U_6(E)/(4\pi).

All other ranks vanish in this resonance prior. The interpolation of
:math:`H_L` between anchors is linear in energy and clamped outside the
anchors; the elementary strength still varies with energy. For
:math:`E\le4` eV,

.. math::

   A_L=(1-w_N)A_L^{\rm low}+w_N A_L^{\rm res},\qquad
   w_N=s((E-1)/0.25).

Thus the resonance multipliers are used in a wider *modeled* interval than
the two measured energies. They should not be described as measured
angle-resolved data across the entire 1.25--4 eV interval.

Gote fractions and the 4--10 eV join
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`g_L(E,\theta)` be the completed Gote percentages divided by 100
and normalized over the reported ranks. The completion retains published
values, represents a ``<1%`` entry by its interval midpoint, and replaces
internally inconsistent angular columns by interpolation between consistent
neighbors. Those choices and the excluded absolute-DCS row are documented
in :ref:`rotation-sources`.

Between tabulated energies :math:`E_a,E_b`, fractions interpolate as

.. math::

   g_L(E,\theta)=(1-t_G)g_L(E_a,\theta)+t_Gg_L(E_b,\theta),\qquad
   t_G=\frac{\log(E/E_a)}{\log(E_b/E_a)}.

Angle interpolation is linear in degrees between 10 and 160 degrees. Below
10 eV the measured-angle fractions are held at their 10 eV values. Outside
the measured angular range, define
:math:`w_0=s(\theta/10^\circ)` and
:math:`w_\pi=s((\theta-160^\circ)/20^\circ)`. The completed shape is

.. math::

   G_L(E,\theta)=
   \begin{cases}
     (1-w_0)S_L(E,\theta)+w_0g_L(E,10^\circ), & 0\le\theta<10^\circ,\\
     g_L(E,\theta), &10^\circ\le\theta\le160^\circ,\\
     (1-w_\pi)g_L(E,160^\circ)+w_\pi S_L(E,\theta), &160^\circ<\theta\le180^\circ.
   \end{cases}

The finite-precision implementation normalizes this result over ranks. The
unmeasured wings use :math:`S_L` at the actual incident energy, including
when extending the 10 eV fractions downward for the join. For 4--10 eV,

.. math::

   A_L=(1-w_G)A_L^{\rm res}+w_G I(E,\theta)G_L(E,\theta),\qquad
   w_G=s\left(\frac{\log(E/4)}{\log(10/4)}\right).

For 10--200 eV, :math:`A_L=I G_L`. The elementary Kutz--Meyer strengths are
therefore used in the low/resonant construction and its join, not imposed
as additional absolute rotational rates throughout 10--200 eV.

O2: complete construction by energy range
-----------------------------------------

O2 nuclear rotation is labeled :math:`N`; the common implementation index
``J`` means this :math:`N`, with electron-spin structure unresolved. Only odd
initial/final levels occur.

.. list-table:: O2 source-kernel regimes
   :header-rows: 1
   :class: rotational-scattering-table
   :widths: 18 42 40

   * - Electron energy
     - Kernel used
     - Evidence and applicability
   * - Up to 1 eV
     - Transition-specific quadrupole/polarization Born kernel
     - Long-range, sub-eV approximation; actual momentum transfer for each transition
   * - 1--20 eV
     - Log-energy blend of Born terms and positive angular completion
     - Modeled bridge; 20 eV integral constraints held while the IAA angular marginal varies with energy
   * - 20--200 eV
     - Positive completion constrained by Bhattacharyya integral and momentum-transfer data
     - Theoretical source, not a complete measured rotational DCS set
   * - 200--1000 eV
     - Equal-momentum-transfer continuation of the 200 eV shape, blended to the spectator prior
     - Modeled continuation
   * - At least 1000 eV
     - Bounded, thermally recoupled spectator outcome spectrum
     - Sudden high-energy approximation; finite-gap reciprocity is approximate

Transition-specific Born kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each allowed :math:`N\to N+2`, let :math:`\Delta_N=B(4N+6)` and, in
atomic units,

.. math::

   k_i=\sqrt{2E/E_h},\quad
   k_f=\sqrt{2\max(E-\Delta_N,0)/E_h},\quad
   q_N=\sqrt{k_i^2+k_f^2-2k_i k_f\mu}.

The reduced differential strength is

.. math::
   :label: rotation-dcs-o2-reduced

   B_N(E,\mu)&=\frac65\frac{(N+1)(N+2)}{(2N+1)(2N+3)}
           \left(\frac{Q}{3}+\frac{\pi\alpha_2q_N}{32}\right)^2a_0^2,\\
   Q&=-0.29,\qquad\alpha_2=4.93.

The origin of the two terms is the anisotropic long-range potential, in
atomic units and the thesis's multipole convention,

.. math::

   V_2(r)=-Q/r^3-\alpha_2/(2r^4).

The first-Born radial matrix element contains
:math:`\int_0^\infty j_2(qr)V_2(r)r^2dr`. Using

.. math::

   \int_0^\infty\frac{j_2(qr)}r dr=\frac13,\qquad
   \int_0^\infty\frac{j_2(qr)}{r^2}dr=\frac{\pi q}{16}

gives the amplitude :math:`-[Q/3+\pi\alpha_2q/32]`. Its square gives
:eq:`rotation-dcs-o2-reduced`; the two potential contributions are added
coherently before squaring, not treated as independent probabilities.
The unpolarized quadrupole angular factor is :math:`(4/5)C_{N,N+2,2}`, equal
to the prefactor :math:`(6/5)(N+1)(N+2)/[(2N+1)(2N+3)]` above. The
outgoing-to-incoming momentum ratio is the channel phase-space/flux factor.
Extending the long-range potential to the full radial integral neglects
short-range and exchange/resonance physics; that is why this formula is not
used as a universal O2 rotational DCS.

The Born DCS is :math:`(k_f/k_i)B_N`; the shared production primitive
supplies the corresponding threshold factor using :math:`P(E)`.
The nonrelativistic momentum in the bracket and the common relativistic
phase convention differ negligibly in the sub-eV region. This evaluates
:math:`q_N` separately for every transition. Rescaling the ground
:math:`1\to3` angular curve to every :math:`N` would not reproduce it.
The source formula is thesis Eq. 11.21b, the quadrupole and anisotropic induced
polarization contribution. There is no permanent-dipole term for O2.

Positive completion of the intermediate-energy constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Denote the paper's three integral cross sections by
:math:`\sigma_{11}^{p},\sigma_{13}^{p},\sigma_{\rm inc}^{p}` and the
corresponding momentum-transfer integrals by
:math:`m_{11}^{p},m_{13}^{p},m_{\rm inc}^{p}`. Values interpolate log--log
in energy. For this completion define :math:`\widehat E=\min(200,\max(20,E))`:
source values and the spectator prior use :math:`\widehat E`, whereas the
normalized IAA angular marginal and its absolute total use the actual
:math:`E`. This distinction matters in the 1--20 eV bridge.

Calculate prior angular moments

.. math::

   S_L^0=\int F(E,\theta)S_L(\widehat E,\theta)d\Omega,\qquad
   M_L^0=\int(1-\mu)F(E,\theta)S_L(\widehat E,\theta)d\Omega.

Let :math:`f_0=f_2=0`, :math:`f_4=5/9`, and :math:`f_L=1` for even
:math:`L\ge6`. The unresolved source residuals are

.. math::

   R_\sigma=\sigma_{\rm inc}^{p}-\sigma_{11}^{p}-\sigma_{13}^{p},\qquad
   R_m=m_{\rm inc}^{p}-m_{11}^{p}-m_{13}^{p}.

For :math:`L\ge4` allocate

.. math::

   a_L=S_L^0\frac{R_\sigma}{\sum_Kf_KS_K^0},\qquad
   m_L=M_L^0\frac{R_m}{\sum_Kf_KM_K^0},

and obtain rank two from recoupling:

.. math::

   a_2=\frac{\sigma_{13}^{p}-(4/9)a_4}{3/5},\qquad
   m_2=\frac{m_{13}^{p}-(4/9)m_4}{3/5}.

These equations follow from
:math:`\sigma_{13}=(3/5)a_2+(4/9)a_4` and the residual identity in
:eq:`rotation-o2-moments`. They provide target rank integrals and first
angular moments, not unique angular functions. The absolute unchanged
:math:`\sigma_{11}^{p}` is not separately imposed on the final
IAA-normalized thermal total.

Use relative-entropy projection to complete the angular functions. Set
:math:`r_2=1`, :math:`r_L=S_L(\widehat E,\theta)` for higher retained ranks,
and a rank-zero prior of one. With
:math:`Z_\theta=1+\sum_{L>0}r_L e^{\lambda_L+\eta_L\mu}`,

.. math::

   q_L=\frac{r_L e^{\lambda_L+\eta_L\mu}}{Z_\theta},\qquad
   q_0=Z_\theta^{-1},\qquad H_L(E,\theta)=I(E,\theta)q_L(E,\theta).

The multipliers minimize

.. math::

   \Phi(\lambda,\eta)=\int F\log Z_\theta\,d\Omega
      -\sum_{L>0}\lambda_L\frac{a_L}{\sigma_{\rm IAA}}
      -\sum_{L>0}\eta_L\frac{a_L-m_L}{\sigma_{\rm IAA}}.

Its gradient enforces the desired :math:`\int Fq_L` and
:math:`\int\mu Fq_L`. The code checks achieved constraints after the solve
and retains only ranks with :math:`a_L>10^{-13}\sigma_{\rm IAA}`.
This positive completion is a stated inference from incomplete theoretical
data. It is not a new measurement or a unique quantum-scattering solution.

The 1--20 eV bridge and unchanged background
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`w_O=s(\log(\max(E,1))/\log20)` for energies in eV, the raw
changing strengths below 20 eV are

.. math::
   :label: rotation-dcs-o2-bridge

   A_{if}^{\rm raw}(E,\theta)=
      w_O\sum_{L>0}C_{ifL}H_L(E,\theta)
      +(1-w_O)\mathbf1_{f=i+2}B_i(E,\theta),\qquad f>i.

Below 1 eV, :math:`w_O=0`: only the transition-specific Born kernel
contributes. At 20 eV, :math:`w_O=1`: the intermediate-energy completion is
fully active. At intermediate energies, :math:`H_L(E)` is recomputed under
the actual-energy elastic marginal with the 20 eV source constraints. It is
not simply the fixed 20 eV angular curve multiplied by a scalar blend.

For 20--200 eV the raw rank strengths are :math:`A_L=H_L`.
The protected-background solve determines the rank-zero unchanged part for
all energies through 200 eV. It also retains the diagonal contributions of
nonzero ranks. Source code's placeholder :math:`A_0=I` below 20 eV must not
be read as an extra unchanged rate to add afterward. The 200--220 eV release
of that background constraint is a normalization transition, distinct from
the 200--1000 eV angular-shape transition.

Common 200--1000 eV angular continuation
----------------------------------------

For N2, let :math:`Q_L^{200}(y)=G_L(200,y)`. For O2, use the completed
200 eV strengths divided by :math:`I(200,y)`, setting the rank-zero remainder
so that :math:`\sum_L Q_L^{200}=1`. Define

.. math::

   r=y\frac{P(E)}{P(200)},\quad y_a=\min(r,1),\quad
   v=s((r-1)/0.05),\quad w=s\left(\frac{\log(E/200)}{\log5}\right),\\
   Q_L^{\rm map}(E,y)=(1-v)Q_L^{200}(y_a)+vS_L(E,y),\\
   A_L(E,y)=I(E,y)[(1-w)Q_L^{\rm map}(E,y)+wS_L(E,y)].

This preserves equal momentum transfer where the anchor covers it. Beyond
that coverage, the additional :math:`1<r<1.05` transition returns to the
spectator prior. N2 evaluates its completed Gote function at the mapped
angle; O2 interpolates its stored anchor linearly in :math:`y`. Neither
uses the 200 eV distribution at the same angle regardless of energy.

From source kernels to final rotational DCS
-------------------------------------------

For :math:`E<1000` eV the raw source primitives are :math:`P(E)A_L` and,
for O2, the transition-specific Born components. The coupled inclusive solve
modifies them into positive primitives :math:`X_L` and :math:`X_i^{B}`.
For an excitation pair,

.. math::

   X_{if}(E,\theta)=\sum_L C_{ifL}X_L(E,\theta)
                      +\mathbf1_{f=i+2}X_i^{B}(E,\theta).

The explicit Born term is absent for N2 and where its O2 bridge weight has
vanished. The final forward, reverse and unchanged DCS are

.. math::
   :label: rotation-dcs-final

   D_{if}(E,\theta)&=\frac{P(E-\Delta_{if})}{P(E)^2}X_{if}(E,\theta),
      \qquad f>i,\ E>\Delta_{if},\\
   D_{fi}(E,\theta)&=\frac{g_i}{g_f}\frac{X_{if}(E+\Delta_{if},\theta)}{P(E)},
      \qquad f>i,\\
   D_{ii}(E,\theta)&=\frac{\sum_L C_{iiL}X_L(E,\theta)}{P(E)}.

Excitation is zero at and below its canonical threshold. The zero-energy
reverse *rate* uses its finite limit rather than evaluating a divergent DCS.
The corrected forward primitive at the shifted energy is also used for its
reverse; no separate post-normalization of incoming state rows is performed.
The precise protected/background and common-factor equations are
:eq:`rotation-quadratic` and :eq:`rotation-primitive`.

Above the 1 meV cold join this construction imposes

.. math::

   \sum_{if}b_iD_{if}(E,\theta)=I(E,\theta).

Consequently the scattering angle is sampled from the inclusive elastic DCS,
and the rotational outcome, conditional on that angle, has probability

.. math::
   :label: rotation-dcs-probability

   \Pr(i,f\mid E,\theta)=\frac{b_iD_{if}(E,\theta)}{I(E,\theta)}.

Summing :math:`i=f` gives the unchanged probability. Summing forward and
reverse pairs gives excitation and de-excitation probabilities. This retains
angle--energy correlations through the DCS; it imposes no additional hard
angular accessibility cutoff. The finite sampling arrays approximate this
joint law through the row mixtures in :ref:`rotation-algorithm`.
Below 1 meV use the cold inclusive sum from :ref:`rotation-continuations`
in the denominator; it need not equal the extrapolated residual IAA DCS.

High-energy rotational DCS and the separate elastic-angle switch
----------------------------------------------------------------

For :math:`E\ge1000` eV, let
:math:`s_L(z)=(2L+1)j_L(z)^2` and retain all ranks needed by the bound final
rotor manifold. Define

.. math::

   H(z)=\sum_i b_i\sum_{f\in\mathcal B}\sum_L C_{ifL}s_L(z),\qquad
   D_{if}^{\rm sudden}(E,\theta)=I(E,\theta)
                   \frac{\sum_L C_{ifL}s_L(z)}{H(z)}.

Thus the thermal conditional probability remains :math:`b_iD_{if}/I` and
the inclusive elastic marginal is preserved. The runtime uses positive
momentum-transfer-bin averages of these probabilities and phase averaging in
the extreme oscillatory tail, as specified in :ref:`rotation-continuations`.
It samples discrete losses and gains, not their mean. Finite-gap phase-space
and shifted-energy corrections are neglected here, so this is an asymptotic
spectator law rather than exact finite-energy differential reciprocity.
All bound rotational excitation gaps are below 1 keV; no closed bound channel
is opened by that approximation in its runtime interval.

The elastic angle itself follows its independent continuation: the IAA shape
through 8 keV, a log-energy smoothstep to screened Rutherford over 8--10 keV,
and the analytic, energy-dependent Rutherford inverse above 10 keV. The
inclusive *integral rate* has its own matched Born join at 6 keV. Those
6, 8 and 10 keV boundaries are not additional changes of rotational state
model. The supported runtime range ends at 1 GeV and is enforced explicitly.
