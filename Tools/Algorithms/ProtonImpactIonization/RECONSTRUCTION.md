# Proton PJG reconstruction: constraints and checked reference formulas

## Status

This is a research record, not a declaration that the production PJG spectrum
has been repaired. The production implementation still uses the provisional
total-only refit and its absolute free-electron cutoff. The independent
kinematics helpers and reference quadratures below are implemented and tested.
The bound-electron spectral response is not yet accepted for production.

The current direction is the smaller, PJG-preserving repair documented in
[PJG_REPAIR.md](PJG_REPAIR.md). It retains the Lorentzian structure and tests
coefficient constraints and a staged refit. The optical/GOS trial recorded
below remains a rejected exploratory calculation, not the preferred design.

The agreed source represents **inclusive electron production** using an
effective electron/singly charged molecular-ion pair. It does not represent
exclusive single ionization or resolve fragmentation and correlated multiple
electrons. The beam remains rigid, and the ion source is thermal. Consequently,
kinematic cross-section bounds do not imply eventwise energy or momentum
conservation by this deliberately non-conservative source operator.

"Printed PJG" means the printed proton formula with the 1977 corrections
applied. Its empirical continuum parameters are not constraints on a new
optical response, and a fit to its Figure 6 is not a substitute for a fit to
experimental data.

## Two different endpoints

Use rest energies in this section: $m=m_ec^2$, $M=M_pc^2$, and $A=M_nc^2$.
Let $E$ be the incident projectile kinetic energy. The stationary-free-electron
endpoint follows from the two-body invariant:

$$
T_{\max}^{\rm free}
=\frac{2mE(E+2M)}{(M+m)^2+2mE}
=\frac{2m\beta^2\gamma^2}{1+2\gamma m/M+(m/M)^2}.
$$

It is fixed, not fitted. Evaluating $\beta^2\gamma^2=u(u+2)$, $u=E/M$,
avoids cancellation for slow massive projectiles. The independent C++ tests
compare the production helper to the first, invariant form, in both particle
precisions.

For a stationary neutral and a specified residual-ion channel of ionization
energy $I$, the residual-ion rest energy is $A_i=A-m+I$. At the maximum
electron energy, the outgoing projectile and residual ion have their minimum
combined invariant rest energy $R=M+A_i$. Thus the three-body endpoint can
be evaluated as a two-body endpoint for an electron and this composite recoil.
Define

$$
s=(M+A)^2+2AE,\qquad
E_{\rm th}=I(1+M/A)+I^2/(2A).
$$

Below $E_{\rm th}$ the reaction is forbidden. Above threshold, use the
rationalized center-of-momentum excess energy

$$
\epsilon=\sqrt{s}-(M+A+I)
=\frac{2A(E-E_{\rm th})}{\sqrt{s}+M+A+I}.
$$

The electron kinetic energy and momentum in that frame are

$$
t_* = \frac{\epsilon(\epsilon+2R)}{2\sqrt{s}}
=\epsilon\left[1-\frac{m+\epsilon/2}{\sqrt{s}}\right],
\qquad cp_* = \sqrt{t_*(t_*+2m)}.
$$

Let $cp=\sqrt{E(E+2M)}$, and evaluate

$$
\gamma_{\rm cm}-1
=\frac{(cp)^2}{\sqrt{s}(E+M+A+\sqrt{s})}.
$$

A forward Lorentz boost then gives

$$
T_{\max}^{\rm mol}
=t_*+(\gamma_{\rm cm}-1)(t_*+m)+\frac{cp}{\sqrt{s}}cp_*.
$$

At the exact reaction threshold all products co-move, so the laboratory
electron kinetic energy need not vanish, although the available phase space
does. The helper returns zero below threshold. Additional residual excitation
must be included in $I$. The formula is evaluated in the neutral rest frame;
it does not itself account for a thermal-frame boost.

The independent test instead solves

$$
(E+M+A)T-cp\sqrt{T(T+2m)}
=E(A-m)-I(M+A-m)-I^2/2
$$

for the upper root and checks the available-energy bound $T\le E-I$.
The molecular endpoint can exceed the free endpoint by orders of magnitude.
This establishes allowed support, not its probability density. A bound tail
must not be created by evaluating the free Bhabha bracket outside its domain.

## Complete hard factor

For a pointlike spin-1/2 projectile, the complete tree-level per-electron result is

$$
\frac{d\sigma_{\rm B}}{dT}
=\frac{C_B}{E_eT^2}
\left[1-\beta^2\frac{T}{T_{\max}^{\rm free}}
      +\frac{T^2}{2(E+M)^2}\right],
\quad C_B=4\pi a_0^2R_y^2,\quad E_e=m\beta^2/2.
$$

For $x=T/T_{\max}^{\rm free}$ in $[0,1]$, evaluate its bracket as

$$
(1-x)+x/\gamma^2+\tfrac12[T/(E+M)]^2>0.
$$

This is checked against the invariant, spin-averaged Dirac trace rather than
another transcription of the same bracket. It is also consistent with the
close-collision expression in [Salvat and Heredia, Eq. (103)](
https://doi.org/10.1016/j.nimb.2023.165157). This spin-1/2 result must not be
identified as an exact spin-independent formula for every bare nucleus.

The $T^{-2}$ term belongs to this **whole** kernel. Appending only its negative
$T^{-1}$ and spin terms to a fitted Lorentzian does not establish this limit.
In particular, positivity must not rely on clipping a negative sum after
integration, or on shifting only one of the cancelling hard terms.

## Optical and finite-momentum response

Write $W=T+I_j$ for a simple resolved ionization channel. In the dipole
approximation its optical oscillator-strength density satisfies

$$
\frac{df_j}{dW}
=\frac{\sigma_{\gamma,j}(W)}{4\pi^2\alpha a_0^2R_y},\qquad
a(T)=\sum_j\frac{1}{T+I_j}\left.\frac{df_j}{dW}\right|_{W=T+I_j}.
$$

Photoabsorption is not automatically photoionization. A total photon cross
section also does not specify how photon energy is partitioned among primary
photoelectrons, shake-off, Auger electrons, and molecular fragments. Inclusive
electron production requires these distinctions in the optical input as well.
More generally, if $Y_e(T\mid W)$ is the conditional electron-number spectrum
per absorption, including branching and multiplicity, the required coefficient is

$$
a(T)=\int dW\,\frac{1}{W}\frac{df_{\rm abs}}{dW}Y_e(T\mid W).
$$

The simple-channel formula follows with $Y_e=\delta(T-W+I_j)$. Absorption
oscillator-strength sum rules apply to the underlying response, not blindly
to a multiplicity-weighted electron-production spectrum. Discrete absorption
followed by autoionization must likewise be distinguished from absorption
that produces no electron.

The relativistic Bethe limit has the structure

$$
\frac{d\sigma}{dT}
\sim\frac{Z^2 C_B}{E_e}
\{a(T)[\ln(E_e\gamma^2/R_y)-\beta^2]+b(T)\}.
$$

The same $a(T)$ multiplies both the logarithm and $-\beta^2$. It is not a
free collision normalization after the optical response is fixed. Optical
data determine $a$, but not the entire non-logarithmic function $b$. See
[Rudd et al. (1992), Eq. (22)](https://doi.org/10.1103/RevModPhys.64.441).

The reference quadrature uses $Q=\sqrt{(cq)^2+m^2}-m$ and
$G(Q,W)=df(Q,W)/dW$. In the common longitudinal/transverse GOS approximation,
[Salvat and Heredia, Eq. (68)](https://doi.org/10.1016/j.nimb.2023.165157) gives

$$
\frac{d^2\sigma}{dQ\,dW}
=\frac{Z^2C_B}{E_e}\left[
\frac{2m}{WQ(Q+2m)}+
\frac{2m\beta^2W\sin^2\theta_r}{[Q(Q+2m)-W^2]^2}
\right]G(Q,W).
$$

Both terms are nonnegative on the allowed momentum interval for nonnegative
$G$. With $p_f c=\sqrt{(E-W)(E-W+2M)}$, the endpoints are
$z_\pm=[c(p\pm p_f)]^2$, where $z=(cq)^2$. The independent numerical
implementation uses

$$
\sin^2\theta_r=\frac{(z-z_-)(z_+-z)}{4(cp)^2z}\ge0,
\qquad y=z-W^2>0.
$$

Integration in $\ln y$ resolves the narrow transverse contribution without
subtracting nearly equal relativistic momenta. The longitudinal integral is
checked analytically; the transverse integral is checked against independent
adaptive integration in $\ln Q$. Its dipole limit recovers
$\ln\gamma^2-\beta^2$, including the subtraction, not just the growing term.

At atomic energy losses the lower recoil is approximately
$Q_-\simeq W^2/(4E_e)$. Integrating a finite-$Q$ response gives the Coulomb
logarithm plus a constant depending on that response. Adding the transverse
term produces the relativistic Bethe form above. This is the origin of the
non-logarithmic matching problem: an arbitrary positive continuation of the
optical function does not uniquely determine the physical $b(T)$.

This kernel alone is not a complete target calculation. It assumes a common
longitudinal/transverse GOS and does not validate an arbitrary response model,
the low-velocity Born approximation, or a hard-scattering replacement. The
free hard kernel remains the independently checked Bhabha expression.

## Reference consistency and rejected trial

The implemented Rudd references use the 1985 recommended total and the 1992
SDCS, Eqs. (41)--(48), Tables I(a) and V. They are recommended model curves;
sampling them densely does not produce new independent experimental points.
The cited differential data cover 5--1700 keV for N2 and 7.5--300 keV for O2.
O2 differential predictions at 1 MeV must not be labeled measurements.

At 5 keV, numerical integration gives:

| Target | Rudd SDCS yield above free Tmax | Integrated SDCS / recommended total |
| --- | ---: | ---: |
| N2 | 0.14073635 | 0.97048504 |
| O2 | 0.22586266 | 1.11972889 |

The O2 5 keV SDCS is a short extrapolation below the cited differential-data
range. At 7.5 keV, where those data start, its integral/total ratio is 1.11118.
These model-evaluation checks are reproducible with `reference.py` and
`test_reference.py`. They are not direct measurements of a tail fraction.

An exploratory three-collision-parameter construction was tested: a complete
hard kernel with a positive bound-state turn-on, plus a positive finite-$Q$
optical response, and the PJG-like low-velocity distortion scale and exponent.
The finite-$Q$ response introduced one matching scale in place of the printed
Bethe constants, rather than a fitted optical amplitude. Even though this
form allows emission above the free endpoint, it was rejected: matching total
yields still left large spectral deficits around the binary-collision region,
including a factor of about eight for N2 near 219 eV at 100 keV incidence.
Its N2 tail fraction at 5 keV was about 0.29, versus 0.14 in the Rudd curve.

These trial numbers are diagnostics, not a new reference parameter set. The
optical inputs themselves were provisional: the N2 optical fit in Rudd's
Table II needs an explicit high-energy continuation, and the coarse O2
branching tables from [NCAR/GLOW](https://github.com/NCAR/GLOW) were used only
as local diagnostic inputs, not redistributed. That diagnostic did not include
a separate core optical channel. O2 discrepancies therefore cannot all be assigned to the
collision closure. No production coefficient is taken from this trial.

## Requirements before production replacement

The remaining work is a bound-electron response that broadens the binary
region, with optical normalization and channel accounting independently
checked. A bound-electron momentum distribution is a physically motivated
input to investigate; simply extending the free Bhabha polynomial above its
endpoint is not. Shell occupancies, binding energies, and available orbital
kinetic energies should be fixed reference inputs before adding collision-fit
parameters. A combined hard/soft response must also be checked for double
counting and the appropriate oscillator-strength sum rules.

The low-velocity distortion can be refitted only after this response is
defined. The empirical PJG delta is not reintroduced into free kinematics.
Total yields, SDCS shapes, tail fractions, energy moments, positivity,
thresholds, and the independent hard and Bethe limits all remain acceptance
criteria. A good total alone is insufficient.

The angular closure must be revised together with the bound tail. The current
free binary polar cosine is outside its physical domain above the free
endpoint. Clamping that cosine is not a bound-electron angular model. The
current bound/free closure also produces an artificial forward point mass
where its unbounded cosine exceeds one; the new closure must use a normalized
distribution with the proper support instead. Neither issue is repaired by
changing an energy-only CDF.

The independent source-weight/energy correlation has been corrected with a
uniformly shifted base-2 radical inverse. A two-parent analytic spectrum
checks conditional uniformity and parent-order independence without assuming
any PJG formula. Thermal ions now use the neutral velocity variance, with
neutral mass inferred from ion mass plus electron mass.

This extra response calculation belongs in host-side table generation. The
existing constant-time inverse-CDF device interface can be retained, so a
more detailed molecular calculation does not require per-event GPU
quadrature. Lookup error, joint sampling bias, memory bandwidth, and
CUDA/HIP/SYCL execution still require separate tests; no accelerator
performance result is implied by the host reference checks.

## Verification on the development host

The CPU library and Python extension built successfully after incorporating
the branch's upstream merges. Two standalone C++ tests and nine Python
reference tests passed. The C++ tests cover both float and double; the Python
checks also test the transverse Bethe subtraction and independent quadrature.
The Python files pass the repository's Ruff checks.

The isolated host quantile microbenchmark, using a release build and the
median of seven runs, measured approximately 0.99 ns per shifted float
quantile versus 0.94 ns for the incorrect ordered baseline, and 0.81 ns versus
0.79 ns in double. These are local host measurements, not accelerator or
end-to-end source speedups.

The rebuilt physics and bounded-creation performance simulations produced
outputs that passed both analysis scripts. Both run stages nevertheless hung
after the final simulation/profile output and reached the 45-second CTest
timeout. The overall end-to-end CTest result is therefore **not a pass**.
CUDA/HIP/SYCL execution and performance have not been measured on this host.
