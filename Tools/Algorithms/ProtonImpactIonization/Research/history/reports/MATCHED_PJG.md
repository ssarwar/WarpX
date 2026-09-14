> Historical research note, recovered 2026-09-14.
> Model names such as 'current' refer to this note's stage, not today's production model.
> Original commands and links retain their historical context; see the
> [archive guide](../../README.md) for status and portable reproduction.

# Matched PJG proton model: earlier formulation and refit

The selected energy-only calibration is now in [FINAL_PJG.md](FINAL_PJG.md).
This document preserves the earlier coefficients and diagnostics for comparison;
they are not the final calibration.

This is the selected research formulation for N2 and O2, not a replacement of
the production tables yet. No production coefficient, collision operator,
input option, or documentation under `Docs/` is changed by this study.
The earlier experiments in `PJG_REPAIR.md` remain historical diagnostics.

The later [source audit](DATASET_AUDIT.md) reassesses this candidate against
additional original and evaluated data. Its provisional trials do not
supersede the coefficients here as an accepted production model.

The target observable is **inclusive electron-production yield**, represented
by an explicitly effective one-electron/one-ion pair source. It is not an
exclusive molecular single-ionization channel cross section. All secondary
energies, including the low-energy peak and the bound-electron tail, belong
to the kinetic source.

## 1. What is retained and what is changed

Start with Porter, Jackman, and Green (1976), including the 1977 erratum.
In particular, the apparent superscript plus in Eq. (16) is ordinary addition,
not a positive-part operation. The baseline comparisons below do not clip
negative printed-model values or silently substitute the corrected cutoff.

Retain PJG's continuum weights and thresholds, a narrow Lorentzian minus a
broad Lorentzian, its projectile-dependent peak position, a Bethe logarithm,
and a low-velocity distortion. Make the following explicit model changes,
which are not claimed to be additional typographical corrections:

1. Put the whole, correctly normalized free Bhabha kernel in the hard term.
2. Use a positive, common-center Lorentzian difference for the soft term.
   Its leading inverse-square tail cancels exactly.
3. Remove low-velocity suppression from the hard-secondary limit, using the
   existing broad width as the interpolation scale.
4. Replace the global free-electron cutoff with a smooth bound-electron
   edge and tail, while retaining the exact free Tmax as its reference.
5. Set the width numerator to zero, retain the printed distortion exponent,
   and replace the fitted logarithm scale with the retained mean binding.

This keeps the PJG architecture rather than replacing it by Rudd's SDCS or
a separately parameterized generalized-oscillator-strength model.

## 2. Kinematics and complete SDCS

Use energy units throughout, with c = 1: E is proton kinetic energy, T is
electron kinetic energy, m is electron rest energy, and M is proton rest
energy. Cross sections below are in cm^2 and energies in eV. Define

\[
\gamma=1+E/M,\qquad
\beta^2=E(E+2M)/(E+M)^2,\qquad E_e=m\beta^2/2,
\]
\[
T_m=\frac{2mE(E+2M)}{(M+m)^2+2mE},\qquad C_B=4\pi a_0^2 R_y^2.
\]

The production kinematic header is already tested independently against
invariant two-body kinematics and a tree-level Dirac trace. The spin term is
the pointlike spin-1/2 result. This is not a claim about proton form factors,
radiative corrections, or the spin term of every bare nuclear species.

With retained PJG fractions f_j and thresholds I_j, set

\[
\bar I=\sum_j f_j I_j,\qquad
L(E)=\ln(4E_e\gamma^2/\bar I+e)-\beta^2,
\]
\[
T_0(E)=T_s-q_t\frac{t_a}{E_e+t_b},\quad
L_n=\frac{1}{(T-T_0)^2+\Gamma^2},\quad
L_b=\frac{1}{(T-T_0)^2+\Gamma^2+\Lambda^2},
\]
\[
D(E)=\frac{1}{1+[J/(mE/M)]^p},\qquad
D_h(E,T)=1-\frac{1-D(E)}{1+(T/\Lambda)^2}.
\]

Use unbounded mE/M in D rather than bounded E_e so that the distortion
actually tends to one. At nonrelativistic energies they agree to leading
order. The selected spectrum is

\[
\boxed{\frac{d\sigma}{dT}=\frac{G(E,T)}{E_e}
\left[D K\Gamma^2 L(E)(L_n-L_b)
+D_h N_e C_B L_b\widetilde F_B(E,T)\right].}
\]

Here N_e is 14 for N2 and 16 for O2. K is the fitted multiplier times the
printed PJG K. The source total is the integral of this same spectrum; there
is no independently normalized total-cross-section table.

### Hard term and the reason for the Lorentzian change

Within the free domain,

\[
F_B(E,T)=1-\beta^2 T/T_m+T^2/[2(E+M)^2].
\]

PJG placed its inverse-square coefficient in the empirical soft term, while
adding only the remaining hard terms. Correcting those remaining terms and
Tmax cannot enforce the coefficient of T^-2: it still depends on the fitted
soft amplitude, width, logarithm, and broad subtraction. A total-only refit
changes that coefficient again.

Instead use the identity

\[
L_n-L_b=
\frac{\Lambda^2}{[(T-T_0)^2+\Gamma^2]
[(T-T_0)^2+\Gamma^2+\Lambda^2]}>0.
\]

Thus L_b is asymptotically T^-2, the soft difference is T^-4, and D_h tends
to one. For T large compared with binding and width scales, but away from
the molecularly broadened binary edge where G differs from one,

\[
\frac{d\sigma}{dT}\longrightarrow
\frac{N_e C_B}{E_e T^2}F_B(E,T).
\]

The finite T^-4 choice is a minimal positive PJG-shaped interpolation, not
an exact molecular photoionization power law at arbitrarily large T.

Above free Tmax, do not evaluate the free polynomial at a forbidden transfer.
For a positive, continuously differentiable continuation use

\[
\widetilde F_B=F_B(E,T_m)
\exp\!\left[\frac{\partial_TF_B(E,T_m)}{F_B(E,T_m)}(T-T_m)\right].
\]

The endpoint value is positive and its slope is nonpositive. This continuation
is an empirical bound-tail prescription, not a QED derivation of that tail.
Holding the endpoint factor constant instead changes the tested total yields
by less than 6.1e-7 relatively and mean secondary energies by less than
2.2e-6 (5, 100, 1000, and 800000 keV). This narrow sensitivity test does not
bound uncertainty in the much more important molecular edge shape itself.

### Molecular edge and exact available support

Use a Rudd-shaped gate separately for each retained PJG continuum:

\[
G=\sum_j f_j\Phi_j\,
\operatorname{expit}\!\left[
\alpha\frac{T_m-2w_j-R_y/4-T}{w_j}\right],\quad
w_j=g(E)\sqrt{E_e I_j}.
\]

The fixed terms are Rudd's binary-edge shift and scale; alpha is refitted.
This introduces one edge-shape freedom while eliminating several PJG
subtraction, width, and cutoff parameters. The gate center is not a new
definition of Tmax, and the removed delta is not hidden inside free kinematics.

For a relativistic projectile, set r=m/M and

\[
g(E)=\frac{\gamma(\gamma+r)(1+r)^3}{(1+2\gamma r+r^2)^2}.
\]

This factor has no fitted coefficient. Its origin is the sensitivity of the
free binary endpoint to a small initial longitudinal electron momentum k.
Let P=E+M, p_p=sqrt(E(E+2M)), epsilon=sqrt(m^2+k^2), and
s(k)=M^2+m^2+2(P epsilon-p_p k). An independent collinear two-body solution is

\[
T(k)=\frac{k^2}{\epsilon+m}
+\frac{2(p_p+k)(p_p\epsilon-Pk)}{s(k)}.
\]

Differentiating at zero gives

\[
T'(0)=-\frac{2p_p(P+m)(M^2-m^2)}{s(0)^2}
=-2\beta\frac{1-r}{(1+r)^2}\,g(E).
\]

Thus g is the exact relativistic/NR sensitivity ratio at the same speed,
and tends to one nonrelativistically. Using it as a width multiplier is
still an approximate impulse-model prescription, not an exact bound-state
momentum convolution.

The factor Phi_j supplies exact three-body support and a distant endpoint
taper. For neutral rest energy A and binding I, the residual composite
minimum mass is M+A-m+I and the reaction threshold is

\[
E_{\rm th}=I(1+M/A)+I^2/(2A).
\]

Let p_p=sqrt(E(E+2M)), p_e=sqrt(T(T+2m)),
C_0=(A-m)(E-E_th)-mI(M+I/2)/A, and C=(E+M+A)T-C_0.
The invariant residual-mass condition gives cos(theta) >= C/(p_p p_e).
Use Phi=clip((1-C/(p_p p_e))/2,0,1), with the appropriate zero-momentum
limit, and zero outside the exact molecular endpoint. This isotropic
allowed-angle fraction is an endpoint taper only, not the angular cross
section. It is effectively one throughout the fitted spectra.

### Soft/Bethe limit: what is and is not enforced

For fast projectiles and fixed low T, the logarithmic coefficient is

\[
a(T)=\frac{K\Gamma^2}{C_B}(L_n-L_b),\qquad
S_{\rm soft}\sim\frac{C_B}{E_e}a(T)
[\ln(4E_e\gamma^2/\bar I)-\beta^2].
\]

The ln(gamma^2) growth and the transverse -beta^2 subtraction have the same
coefficient, as required by the relativistic dipole limit. The optical
diagnostic compares a(T) with sum_j[(df_j/dW)/W] at W=T+I_j. The non-logarithmic
target response is still approximated by the PJG shape and logarithm scale;
this is not a full relativistic PWBA calculation with an exact molecular
response. Also, beta approaching zero and T approaching zero are different
limits: Bethe theory does not supply low-projectile-energy collision dynamics.

The spectrum is finite and nonnegative at T=0. It does **not** reproduce the
pointwise threshold optical coefficient exactly. On the 40-point 2--100 eV
fit grid, the model/reference optical ratios span 0.867--1.425 for N2 and
0.652--1.387 for O2. At T=0 they
are 1.80 and 3.88, respectively; the coarse O2 bins are not a precision
threshold reference, but this discrepancy must not be concealed. Enforcing
the whole optical threshold structure would need more molecular response
detail or additional shape flexibility, contrary to the current parsimonious
PJG-preserving choice.

The effective continuum-strength integrals int (T+mean(I))*a(T) dT are
10.897 (N2) and 12.472 (O2), below electron counts 14 and 16. This is a useful
budget check, not a proof of exact oscillator-strength sum rules for an
inclusive decay spectrum. Ionization energy moments at 2--1000 MeV remain
below a simple Bethe total-stopping estimate (ratios about 0.896--0.971).
That estimate excludes density, shell, and higher-order corrections; it
does not establish the excitation budget or constitute a PSTAR validation.

The subsequent [optical and stopping audit](OPTICAL_STOPPING_AUDIT.md)
compares actual NIST PSTAR tables. It supersedes any inference from that
simple Bethe estimate that the excitation budget is established: baseline
O2 reaches 1.00573 of NIST electronic stopping at 1.75 MeV. The audit also
records denser optical sampling, parameter-count-neutral optical trials,
and independent first/second-moment tests. The baseline coefficients below
remain unchanged; the O2 trial is not promoted.

## 3. Refit, data, and parameter reduction

The selected fit has six adjustable quantities per gas:

| Quantity | N2 | O2 |
|---|---:|---:|
| J (eV) | 14.20737 | 6.960482 |
| K / printed K | 1.279760 | 1.078892 |
| Gamma (eV) | 11.31184 | 14.15467 |
| Lambda (eV) | 99.98727 | 155.74550 |
| q_t, multiplier of printed t_a | 0.7104965 | 0.8406507 |
| alpha | 0.7871437 | 0.6825391 |

Keep p=1+nu at 0.807 and 1.314, respectively. Keep printed T_s, t_a, t_b,
f_j, I_j, and electron counts. Set gamma_1=0, so gamma_2 is unused. Remove
independent roles for B0, B1, E0, the broad center T1, and delta. Fix the
logarithm with mean(I)=19.59248 and 19.945 eV instead of a fitted C scale.
The apparent collection of printed C_j is already almost proportional to I_j;
it should not be counted as many genuinely independent fit parameters.

Objective blocks are logarithmic residuals, each divided by sqrt(number of
samples): recommended totals weight 1; recommended SDCS weight 1; optical
coefficient weight 0.25; original N2 1979 totals weight 0.4; original N2 1979
mean electron energies weight 0.5. These are engineering compromise weights,
not inverse uncertainties. The dense recommended curves and old measurements
are correlated; the objective is not a chi-squared statistic.

- Rudd (1985) recommended electron-production totals: 31 logarithmic proton
  energies from 5 to 4000 keV.
- Rudd (1992) recommended SDCS: N2 at 5, 10, 30, 50, 100, 300, 1000 keV;
  O2 at 7.5, 10, 30, 50, 100, 300 keV. Up to 48 logarithmic secondary
  energies per proton energy, starting at 2 eV and excluding values below
  1e-4 of that spectrum's value at 2 eV. These are reference-curve samples,
  not newly digitized independent experimental points.
- Original Crooks/Rudd (1971) Table I totals are reported as a separate
  check, not fitted as independent calibrations.
- Original Rudd (1979) N2 Table I totals and mean energies are retained
  without renormalizing them to the later recommendation.
- N2 optical shape: Rudd (1992), Table II; O2: NCAR/GLOW's coarse partial
  photoionization table and its stated thresholds. Both use 40 logarithmic
  secondary energies from 2 to 100 eV, not an exact threshold constraint.

Parameter-removal trials support the selected simplification. Freeing a
seventh logarithm-scale parameter gives objective norms 0.155743/0.100738
(N2/O2), compared with 0.155913/0.100767 with the mean-binding scale fixed.
Allowing an energy-dependent width made only small changes; freeing the
printed exponent was not compelling. Fixing q_t=1 worsened the objectives
to 0.167261/0.113114. Four independent multiplicative starts, each differing
by up to exp(0.7), reached the same selected minimum; maximum parameter
differences were 6.5e-6 relatively for N2 and 8.5e-7 for O2. These are
numerical reproducibility checks, not parameter confidence intervals.

## 4. Quantitative comparisons

All printed comparisons include the 1977 erratum. "Corrected" changes free
Tmax, repairs the hard remainder algebra, and sets delta=0, without refitting.
"Previous refit" additionally uses J=59.65 for N2, J=19.28 and
K=4.581e-16 cm^2 for O2; it retains the old unnormalized soft tail.

Total / Rudd recommended total over 201 incident energies, 5--4000 keV:

| Model | N2 minimum--maximum | O2 minimum--maximum |
|---|---:|---:|
| Printed + erratum | 1.081--2.269 | 0.118--1.435 |
| Corrected, no refit | 1.099--6.915 | 0.451--1.454 |
| Previous corrected refit | 0.813--1.190 | 0.795--1.110 |
| New matched PJG | 0.925--1.054 | 0.999--1.089 |

For 800-MeV protons, SDCS / complete Bhabha at secondary energies
100 keV, 500 keV, and 2 MeV:

| Model | N2 | O2 |
|---|---|---|
| Printed + erratum | 0.7472, 0.1516, outside printed support | 1.0633, 0.5100, outside printed support |
| Corrected, no refit | 0.8712, 0.8544, 0.7091 | 1.1873, 1.2128, 1.4257 |
| Previous corrected refit | 0.8700, 0.8532, 0.7081 | 0.8216, 0.7981, 0.5965 |
| New matched PJG | 1.000079, 1.000016, 1.000004 | 1.000128, 1.000025, 1.000006 |

The new spectrum is deliberately broadened near the binary endpoint; exact
agreement with a stationary free electron is not claimed there. Its high-energy
agreement is a theoretical matching check, not a comparison with relativistic
N2/O2 SDCS measurements.

On the fitted SDCS grid, the central 80% of model/reference ratios are
0.874--1.203 for N2 and 0.913--1.043 for O2. The full ranges are
0.767--1.373 and 0.692--1.207. A denser independent grid of 51 incident
energies, extending to 1700 keV for N2 and 300 keV for O2, gives full ranges
0.738--1.377 and 0.692--1.209. Thus this is a useful smooth approximation,
not a uniformly few-percent SDCS fit. On 0--2 eV the dense reference ratios
span 0.613--1.223 and 0.538--0.988; the lowest-energy electron shape remains
an approximation, although those electrons are included in the integral
and must not be removed from the kinetic source.

The previous corrected refit's strict free cutoff discarded nonzero Rudd
SDCS at 13.75% of the sampled N2 points and 15.66% of O2 points. Those are
fractions of comparison-grid points, not fractions of the total yield.
The new above-free-Tmax yield fractions at 5 keV are 17.43% and 23.34%;
at 100 keV, 0.290% and 0.433%. O2's differential curve at 5 keV is an
extrapolation below its 7.5-keV experimental coverage.

New / original 1971 total ratios span 0.968--1.066 for N2 and 0.949--1.031
for O2. The 1979 N2 5-keV point remains the largest original-total mismatch:
new total 1.842e-16 cm^2 versus 2.460e-16 cm^2, and new mean secondary
energy 6.403 eV versus 5.420 eV. At 10 keV the corresponding values are
2.935e-16 versus 3.220e-16 cm^2 and 8.474 versus 7.740 eV. Not every
individual mean-energy datum improves over the previous total-only refit.
The original paper has substantial low-energy uncertainties and a shared
normalization with the 1971 data. Its 5-keV total is already 23.8% above
the later recommended total. This later-data discrepancy must not be
misrepresented as an explanation of the earlier 1976 PJG Figure 6 inconsistency.

## 5. Numerical checks, angular closure, and production boundaries

The Python suite has 33 passing tests, including ten new matched-model tests:
the complete Bhabha bracket, positive C1 continuation, independent moving-target
derivative, realistic hard-secondary matching, relativistic logarithmic
coefficient, positivity from 20 eV to 1 PeV incident energy, exact molecular
support, independent CM endpoint, quadrature, and conditioned angular closure.
This wide arithmetic stress test does not establish physical accuracy over
that whole incident-energy range. The fitted-range quadrature differs by
less than 7e-11 in totals and means when increasing 1025 to 8193 points.
At relativistic energies the narrow binary edge requires a finer integration
grid; 8193-point quadrature agrees with independent adaptive integration to
the tested 2e-7 relative tolerance. Test tolerances were not loosened.

The independent C++ kinematics and mixed-parent sampling tests also pass.
Their host quantile microbenchmark remains about 0.8 ns/quantile on this
machine, but that is neither a complete collision timing nor a GPU result.
No CUDA, HIP, or SYCL performance claim follows from host tests.

The research angular closure retains practical IAA's bound/free interpolation.
Let t_b=min(T,T_m) for this paragraph only, and

\[
\mu_f=\sqrt{\frac{t_b(T_m+2m)}{T_m(t_b+2m)}},\quad
a=\mu_f\frac{T+I/2}{T+I},\quad b=\frac{I}{T+I}.
\]

Sample cos(theta) uniformly on [a-b,a+b] intersected with the physical range
and the allowed molecular recoil cap. Condition the interval rather than
clipping samples, which would generate artificial endpoint point masses.
It tends to isotropy at T=0 and to heavy-projectile binary kinematics for
hard electrons. Above free Tmax, the reference direction is forward with a
finite binding width; this is an empirical closure, not a forbidden free
collision. A near-threshold nonoverlap falls back to the allowed angular cap.
The closure supports a general projectile mass but has not been fitted to
N2/O2 DDCS. The SDCS coefficients reported here are proton fits, not new
bare-ion cross-section measurements or a derivation of arbitrary-ion spin terms.

For production, these formulas belong in host-side table construction. The
GPU sampler can retain fixed-cost inverse-CDF interpolation and randomized
conditional stratification; there is no need to evaluate exponentials,
sum molecular continua, or perform rejection loops per emitted electron.
The conditional mean binding can be pretabulated from the same gate weights
for the effective angular closure. Table accuracy, full operator cost, and
accelerator execution must be verified during production integration.

The prescribed beam remains rigid/ballistic and the effective ion is drawn
from the neutral thermal distribution. Exact molecular support specifies
kinematically possible electrons; it does not make that simplified source
event-by-event energy/momentum conserving. Recoil and projectile depletion
are deliberately not deposited into the tracked beam/thermal ion. These
approximations must remain explicit in production documentation.

## 6. Reproduction and sources

Run from the repository root with NumPy and SciPy available:

```sh
python Tools/Algorithms/ProtonImpactIonization/fit_pjg_matched.py \
  --target N2 --free-edge --binding-log
python Tools/Algorithms/ProtonImpactIonization/fit_pjg_matched.py \
  --target O2 --free-edge --binding-log --o2-optical /path/to/ephoto_xo2.dat
python Tools/Algorithms/ProtonImpactIonization/audit_pjg_matched.py \
  --figures tmp/pjg-matched-results
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
```

The optional figures also need Matplotlib. The audit and tests do not need
the external GLOW table; reproducing the O2 refit does. GLOW's table is not
needed to evaluate the selected O2 model. Add `--multistart 4` to a refit
command to repeat the deterministic starting-point check. GLOW's table is not
redistributed under WarpX's license. The source revision is
`6cb880d0e112810dbe9a79203f0183f9b07e3473`, with SHA-256
`9d803731e45ce7906d60d190681968c21fdb0ae9564c24fe98cfac58cd1b1eb2` for
`data/ephoto_xo2.dat`. Its ionization column is used, not total absorption;
rounded branch fractions are normalized and evaluated at T+I_j. GLOW's
four O2 thresholds are 12.07, 16.10, 18.20, 20.00 eV. The finite optical
audit window is below the Auger region; it is not an inclusive Auger response.

- Porter, Jackman, Green (1976): https://doi.org/10.1063/1.432812.
- Erratum (1977): https://doi.org/10.1063/1.323427.
- Rudd et al. (1985), recommended total yields:
  https://doi.org/10.1103/RevModPhys.57.965.
- Rudd et al. (1992), especially Eqs. (32)--(33), (41)--(48), Tables I, II, V:
  https://doi.org/10.1103/RevModPhys.64.441.
- Crooks and Rudd (1971), Table I: https://doi.org/10.1103/PhysRevA.3.1628.
- Rudd (1979), Table I: https://doi.org/10.1103/PhysRevA.20.787.
- NCAR/GLOW optical input:
  https://github.com/NCAR/GLOW/blob/6cb880d0e112810dbe9a79203f0183f9b07e3473/data/ephoto_xo2.dat.
- NCAR/GLOW thresholds, ionization conversion, and Auger treatment:
  https://github.com/NCAR/GLOW/blob/6cb880d0e112810dbe9a79203f0183f9b07e3473/ephoto.f90.
