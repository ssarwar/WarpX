> Historical research note, recovered 2026-09-14.
> Model names such as 'current' refer to this note's stage, not today's production model.
> Original commands and links retain their historical context; see the
> [archive guide](../../README.md) for status and portable reproduction.

# Calibrated PJG proton energy model for N2 and O2

The energy-model formulation and joint calibration are complete. The selected
reference is `calibrated_pjg.py`: six adjustable coefficients for N2 and five
for O2. The angular model is outside this decision. No production collision
operator or energy table has been replaced, and no commit or push is made.

The observable is inclusive electron-production yield, represented by one
effective electron/ion pair. It is not exclusive single molecular ionization.
All secondary energies from zero through the molecular endpoint are included.
The prescribed projectile remains rigid; the energy-loss moments below are
diagnostics, not an instruction to apply collisional slowing to the beam.

## 1. Formulation

Use eV for energies, cm2 per molecule for total cross sections, and cm2/eV
per molecule for the SDCS. Let E be proton kinetic energy, T electron kinetic
energy, m and M the electron and proton rest energies, and

```math
\gamma=1+E/M,\qquad E_e=m\beta^2/2,\qquad K_{\rm nr}=mE/M,
\qquad C_B=4\pi a_0^2R_y^2.
```

The exact free stationary-electron maximum is

```math
T_m=\frac{2mE(E+2M)}{(M+m)^2+2mE}
   =\frac{2m\beta^2\gamma^2}{1+2\gamma m/M+(m/M)^2}.
```

The selected SDCS is

```math
S(E,T)=\frac{G(E,T)}{E_e}\left[
 D(E)K\Gamma^2 L(E)(L_n-L_b)
 +D_h(E,T)N_e C_B L_b\widetilde F_B(E,T)\right],
```

where

```math
T_0=T_s-\frac{q_t t_a}{E_e+t_b},\qquad
L_n=\frac{1}{(T-T_0)^2+\Gamma^2},\qquad
L_b=\frac{1}{(T-T_0)^2+\Gamma^2+\Lambda^2},
```

```math
D=\frac{1}{1+(J/K_{\rm nr})^p},\qquad
D_h=1-\frac{1-D}{1+(T/\Lambda)^2},\qquad
L(E)=\ln\left(4E_e\gamma^2/\bar I+e\right)-\beta^2,
\qquad \bar I=\sum_j f_j I_j.
```

Here the `e` inside the logarithm is Euler's number, not electron charge.
The mean binding scale is not the stopping-power mean excitation energy
(82/95 eV in PSTAR), and is not an exact energy-loss-dependent Bethe constant.
The fixed PJG continuum allocation is

| Target | I_j (eV) | f_j | N_e | Mean I_j (eV) |
|---|---|---|---:|---:|
| N2 | 15.58, 16.73, 18.75, 22, 23.6, 40 | .456, .2, .104, .07, .07, .1 | 14 | 19.59248 |
| O2 | 12.1, 16.1, 16.9, 18.2, 20.3, 23, 37 | .08, .19, .19, .17, .11, .16, .1 | 16 | 19.945 |

### Hard term and positive bound tail

Inside the free domain, the complete pointlike spin-1/2 tree-level factor is

```math
F_B(E,T)=1-\beta^2 T/T_m+\frac{T^2}{2(E+M)^2}.
```

This factor and its normalization are not fitted. For T greater than T_m,
continue it as

```math
\widetilde F_B(E,T)=F_B(E,T_m)
 \exp\left[\frac{\partial_T F_B(E,T_m)}{F_B(E,T_m)}(T-T_m)\right].
```

The continuation matches value and derivative, stays positive, and introduces
no coefficient. It is an empirical bound-electron continuation, not a claim
that forbidden free scattering occurs. The earlier exponential/constant
continuation sensitivity checks are recorded in `MATCHED_PJG.md`.

The smooth molecular gate is

```math
G(E,T)=\sum_j f_j\Phi_j(E,T)
 \operatorname{expit}\left[\alpha\frac{T_m-2w_j-R_y/4-T}{w_j}\right],
\qquad w_j=g(E)\sqrt{E_e I_j},
```

```math
g(E)=\frac{\gamma(\gamma+r)(1+r)^3}{(1+2\gamma r+r^2)^2},\qquad r=m/M.
```

This retains the binary-edge form used in the preceding matched-PJG study.
The broadening factor follows from differentiating the exact binary endpoint
with respect to initial longitudinal electron momentum and dividing by its
nonrelativistic value. It is not an exact convolution with a molecular
momentum distribution. Alpha is 0.70 times the fitted edge multiplier for N2
and fixed at 0.59 for O2, using the Rudd edge scales [6].

Phi_j is only a parameter-free molecular-support taper. With neutral rest
energy A, residual minimum mass R=M+A-m+I_j, projectile momentum p_E, and
electron momentum p_T, the recoil condition gives

```math
c_{\min}=\frac{R^2-[(E+M+A-m-T)^2-p_E^2-p_T^2]}{2p_Ep_T},
\qquad \Phi_j=\operatorname{clip}[(1-c_{\min})/2,0,1].
```

It is set to zero below the exact channel threshold
I_j(1+M/A)+I_j^2/(2A) and at or above its exact three-body electron endpoint.
The stable endpoint and zero-momentum evaluations are in `pjg_matched.py`;
independent invariant and recoil-quadratic tests cover them. This phase-space
taper is effectively one over the fitted spectra, not an angular closure.

### Why the limits improve

The identity

```math
L_n-L_b=\frac{\Lambda^2}
 {[(T-T_0)^2+\Gamma^2][(T-T_0)^2+\Gamma^2+\Lambda^2]}>0
```

makes the soft term positive and asymptotically proportional to T^-4.
The broad hard term is proportional to T^-2. Thus in the overlap region
where binding and widths are small compared with T, but the molecular edge
is still distant, G and D_h approach one and

```math
S(E,T)\longrightarrow\frac{N_e C_B}{E_e T^2}F_B(E,T).
```

This is an overlap limit, not T going to infinity at fixed finite projectile
energy. Near the binary edge, molecular broadening is deliberately retained.
The exact factor itself is positive because it can be written as
(1-x)+x/gamma^2+T^2/[2(E+M)^2], with x=T/T_m between zero and one.

At T=0 the model is finite and positive, without a low-energy electron cut.
The soft logarithm contains the relativistic dipole enhancement ln(gamma^2)
and transverse subtraction -beta^2. Its spectral coefficient is an optical-
constrained PJG approximation, not an exact molecular PWBA/GOS response.
Narrow resonances, detailed thresholds, and Auger peaks are not resolved.
Keeping that distinction is important: this refit does not prove exact Bethe
or oscillator-strength sum-rule agreement at every secondary energy.

## 2. Original and final parameters

"Original" always includes the 1977 corrections in Garvey, Porter and Green,
Ref. 26 [2]. In particular the apparent superscript plus in PJG Eq. (16)
means ordinary addition. It is not positive-part clipping. The Eq. (14) and
Table II corrections belong to the GOS construction, and the negative O2
E_Gamma correction belongs to the electron branch, not the proton rows.

The printed proton model put its inverse-square coefficient inside the
empirical Lorentzian combination and added only a hard remainder. Fixing
that remainder and T_m alone did not protect the Bhabha coefficient. The
present positive Lorentzian difference removes that overlap and fixes the
complete hard term separately.

| Quantity | Original N2 | Final N2 | Original O2 | Final O2 |
|---|---:|---:|---:|---:|
| J (eV) | 3.39 | 13.82353 | 40.3 | 8.76636 |
| p=1+lambda | .807 | .807, fixed | 1.314 | 1.314, fixed |
| K (cm2) | 7.58e-16 | 6.67952e-16 | 6.55e-16 | 5.11905e-16 |
| Gamma_s (eV) | 11.1 | 12.50323 | 13.1 | 16.75750 |
| Width numerator gamma_1 (eV2) | 1.27e4 | 0 | 5e5 | 0 |
| Width denominator gamma_2 (eV) | 1810 | removed | 76000 | removed |
| Peak numerator t_a (eV2) | 20300 | 3635.677 | 2520 | 943.716 |
| T_s (eV) | 4 | 4, fixed | 6.34 | 6.34, fixed |
| t_b (eV) | 1970 | 1970, fixed | 128 | 128, fixed |
| Separate broad center T_1 (eV) | 53.3 | use T_0 | 68.3 | use T_0 |
| Original broad width (eV) | 115 | replaced by Lambda=88.54008 | 189.1 | replaced by Lambda=159.26980 |
| delta (eV) | 84 | removed | 132.1 | removed |
| Empirical B(E_e) | fitted | removed | fitted | removed |
| Edge multiplier | no smooth molecular edge | 1.191140 | no smooth molecular edge | 1, fixed |

Lambda is a width excess: the actual broad width is sqrt(Gamma^2+Lambda^2).
It is not numerically identical to the printed broad-width parameter.
The peak multipliers q_t are .1790974103 and .3744906480; they replace the
existing t_a values and do not add a second peak-numerator parameter.
The printed C_j constants are replaced by the one fixed mean-binding scale
in L(E). The gamma_1=0 reduction removes the entire energy-dependent width
numerator/denominator pair. The original units printed beside numerator
parameters do not override the eV2 dimensions required by their equations.

The final free coefficients are J, K, Gamma, Lambda, q_t, and N2's edge
multiplier. O2's unrestricted sixth coefficient was 1.0007327: fixing it to
one changes the objective norm from .26689724 to .26689843. It is unnecessary.
For N2, fixing the edge to one worsens the norm from .24739024 to .26614192.
Setting the peak numerator to zero instead gives .25321057 and worsens the
maximum measured-mean ratio from 1.279 to 1.335. Retain six for N2.
These norms are model-selection diagnostics, not statistical significance.

## 3. Data and weighting

The source search and full-text audit are in `DATASET_AUDIT.md`. A recent
publication date is not sufficient to displace a better normalized evaluated
spectrum. The final fit uses these sources with distinct roles:

| Source | Role |
|---|---|
| Crooks and Rudd 1971, Table I [3] | Measured N2/O2 totals, six energies each |
| Rudd 1979, Table I [4] | Measured N2 totals and mean electron energies, 5--70 keV |
| Rudd et al. 1983, Table V [5] | Authors' fit to measured electron yields, 18 rows 5--3000 keV; not raw points |
| Crooks/Rudd and Toburen markers reproduced in PJG Fig. 5 [1,3,7] | 49 N2 measured SDCS markers, 50/100/300/1000 keV |
| Cheng et al. 1989, Fig. 1 [8] | 30 O2 measured SDCS markers, 7.5/50/150 keV |
| Rudd et al. 1992 [6] | Weak analytic SDCS interpolation prior, not measured markers |
| NIFS-DATA-109 and Leiden [9,10] | N2 and O2 optical constraints on W=25--100 eV, respectively |
| Rudd 1979, Fig. 7 [4] | 18 N2 held-out spectral markers at 5/20/70 keV; no additional fit weight |
| NIST PSTAR [11] | Independent stopping-budget comparison; no fitting penalty |

The low-energy N2 held-out markers are not an entirely independent experiment:
the same paper supplies fitted totals and means. Cheng's displayed SDCS are
rescaled to recommended totals. Crooks' absolute uncertainty is about 17%;
Toburen quotes about 25% for most secondary energies. The Rudd-1983 estimated
uncertainties decrease from 25% at 5 keV to 8% at and above 500 keV. Shared
normalizations are not counted as independent precise point constraints.
Rudd's lowest-energy electron measurements are specifically less reliable.

Fit log residuals. Normalize each block by the square root of its point
count. The Rudd-1983 total block has coefficient 0.15 divided by the published
fractional uncertainty. Crooks totals have weight .5, Rudd-1979 totals .4,
and Rudd-1992 SDCS .25. The N2 means have weight one. Multiply the 1979
total/mean weights by .5 at 5 keV, 2/3 at 10 keV, and one at 30 keV and above,
interpolating in log energy. The optical block weights are one for N2 and two
for O2; O2 has an ionization column, while N2 uses an absorption proxy.

For each measured SDCS spectrum, profile out its shared log normalization
and fit its shape. Half-weight T below 10 eV. The final absolute SDCS is not
rescaled after fitting; its normalization comes from the coupled total and
spectral model. These are explicit judgment weights, not fabricated error
bars, confidence intervals, or chi-squared. Small digitization uncertainties
are additional to experimental uncertainties.

The optical forward diagnostic is

```math
(df/dW)_{\rm model}=W\sum_j f_j a(W-I_j)\Theta(W-I_j),\qquad
a(T)=K\Gamma^2(L_n-L_b)/C_B
```

using the high-energy limiting peak. Its channel-opening steps near 40 eV
for N2 and 37 eV for O2 are approximation artifacts, not measured resonances.
The evaluated photoionization-to-oscillator conversion is 109.76097 Mb per
eV^-1. Optical W is energy loss/photon energy, not electron kinetic T.

Gallagher et al. (1988), Mahla and Mehnen (2025), and the N2/O2 finite-q
studies were checked, but are not additional numerical fit blocks [12--15].
The 2025 O2 calculation has shifted raw thresholds and incompletely closing
tabulated partials; it is a comparison, not the absolute anchor. The finite-q
studies constrain the interpretation but do not determine a unique complete
GOS from this already momentum-integrated SDCS. No quantitative GOS fit is
claimed. Lopez-Patino et al. (2016) remains excluded because full text and
inclusive-yield channel semantics were not available [16]; this does not
block the present calibration.

## 4. Comparison results

All ranges below refer to the named comparison grid, not a uniform error
bound over all energies or a fitted uncertainty interval.

| Model/1983 authors' total fit | Original + 1977 | Earlier corrected total-only refit | Final |
|---|---:|---:|---:|
| N2 | 1.079--2.007 | .796--1.093 | .834--1.096 |
| O2 | .090--1.424 | .605--1.012 | .879--1.023 |

Final total log-RMS differences are .0908 (N2) and .0695 (O2). Ratios to
Crooks' measured totals are 1.016--1.174 and .976--1.177. Some differences
remain outside an individual quoted nominal band, particularly high-energy
total normalization. The model is a joint compromise, not the closest
possible total-only fit. The old N2 curve printed in PJG Fig. 6 remains
inconsistent with evaluating its corrected printed formula and parameters;
the comparison does not silently substitute that plotted curve.

| Measured spectral subset | Final/data range | Shape log-RMS after shared normalization |
|---|---:|---:|
| N2, 50 keV, PJG markers | 1.009--1.519 | .128 |
| N2, 100 keV, PJG markers | .919--1.509 | .151 |
| N2, 300 keV, PJG markers | .784--1.569 | .184 |
| N2, 1 MeV, PJG markers | .939--1.244 | .087 |
| N2, 5 keV, held-out 1979 markers | .857--1.131 | .098 |
| N2, 20 keV, held-out 1979 markers | .571--1.191 | .285 |
| N2, 70 keV, held-out 1979 markers | .962--1.336 | .124 |
| O2, 7.5 keV, Cheng markers | 1.101--1.345 | .073 |
| O2, 50 keV, Cheng markers | .873--1.481 | .158 |
| O2, 150 keV, Cheng markers | .750--1.278 | .151 |

The shape diagnostic removes one mean log ratio per displayed energy; the
300-keV row pools the two measured sources. The fit itself profiles those
sources separately. The largest held-out N2 discrepancy is the last 20-keV spectral point near
T=160 eV. The O2 residuals include a roughly 48% excess near T=202 eV at
50 keV. These are not hidden by renormalizing plotted curves.

For O2, compare means only over the observed windows, without extrapolating
digitized markers into an invented complete experimental spectrum:

| Proton energy | Measured window (eV) | Digitized-window mean T (eV) | Model, same window (eV) |
|---|---:|---:|---:|
| 7.5 keV | 3.04--40.31 | 11.779 | 12.040 |
| 50 keV | 2.95--294.71 | 25.890 | 25.037 |
| 150 keV | 2.90--390.69 | 42.494 | 39.081 |

These are means of piecewise log-log interpolated digitizations, not the
authors' full-spectrum means. The corresponding window-integral ratios are
1.265, 1.123, and 1.067. Full-spectrum model means follow below.

### Mean energies, tails, and stopping

Define sigma=integral S dT, M1=integral T S dT, and M2=integral T^2 S dT.
The effective-pair loss cross section is integral [T+I_eff(E,T)] S dT, where
I_eff is the gate-weighted conditional PJG threshold. Convert eV cm2 per
molecule to MeV cm2/g by multiplying N_A divided by molar mass and by 1e-6.
Do not multiply again by gas density. This is not the cascade-inclusive W
value and is not an exact projectile stopping reconstructed from inclusive
electron yield. Excitation, detailed decay multiplicities and recoil losses
are absent.

| Proton energy | N2 sigma (cm2) | O2 sigma (cm2) | N2 mean T (eV) | O2 mean T (eV) |
|---|---:|---:|---:|---:|
| 5 keV | 1.75922e-16 | 1.53538e-16 | 6.935 | 9.035 |
| 10 keV | 2.98413e-16 | 2.84783e-16 | 9.109 | 11.609 |
| 50 keV | 6.48962e-16 | 6.09444e-16 | 19.827 | 21.997 |
| 100 keV | 5.91061e-16 | 5.70481e-16 | 28.847 | 29.986 |
| 1 MeV | 1.35471e-16 | 1.43274e-16 | 52.848 | 57.424 |
| 800 MeV, extrapolation | 9.32765e-19 | 9.79091e-19 | 68.491 | 78.468 |

N2's measured means at 5/10/30/50/70 keV are 5.42/7.74/15.2/20.4/24.4 eV;
the model gives 6.935/9.109/15.039/19.827/23.907 eV. The 5-keV mean is still
28% high. Forcing it exactly would overstate the low-energy source accuracy.

| Proton energy | N2 pair loss / PSTAR electronic | O2 pair loss / PSTAR electronic |
|---|---:|---:|
| 5 keV | 100.54 / 283.0 | 84.17 / 222.4 |
| 50 keV | 541.98 / 722.9 | 476.05 / 584.4 |
| 100 keV | 611.61 / 759.4 | 532.78 / 643.3 |
| 1 MeV | 210.96 / 225.9 | 208.61 / 216.1 |
| 800 MeV, extrapolation | 1.7663 / 2.081 | 1.8134 / 2.050 |
| 10 GeV, extrapolation | 1.7153 / 2.062 | 1.7799 / 2.040 |

All stopping entries are MeV cm2/g. The maximum ratios over PSTAR grid
energies from 5 keV through 10 GeV are .93850 (N2 at .7 MeV) and .96534
(O2 at 1 MeV). In particular the relativistic first moment does not run above
the reference electronic budget. O2 still leaves only 3.5% at its maximum
for omitted electronic losses under this pair-cost convention; this is a
limitation, not a complete stopping validation or a reason to force arbitrary
agreement with PSTAR.

At 5 keV, 19.73% of N2 and 29.49% of O2 electrons are above free T_m. They
carry 50.38% and 64.62% of the kinetic first moments. At 50 keV the count
fractions fall to .825% and 1.582%. Discarding those tails would bias both
the fitted spectrum and its mean. The molecular energy bound is preserved.

The median/90th/99th percentiles at 50 keV are 12.07/48.59/104.45 eV for N2
and 13.54/52.27/123.18 eV for O2. At 1 MeV they are 14.42/108.37/755.81
and 18.33/114.37/786.65 eV. A finite mean does not imply low sampling noise:
Var(T)/mean(T)^2 at 800 MeV is approximately 1839/1525, because rare hard
electrons dominate high moments. A future GPU implementation should retain
the full distribution and use stratified quantiles, not remove hard electrons.

### Optical and hard-limit checks

Over W=25--100 eV, optical model/data ratios are .842--1.320 for N2 and
.781--1.185 for O2, with log-RMS .140 and .097. The earlier matched candidate
had ratios 1.106--1.749 and .925--1.458. The integrated effective soft
oscillator strengths are 8.076 and 10.940, not the full electron counts 14/16.
Do not report these as complete TRK sum-rule matches: excitation and inner-
shell/channel physics cannot be reconstructed by rescaling the soft fit.

At 800 MeV, full SDCS/free-Bhabha ratios at T=1/10/100 keV are
1.00994/1.000816/1.000080 for N2 and 1.02571/1.001399/1.000128 for O2.
The earlier corrected total-only model gives .845/.871/.870 and
.801/.823/.822. The hard-secondary normalization defect has therefore been
removed without using an empirical normalization on the hard coefficient.
At 10 GeV, M2 divided by the analytic free Bhabha second moment is .999900
for N2 and .999881 for O2. These are checks of the specified pointlike
tree-level model, not proton form factors or higher-order radiative physics.

## 5. Numerical verification and reproduction

The canonical wrapper guards the audited incident-energy range 5 keV--10 GeV.
Measured-total calibration is 5--4000 keV; measured N2 SDCS used here extend
to 1 MeV and measured O2 SDCS to 150 keV. Relativistic checks are extrapolation
tests, not experimental validation at GeV energies. Near-threshold incident
energies below 5 keV require a separate validation even though all secondary
energies down to zero are retained for allowed incident energies.

Tests cover nonnegativity, finite T=0 behavior, exact support, a nonzero bound
tail, total/SDCS identity against independent adaptive quadrature, moment
inequalities, hard-tail normalization, preserved PSTAR units, parameter
reductions, normalization-invariant shape residuals, and digitization axes.
The earlier independent invariant Bhabha/PWBA/endpoint tests remain applicable.
All 67 Python checks and both host CTests pass. These host results do not
establish CUDA/HIP/SYCL execution or performance.
The 2049-to-8193-point moment changes are below 1.2e-9 for the audited energies;
quantile changes are below 5.5e-5. Alternate optimizer starts change the fitted
spectra by at most 1.2e-6 for N2 and 7.4e-8 for O2, not an uncertainty estimate.

Use NumPy, SciPy, and Matplotlib in a host Python environment:

```sh
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
python Tools/Algorithms/ProtonImpactIonization/finalize_pjg.py \
  --nifs tmp/pdfs/dataset-refresh-20260909/nifs-optical.json \
  --o2-leiden tmp/pdfs/dataset-refresh-20260909/leiden-o2.txt \
  --output output/pjg-final/fit-results.json \
  --figures output/pjg-final/figures
```

Add `--read-results` to regenerate plots without optimization. The optional
`--initial-fits` selects the earlier source-refit starting values; omitting it
starts from the frozen calibration. The output includes full coefficients,
residuals, data ratios, moments, quantiles, stopping rows, and input SHA-256.
Input extraction/download instructions are in `DATASET_AUDIT.md`.

N2 PJG Fig. 5 and Rudd Fig. 7 marker coordinates, crop/axis calibration, PDF
digests, and unit conversions are preserved in `source_datasets.py`. The
original printed comparison preserves signed SDCS values and its printed
cutoff. Negative/zero values cannot be plotted on log axes; no negative
value is clipped before calculating a total.

The final figures are `total_cross_sections.png`, `n2_sdcs.png`, `o2_sdcs.png`,
and `physics_validation.png` in the chosen figures directory. Gray means
original plus 1977, magenta the previous corrected total-only fit, and blue
the final calibration. Measured markers are distinct from the green Rudd-1983
authors' fitted total and uncertainty envelope. O2's measured markers retain
the authors' normalization; they are not renormalized to the new curve.

This Python implementation performs offline fitting and integration. It is
not intended for per-collision GPU execution. Production table generation,
interpolation error, CUDA/HIP/SYCL execution and performance remain separate
implementation/verification steps after the formulation review.

## 6. References

1. Porter, Jackman and Green (1976), JCP 65, 154,
   [PJG](https://doi.org/10.1063/1.432812), Eq. (16), Table III, Figs. 5--6.
2. Garvey, Porter and Green (1977), JAP 48, 4353,
   [1977 corrections, Ref. 26](https://doi.org/10.1063/1.323427).
3. Crooks and Rudd (1971), PRA 3, 1628,
   [measured proton totals and spectra](https://doi.org/10.1103/PhysRevA.3.1628).
4. Rudd (1979), PRA 20, 787,
   [N2 totals, spectra and means](https://doi.org/10.1103/PhysRevA.20.787).
5. Rudd et al. (1983), PRA 28, 3244,
   [proton electron-production and capture](https://doi.org/10.1103/PhysRevA.28.3244).
6. Rudd et al. (1992), RMP 64, 441,
   [SDCS recommendations](https://doi.org/10.1103/RevModPhys.64.441).
7. Toburen (1971), PRA 3, 216,
   [N2 measured spectra](https://doi.org/10.1103/PhysRevA.3.216).
8. Cheng, Rudd and Hsu (1989), PRA 40, 3599,
   [O2 measured spectra](https://doi.org/10.1103/PhysRevA.40.3599).
9. Sakamoto et al. (2010), NIFS-DATA-109,
   [evaluated oscillator strengths](https://nifs-repository.repo.nii.ac.jp/records/11706).
10. Leiden [O2 numerical ionization file](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/O2/O2.txt),
    [Heays et al. (2017)](https://doi.org/10.1051/0004-6361/201628742),
    [Hrodmarsson and van Dishoeck (2023)](https://doi.org/10.1051/0004-6361/202346645).
    O2 continuum sources are Brion et al. (1979) and Holland et al. (1993),
    detailed in `DATASET_AUDIT.md`; these are not new 2026 measurements.
11. [NIST PSTAR electronic stopping](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html),
    nitrogen material 007 and oxygen 008; preserved in `pstar_reference.json`.
12. Gallagher et al. (1988), JPCRD 17, 9,
    [optical channel review](https://doi.org/10.1063/1.555821).
13. Mahla and Mehnen (2025), JCP 163, 234312,
    [O2 photoionization](https://doi.org/10.1063/5.0307721) and
    [numerical supplement](https://doi.org/10.60893/figshare.jcp.30757583.v1).
14. Sun et al. (2005), Chinese Phys. 14, 1378,
    [N2 finite-q densities](https://doi.org/10.1088/1009-1963/14/7/019).
15. Lin et al. (2013), CPB 22, 023404,
    [O2 finite-q response](https://doi.org/10.1088/1674-1056/22/2/023404).
16. Lopez-Patino et al. (2016), IJMS 405, 59,
    [low-energy ion channels, not fitted](https://doi.org/10.1016/j.ijms.2016.05.014).
