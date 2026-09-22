> Historical research note, recovered 2026-09-14.
> Model names such as 'current' refer to this note's stage, not today's production model.
> Original commands and links retain their historical context; see the
> [archive guide](../../README.md) for status and portable reproduction.

# Optical and stopping audit of the matched PJG reference

This audit retains the N2/O2 coefficients in `SELECTED_PARAMETERS` as the
baseline. The additional fits are research trials, not promoted production
coefficients. No commits or pushes are part of this investigation.

## Findings

Increasing the optical objective weight improves N2 substantially without
changing the SDCS formula or its six fitted quantities. O2 does not obtain
a comparable improvement from reweighting alone. Exchanging the fitted
multiplier of PJG's `t_a` for a fit of its existing `T_s` improves O2 over
2--100 eV, but compromises the recommended SDCS comparison and does not
solve the threshold optical structure. An unguarded O2 optical fit also
increases the effective ionization energy loss above NIST electronic stopping.

The baseline has no runaway relativistic energy moment on the audited
5 keV--10 GeV proton range. Nevertheless, O2 nearly consumes the entire
PSTAR electronic-stopping budget around 1--2 MeV and peaks 0.573% above it.
That small excess is within reference/model uncertainties, but it is not
evidence of a validated ionization-plus-excitation energy budget. The former
comparison against an uncorrected Bethe estimate was insufficient to find
this issue.

## 1. Optical fit without additional independent parameters

The asymptotic logarithmic coefficient remains

\[
a(T)=\frac{K\Gamma^2}{C_B}
\left[\frac{1}{(T-T_0)^2+\Gamma^2}
-\frac{1}{(T-T_0)^2+\Gamma^2+\Lambda^2}\right],
\qquad T_0=T_s-q_t\frac{t_a}{E_e+t_b}.
\]

Here `K` includes the fitted multiplier of printed PJG K. The relativistic
optical limit uses `E_e -> m_e c^2/2`, not `E_e -> infinity`. The positive
Lorentzian difference falls as `T^-4`; changing its amplitude, widths, or
finite center does not change the hard term's `N_e C_B/(E_e T^2)` coefficient
or the fixed Bhabha bracket. The same coefficient multiplies the soft
`ln(gamma^2)` and transverse `-beta^2` terms.

The baseline six fitted quantities are J, K, Gamma, Lambda, q_t, and alpha.
The peak-swap trial fixes q_t=1 and fits T_s instead, retaining six fitted
quantities and the same analytic form. A negative fitted T_s is a Lorentzian
center below the physical integration domain, not a negative emitted energy.
The logarithm scale, printed distortion exponent, free Tmax, and Bhabha
factors remain fixed. No energy-dependent normalization is applied afterward.

### Reference distinction

- The N2 optical reference is Rudd et al. (1992), Eqs. (32)--(33), Table II:
  a parameterization derived from photoionization information. It is not
  the Rudd proton SDCS formula.
- The O2 optical reference is the coarse NCAR/GLOW partial-photoionization
  table documented in `MATCHED_PJG.md`. The channel photon energy is T+I_j;
  using the total absorption cross section at photon energy T would be wrong.
  Its bins and effective dissociation thresholds cannot establish the exact
  pointwise threshold electron spectrum.
- The proton SDCS fit comparison is against Rudd's recommended SDCS formula,
  not a new direct fit to raw measured differential points. The total fit
  and original N2 total/mean-energy measurements are as specified in
  `MATCHED_PJG.md`. These correlated references do not define a chi-squared.

The optical objective has 40 logarithmically spaced samples on 2--100 eV,
with each residual block normalized by the square root of its sample count.
Its baseline weight is 0.25, compared with weight 1 for recommended total
and SDCS blocks. Increasing that weight changes the calibration objective,
not the number of SDCS parameters or runtime operations.

| Gas / trial | Optical log-RMS, 40 fit samples | Optical ratio, 40 samples | Recommended total ratio, 31 energies | SDCS log-RMS | Maximum pair loss / PSTAR |
|---|---:|---:|---:|---:|---:|
| N2 baseline | 0.19080 | 0.8666--1.4255 | 0.9249--1.0543 | 0.12086 | 0.96717 |
| N2 optical weight 1 | 0.07892 | 0.8391--1.1755 | 0.9063--1.0528 | 0.13374 | 0.91255 |
| N2 optical weight 2 | 0.06988 | 0.8666--1.1451 | 0.9001--1.0436 | 0.14137 | 0.90195 |
| O2 baseline | 0.22054 | 0.6522--1.3866 | 0.9987--1.0885 | 0.07200 | 1.00573 |
| O2 optical weight 1 | 0.21230 | 0.6759--1.4369 | 0.9993--1.0843 | 0.07191 | 1.02069 |
| O2 peak swap, weight 2 | 0.10504 | 0.7983--1.2783 | 0.9599--1.0788 | 0.16532 | 1.04040 |
| O2 peak swap, weight 2, PSTAR ceiling | 0.11781 | 0.7883--1.2705 | 0.8836--1.0826 | 0.16763 | 1.00000 |

The ceiling is a one-sided inequality on the original PSTAR table energies:
effective-pair ionization loss must not exceed tabulated electronic stopping.
It introduces no SDCS parameter and does not guess an excitation fraction.
Touching that ceiling is not a complete physical validation. Since PSTAR is
uncertain and has contributions outside this model, it must not be fitted
as an equality to ionization or used to rescale the hard Bhabha term.

A denser 4001-point audit exposes features between the 40 fit samples:

| Gas / trial | Optical ratio over the dense 2--100 eV grid |
|---|---:|
| N2 baseline | 0.86656--1.42680 |
| N2 optical weight 2 | 0.86659--1.14655 |
| O2 baseline | 0.60542--1.39816 |
| O2 peak swap, weight 2, PSTAR ceiling | 0.74200--1.27536 |

Do not extrapolate this improvement to T=0. The N2 optical ratio there changes
from 1.80 to 1.46; the O2 ratio against the coarse proxy changes from 3.88 to
6.23. For O2 this is another reason not to promote the trial solely because
its 2--100 eV objective improved. Smooth single-center Lorentzians cannot
reproduce all channel thresholds and resonances. More accurate optical
reference data need not increase model complexity, but are required before
interpreting a threshold-fit residual as a precise physical error.

Relevant measured O2 partial-channel data include "Total and partial
photoionization cross-sections of O2 from 100 to 800 Angstrom" (1977), DOI
10.1016/0368-2048(77)85079-2, and Brion et al. (1979),
J. Electron Spectrosc. Relat. Phenom. 17, 101--119,
DOI 10.1016/0368-2048(79)85032-X. These were located, not digitized or used
in this trial. An attempted download of the NIST-hosted Gallagher et al.
(1988) compilation, DOI 10.1063/1.555821, was unavailable. Total absorption
data alone cannot replace the missing partial electron-energy mapping.

Subsequent access update: the user has supplied the Gallagher full text.
The new [source audit](DATASET_AUDIT.md) records its channel definitions,
available figures, and the outstanding numerical companion. The preceding
access statement describes the earlier trial, not current availability.

The two plotted trial parameter sets are:

| Quantity | N2 optical weight 2 | O2 peak swap + PSTAR ceiling |
|---|---:|---:|
| J (eV) | 13.6343724083 | 4.41380903583 |
| K / printed K | 1.02377674678 | 1.47184534759 |
| Gamma (eV) | 11.9158741942 | 15.0261873902 |
| Lambda (eV) | 114.330266379 | 140.294232823 |
| q_t | 0.443381243348 | 1, fixed |
| T_s (eV) | 4, fixed | -1.90726330950 |
| alpha / printed alpha | 1.11184680503 | 1.26356175078 |

Fix p=0.807/1.314, the binding logarithm scale, and zero width numerator as
in the baseline. Reallocating center freedom for N2 is less favorable than
simple reweighting: it degrades the recommended low-energy total comparison.

## 2. What stopping can be inferred from this SDCS?

Let S_j(E,T) be the model's inclusive effective-continuum contribution, in
cm^2/eV per molecule, and S=sum_j S_j. Define the cross-section moments

\[
\mathcal L_{\rm kin}(E)=\int T S(E,T)\,dT,
\qquad
\mathcal L_{\rm pair}(E)=\sum_j\int (T+I_j)S_j(E,T)\,dT.
\]

The first follows directly from the emitted-electron spectrum. The second
is the implied effective one-electron/one-ion source cost. Its binding term
uses the conditional I_j weights of the actual continuum gates, not a fixed
global mean binding energy. Integration includes the allowed bound-electron
tail through the molecular endpoint, not just free Tmax.

Inclusive electron yield alone does not uniquely determine projectile loss:
multiple ionization, inner-shell relaxation, and dissociation can give several
electrons per parent event. Assigning a binding cost to each effective pair
is therefore a model convention, not an exact channel-resolved stopping
reconstruction. Neither expression includes discrete excitation or nuclear
stopping. Neutral-thermal ion birth also does not supply the true recoil
energy of an individual ionizing event. The rigid production beam remains
rigid: this is an energy-budget diagnostic, not an energy-loss update.

The mass stopping-power conversion is

\[
\frac{S_{\rm pair}}{\rho}
\;[\mathrm{MeV\,cm^2/g}]
=10^{-6}\frac{N_A}{M_{\rm mol}}
\mathcal L_{\rm pair}\;[\mathrm{eV\,cm^2}],
\]

with N_A=6.02214076e23, M_N2=28.0134 g/mol, and M_O2=31.9988 g/mol.
Using atomic rather than molecular mass would introduce a factor-of-two
error. Density is not an extra multiplier in this mass-stopping conversion;
linear stopping is the molecular number density times the loss cross section.

### Actual NIST reference, not an uncorrected Bethe surrogate

`pstar_reference.json` preserves the original rounded seven-column default
PSTAR rows retrieved on 2026-09-07 for materials 007 (nitrogen) and 008
(oxygen). There are 132 rows per material, from 0.001 to 10000 MeV. This
audit uses the rows at and above 0.005 MeV and the electronic column,
not electronic plus nuclear stopping. The parser checks material, units,
monotonic energy, the complete range, and rounded electronic+nuclear=total.

The listed gas densities are 1.16528e-3 and 1.33151e-3 g/cm^3. NIST's
mean excitation energies are 82 and 95 eV, respectively. These are Bethe
logarithmic mean excitation energies, not PJG continuum thresholds or the
energy cost to create one effective pair.

| Proton energy | N2 pair model | N2 NIST electronic | Ratio | O2 pair model | O2 NIST electronic | Ratio |
|---|---:|---:|---:|---:|---:|---:|
| 5 keV | 102.950 | 283.0 | 0.3638 | 66.455 | 222.4 | 0.2988 |
| 100 keV | 547.280 | 759.4 | 0.7207 | 479.336 | 643.3 | 0.7451 |
| 1 MeV | 206.675 | 225.9 | 0.9149 | 215.678 | 216.1 | 0.9980 |
| 2 MeV | 128.410 | 138.8 | 0.9251 | 134.069 | 133.4 | 1.0050 |
| 10 MeV | 38.594 | 40.43 | 0.9546 | 38.492 | 39.32 | 0.9789 |
| 800 MeV | 1.98810 | 2.081 | 0.9554 | 1.93095 | 2.050 | 0.9419 |
| 10 GeV | 1.93252 | 2.062 | 0.9372 | 1.89148 | 2.040 | 0.9272 |

Stopping columns are MeV cm^2/g. Baseline maximum pair/NIST ratios are
0.967173 at 50 MeV for N2 and 1.005731 at 1.75 MeV for O2. Emitted-electron
kinetic energy alone remains below 0.708 and 0.737 of NIST, respectively.
At 800 MeV those kinetic-only values are 1.46849 and 1.48681 MeV cm^2/g.

The N2 optical-weight-2 trial reduces its largest pair/NIST ratio to 0.90195;
its 800 MeV value is 1.85124 MeV cm^2/g. O2 optical reweighting alone makes
the budget worse. The O2 constrained trial touches NIST at 2.25 MeV and
has an 800 MeV value of 1.94250 MeV cm^2/g. It has not established a
complete excitation budget merely by eliminating the excess.

NIST includes shell, Barkas, Bloch, and density corrections at high energy
and data-based low-energy stopping. It quotes uncertainties around 2--5%
at 1 MeV and 5--10% at 100 keV. Thus O2's 0.573% excess is not a
statistically resolved failure by itself; the missing-process budget and
effective-pair convention remain the more important concerns. N2 also
requires excitation accounting before treating its small high-energy gap
as adequate. No universal fixed 80--90% ionization fraction is imposed.

### Nonrelativistic Rudd comparison

Integrating the original recommended Rudd orbital spectra with their own
nominal I_j gives the following pair-loss/NIST ratios:

| Proton energy | N2 Rudd | N2 matched PJG | O2 Rudd | O2 matched PJG |
|---|---:|---:|---:|---:|
| 100 keV | 0.7032 | 0.7207 | 0.7597 | 0.7451 |
| 1 MeV | 0.8699 | 0.9149 | 0.9537 | 0.9980 |
| 2 MeV | 0.8899 | 0.9251 | 0.9603 | 1.0050 |

O2's differential reference above 300 keV and N2's above 1700 keV are
extrapolations, not additional measured spectra. This diagnostic is also
correlated with the SDCS fitting reference and retains the inclusive/nominal
binding caveat. It shows that the O2 MeV budget issue can arise from a
modest energy-weighted spectral difference, without a catastrophic
relativistic tail. It does not diagnose an unspecified relativistic Rudd
variant or justify extending the nonrelativistic Rudd formula to GeV energy.

## 3. Additional physical and numerical checks

### Free Bhabha second moment

The electron kinetic-energy second moment strongly weights hard secondaries.
For the complete free kernel its integral is analytic:

\[
\mathcal M_{2,B}=\frac{N_eC_B}{E_e}
\int_0^{T_m}\left[1-\beta^2\frac{T}{T_m}
+\frac{T^2}{2(E+M)^2}\right]dT
=\frac{N_eC_B}{E_e}
\left[T_m(1-\beta^2/2)+\frac{T_m^3}{6(E+M)^2}\right].
\]

This is the free kinetic-transfer contribution to straggling, not a
complete molecular energy-loss variance including all binding channels.
For inclusive emission, single-electron second moments also omit the
correlations among electrons emitted in the same parent event.

| Proton energy | N2 model M2 / free M2 | O2 model M2 / free M2 |
|---|---:|---:|
| 1 MeV | 0.86488 | 0.86911 |
| 10 MeV | 0.97039 | 0.97224 |
| 800 MeV | 0.997743 | 0.997770 |
| 10 GeV | 0.999898 | 0.999891 |

No equality is expected at low projectile energy, where binding and
distortion are important. At 800 MeV the original baseline pointwise
SDCS/Bhabha errors in the 10 keV--2 MeV secondary window remain below
0.084% (N2) and 0.139% (O2). The plotted optical trials retain that limit;
their corresponding maximum errors are 0.084% and 0.015%.

### Above-free tail does not conceal a large relativistic energy loss

At 800 MeV the part above free Tmax contributes 0.001303% (N2) and
0.001810% (O2) of the emitted-electron kinetic-energy moment. At 10 GeV
those fractions fall to 0.00001465% and 0.00001929%. At 5 keV the fractions
are much larger, 49.1% and 58.4%, because a free stationary-electron endpoint
is not a valid molecular cutoff there. Both free Tmax and the molecular
energy-conservation support remain unchanged.

### Soft strength and threshold limitations

Direct integration and an independent analytic Lorentzian primitive agree
for the effective continuum strength
`F_eff = integral (T+mean(I_j))*a(T) dT`: 10.8971 for baseline N2 and
12.4718 for baseline O2, compared with electron counts 14 and 16. The
plotted trials give 9.5168 and 12.5811. This is a coarse magnitude check,
not a rigorous sum-rule bound on inclusive decay yield. The deficit is not automatically discrete excitation;
inner-shell response and inclusive relaxation also matter. Neither a
complete Thomas--Reiche--Kuhn sum nor NIST's mean excitation energy has
been reconstructed from this effective spectrum.

Measured nonionizing optical strength is not zero: Huebner et al. (1975)
report an O2 oscillator-strength sum 0.198 below its first ionization
potential. Optical strengths alone do not fix the full finite-momentum,
finite-projectile-energy excitation stopping. Likewise, a gas W-value
includes the subsequent electron cascade and must not be compared directly
to the primary moment `L_pair / sigma`.

The existing tests also check finite, nonnegative T=0 spectra, vanishing
subthreshold yield, exact molecular support, continuity of the bound-tail
join, the relativistic dipole logarithmic coefficient, and the free hard
coefficient. Finite and positive is not the same as an exact threshold
optical spectrum or a complete relativistic PWBA calculation.

### Numerical audit

The stopping quadrature resolves the low-energy peak in log(1+T), both
sides of free Tmax on separate intervals, and the distant molecular
endpoint. Each interval retains its exact Jacobian. Keeping a logarithmic
soft interval even when the edge width exceeds Tmax avoids under-resolving
the peak with a long linear interval.

At 5 keV, 100 keV, 1.75 MeV, 800 MeV, and 10 GeV, 1025 versus 4097 points
per segment differ by less than 9.3e-9 in total/first/binding/second moments
and 5e-8 in the above-free moments. Independent adaptive quadrature agrees
to the 1e-7 test criterion. Other checks verify molecular/atomic mass units,
the independently integrated Bhabha second moment, positive moment bounds,
`M2*sigma >= M1^2`, and the available projectile-energy budget.

These checks detect numerical and asymptotic errors; they do not establish
experimental accuracy below 5 keV, density-effect accuracy at arbitrary
gas density, or reliability beyond the 10 GeV PSTAR comparison. The model
has no explicit dielectric density effect. It must not be certified for
arbitrarily high energy or condensed matter from this audit.

## 4. Reproduction and implementation status

```sh
python Tools/Algorithms/ProtonImpactIonization/stopping_pjg_matched.py
python Tools/Algorithms/ProtonImpactIonization/study_optical_stopping.py \
  --o2-optical /path/to/ephoto_xo2.dat \
  --output tmp/pjg-optical-stopping-trials.json \
  --figures tmp/pjg-stopping-figures
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
ctest --test-dir build/proton-impact-physics --output-on-failure
```

The stored PSTAR rows make stopping checks independent of live network
availability. `--n2-pstar` and `--o2-pstar` accept fresh text-only HTML
tables; request `prog=PSTAR`, `matno=007` or `008`, `ShowDefault=on`, and
empty `Energies` at the POST endpoint below. No authentication is required.
`--read-results` with `--figures` regenerates figures without repeating fits.

The 40 Python checks and two C++ host tests pass. All production tables,
source sampling, beam motion, and neutral-thermal ion birth are unchanged.
The trial center override is confined to host-side research code. Changing
the selected calibration would not require an additional GPU lookup or
shape function, but no CUDA/HIP/SYCL execution or timing claim follows from
these host calculations.

## Sources

- NIST PSTAR text interface:
  https://physics.nist.gov/PhysRefData/Star/Text/PSTAR-t.html.
- NIST table endpoint:
  https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl.
- NIST methods, uncertainties, and corrections:
  https://physics.nist.gov/PhysRefData/Star/Text/programs.html.
- NIST stopping definitions:
  https://physics.nist.gov/PhysRefData/Star/Text/appendix.html.
- Material data:
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=007 and
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=008.
- Rudd et al., Rev. Mod. Phys. 64, 441 (1992), especially Sec. V F and
  Table II: https://doi.org/10.1103/RevModPhys.64.441.
- O2 partial optical input:
  https://github.com/NCAR/GLOW/blob/6cb880d0e112810dbe9a79203f0183f9b07e3473/data/ephoto_xo2.dat.
- Huebner et al. (1975), measured O2 excitation oscillator strengths:
  https://www.nist.gov/publications/apparent-oscillator-strengths-molecular-oxygen-derived-electron-energy-loss-0.
