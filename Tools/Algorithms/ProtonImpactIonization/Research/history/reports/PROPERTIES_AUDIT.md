> Historical research note, recovered 2026-09-14.
> Model names such as 'current' refer to this note's stage, not today's production model.
> Original commands and links retain their historical context; see the
> [archive guide](../../README.md) for status and portable reproduction.

# Energy-distribution and stopping audit of the PJG source trials

Date: 2026-09-09. Both N2 and O2 are included. These are research results,
not promoted coefficients, a production model release, or GPU validation.
No commits or pushes were made for this study.

Unless explicitly stated otherwise, **source trial** below means
`six_optical_weight_1` from `source-refits.json`, not `SELECTED_PARAMETERS`
in the Python module. The latter remains the earlier matched baseline.
The new mean-energy reweightings are separate sensitivity trials.

## Principal findings

1. The source trials have positive, numerically converged energy moments and
   recover the free hard Bhabha spectrum substantially better than the
   printed, algebraically corrected, or earlier J/amplitude-refitted PJG.
   This does not validate their bound-electron tail or full molecular PWBA
   response from first principles.
2. Total electron-production cross sections are close to Rudd's 1983
   measurement-derived fit: model/reference ranges are 0.890--1.042 for N2
   and 0.874--1.001 for O2. The table is the authors' fit, not raw data.
3. N2's original measured mean energies expose a remaining shape problem:
   +28.1% at 5 keV and approximately -11.8% at 50--70 keV. Reweighting these
   measurements, with the same six SDCS coefficients, improves the means
   but worsens other comparisons. There is no simultaneously dominant fit.
4. O2 means over the digitized Cheng measurement windows agree within
   +1.0%, -1.8%, and -6.3% at 7.5, 50, and 150 keV. These are not full-spectrum
   mean-energy measurements, and the digitized points have shared systematics.
5. Source-trial effective-pair stopping does not exceed NIST electronic
   stopping on the audited 5 keV--10 GeV grid. O2 nevertheless reaches
   98.53% around 1 MeV. The remaining budget is small, not a proven failure
   or a validated excitation budget. Increasing the optical weight to 2
   lowers that maximum to 96.14%, without adding an SDCS parameter.
6. The current practical IAA angular closure fails a direct qualitative
   comparison with Cheng Fig. 5: it excludes measured backward electrons.
   This is a dynamics/closure failure, not a molecular-kinematics restriction.

## Definitions and normalization

Let T denote the ejected electron's kinetic energy, E the proton kinetic
energy, and S(E,T) the inclusive electron-production SDCS per molecule.
Every integral includes the bound-electron tail through the molecular
endpoint; free Tmax is not used as the integration cutoff.

\[
M_n(E)=\int_0^{T_{\rm mol}}T^n S(E,T)\,dT,\qquad
\sigma_e=M_0,\qquad \langle T\rangle=M_1/M_0,
\]

\[
\operatorname{Var}(T)=M_2/M_0-(M_1/M_0)^2.
\]

The kinetic-energy moment M1 is determined by the SDCS. For the explicitly
effective one-electron/one-ion source, the additional binding cost is

\[
B(E)=\int I_{\rm eff}(E,T)S(E,T)\,dT,\qquad
I_{\rm eff}=\frac{\sum_j I_j G_j(E,T)}{\sum_j G_j(E,T)},\qquad
L_{\rm pair}=M_1+B.
\]

G_j are the actual channel gates, including their endpoint factors. A fixed
global average binding energy is not substituted at low proton energy.
The resulting effective binding costs are about 19--20 eV, not the NIST
logarithmic mean excitation energies of 82 and 95 eV.

Inclusive yield does not uniquely identify parent events: multiple ionization,
fragmentation, and relaxation can create several electrons per event. Thus
L_pair is the effective source's implied cost, not an exact reconstruction
of measured projectile energy loss. It omits discrete excitation, capture,
true ion recoil and nuclear stopping. Neutral-thermal ion birth does not
reconstruct the physical recoil of an individual ionizing collision.

Mass stopping is obtained using molecular, not atomic, molar mass:

\[
S_{\rm pair}/\rho=10^{-6}(N_A/M_{\rm mol})L_{\rm pair}
\quad[\mathrm{MeV\,cm^2/g}],
\]

where cross-section moments are in eV cm2, M_N2=28.0134 g/mol,
M_O2=31.9988 g/mol, and N_A=6.02214076e23/mol. No extra gas-density multiplier
belongs in this mass-stopping conversion. For molecular density n, the
collision rate is n*v*sigma_e and the primary mean free path is 1/(n*sigma_e).

Neither L_pair/sigma_e nor NIST electronic stopping divided by sigma_e is
the usual gas W value. That W value counts the eventual ion-pair cascade,
including secondary-electron ionizations. These quantities count primary
emitted electrons. Rudd (1992), Sec. V F explicitly distinguishes them.

The production beam remains externally prescribed/ballistic. None of these
diagnostics subtracts energy or momentum from it. A full-system energy
conservation assertion would be inappropriate for that rigid-source choice.

## 1. Total electron-production yield and full-spectrum mean

Cross sections in this table are in 1e-16 cm2 per molecule; means are in eV.
The reference columns are Rudd et al. (1983), Table V, computed by those
authors from their Eq. (17). Do not label them individual raw measurements.
The experimental range ends at 4000 keV; their printed 5000-keV extrapolation
is deliberately excluded. The mean columns here are model predictions.

| Proton keV | N2 Rudd fit | N2 source | O2 Rudd fit | O2 source | N2 mean T | O2 mean T |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 2.110 | 1.919 | 1.670 | 1.460 | 6.94 | 8.72 |
| 7 | 2.720 | 2.422 | 2.260 | 2.000 | 7.82 | 9.89 |
| 10 | 3.450 | 3.076 | 3.010 | 2.712 | 8.89 | 11.32 |
| 15 | 4.320 | 3.962 | 3.960 | 3.657 | 10.37 | 13.24 |
| 20 | 4.900 | 4.642 | 4.630 | 4.355 | 11.67 | 14.85 |
| 30 | 5.570 | 5.529 | 5.440 | 5.237 | 13.97 | 17.59 |
| 50 | 5.960 | 6.157 | 5.980 | 5.899 | 18.00 | 22.06 |
| 70 | 5.860 | 6.107 | 5.950 | 5.931 | 21.51 | 25.75 |
| 100 | 5.470 | 5.664 | 5.600 | 5.606 | 25.94 | 30.27 |
| 150 | 4.780 | 4.841 | 4.930 | 4.904 | 31.45 | 35.90 |
| 200 | 4.190 | 4.175 | 4.360 | 4.296 | 35.29 | 39.91 |
| 300 | 3.360 | 3.266 | 3.530 | 3.426 | 40.15 | 45.14 |
| 500 | 2.420 | 2.302 | 2.580 | 2.460 | 45.07 | 50.57 |
| 700 | 1.910 | 1.799 | 2.050 | 1.940 | 47.65 | 53.45 |
| 1000 | 1.460 | 1.371 | 1.580 | 1.489 | 49.94 | 56.02 |
| 1500 | 1.060 | 0.999 | 1.170 | 1.089 | 52.11 | 58.50 |
| 2000 | 0.842 | 0.794 | 0.933 | 0.867 | 53.43 | 60.03 |
| 3000 | 0.603 | 0.571 | 0.675 | 0.624 | 55.04 | 61.94 |

Rudd's quoted uncertainty of the fitted totals is 25/20/15/10/8% at
5/10/25/100/500 keV, and 8% above 500 keV. These are not independent Gaussian
errors. Against Crooks and Rudd (1971) Table I measurements at 50--300 keV,
the source ratios span 1.008--1.113 (N2) and 0.985--1.139 (O2). That experiment
has approximately 17% shared normalization uncertainty.

The difference between source references is material. At 50 keV, N2 has
5.33 in Rudd (1979), 5.53 in Crooks (1971), and 5.96 in the 1983 fit. The
source prediction is 6.157. A close match to one compilation is not a
simultaneous percent-level match to all experiments.

### Original measured N2 mean electron energies

The following means were checked visually in Rudd (1979), Table I, p. 791.
The nitrogen column does not include the hydrogen-only 100-keV entry.
The 1979 absolute scale was tied to Crooks' 50-keV data; do not count that
scale twice. The paper does not provide pointwise covariance for these means.

| Proton keV | Measured mean, eV | Previous matched | Source trial | Source error |
|---:|---:|---:|---:|---:|
| 5 | 5.42 | 6.403 | 6.942 | +28.1% |
| 7 | 6.67 | 7.319 | 7.816 | +17.2% |
| 10 | 7.74 | 8.474 | 8.891 | +14.9% |
| 15 | 9.99 | 10.107 | 10.374 | +3.8% |
| 20 | 11.90 | 11.556 | 11.665 | -2.0% |
| 30 | 15.20 | 14.180 | 13.968 | -8.1% |
| 50 | 20.40 | 18.792 | 18.001 | -11.8% |
| 70 | 24.40 | 22.746 | 21.514 | -11.8% |

### Same-six-parameter mean-energy reweighting

The baseline source objective assigns the measured-mean block weight 0.5.
Increasing it changes the fitting priorities, not the SDCS formula or number
of coefficients. Optical weight remains 1 in the following rows. Weights
are not inverse standard errors, and objective norms across weights are
not statistically comparable.

| Mean-block weight | Model/measured mean range | Mean log-RMS | Model/Rudd83 total range | Optical W-response range, 25--100 eV | Max pair/NIST |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.882--1.281 | 0.1351 | 0.890--1.042 | 0.869--1.329 | 0.9076 |
| 1 | 0.897--1.212 | 0.1063 | 0.875--1.060 | 0.852--1.344 | 0.9086 |
| 2 | 0.928--1.140 | 0.0695 | 0.834--1.096 | 0.806--1.362 | 0.9121 |
| 4 | 0.948--1.089 | 0.0435 | 0.786--1.143 | 0.770--1.375 | 0.9155 |

At weight 4 the means at 5/50/70 keV become 5.901/19.901/24.551 eV. However,
the Crooks total ratio reaches 1.232 and the Rudd-formula SDCS log-RMS
increases to 0.2371, from 0.1920 in the source trial. This does not establish
a globally better fit. Optical weight 0.5 with mean weight 2 gives a
0.931--1.123 mean ratio, but expands the optical range to 0.742--1.505.
All four additional fits converged in 14--21 function evaluations.

Changing J alone cannot repair a mean-energy discrepancy when it is just
a common multiplicative factor; in the matched formulation the hard fading
factor gives it some shape influence, but the soft-center/width and edge
trade-off remains. None of these trials changes free Tmax or Bhabha factors.

## 2. Measured O2 SDCS and restricted means

Cheng, Rudd and Hsu (1989), Fig. 1 supplies three spectra after adjustment
to recommended total cross sections. Thirty open-circle points were
digitized with preserved pixel coordinates. The absolute normalization is
therefore not independent of total-cross-section recommendations.

For each actual digitized window [a,b], define

\[
\langle T\rangle_{[a,b]}=
\frac{\int_a^b T S(E,T)\,dT}{\int_a^b S(E,T)\,dT}.
\]

Positive measured points are interpolated piecewise as power laws in T;
these integrals are evaluated analytically on each cell, with the logarithmic
limit retained. Models are integrated independently over exactly the same
window. No values below a or above b are inferred from measurements.
The mean cancels a constant normalization shift, but not energy-dependent
detector or digitization errors. Sparse data do not resolve every feature.

| Proton keV | Digitized window, eV | Data-window mean, eV | Source-window mean | Previous matched | Rudd SDCS formula |
|---:|---:|---:|---:|---:|---:|
| 7.5 | 3.041--40.313 | 11.779 | 11.892 | 11.263 | 11.209 |
| 50 | 2.947--294.705 | 25.890 | 25.417 | 24.685 | 24.967 |
| 150 | 2.901--390.694 | 42.494 | 39.799 | 39.558 | 39.701 |

For example, the source's **full-spectrum** mean at 50 keV is only 22.064 eV:
the restricted comparison excludes many electrons below 2.947 eV. The two
numbers are not inconsistent. No full-spectrum O2 mean table has been
verified in the currently accessible Cheng/Crooks source set.

At the digitized points, source/data SDCS ratios span 1.054--1.209,
0.752--1.382, and 0.749--1.210 at 7.5/50/150 keV. Shape log-RMS values after
removing one constant log normalization are 0.0443/0.1664/0.1390. Thus good
restricted means do not imply an equally good pointwise spectrum.

The 50-keV last measured point is near 295 eV, above the free Tmax of
108.8 eV. The Rudd-formula/data ratio there is 0.438; shifting the digitization
by two pixels in each coordinate gives 0.330--0.580. That is a reading
sensitivity, not an experimental confidence interval. It illustrates why
agreement with Rudd's formula must not be presented as agreement with data.

## 3. Stopping compared with NIST

All entries below are mass stopping powers in MeV cm2/g. The comparison
uses the NIST electronic column, not electronic plus nuclear stopping.
The stored tables are gas nitrogen/oxygen, materials 007/008, 132 rows each,
retrieved 2026-09-07. The complete 5-keV--10-GeV comparison is in the JSON.

| Gas | Proton energy | Electron kinetic | Effective binding | Pair sum | NIST electronic | Pair/NIST |
|---|---:|---:|---:|---:|---:|---:|
| N2 | 5 keV | 28.639 | 81.023 | 109.662 | 283.0 | 0.3875 |
| N2 | 10 keV | 58.786 | 127.270 | 186.056 | 400.2 | 0.4649 |
| N2 | 100 keV | 315.887 | 234.732 | 550.619 | 759.4 | 0.7251 |
| N2 | 1 MeV | 147.252 | 57.762 | 205.014 | 225.9 | 0.9075 |
| N2 | 2 MeV | 91.148 | 33.422 | 124.570 | 138.8 | 0.8975 |
| N2 | 10 MeV | 26.414 | 8.820 | 35.234 | 40.43 | 0.8715 |
| N2 | 800 MeV | 1.356 | 0.403 | 1.759 | 2.081 | 0.8452 |
| N2 | 10 GeV | 1.340 | 0.378 | 1.718 | 2.062 | 0.8330 |
| O2 | 5 keV | 23.952 | 55.169 | 79.120 | 222.4 | 0.3558 |
| O2 | 10 keV | 57.789 | 101.058 | 158.846 | 314.5 | 0.5051 |
| O2 | 100 keV | 319.369 | 207.225 | 526.595 | 643.3 | 0.8186 |
| O2 | 1 MeV | 157.021 | 55.895 | 212.916 | 216.1 | 0.9853 |
| O2 | 2 MeV | 97.938 | 32.541 | 130.479 | 133.4 | 0.9781 |
| O2 | 10 MeV | 28.454 | 8.541 | 36.995 | 39.32 | 0.9409 |
| O2 | 800 MeV | 1.461 | 0.388 | 1.849 | 2.050 | 0.9022 |
| O2 | 10 GeV | 1.446 | 0.365 | 1.811 | 2.040 | 0.8878 |

For comparison, nominal nonrelativistic Rudd orbital-binding loss divided
by PSTAR is 0.334/0.703/0.870 for N2 and 0.286/0.760/0.954 for O2 at
5 keV/100 keV/1 MeV. This calculation uses Rudd's own orbital thresholds.
It is a correlated model comparison, not independent experimental stopping.
The NR formula is not extended past 4 MeV here; O2 SDCS above 300 keV and
N2 above 1700 keV are already reference-curve extrapolations.

N2's maximum source pair/NIST ratio is 0.90760 at 0.85 MeV. O2's maximum
is 0.98527 at 1 MeV, leaving 3.184 MeV cm2/g, or 1.47%, below NIST. The
optical-weight-2 O2 trial has pair stopping 207.760 at 1 MeV, a 0.96141 ratio.
Its mean T is slightly **higher**, 57.243 versus 56.025 eV: stopping depends
on the yield and binding term as well as mean electron energy.

One must not force ionization-only loss to equal PSTAR or subtract a guessed
universal excitation fraction. The low-energy PSTAR treatment includes
experimental fitting and charge-changing effects outside a permanently
bare, ionization-only source. Its high-energy calculation includes shell,
Barkas, Bloch and density corrections. Reference uncertainty and the
inclusive-pair convention preclude interpreting a small residual as a
precise excitation measurement. The molecular optical sums also do not
justify a universal 10--20% discrete-excitation fraction for both gases.

The model has no explicit density-effect response. Results at hundreds of
MeV and GeV are asymptotic consistency/extrapolation tests, not a validation
at arbitrary material density. Remaining below NIST is necessary as a
coarse budget check, but is not sufficient to establish correct stopping.
No ionization-only range is reported as a NIST CSDA range.

## 4. Distribution widths, tails and sampling noise

The following values describe the entire source-trial distribution.
Quantiles are obtained from a positive mapped CDF, not from an energy cutoff.

| Gas | Proton E | Mean T, eV | Median, eV | 90th, eV | 99th, eV | Std. dev., eV |
|---|---:|---:|---:|---:|---:|---:|
| N2 | 5 keV | 6.942 | 5.133 | 14.979 | 31.784 | 6.729 |
| N2 | 100 keV | 25.942 | 12.364 | 69.505 | 172.866 | 35.882 |
| N2 | 1 MeV | 49.944 | 12.705 | 99.766 | 748.215 | 143.187 |
| N2 | 800 MeV | 65.913 | 11.557 | 58.064 | 530.232 | 2899.564 |
| N2 | 10 GeV | 69.433 | 11.310 | 52.417 | 409.250 | 16575.127 |
| O2 | 5 keV | 8.719 | 6.002 | 19.711 | 42.912 | 9.092 |
| O2 | 100 keV | 30.272 | 15.842 | 77.767 | 189.385 | 39.387 |
| O2 | 1 MeV | 56.025 | 17.716 | 111.620 | 769.520 | 146.309 |
| O2 | 800 MeV | 75.038 | 16.855 | 75.059 | 568.388 | 2981.011 |
| O2 | 10 GeV | 78.933 | 16.566 | 69.598 | 444.173 | 17014.900 |

These distributions are not narrow or approximately Gaussian. The median
can stay nearly fixed while a very rare high-energy tail increases M1 and,
especially, M2. An unchanged or decreasing 99th percentile does not imply
that the highest-energy electrons are unimportant.

The fractions **above free Tmax** are:

| Gas | Proton E | Electron count | Kinetic-energy moment | Second moment |
|---|---:|---:|---:|---:|
| N2 | 5 keV | 19.523% | 49.294% | 76.066% |
| O2 | 5 keV | 27.852% | 63.668% | 87.111% |
| N2 | 100 keV | 0.2939% | 2.8689% | 9.8092% |
| O2 | 100 keV | 0.5285% | 4.5631% | 15.0538% |
| N2 | 1 MeV | 0.00588% | 0.2716% | 1.3644% |
| O2 | 1 MeV | 0.01015% | 0.4237% | 2.2724% |

At 5 keV free Tmax=10.8805 eV. Truncating there removes about half of N2's
and nearly two thirds of O2's predicted electron kinetic energy, not merely
a small number of tail particles. These tail probabilities remain empirical
model predictions; the molecular endpoint permits them but does not determine
their magnitudes. At high proton energy most large energy transfers are
below free Tmax; a small above-free fraction does not mean a small hard tail.

All low-energy electrons remain represented. For example, at 100 keV the
fractions below 2 eV are 9.07% (N2) and 8.21% (O2); the fractions below the
first ionization threshold are 57.63% and 41.59%. Being unable to cause
further ionization is not a reason to remove a primary electron kinetically.

The explicit hard component accounts for 17.7%/15.8% of the N2/O2 electron
count at 100 keV, but 48.8%/40.3% of electron kinetic energy. At 800 MeV
these fractions are 7.31%/6.15% of count and 74.4%/67.5% of kinetic energy.
Component fractions depend on the model's matching convention; they are
not experimentally separated soft and hard event counts.

For N independent equal-weight emitted-energy samples,

\[
\frac{\operatorname{Var}(\widehat{\langle T\rangle})}{\langle T\rangle^2}
=\frac{1}{N}\left[\frac{M_2M_0}{M_1^2}-1\right].
\]

The bracket is 1.913/1.693 at 100 keV, 8.219/6.820 at 1 MeV,
1935/1578 at 800 MeV, and 56988/46467 at 10 GeV for N2/O2.
Naive iid sampling would consequently require about 19.4/15.8 million
samples for 1% relative error of mean energy at 800 MeV. This concerns a
Monte Carlo estimator, not the physical number of electrons in a beam.

If the event count itself is Poisson, the variance of deposited electron
kinetic energy over path length L is n*L*M2, not n*L*M0*Var(T). For the
relative variance at expected count lambda the bracket becomes
M2*M0/M1^2, without the subtraction of one. This distinction is tested.

These results favor testing stratified/importance-aware energy sampling
with correct particle weights and parent conditioning. They do not justify
discarding energetic electrons or changing their physical cross section.
An unbiased mean-energy estimator alone also does not establish low noise
for charge/current deposition or other PIC observables.

## 5. Hard-tail comparison with the old models

At E=800 MeV, each entry below is SDCS divided by the **complete pointlike
spin-1/2 free Bhabha SDCS**, evaluated at the same secondary energy. These
energies are well below free Tmax=2.481 MeV. The ratios are not renormalized.
"Printed" always includes the 1977 erratum.

| Gas / model | T=1 keV | T=10 keV | T=100 keV |
|---|---:|---:|---:|
| N2 printed | 0.8451 | 0.8597 | 0.7472 |
| N2 corrected, no refit | 0.8463 | 0.8717 | 0.8712 |
| N2 earlier J refit | 0.8451 | 0.8705 | 0.8700 |
| N2 source trial | 1.0091 | 1.00081 | 1.00008 |
| O2 printed | 1.1446 | 1.1658 | 1.0633 |
| O2 corrected, no refit | 1.1457 | 1.1778 | 1.1873 |
| O2 earlier J/amplitude refit | 0.8012 | 0.8229 | 0.8216 |
| O2 source trial | 1.0244 | 1.00138 | 1.00013 |

The cancellation of the soft term's T^-2 tail and restoration of its
coefficient in the fixed hard term matter: fitting only J or a common
amplitude does not preserve that coefficient. The source Lorentzian
difference falls as T^-4 and its hard part has the required T^-2 coefficient.

With E_e=m_e*c^2*beta^2/2, C_B=4*pi*a_0^2*Ry^2 and M=m_p*c^2,
the independently integrated free hard second moment is

\[
M_{2,\rm free}=\frac{N_e C_B}{E_e}
\left[T_{\max}(1-\beta^2/2)
+\frac{T_{\max}^3}{6(E+M)^2}\right].
\]

At 10 GeV the source/free ratios are 0.999895 for N2 and 0.999883 for O2.
This is a useful high-tail straggling check and shows no artificial exploding
second moment over the tested range. It is not an experimental validation
of all straggling processes.

Legacy spectra are integrated with their signed values, never silently
clipped positive. A negative SDCS makes the corresponding result an invalid
probability distribution even if its integrated total or mean looks plausible.
For example, the sampled printed N2 curve becomes negative at some secondary
energies for E=100 keV. The JSON flags such sampled violations; an unflagged
finite grid is not a proof of global positivity.

## 6. Angular closure: failed experimental-support check

Cheng Fig. 5, p. 3602, shows nonzero O2 DDCS at backward angles for
50-keV protons and T=10/50/100/150 eV. Its solid curves are the present
measurements; the dashed/dotted curves refer to earlier experiments.

The current practical IAA closure samples uniformly in a finite cosine
interval. After allowing every retained binding channel, its largest
possible angles are:

| T, eV | N2 maximum angle | O2 maximum angle | O2 backward probability |
|---:|---:|---:|---:|
| 10 | 128.18 deg | 127.11 deg | 0.3408 |
| 50 | 85.25 deg | 83.77 deg | 0 |
| 100 | 57.59 deg | 56.00 deg | 0 |
| 150 | 46.83 deg | 45.32 deg | 0 |

Molecular kinematics allows the full solid angle at every listed point for
every retained channel. The closure's hard angular edges are therefore
not required by energy-momentum conservation. Even at T=10 eV its maximum
angle excludes part of the observed far-backward range.

Passing normalization, no-endpoint-point-mass, and kinematic tests cannot
overrule this source evidence. The current closure must not be promoted as
a validated proton DDCS model. A next candidate should retain a forward
binary feature with nonzero molecular-scattering wings, conserve the chosen
SDCS normalization, and be compared with Cheng and original angular data.
A guessed isotropic mixture coefficient is not introduced by this audit.

## 7. Numerical verification and implementation status

The moment integration resolves the soft peak and both sides of the free
edge using four mapped segments, including a separate above-free integral.
Simpson quadrature checks M0, M1, binding and M2. A separate positive
trapezoidal CDF provides quantiles. Its normalization differs from the
4097-point Simpson total by at most 3.09e-6 (N2) and 2.11e-6 (O2) over the
25 reported energies; this is recorded rather than conflated with moment
accuracy. Independent adaptive quadrature checks the reported quantiles.

Existing source-trial 1025-to-4097 moment convergence errors are below
1.77e-8 (N2) and 6.73e-9 (O2). The new 2049-to-8193 check changes moments by
at most 1.31e-9 and quantiles by at most 3.47e-5 relative, over five energies
from 5 keV through 10 GeV. Quantiles are not assigned the smaller moment
integration error. Tests verify positive density, monotone CDF,
M2*M0 >= M1^2, M2 <= T_mol*M1, binding bounds, molecular energy support,
analytic power-law window integrals and both sampling-variance identities.

The complete host Python suite passes 58 tests, including six added here.
The two independent C++ host tests also pass. Their kinematic and Bhabha
trace tests use independent expressions, not copies of the final SDCS fit.
The existing host quantile microbenchmark gives approximately 1.03 ns for
float and 0.81 ns for double for shifted quantile generation on this run.
These are not full collision timings, not SDCS-evaluation timings, and not
CUDA/HIP/SYCL benchmarks. No accelerator speedup is claimed.

The production tables, source implementation and coefficients are unchanged.
Remaining physics work includes the angular support failure, the N2
mean/total/optical trade-off, a channel-consistent O2 energy budget, and a
target-response/GOS calibration. A one-dimensional empirical SDCS does not
uniquely specify the full generalized oscillator strength F(q,W).

## Reproduction and sources

Use NumPy/SciPy/Matplotlib in the `warpx-cpu-mpich-dev` environment. Source
extraction and source-refit commands are recorded in [DATASET_AUDIT.md](DATASET_AUDIT.md).

```sh
python Tools/Algorithms/ProtonImpactIonization/audit_pjg_properties.py \
  --fits tmp/pdfs/dataset-refresh-20260909/source-refits.json \
  --output tmp/pdfs/dataset-refresh-20260909/properties.json \
  --figures tmp/pjg-properties-figures
python Tools/Algorithms/ProtonImpactIonization/study_mean_energy.py \
  --fits tmp/pdfs/dataset-refresh-20260909/source-refits.json \
  --nifs tmp/pdfs/dataset-refresh-20260909/nifs-optical.json \
  --output tmp/pdfs/dataset-refresh-20260909/mean-sensitivity.json \
  --figures tmp/pjg-properties-figures
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
ctest --test-dir build/proton-impact-physics --output-on-failure
```

The source-refit digest is preserved in `properties.json`; source-table
digests and the full precision coefficients are in its input JSON. The
additional sensitivity coefficients are preserved in `mean-sensitivity.json`.
Raw source PDFs and copyrighted numerical supplements are not copied into
the WarpX source tree by these scripts.

- PJG (1976), model and tables: [10.1063/1.432812](https://doi.org/10.1063/1.432812).
- PJG erratum (1977): [10.1063/1.323427](https://doi.org/10.1063/1.323427).
- Crooks and Rudd (1971), Table I totals and measured spectra:
  [10.1103/PhysRevA.3.1628](https://doi.org/10.1103/PhysRevA.3.1628).
- Rudd (1979), Table I N2 totals/means:
  [10.1103/PhysRevA.20.787](https://doi.org/10.1103/PhysRevA.20.787).
- Rudd et al. (1983), Table V measurement-derived total fits:
  [10.1103/PhysRevA.28.3244](https://doi.org/10.1103/PhysRevA.28.3244).
- Cheng, Rudd and Hsu (1989), Figs. 1 and 5:
  [10.1103/PhysRevA.40.3599](https://doi.org/10.1103/PhysRevA.40.3599).
- Rudd et al. (1992), SDCS formula and stopping/W-value discussion:
  [10.1103/RevModPhys.64.441](https://doi.org/10.1103/RevModPhys.64.441).
- NIST [PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) and
  [methods and uncertainties](https://physics.nist.gov/PhysRefData/Star/Text/programs.html).
- Optical anchors: [NIFS-DATA-109](https://nifs-repository.repo.nii.ac.jp/records/11706)
  for the N2 evaluated response, and the [Leiden O2 numerical photoionization table](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/O2/O2.txt)
  using the experimental/evaluated lineage detailed in DATASET_AUDIT.md.
  Photon energy W is not secondary energy T. Mahla (2025) and the measured
  GOS sources are audited comparisons, not new fitted inputs in this turn.
