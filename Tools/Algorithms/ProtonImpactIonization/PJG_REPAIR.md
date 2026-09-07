# A smaller, PJG-preserving repair

## Status and scope

This is a research checkpoint, not the production SDCS. The current direction
is to retain PJG's empirical structure and repair its limiting behavior,
rather than replace it by a new generalized-oscillator-strength model.
`pjg_repair.py` contains two deliberately separate experiments:

- A relativistic, free-domain core for testing the Lorentzian decomposition,
  positivity, and the hard and Bethe limits. It is not a molecular cutoff model.
- A nonrelativistic data-fit diagnostic, including an empirical bound tail,
  restricted to 5--4000 keV. It cannot be extrapolated to relativistic energies.

Neither experiment is used to generate production tables. The rigid beam,
inclusive effective pair source, thermal ion source, and all-kinetic electron
representation remain unchanged.

## Which constraints should be enforced?

Exact free two-body kinematics and the complete pointlike spin-1/2 Bhabha
factor are fixed inputs. Positivity, finite integrals, correct support, and
using the same SDCS for the total and energy CDF are mandatory numerical
properties. The molecular endpoint is distinct from the free endpoint.

The Bethe constraint is a fast-projectile limit. Rudd et al. (1992), Sec. V,
pp. 453--454, discuss its use for protons above a few hundred keV and require
the equivalent projectile energy to exceed the relevant orbital kinetic
energies. It is not a low-incident-energy law. Their discussion on p. 462
explicitly excludes applying Kim's dipole-based fit below 200 keV.
See [Rudd et al.](https://doi.org/10.1103/RevModPhys.64.441).

At fixed secondary energy, retain the relativistic logarithm and its
associated `-beta^2` term with the same coefficient. The coefficient should
be checked against reliable optical electron-production information, with
the channel and inclusive-yield distinctions in `RECONSTRUCTION.md`.
An arbitrary finite-momentum continuation of optical data is not an exact
constraint. Neither is pointwise agreement with a coarse, incomplete optical
table. Conversely, retaining the logarithm alone does not validate a badly
wrong optical coefficient at high incident energies.

PJG itself notes on p. 160 that measured spectra extend beyond its classical
cutoff; the authors adopted a sharp cutoff for practical purposes. A bound
tail therefore repairs a stated PJG approximation, rather than contradicting
the model's intended physics. See [PJG](https://doi.org/10.1063/1.432812).
The [1977 corrections](https://doi.org/10.1063/1.323427) remain applied; in
particular, Eq. (16)'s apparent superscript plus denotes ordinary addition.

## Derivation within the two-Lorentzian structure

Let $m=m_ec^2$, $M=M_pc^2$, $E_e=m\beta^2/2$, and
$C_B=4\pi a_0^2R_y^2$. Retain PJG's width and center functions,

$$
\Gamma(E)=\Gamma_s+\frac{\gamma_1}{E_e+\gamma_2},\qquad
T_0(E)=T_s-\frac{t_a}{E_e+t_b},
$$

When energies are in eV, the numerators $\gamma_1$ and $t_a$ have units of
eV squared, despite the eV labels in printed Table III. The numerical values
are unchanged; the dimensions follow from the denominators in these equations.

Also retain its continuum fractions and logarithms,

$$
L_j(E)=\ln(4E_eC_j\gamma^2/I_j+e)-\beta^2.
$$

The printed $C_j/I_j$ values are nearly constant for each molecule. Their
individual rounding is retained; they are not refitted independently.

Use the two Lorentzians

$$
\ell_n=\frac{1}{(T-T_0)^2+\Gamma^2},\qquad
\ell_b=\frac{1}{(T-T_0)^2+\Gamma^2+\Lambda^2}.
$$

This ties the broad center to $T_0$ and parametrizes its squared width as
$\Gamma_b^2=\Gamma^2+\Lambda^2$. It is a structural change, not a printed
misprint correction. It retains one broad-width degree of freedom while
guaranteeing

$$
\ell_n-\ell_b
=\frac{\Lambda^2}
 {[(T-T_0)^2+\Gamma^2][(T-T_0)^2+\Gamma^2+\Lambda^2]}>0.
$$

The common center is important: merely setting the subtraction coefficient
to one while retaining distinct centers leaves a $T^{-3}$ term whose sign
can be negative. With the common center, the difference falls as $T^{-4}$,
whereas $\ell_b$ falls as $T^{-2}$. Evaluate the difference in factored form;
direct subtraction loses relative accuracy in the tail.
This is a simple sufficient positivity construction, not a claim that every
physical optical spectrum must have a $T^{-4}$ asymptote. Other constrained
centers can produce positive differences too; the printed unconstrained
centers do not provide that guarantee.

To see how this repairs PJG's amplitude, put $A_j=K\Gamma^2L_j$ and first
consider the nonrelativistic hard limit. The leading coefficient of
$A_j(\ell_n-B_j\ell_b)$ is $A_j(1-B_j)$. Fixing it to $N_eC_B$ gives

$$
B_j=1-\frac{N_eC_B}{A_j},\qquad
A_j(\ell_n-B_j\ell_b)
=A_j(\ell_n-\ell_b)+N_eC_B\ell_b.
$$

Thus the empirical subtraction is replaced by a coefficient constraint,
without introducing a finite-momentum target-response model. The factored
form remains positive even where the equivalent $B_j$ is negative.
In the relativistic free-domain test, replace the complete hard term by

$$
S_{\rm core}(E,T)=\frac{D(E)}{E_e}\sum_j f_j
\left\{K\Gamma^2L_j(\ell_n-\ell_b)
       +N_eC_B\ell_b F_B(E,T)\right\},
\quad 0\le T\le T_{\max}^{\rm free},
$$

with the fixed exact endpoint and complete Bhabha factor from
`RECONSTRUCTION.md`. The same PJG-shaped distortion multiplies both terms,
but uses an unbounded energy argument,

$$
D(E)=\frac{K_{\rm nr}^{p}}{J^{p}+K_{\rm nr}^{p}},\qquad
K_{\rm nr}=mE/M,\qquad p=1+\nu>0.
$$

Consequently $D\to1$, $T^2\ell_b\to1$, and the soft difference is subleading
in the joint fast-projectile/hard-secondary limit. At finite incident energy,
$D<1$ remains an empirical bound-target distortion: the molecular candidate
is not exactly a free-electron cross section there. Tests at extremely large
energies check this algebraic limit, not the physical validity of pointlike
protons at those energies.

The construction removes the independently fitted $B_0$, $B_1$, $E_0$, and
$T_1$, and keeps the PJG cutoff shift $\delta$ removed. It preserves the
functional forms of $F_j$, $\Gamma$, and $T_0$. It does not prove an accurate
optical coefficient, shell decomposition, or binary-endpoint broadening.

## Nonrelativistic shape-fit experiment

To test the data fit separately from relativistic endpoint broadening, use
$E_e=K_{\rm nr}$, omit the relativistic terms, and multiply each continuum
by the published Rudd bound-state factor

$$
P_j=\left[1+\exp\left(\frac{\alpha(T-T_{c,j})}{\sqrt{K_{\rm nr}I_j}}\right)\right]^{-1},
\quad T_{c,j}=4K_{\rm nr}-2\sqrt{K_{\rm nr}I_j}-R_y/4.
$$

The trial is $D/K_{\rm nr}$ times the continuum sum of
$f_jP_j\{K\Gamma^2\ln(4K_{\rm nr}C_j/I_j+e)(\ell_n-\ell_b)+N_eC_B\ell_b\}$.
Use stable logistic evaluation, and exclude $T+I_j>E$. This last condition
is necessary energy accounting, not an exact three-body threshold closure.
The cutoff center is empirical and is not a modification of free kinematics.
$\alpha=0.70$ for N2 and $0.59$ for O2 are fixed Rudd reference values.
Borrowing this factor is an explicit model modification, not a PJG erratum.

`fit_pjg_repair.py --refit` first varies J alone, then J/p/K, then
J/p/K/Gamma_s/Lambda. The objective gives equal aggregate weight to logarithmic
total and SDCS residuals. Totals sample the 1985 recommendation at 31
logarithmic energies from 5 to 4000 keV. The SDCS grid samples the 1992
recommendation at N2 energies 5, 10, 30, 50, 100, 300, and 1000 keV, and O2
energies 7.5, 10, 30, 50, 100, and 300 keV. Each spectrum uses 48 logarithmic
secondary energies from 2 eV to min(8 K_nr, 3000 eV), retaining points above
10^-4 of its 2-eV value. These are model comparisons, not independent
experimental points or uncertainty-weighted chi-squared fits.

| Five-parameter diagnostic | N2 | O2 |
| --- | ---: | ---: |
| J (eV) | 30.8026 | 6.89268 |
| p | 0.737826 | 1.35140 |
| K / printed K | 1.58220 | 0.884700 |
| Gamma_s (eV) | 4.48891 | 8.47138 |
| Lambda (eV) | 80.4867 | 319.984 |
| Total / recommended total, min--max | 0.905--1.058 | 0.966--1.093 |
| SDCS ratio, 10th--90th percentile | 0.856--1.197 | 0.868--1.095 |
| SDCS ratio, full sampled range | 0.757--1.436 | 0.709--1.208 |

The raw 1971 Table I totals were not used in this objective. Predicted/original
ratios are 1.017--1.069 for N2 and 0.944--1.016 for O2 at 50--300 keV.
These comparisons are useful but not statistically independent validation
against the evidence underlying Rudd's recommendations.

For the original 1979 N2 totals at 5--70 keV the ratios are 0.731--1.115.
The predicted mean energies, also not fitted, differ from the tabulated means
by -13.5% to +26.6%. At 5 keV the mean is 6.86 eV versus 5.42 eV, and the
tail above the nonrelativistic free endpoint is 0.199, versus about 0.141 for
the Rudd recommendation. These differences remain visible; a good total fit
does not excuse them. An eight-parameter exploratory fit drove K to its
upper search bound and was not retained as an identifiable improvement.

## Remaining checks and implementation boundary

The implied fast-projectile coefficient is
$a(T)=K\Gamma_\infty^2(\ell_{n,\infty}-\ell_{b,\infty})/C_B$.
It is finite and positive and gives a finite integral of $(T+I_j)a(T)$.
That is a mathematical improvement over assigning an unrestricted
$T^{-2}$ dipole coefficient, but is not a measurement of oscillator strength.
For the recorded N2 fit, the formal nonrelativistic coefficient at T=100 eV
is only about 0.14 of the Kim optical parameterization in Rudd's Table II.
That large discrepancy prevents claiming an optically validated asymptote.
The limiting relativistic PJG width uses E_e=m/2, not E_e=infinity; both
versions are exposed explicitly in the reference function.

### Constraining the existing N2 shape with optical information

`fit_pjg_optical.py` tests whether the same PJG shape can reproduce the
optical coefficient without a new target-response function. It transcribes
Rudd's Table II, Eqs. (32)--(33), for the N2 L_A, L_B, and K groups, with
thresholds 15.59, 28.8, and 410 eV. The tabulated quantity is $W\,df/dW$,
so the photoelectron coefficient is the sum of that quantity divided by
$W^2$, evaluated at $W=T+I_j$. Orbital populations are already included.

Fit K/Gamma_s/Lambda to 40 logarithmic secondary energies from 2 to 100 eV,
using logarithmic residuals. This chosen audit window avoids appending an
unverified high-energy continuation to the optical fit. It is not a claim
of experimental accuracy throughout that interval. The optical inputs are
photoelectron parameterizations, not a complete inclusive response including
all subsequent decay electrons.

The fitted values are K / printed K = 1.01808, Gamma_s = 11.5447 eV, and
Lambda = 138.389 eV. Thus the amplitude and narrow width remain close to
printed PJG. The predicted/reference optical coefficient spans 0.867--1.143
over the fitting grid, in place of the large discrepancy of the proton-only
fit. Freeze these three values before refitting proton-energy dependence.

Then vary only J, p, and a common multiplier of the printed C_j, using the
same proton objective as above. A common C_j multiplier changes the
non-logarithmic Bethe contribution, not the coefficient of the asymptotic
logarithm or its associated -beta^2 term. Compare retaining the printed
width numerator with eliminating it:

| Optically constrained N2 diagnostic | Retain gamma_1 | Set gamma_1 = 0 |
| --- | ---: | ---: |
| J (eV) | 27.9340 | 19.8814 |
| p | 0.648530 | 0.828636 |
| Common C_j multiplier | 0.669456 | 3.37170 |
| Total / recommended total, min--max | 0.903--1.025 | 0.880--1.076 |
| SDCS ratio, full sampled range | 0.627--1.799 | 0.762--1.299 |
| SDCS ratio, 10th--90th percentile | 0.791--1.355 | 0.824--1.179 |
| Mean-energy ratio at 5 keV, to original 1979 value | 1.416 | 1.216 |
| Fraction above NR free endpoint at 5 keV | 0.239 | 0.186 |

Removing gamma_1 also eliminates its unused denominator gamma_2. This is
a justified simplification to investigate further: with the optical shape
held fixed, it improves the proton spectral comparison substantially. Six
existing PJG quantities are adjusted across these two fits (three from the
optical reference, three from proton curves); no new independent shape
parameter is introduced. The original C_j are not fitted separately.

It is still not an accepted production fit. With gamma_1=0, predicted/original
1971 totals span 0.969--1.085, but the original 5-keV 1979 total is
underestimated by 28.9% and its mean energy is overestimated by 21.6%.
The comparisons inherit the reference normalization and channel limitations
already discussed. Rudd's Table II does not provide an O2 optical fit, so
this particular optical-constrained experiment is N2 only. The relativistic
coefficient uses the saturated center function rather than T_s. A separate
evaluation with gamma_1=0 gives optical ratios 0.869--1.141 on the same
2--100 eV grid. Neither comparison establishes the complete inclusive optical
response or behavior outside that finite audit window. At T=0 the
predicted/reference coefficient is 1.445 in the nonrelativistic limit and
1.451 in the relativistic limit. The spectrum is finite there, but its
threshold normalization is not established by fitting points above 2 eV.
This also needs a channel-resolved optical check; the grouped Table II
parameterization must not silently become an exact threshold constraint.

The next work is the joint optical/spectral validation for both targets and
combining the exact hard kernel with a positive molecular-tail treatment.
Do not append the free Bhabha polynomial outside its domain or hide a failed
comparison by clipping the SDCS. The endpoint-region angular closure remains
separate unfinished work.

Twenty-three Python checks cover the original tables, reference quadratures, the
stable Lorentzian identity, finite optical-strength integral, positive
spectra, numerical integration, distortion limit, and formal Bhabha/Bethe
limits. Direct-subtraction checks use a 60-digit independent reference to
avoid validating a stable expression against an inaccurate double-precision
subtraction. Optical checks reproduce the stated finite-window approximation
and verify that changing the Bethe constant does not rescale its coefficient.
These tests validate algebra and implementation, not the fit's
physical acceptance. All added work is host-side reference work; the constant-
time inverse-CDF GPU interface is unchanged, and no new CUDA/HIP/SYCL
performance claim is made.
