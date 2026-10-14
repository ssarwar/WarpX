> Historical research note, recovered 2026-09-14.
> Model names such as 'current' refer to this note's stage, not today's production model.
> Original commands and links retain their historical context; see the
> [archive guide](../../README.md) for status and portable reproduction.

# N2/O2 source audit and refit sensitivity

Research status, 2026-09-09: source-updated sensitivity fits have been run for
both gases. These are not final production coefficients. No commit, push,
production table, collision operator, or `Docs/` file is changed by this
study. The [matched-PJG candidate](MATCHED_PJG.md) is the baseline; the
[earlier optical/stopping trials](OPTICAL_STOPPING_AUDIT.md) remain
historical comparisons, not accepted replacements.

The exact free-electron maximum transfer, complete pointlike spin-1/2
tree-level Bhabha factor, inclusive effective-pair observable, and allowed
molecular tail remain unchanged.

## 1. Source selection and observable definitions

Here E is proton kinetic energy, T emitted-electron kinetic energy, and W
photon/energy-loss energy. For a resolved direct-ionization channel j,
W = T + I_j. A total photon spectrum cannot be converted to an electron
spectrum with one average binding energy. A fragment-ion yield is not
automatically an electron-production cross section.

| Source | Data and range | Role and qualifications |
|---|---|---|
| Crooks and Rudd (1971), Table I [1] | Original N2/O2 electron-production totals, 50--300 keV | Measured constraints; approximately 17% shared normalization uncertainty |
| Rudd (1979), Table I [2] | Original N2 totals and mean T, 5--70 keV | Total and first-moment constraints; absolute scale tied to the 1971 measurements |
| Rudd et al. (1983), Table V [3] | Authors' fit to absolute electron-production measurements, experiment 5--4000 keV | Measurement-derived total constraint, not raw point data |
| Rudd et al. (1985, 1992) [4,5] | Recommended total and analytic SDCS curves | The 1992 SDCS is a weakened interpolation prior; 1985 totals are a comparison |
| Cheng, Rudd and Hsu (1989), Fig. 1 [6] | O2 measured SDCS shapes, displayed at 7.5, 50, 150 keV | Thirty digitized open circles; separately normalized shape constraints |
| NIFS-DATA-109 (2010) [7] | Evaluated photoabsorption, including lines and core absorption | N2 25--100 eV optical diagnostic and N2/O2 integral-moment checks |
| Leiden N2/O2 files [8] | Evaluated absorption, dissociation and photoionization | O2 25--100 eV photoionization constraint; alternative N2 comparison |
| Gallagher et al. (1988) [9] | Evaluated state-resolved and fragment-resolved optical measurements | Channel definitions and source reliability; not yet numerical partial-channel fit inputs |
| Mahla and Mehnen (2025), paper and supplement [10] | O2 R-matrix total and three grouped partial spectra | Quantitative comparison only; raw theoretical thresholds are not adopted |
| Sun et al. (2005), Lin et al. (2013) [11,12] | N2/O2 optical and finite-q energy-loss spectra | Normalization and finite-q checks, not a claimed complete GOS fit |
| NIST PSTAR [13] | Electronic mass stopping, gas N/O, 1 keV--10 GeV | Independent energy-budget comparison at tabulated energies at or above 5 keV |

The local `elmolcs/Data/oos` N2/O2 files contain a few neutral-excitation
oscillator strengths, not photoionization continua. They cannot replace
the continuum inputs.

### Newly uploaded proton papers

Rudd (1983) explicitly identifies Table V as calculated from Eq. (17).
The negative charge yield is inclusive electron production; positive-ion
yield also includes capture. The code transcribes only the former, using
18 rows through 3 MeV and excluding the extrapolated 5 MeV row. The authors'
estimated fit uncertainties are 25%, 20%, 15%, 10%, and 8% at 5, 10, 25,
100, and 500 keV, with 8% retained above 500 keV. These are not independent
standard errors. Their discussion documents systematic differences between
accelerator datasets; Table V must not be called newly recovered raw data. [3]

Cheng's displayed SDCS were adjusted to recommended total cross sections.
Their shapes are additional evidence, but their absolute normalization
must not be counted independently. The reduced ordinate is

    Y = K_NR (T + 13.1 eV)^2 (d sigma/dT) / C_B,
    K_NR = (m_e/M_p) E,   C_B = 4 pi a_0^2 Ry^2.

The 13.1 eV binding is from the paper's Table I, not the 12.0697 eV
adiabatic O2 threshold. `source_datasets.py` preserves the raw marker
coordinates and calibration. This is approximate manual digitization,
not an author-supplied numerical table. [6]

At E = 50 keV and T approximately 295 eV, Rudd's 1992 formula is 0.438
times the selected Cheng marker. Shifting both reading coordinates by
plus/minus two pixels gives 0.330--0.580. This is reading sensitivity,
not an experimental confidence interval. The baseline PJG ratio is about
0.355; the weight-one source trial is about 0.752. Thus poorer agreement
with the Rudd formula can accompany better agreement with measured shapes.

## 2. Optical audit

The NIFS report was published in April 2010; its later repository deposit
is not a new measurement. Extraction retains 239 N2 and 729 O2 continuum
rows, discrete strengths, and O2 integrated bands. Redundant wavelength,
cross-section and mass-attenuation columns check transcription and units.
Its source selection already used sum rules and polarizability, so these
checks are not independent validation. The source notes document the
experimental/high-energy evaluations, primarily through Berkowitz (2002). [7]

Leiden's O2 continuum uses Brion et al. (1979) at 4--49 nm and Holland
et al. (1993) at 49--103 nm, with the latter's ionization efficiency. The
file was compiled in March 2020; website maintenance in 2026 does not make
these 2026 measurements. N2's comparison file cites Shaw et al. (1992),
Samson et al. (1987), and Heays et al. (2014), the last for neutral bands. [8]

Only O2's ionization column is fitted. N2 NIFS absorption is an approximate
smooth-continuum ionization diagnostic on 25--100 eV, not exact inclusive
electron yield. On that window, Leiden/NIFS spans 0.898--1.117 for N2;
NIFS absorption/Leiden ionization spans 0.927--1.028 for O2. Shared
underlying measurements preclude treating these as independent replicates.

Two qualifications are enforced in the readers:

- Tiny negative O2 dissociation values, of order 1e-32 cm^2 where absorption
  is of order 1e-17 cm^2, are subtraction roundoff. Only roundoff relative
  to local absorption is admitted, not physical negativity. The source is
  unchanged.
- N2's 0.1 nm resampling fails pointwise additive closure across some narrow
  sub-ionization neutral bands. Those are excluded. Absorption equals
  ionization on the fitted 25--100 eV interval. Interpolation also never
  bridges NIFS's separately tabulated O2 band interval, 9.75--12.07 eV.

The conversion is sigma_photo [Mb] = 109.76097 (df/dW) [eV^-1]. Changing
the cross section's abscissa from wavelength to energy introduces no
Jacobian into sigma itself; changing an integration variable does.
Integrated line/band strengths are not density samples.

### Integral moments of the full evaluated response

These calculations include the tabulated discrete/band strengths, stop at
100 keV, and add no invented high-energy tail:

| Quantity | N2, linear interpolation | O2, linear interpolation |
|---|---:|---:|
| Integrated oscillator strength | 14.0307 | 15.9725 |
| Electron count | 14 | 16 |
| Logarithmic mean excitation energy, eV | 83.8438 | 97.5802 |
| Polarizability, a0^3 | 11.9571 | 10.6479 |
| Strength below first ionization threshold | 1.22180 | 0.189684 |
| Continuum strength on 25--100 eV | 5.54841 | 7.18843 |

Log-log interpolation instead gives strengths 13.8994/15.8098 and mean
excitation energies 83.1648/96.9754 eV. This is interpolation sensitivity,
not experimental uncertainty. PSTAR uses 82/95 eV; those are stopping-power
mean excitation energies, not molecular binding thresholds. [7,13]

Itikawa et al. (1989), Table 2.10, provides a historical O2 check: strength
16, mean excitation energy 95.07 eV, and Ry-scaled inverse-square moment
2.648, corresponding to polarizability 4 * 2.648 = 10.592 a0^3. Its
partial optical spectra are graphical, not additional numerical tables. [14]

### The new O2 calculation needs qualifications

Mahla's full paper and public numerical supplement were examined. The
calculated ground-channel onset is 10.993 eV; NIST's adiabatic ionization
energy is 12.0697 eV. The calculation applies no term-energy correction.
The start of a measured spectral window must not be interpreted as the
physical ionization threshold. [10,15]

The supplement has one total and three grouped partial spectra, not seven
separate electronic-channel tables. At W = 23.313 eV their sum is 81.593 Mb
versus total 92.81 Mb, a ratio of 0.87914. This may reflect unreported
channels or tabulation mismatch; the cause is unresolved and the partials
are not silently renormalized. The raw subthreshold integral is 11.6712
Mb eV, or 2.069% of its integral from 10.993 to 40 eV.

Ratios of integrated total cross section to Leiden photoionization are
0.931/0.864/1.031/1.172/1.041 on W intervals 12.1--16/16--20/20--25/
25--30/30--40 eV. Each native grid is integrated separately, avoiding
comparisons of unresolved peak heights. The calculation is retained as a
comparison, not the absolute or threshold anchor.

## 3. Finite-q information

Sun's N2 measurements use 2500 eV electrons at 0, 2, 4, and 6 degrees.
The optical scale is normalized to Chan (1993) at 50 eV; the finite-q scale
uses a valence sum rule, Pauli correction, and tail extrapolation. These
are not new independent absolute calibrations. The approximately 23 eV
and 31.4 eV structures have different momentum dependence: the latter is
dipole-forbidden but quadrupole-allowed. [11]

Lin's O2 measurements use 2500 eV electrons at 0, 2, and 4 degrees and
electron-ion coincidences. The approximately 0.23 and 0.91 atomic-unit
values denote q^2, not q. The normalization uses Fan et al. (2005) at
24--26 eV; quoted uncertainty is 20--40%. Normalizing ion yield above
20 eV is not an independent demonstration of unity ionization efficiency. [12]

These data preclude identifying every finite-q feature with an optical
dipole feature. The matched SDCS is already integrated over momentum
transfer and does not define a unique F(q,W). A quantitative GOS fit is
therefore not claimed. The separately tested relativistic PWBA kernel in
`reference.py` still requires a target response.

## 4. Same-parameter-count refits

The formula, Bhabha factor, free maximum transfer, gate, printed channel
weights/thresholds, printed p, and limiting peak position are retained.
Six existing quantities are fitted: J, K/printed K, Gamma, Lambda,
peak-numerator multiplier q_t, and edge steepness. A five-parameter trial
sets q_t exactly to zero. No new independent shape parameter is added.

The optical forward diagnostic assumes the existing PJG allocation:

    A_j(T) = f_j A(T),
    (df/dW)_model = W sum_j f_j A(W-I_j) Theta(W-I_j).

This is not a measured photoelectron spectrum. Agreement in W does not
establish agreement in T. Features near W = 40 eV (N2) and 37 eV (O2)
in the plots are fixed model channel openings, not data discontinuities.

Residual blocks are divided by the square root of their sample counts.
The Rudd-1983 total block uses log ratios weighted by 0.15 divided by the
published fractional uncertainty. Other weights are 0.5 for the Rudd-1992
SDCS prior, 0.5 for Crooks totals, and 0.4/0.5 for N2's 1979 totals/mean T.
Each Cheng spectrum has its mean log normalization removed. Optical
weights 0.5, 1, and 2 are tested at 60 logarithmically spaced W points on
25--100 eV. This is not chi-squared or an uncertainty estimate. Stopping
power is checked, not fitted as a penalty in these source-refresh trials.

| Check | N2 baseline | N2 weight-one trial | O2 baseline | O2 weight-one trial |
|---|---:|---:|---:|---:|
| Optical W = 25--100 eV ratio range | 1.106--1.749 | 0.869--1.329 | 0.925--1.458 | 0.843--1.260 |
| Optical log-RMS residual | 0.337 | 0.135 | 0.218 | 0.123 |
| Total/Rudd-1983-fit range | 0.851--1.029 | 0.890--1.042 | 0.764--1.032 | 0.874--1.001 |
| SDCS log-RMS against Rudd-1992 formula | 0.121 | 0.192 | 0.072 | 0.178 |
| Maximum pair stopping/PSTAR | 0.9672 | 0.9076 | 1.0057 | 0.9853 |

Ranges refer to the stated comparison grids, not all energies. The optical
window is not the earlier 2--100 eV secondary-energy comparison.

The following trial values identify the calculation, not a production
recommendation. Full precision is preserved in the output JSON.

| Parameter | N2 | O2 |
|---|---:|---:|
| J, eV | 15.0675338 | 8.44061243 |
| K/printed K | 0.991878236 | 0.847120840 |
| Gamma, eV | 10.9414922 | 16.3277099 |
| Lambda, eV | 115.010246 | 144.971060 |
| q_t | 0.103638681 | 0.549686525 |
| Edge steepness | 1.09065475 | 1.03023837 |

J-only cannot improve the asymptotic optical coefficient; it changes the
projectile distortion instead. J-only fits give about 12.42/5.707 eV for
N2/O2. Eliminating q_t changes the weight-one N2 objective norm only from
0.20395 to 0.20626, but worsens O2 from 0.21638 to 0.29763. Objective norms
with different optical weights are not directly comparable.

Important unresolved tradeoffs:

- N2 mean-T ratios are 0.882--1.281 against the original 1979 values;
  a smaller optical residual alone does not justify accepting the trial.
- O2 Cheng shape log-RMS residuals improve from 0.127/0.299/0.179 to
  0.044/0.166/0.139 at 7.5/50/150 keV. Absolute ratios still span roughly
  0.75--1.38 at 50 keV, with non-independent normalization.
- O2's pair/PSTAR maximum is below one, but leaves only about 1.5% for
  other electronic loss near 1 MeV. This is a necessary energy-budget
  check, not complete stopping validation. Optical weight two leaves
  about 3.9% at its maximum, with different spectral tradeoffs.
- Integrated model soft oscillator strengths change from 10.897 to 8.439
  for N2 and 12.472 to 11.249 for O2. A finite-window fit does not supply
  missing inner-shell response or settle inclusive decay multiplicities.
  Missing excitation must not be concealed by an arbitrary rescaling.

The hard constraint remains intact. At 800 MeV incident energy, full
SDCS/free-Bhabha ratios at T = 1/10/100 keV are 1.00909/1.000809/1.000080
for N2 and 1.02443/1.001385/1.000128 for O2. At 10 GeV, second kinetic
moments divided by the analytic free Bhabha moment are 0.999895/0.999883.
These check the specified pointlike model, not finite-size or higher-order
QED effects.

## 5. Numerical verification

The new reader tests cover optical units, abscissa conversion, additive
channels, negativity, interpolation support, analytic moments, lines/bands,
gaps, source-table semantics, and digitized spectrum grouping. Together
with existing physics checks, 52 Python tests pass; both independent host
C++ kinematics/sampling tests also pass. No GPU benchmark is claimed for
these host research fits.

Refining stopping segments from 1025 to 4097 points changes checked total,
first, and second moments by less than 1.8e-8 relative at 5 keV, 50 keV,
1 MeV, 800 MeV, and 10 GeV. The fitting total grid agrees with a differently
segmented grid within 3.5e-11 on the total-fit energies. Two additional
optimizer starts per gas return spectra within 2.6e-6 relative on the
SDCS prior grid. Numerical stability does not establish global
identifiability or experimental accuracy.

## 6. Access gaps and newer work screened

The uploaded Rudd, Cheng, Gallagher, and Mahla papers are available.
Mahla's public publisher-linked numerical supplement was downloaded;
no upload is needed for it. The most useful outstanding items are:

1. Lopez-Patino et al. (2016), *Low energy ionization and fragmentation
   cross sections for H+ impact on N2 and O2*, IJMS 405, 59--63,
   [DOI 10.1016/j.ijms.2016.05.014](https://doi.org/10.1016/j.ijms.2016.05.014).
   Only abstract/highlights were accessible. Its 2--10 keV target-ion
   channels need a full-text capture/multiplicity check before use as
   inclusive electron yield. No numerical values were fitted.
2. Gallagher/Brion/Samson/Langhoff, **JILA Data Center Report No. 32**.
   The 1988 review identifies this as its numerical figure tabulation
   (printed p. 20). It was not located in accessible full text. The review's
   figures remain available for approximate digitization. [9]

Additional sources screened, not numerical fit constraints:

- Zammit et al. (2025), [O2 photoionization theory](https://doi.org/10.1088/1361-6455/add43c):
  accessible theory, not selected as a more accurate absolute anchor.
- [N2 Fano-resonance study (2025)](https://doi.org/10.1002/jcc.70067):
  narrow near-threshold and multiphoton scope; no numerical input extracted.
- Liu et al. (2019), [low-lying GOS transitions](https://doi.org/10.1063/1.5087603):
  discrete transitions, not a complete continuum dataset.
- Knudsen et al. (1995), [proton/antiproton molecular ionization](https://doi.org/10.1088/0953-4075/28/16/011):
  possible additional N2 ion-channel check, not O2; no numerical use here.
- Song et al., [new oxygen electron-collision review, accepted manuscript](https://discovery.ucl.ac.uk/id/eprint/10222954/1/JPCRD25-AR-00038-revised.pdf):
  revised November 2025, literature surveyed through 2024. Its evaluated
  electron-impact cross sections are not new proton-impact measurements.

The search does not establish that no newer absolute proton SDCS exist.
The pending optical-channel and low-E issues prevent promoting these trials
as a completed physical calibration.

## 7. References for numerical inputs and checks

1. Crooks and Rudd, PRA 3, 1628 (1971), Table I:
   [original proton totals](https://doi.org/10.1103/PhysRevA.3.1628).
2. Rudd, PRA 20, 787 (1979), Table I:
   [N2 totals and mean electron energies](https://doi.org/10.1103/PhysRevA.20.787).
3. Rudd, DuBois, Toburen, Ratcliffe and Goffe, PRA 28, 3244 (1983):
   [proton ionization and capture](https://doi.org/10.1103/PhysRevA.28.3244).
4. Rudd et al., RMP 57, 965 (1985):
   [recommended electron-production totals](https://doi.org/10.1103/RevModPhys.57.965).
5. Rudd et al., RMP 64, 441 (1992):
   [differential cross sections and optical fits](https://doi.org/10.1103/RevModPhys.64.441).
6. Cheng, Rudd and Hsu, PRA 40, 3599 (1989):
   [7.5--150 keV O2/CO2 electron distributions](https://doi.org/10.1103/PhysRevA.40.3599).
7. Sakamoto et al., *Oscillator strength spectra and related quantities of
   9 atoms and 23 molecules over the entire energy region*, NIFS-DATA-109
   (2010): [evaluated report and numerical tables](https://nifs-repository.repo.nii.ac.jp/records/11706).
   Its N2/O2 source notes identify Chan (1993), Holland (1993), Berkowitz
   (2002), Henke, Zhadenov (1987), Barrus (1979), and Chantler (1995) for
   the respective intervals. These are provenance through the evaluation,
   not separately refitted independent datasets.
8. Leiden: [N2 file](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/N2/N2_0.1nm.txt),
   [O2 file](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/O2/O2.txt),
   [source bibliography](https://home.strw.leidenuniv.nl/~ewine/photo/references.html).
   Cite [Heays, Bosman and van Dishoeck (2017)](https://doi.org/10.1051/0004-6361/201628742)
   and [Hrodmarsson and van Dishoeck (2023)](https://doi.org/10.1051/0004-6361/202346645).
   Relevant continuum sources:
   [Brion (1979), O2](https://doi.org/10.1016/0368-2048(79)85032-X),
   [Holland (1993), O2](https://doi.org/10.1016/0301-0104(93)80148-3),
   [Shaw (1992), N2](https://doi.org/10.1016/0301-0104(92)80097-F),
   [Samson (1987), N2](https://doi.org/10.1063/1.452452).
   [Heays (2014)](https://doi.org/10.1051/0004-6361/201322832) supplies
   neutral-band physics, not a proton continuum measurement.
9. Gallagher, Brion, Samson and Langhoff, JPCRD 17, 9--153 (1988):
   [optical compilation](https://doi.org/10.1063/1.555821).
10. Mahla and Mehnen, JCP 163, 234312 (2025):
    [O2/O3 photoionization](https://doi.org/10.1063/5.0307721);
    [numerical supplement v1, CC BY 4.0](https://doi.org/10.60893/figshare.jcp.30757583.v1).
11. Sun et al., Chinese Phys. 14, 1378 (2005):
    [N2 oscillator-strength densities below 100 eV](https://doi.org/10.1088/1009-1963/14/7/019).
12. Lin et al., CPB 22, 023404 (2013):
    [momentum dependence of O2 ionization/dissociation](https://doi.org/10.1088/1674-1056/22/2/023404).
    Normalization source: [Fan et al., PRA 71, 032704 (2005)](https://doi.org/10.1103/PhysRevA.71.032704).
13. [NIST PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html),
    nitrogen material 007 and oxygen 008. Original rows and retrieval metadata
    are in `pstar_reference.json`; use electronic, not total, mass stopping.
14. Itikawa et al., JPCRD 18, 23--42 (1989), Table 2.10 and Figs. 3.1--3.2:
    [oxygen electron/photon review](https://doi.org/10.1063/1.555841).
15. [NIST WebBook: O2 ion energetics](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=20).

## 8. Reproduction

Scripts require NumPy/SciPy, with `pypdf` for extraction and Matplotlib
for plots. They do not download data, initialize WarpX, or modify production
files. The current external-input directory is
`tmp/pdfs/dataset-refresh-20260909` under the worktree root.

```sh
python Tools/Algorithms/ProtonImpactIonization/extract_nifs_oscillators.py \
  tmp/pdfs/dataset-refresh-20260909/nifs-data-109.pdf \
  --output tmp/pdfs/dataset-refresh-20260909/nifs-optical.json
python Tools/Algorithms/ProtonImpactIonization/optical_dataset_audit.py \
  tmp/pdfs/dataset-refresh-20260909/nifs-optical.json \
  --output tmp/pdfs/dataset-refresh-20260909/optical-moments.json
python Tools/Algorithms/ProtonImpactIonization/study_source_refresh.py \
  --nifs tmp/pdfs/dataset-refresh-20260909/nifs-optical.json \
  --o2-leiden tmp/pdfs/dataset-refresh-20260909/leiden-o2.txt \
  --n2-leiden tmp/pdfs/dataset-refresh-20260909/leiden-n2-0.1nm.txt \
  --mahla tmp/pdfs/dataset-refresh-20260909/mahla-o2-o3-supplement.zip \
  --output tmp/pdfs/dataset-refresh-20260909/source-refits.json
python Tools/Algorithms/ProtonImpactIonization/audit_source_refresh.py \
  --results tmp/pdfs/dataset-refresh-20260909/source-refits.json \
  --nifs tmp/pdfs/dataset-refresh-20260909/nifs-optical.json \
  --o2-leiden tmp/pdfs/dataset-refresh-20260909/leiden-o2.txt \
  --output tmp/pdfs/dataset-refresh-20260909/source-checks.json \
  --figures tmp/pjg-source-refresh-figures --multistart
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
```

Cheng's digitization refers to this 510x900 Poppler crop; raw coordinates
remain in `source_datasets.py`:

```sh
pdftoppm -f 2 -l 2 -scale-to 3000 -x 1440 -y 300 -W 510 -H 900 -png \
  /Users/ssarwar/Research/warpx-with-resources/resources/PhysRevA.40.3599.pdf \
  tmp/pdfs/dataset-refresh-20260909/cheng-figure1
```

Source SHA-256 digests (input files are unchanged):

| Input | SHA-256 |
|---|---|
| Rudd 1983 PDF | `51d89fa39e0583633e8d22ca6d766e94cfa9d9535999ca1d4703cbea4a86d482` |
| Cheng 1989 PDF | `9f7850bc3a58a880bf6c0e133b3d8e6555c57f7cb65bb57e015e46c4ff658189` |
| Gallagher 1988 PDF | `4bc711dd7636cbd6622dc33a78abd77d51ceb33723144146581b20a85e4fe09a` |
| Mahla 2025 PDF | `daa5d6be44241ef1aa78f64c14f3ba57728e734787f4f5862b4cb5a7eff4b216` |
| Mahla supplement | `14555867f779ca461c01d141cab000b4fbdf6997f9f3c38afe94c90ccc8d7310` |
| NIFS report | `432332e8835a00c1e0fcee32d1cc8a85d4a983f8841c799abd832f9e19af76d7` |
| Leiden O2 | `3787a1dd6f474f6c0bc691dfb7dfbb3a76901afdde892826c36b4b098a15d72c` |
| Leiden N2 0.1 nm | `075651a0beeb895f3492f7ad2ff0d2d2fdebafac5480e924c11358ab1c349e44` |
| Sun 2005 | `352a0807b51f0d2a49d52556ae784c063d23308436a6553bba3264644ea50ebf` |
| Lin 2013 | `061d37222c8c6cbd1cb88ec45045aa3cc06097aca7b618156076d013bd735d50` |
| Itikawa 1989 | `6671014fa1b4718219686eff3b1fada9f905d63fda1b7eb17070f3bfa80e40e5` |

Output JSON preserves input hashes and full precision. Plots are
`source_optical_stopping.png` and `cheng_1989_sdcs.png` in the requested
figure directory. Neither output is an accepted production table.
