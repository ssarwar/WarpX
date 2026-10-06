.. _mcc-iaa-sources:

Electron--air collision source audit
===================================

The IAA models refer to A. Schmalzried's 2023 thesis, especially
Sections 11.2 and 11.6, Table 11.12, and Eqs. (2.66), (11.119)--(11.122).
The alternative elmolcs model refers to the supplied archive with SHA256
``b864381086120fff37eb8c658dfcb0f150b80969cf9add626a38d89bab0a3c38``.
This identifier describes a fixed source snapshot, not an installed package version.

RBEQ parameterizations
---------------------

The original WarpX N2 and O2 shell parameters agree with thesis Table 11.12.
The fitted elmolcs ``iaa*`` N2 block instead splits the thesis 18.72 eV shell
into 18.746 and 23.6 eV contributions, changes 37.30 to 37.8 eV, and changes
16.74 to 16.716 eV. Its O2 inner binding is 531 eV instead of 543.8 eV.
The package's default O2 orbital data are a third, eight-contribution model
and must not substitute for either fitted model.

The thesis uses the asymptotic dipole correction in Eq. (11.121).
The package evaluates its logarithmic moment at 1 MeV instead.
Furthermore, ``genRBEQ`` uses ``1+(E+B)/mc2`` in its exchange and squared-binding
terms, whereas ``_RBEQ`` uses ``1+(E+U+B)/mc2``. The latter is the common
kernel underlying the differential and cumulative package helpers and the
selected snapshot model. This small additional helper inconsistency must
be distinguished from the fitted-parameter differences.
The ``lnBm`` antiderivative also differs from numerical integration of its
stated logarithmic moment: its last logarithmic coefficient is 1/4 instead
of 1/2. The snapshot option preserves that formula; it is not a correction
of the underlying moment integral. Independent quadrature checks include
this explicit difference.

The numerical N2 ``iaa`` table is close to the thesis parameterization;
the analytical ``iaa*`` fit is a different model. O2's numerical table is
close to its package fit. Matching production tables must be regenerated.
Negative shell totals are replaced by zero before summing production rates.
Where the differential formula is negative, the conditional sharing is
uniform until the whole differential distribution is nonnegative, including
the narrow continuation within 0.1 percent of a shell's binding energy.
Signed source totals are retained only for reference comparisons.

IAA ionization uses Eq. (2.66) for the secondary at all energies and
Eq. (2.60) for the primary. For 2.5 MeV impact, binding 15.58 eV and equal
sharing, the secondary mean cosine is approximately 0.8804. The ion recoil
is enforced separately through the three-product kinematics solve.
The omission of Coulomb interference is intentional.

Rotation and inclusive elastic scattering
----------------------------------------

The residual elastic integral in Eq. (11.10) includes rotation. It must not
be used unchanged alongside additional rotational channels. Its normalization
is intentionally independent of the elastic DCS integral. The DCS continuation
above 10 keV is the screened Rutherford model.

The generated numerical N2 rotational tables are quarantined: their 0-to-2
and 0-to-4 resonances occur near 0.800 eV, while the elementary tables peak
near 2.3 eV. The supplied generator is not implemented, so an energy-axis
repair cannot be justified. The elementary 0-to-6 threshold header is also
incorrect: the rigid rotor gives 42 B, approximately 0.010404 eV.
All thresholds must be computed from canonical levels, not rounded labels.

The canonical rotational constants are 1.998 and 1.438 inverse centimeters
for N2 and O2. Nuclear-spin weights are 6:3 for even:odd N2 levels and 0:1
for O2. At 1000 K, the short lists of initial states in the package omit
approximately 19.5 and 18.6 percent of the respective equilibrium populations.
A thermal sum must extend the state set, rather than renormalize those lists.

IAA Eq. (11.24) constructs transitions from elementary integral rates using
the sudden approximation. The V3/V4 rate models use the canonical internal-energy
threshold :math:`\Delta=E_{J'}-E_J`, without molecular recoil shifts. For
relativistic electron momentum :math:`p(E)=\sqrt{E(E+2m_ec^2)}/c`, reverse rates
satisfy the heavy-target integral relation

.. math::

  g_Jp(E+\Delta)^2\sigma_{J\to J'}(E+\Delta)
  =g_{J'}p(E)^2\sigma_{J'\to J}(E).

For the V3/V4 rate models, the runtime draws a discrete rotational change
independently of the angle from the existing inclusive elastic DCS.
At exactly zero relative momentum there
is no incident axis, so emitted electrons use the isotropic limit. A selected
excitation below its exact level spacing, possible through energy-grid
rounding, becomes an unchanged event. No other rotational rate is increased.
There is no angular conditioning, transition search, or rejection in the
alias sampler, and it requires no cumulative table.

The spectator relation in Eq. (11.29), :math:`J_R=kR\sin(\theta/2)`, is not a
hard quantum selection rule. The V3/V4 models use neither this relation nor
a separate rotational DCS.
The elastic DCS supplies the vibrationally elastic angular marginal, including
unresolved rotations; see thesis Eq. (8.129), the end of Section 11.1 and
Section 12.1. Its integral and the residual rate have separate normalization.

Recoil uses the sampled laboratory angle and relativistic transformations.
Exact two-body kinematics applies outside the small energy-only band
:math:`0\le E-\Delta\le 2(m_e/M)E` for excitation. Inside that band, the electron
retains :math:`E'=E-\Delta`; its virtual neutral receives the momentum
difference, and its recoil energy is neglected in the electron update.
The resulting energy defect is bounded by :math:`2(m_e/M)E` for N2/O2.
This explicit heavy-target continuation retains the nominal threshold and
independent angle without rejecting outcomes or projecting angles. It does
not claim exact four-momentum conservation in the continuation band.
Unchanged events retain the ordinary elastic recoil implementation.

Integral detailed balance in the V3/V4 reference rates does not imply differential
detailed balance for an energy-dependent elastic DCS used independently of
rotational outcomes. This is the selected approximation. Independent tests
cover level thresholds, unchanged/loss/gain probabilities, angular independence,
signed recoil, the continuation bound and thermal power balance.

Physical bundles are generated offline and kept in ``warpx-data``. Runtime
initialization only validates loaded data, weights equilibrium populations,
prepares shared in-memory samplers and uploads them to the device. It neither
evaluates elmolcs source models nor creates cross-section files.

The elementary N2 tables expose a normalization conflict near the resonance:
at 2.3 eV their excitation sum is :math:`2.40\times10^{-19}\,\mathrm{m^2}`,
larger than the residual integral :math:`1.73\times10^{-19}\,\mathrm{m^2}`.
An earlier coarse candidate audit missed this interval. A nonnegative
unchanged rate cannot be obtained by subtracting these particular source sets.
The separate Kutz--Meyer :math:`J=0\to0` table is a fixed-initial-state
theoretical result, not the measured residual minus rotation and not a
thermal average over :math:`J\to J` channels. Replacing the residual by its
sum with the excitation tables changes the physical rate substantially.
In the sudden approximation the unchanged channel for an arbitrary initial
state is

.. math::

  \sigma_{J\to J}(E)=A_0(E)+\sum_{\lambda=2,4,\ldots}
  |C_{J0,\lambda0}^{J0}|^2 A_\lambda(E),
  \qquad A_\lambda=\sigma_{0\to\lambda}p_\mathrm{in}/p_\mathrm{out}.

Simply reusing :math:`\sigma_{0\to0}` for every initial state misses the
second term. V4 bundles carry these state-resolved unchanged rates, which
are thermally averaged and collapsed to one outcome during initialization.
At 300 K the source-based N2 reconstruction at 2.3 eV still totals about
:math:`5.09\times10^{-19}\,\mathrm{m^2}`; Boltzmann averaging does not cure
the normalization conflict. States through J=48 for N2 and odd J=57 for O2
leave less than :math:`10^{-10}` population outside a 300 K partition sum.

Kinetic spectator closure
~~~~~~~~~~~~~~~~~~~~~~~~~

The separate ``iaa_spectator`` model implements discrete outcomes in place of
the mean rotational loss in Eq. (2.48). It retains the IAA inclusive elastic
rate and samples its angle first. For each initial rotational level, Eqs.
(11.24) and (11.35) give nonnegative differential weights

.. math::

   S_{if}(E,\theta) = p(E-\Delta_{if})
     \sum_{L=0,2,4,6}|C_{i0,L0}^{f0}|^2 A_L(E)h_L(E,\theta),
   \qquad
   P_{if}(E,\theta) = \frac{S_{if}(E,\theta)}{\sum_k S_{ik}(E,\theta)}.

Here :math:`p(E)=\sqrt{E(E+2m_ec^2)}/c`, closed excitation channels have zero
weight, and :math:`h_L` is the angularly normalized squared spherical Bessel
function from Eq. (11.35). The elementary amplitudes are
:math:`A_0(E)=\sigma_{00}(E)` and
:math:`A_L(E)=\sigma_{0L}(E)p(E)/p(E-\Delta_{0L})` for positive ranks,
The reduced amplitude is held at its first positive source value between
the canonical threshold and the first datum, as in the offline source exporter.
The common incident-momentum factor cancels in the conditional normalization.
The thermal outcome probability is
:math:`p_i(T_\mathrm{rot})P_{if}(E,\theta)`. The rate of each outcome is therefore
determined by its angular average under the inclusive IAA DCS. It is not the
unmodified elementary rotational cross section multiplied by a population.
No angular balancing step restores the original integral channel rates.
If :math:`F_{\mathrm{el}}(E,\Omega)` denotes the normalized elastic angular
distribution, the effective differential cross sections are

.. math::

   \frac{d\sigma^{\mathrm{eff}}_{if}}{d\Omega}
   = \sigma_{\mathrm{el}}(E) F_{\mathrm{el}}(E,\Omega)
     p_i(T_\mathrm{rot})P_{if}(E,\theta).

Their sum recovers the inclusive elastic rate and angular marginal, without
subtracting incompatible absolute rotational differential cross sections.

The average sampled loss is Eq. (2.48) evaluated with these probabilities,
while the discrete draw also retains energy diffusion and angle--energy
correlations. The source set uses the supplied N2 elementary tables, including
the low-energy Itikawa recommendation, corrected level thresholds and N2
nuclear-spin weights. The spectator geometry uses the Kutz--Meyer separation
:math:`R=2.068a_0`. De-excitation uses the same incident-energy sudden
weights as excitation. Normalizing these weights does not enforce detailed
balance, even if a preceding differential source model were reciprocal.
Equilibrium heating/cooling must be measured as a model property; the closure
must not be described as preserving a Maxwellian electron bath.
The optional ``--thermal-balance`` source check reports this defect separately
from tabulation accuracy; it does not assert that the defect vanishes.

``Tools/CrossSections/thermiaa_spectator.py`` tabulates the angular basis and
finite-threshold transition coefficients offline in a V5 bundle. Factoring
the differential weights avoids a large table of normalized probabilities.
The GPU selects a virtual initial level with a small precomputed Boltzmann
alias and evaluates only seven candidate final levels from loaded weights.
There is no loop over all populated states, no Bessel or Clebsch--Gordan
evaluation and no neutral-state evolution. Canonical losses use double
precision; nonnegative interpolation weights are stored in binary32.
Energy knots remain in binary64, and the interval immediately above the
first excitation threshold uses square-root interpolation of row probabilities.
This retains the threshold law in the narrow interval where further energy
grid refinement would no longer give distinct binary32 knots.
The independent source check uses triple-Legendre quadrature for the angular
momentum coefficients and direct quadrature of Eq. (2.48).

The supplied reference set stops at rank six and 1000 eV. Its availability
does not establish convergence of higher-rank rotational rainbows. The
spectator angular approximation is known to be inaccurate in the thermal
and resonance regimes (Fig. 11.15); using it there is an explicit part of
this Thetaermiaa-style closure. It is not a hard angular accessibility rule.
The finite elastic rate also differs from a superelastic rate coefficient
that remains nonzero as the incident energy vanishes. These physical limits
are distinct from tabulation and sampling errors.

``Tools/CrossSections/export_elmolcs.py`` exports the source elastic,
DCS and elementary rotational cross sections offline. It preserves positive
source knots and refines for the ordinary linear cross-section reader.
The residual elastic continuation above 6 keV is the package's fitted Born
model, matched at the last datum. DCS rows retain the source resolution.
Rotational tables use canonical thresholds and interpolate
:math:`\sigma p_\mathrm{in}/p_\mathrm{out}`, with a constant reduced amplitude
between threshold and the first positive datum. No rotational tail is
extrapolated beyond the source endpoint.

O2 export uses the actual numerical integral table attributed to
Takayanagi--Itikawa, rather than the previous replacement polynomial for
Eq. (11.21b). N2 elementary coverage ends at 1000 eV and O2 at 20 eV.
These endpoints do not establish physical validity of the approximations
throughout their energy ranges, or bound omitted higher-rank transitions.
The optional ``--test-bundles`` output belongs in a build directory and
exercises the software with these real sources. It uses N2's separate
:math:`0\to0` table and a ground-rotor reference decomposition for O2;
it is not a resolution of the production source-consistency questions.
Neither clipping a negative remainder nor publishing synthetic fixtures as
physical data is an acceptable resolution.

``Tools/CrossSections/test_elmolcs.py`` compares exported tables to source
knots, checks interpolation at independent probes, and checks equilibrium
power at 100, 300 and 1000 K. Its device references use independent angular
quadrature. The sampler tests, full MCC inputs and benchmark driver accept
prepared data; no source model or cross-section file is generated at runtime.
Finite inclusive cross sections imply :math:`v\sigma\to0`, whereas thermal
superelastic rates can remain finite. This low-energy distinction also needs
an explicit reference-population model in a production combined family.

Reciprocal hybrid rotational model
---------------------------------

``rotation_model = reciprocal_hybrid`` selects the V6 production model for
N2 and O2. The rotational bath has a fixed temperature; the 300 K bundles
are distributed in ``warpx-data/MCC_cross_sections/{N2,O2}/IAA/reciprocal_hybrid_300K``.
Other temperatures require an offline export. The inclusive elastic input
remains a consistency check; the bundle supplies the complete collision
family, including its rate and angular tables.

The N2 low-energy model retains the elementary elmolcs rates and an isotropic
leading rotational angular approximation below 1 eV. Between 1 and 1.25 eV
it joins the Read resonance tensors with Kutz--Meyer strengths and a direct
background constrained by Jung's vibrationally elastic data. Between 4 and
10 eV it joins the completed Gote transfer-rank distributions, which constrain
the model through 200 eV. O2 uses the transition-specific Born construction
of thesis Eq. (11.21b) below 1 eV, a positive bridge to 20 eV, and
Bhattacharyya's integral and momentum-transfer constraints through 200 eV.
Unresolved O2 higher ranks are completed with a positive relative-entropy
fit to those constraints and the inclusive IAA angular marginal.
Both gases transition from 200 to 1000 eV at equal momentum transfer to a
two-centre spectator model restricted to bound rigid-rotor levels.
``Tools/CrossSections/reciprocal_rotation/SOURCES.rst`` and
``IMPLEMENTATION.rst`` specify the source assumptions and completions.

For an excitation :math:`i\to f` with level gap :math:`\Delta`, a common
positive primitive :math:`X_{if}` defines both directions:

.. math::

   D_{if}(E,\theta) &= \frac{p(E-\Delta)}{p(E)^2}X_{if}(E,\theta),\\
   D_{fi}(E-\Delta,\theta) &= \frac{g_i}{g_f p(E-\Delta)}X_{if}(E,\theta).

This enforces :math:`g_i p(E)^2 D_{if}=g_f p(E-\Delta)^2 D_{fi}` in the
heavy-target reference. The inclusive normalization couples forward kernels
at :math:`E` to reverse partners evaluated at :math:`E+\Delta`. It does not
independently rescale the two incoming-energy rows after constructing reverse
rates. The low-energy changing rates are protected while a nonnegative
unchanged background satisfies the inclusive constraint. A negative physical
remainder is an export failure. The continuous reference is reciprocal;
finite mixture and alias tables approximate it to the measured accuracy.

The initial-state Boltzmann sum is extended until the omitted population is
below :math:`10^{-10}` and transfer moments converge. Canonical double-precision
level differences remain discrete. N2 nuclear-spin weights are 6:3 for even:odd
levels; O2 uses odd rotational levels with unresolved electronic spin.
The high-energy closure excludes levels at or above the model dissociation
limit in both directions before normalization. This bounded rigid-rotor
continuation does not resolve centrifugal distortion or vibration at high
angular momentum. The O2 bridge below 20 eV is an interpolation model and
does not establish a measured low-energy resonance DCS.

Offline energy continuations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The supported relative electron energy is 0--1 GeV. The inclusive residual
elastic source extends through 6 keV and then uses the continuously matched
elmolcs Born continuation. The angular distribution blends to screened
Rutherford between 8 and 10 keV, retaining its energy dependence above the
join. The short elementary rotational endpoints (1 keV for N2, 20 eV for O2)
are source limits, not constant rotational tails. Additional upper-energy
support for reverse evaluations is internal to the exporter.

Below 1 meV a cold continuation joins the reciprocal changing kernels to an
unchanged s-wave background. A finite zero-energy superelastic rate requires
storing :math:`K=v\sigma` directly. It cannot be represented by a finite
cross section at :math:`v=0`. The zero-momentum emission is isotropic.
The inclusive rate equals the IAA rate above this join. Excitations remain
exactly forbidden below their canonical level gaps.

Queries outside the supported interval fail. The relative endpoint allowance
is eight machine epsilons of particle precision, solely for roundoff. Marked
ordinary tables also enforce their declared domains before either cached
selection or direct lookup. The MCC block checks the intersection once per
candidate collision; its majorant uses only that supported interval.
Unmarked legacy tables retain endpoint clamping and the associated
infinite-energy majorant bound.

Sampling and storage
~~~~~~~~~~~~~~~~~~~

Below 1 keV, the collision kernel locates the energy interval, interpolates
its rate, and chooses an endpoint row with probability proportional to its
contribution to that rate. It samples the angle from that same row's inverse
inclusive CDF. At the resulting angular quantile it draws unchanged/change
and, for a change, a discrete signed gap from a prepared alias distribution.
Mixtures interpolate probability distributions, not level energies. Above
1 keV, the actual energy and deflection select a compact momentum-transfer
distribution. Every unchanged event still scatters and recoils.

The sampler carries :math:`1-\cos\theta` in double precision through the
momentum-transfer and recoil calculations. The signed-loss kinematics retain
the documented nominal-threshold continuation in the narrow recoil band;
there is no new angular veto. Sufficient-resolution uniforms are used even
when both field and particle precision are single. Single-precision particle
storage can nevertheless round away meV changes in MeV electrons; table and
sampling accuracy must be distinguished from the accuracy of stored velocities.

All source evaluation, population sums, inverse CDFs, and alias construction
are offline. Runtime initialization reads, validates, shares, and uploads
immutable arrays, releasing bulk host staging storage. Alias entries use
8 bytes (binary32 cutoff and two checked uint16 indices); energies remain
binary64. The 300 K N2 and O2 payloads total 751,195,072 bytes, about 716 MiB,
including auxiliary arrays and exact-search indices. The reader enforces a
1 GiB combined per-process table budget. One MPI process per GPU therefore meets the table budget;
multiple MPI processes on one GPU each own a copy.

Independent intermediate-source checks of the 300 K low-energy tables found
maximum rate/transfer-moment interpolation errors of 0.120% for N2 and 0.100%
for O2, including angular weighting. The maximum equilibrium power imbalance
was 0.0065% of heating plus cooling. High-energy rate/angular interpolation
errors were below 0.103%, with the separately bounded momentum-transfer
approximation below 0.038%. These quantify discretization of the prescribed
model, not experimental accuracy.

``Tools/CrossSections/reciprocal_rotation/README.rst`` gives the offline
export and validation commands. The prepared cumulative reference decodes
the alias probabilities offline and samples identical distributions. The
particle tests exercise heating, cooling, equilibrium, endpoints, recoil,
precision, and restart. ``benchmark_reciprocal_rotation.py`` measures warmed
complete MCC operators and complete PIC timesteps, with repeated seeds and
variance times computational cost. Performance conclusions require those
complete measurements, not only isolated table lookup timings.
Measured workloads, precision limits and performance tradeoffs are documented
in ``Tools/CrossSections/reciprocal_rotation/VALIDATION.rst``.

Other table conversion issues
-----------------------------

* O2's longest-band threshold is 9.97 eV, not the 9.7 eV header.
* The two numerical O2 SR tables have an artificial zero at their 10 eV join.
  Regeneration must use the continuous analytical rate.
* Fixed excitation losses do not represent the vibronic and continuum loss
  distributions in the thesis.
* Source interpolation and continuation must be explicit. Ordinary WarpX
  tables use linear interpolation and endpoint clamping; nonzero source
  endpoints cannot be replaced by zero without a physical continuation.
* elmolcs includes dissociative O2 attachment but no trusted three-body model.
  The supplied attachment data replace, rather than supplement, the overlapping
  dissociative channel. Three-body m5 data require the appropriate stabilizer.

The original PDFs and archive remain provenance inputs, not redistributed
runtime dependencies. Full electron-air export follows the source and physics
checks; an unvalidated rotational bundle is not a beam-production data set.
