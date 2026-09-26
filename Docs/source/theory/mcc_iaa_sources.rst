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
the sudden approximation. The V3 rate model uses the canonical internal-energy
threshold :math:`\Delta=E_{J'}-E_J`, without molecular recoil shifts. For
relativistic electron momentum :math:`p(E)=\sqrt{E(E+2m_ec^2)}/c`, reverse rates
satisfy the heavy-target integral relation

.. math::

  g_Jp(E+\Delta)^2\sigma_{J\to J'}(E+\Delta)
  =g_{J'}p(E)^2\sigma_{J'\to J}(E).

The runtime draws a discrete rotational change independently of the angle from
the existing inclusive elastic DCS. At exactly zero relative momentum there
is no incident axis, so emitted electrons use the isotropic limit. A selected
excitation below its exact level spacing, possible through energy-grid
rounding, becomes an unchanged event. No other rotational rate is increased.
There is no angular conditioning, transition search, or rejection in the
alias sampler, and it requires no cumulative table.

The spectator relation in Eq. (11.29), :math:`J_R=kR\sin(\theta/2)`, is not a
hard quantum selection rule. Neither it nor a separate rotational DCS is used.
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

Integral detailed balance in the reference rates does not imply differential
detailed balance for an energy-dependent elastic DCS used independently of
rotational outcomes. This is the selected approximation. Independent tests
cover level thresholds, unchanged/loss/gain probabilities, angular independence,
signed recoil, the continuation bound and thermal power balance.

Physical bundles are generated offline and kept in ``warpx-data``. Runtime
initialization only validates loaded data, weights equilibrium populations,
prepares shared in-memory samplers and uploads them to the device. It neither
evaluates elmolcs source models nor creates cross-section files.

The zero-reference-temperature candidate audit has nonnegative unchanged
integral rates. This is an explicit reference-population assumption, not a
demonstrated temperature for all source measurements. A 300 K diagnostic
exposes negative low-energy residuals with the current extrapolations. Finite
inclusive cross sections imply :math:`v\sigma\to0`, which cannot contain the
finite superelastic contribution of a thermally populated rotor. Validated
thermal low-energy continuations, omitted-rank convergence and coverage through
beam energies remain production requirements. O2's angular-integrated Born
rate approximation retains its limited validity; the elastic angular sampler
does not extend it. Synthetic V3 verification bundles establish software
correctness, not N2/O2 beam-physics validity.

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
