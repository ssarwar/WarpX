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

The secondary-angle prescriptions in Eqs. (11.132) and (2.66) are different.
For 2.5 MeV impact, binding 15.58 eV and equal sharing, their mean cosines
are approximately 0.7071 and 0.8804. The primary retains Eq. (2.60).
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

IAA Eq. (11.24) constructs transitions from elementary rates using the sudden
approximation. Eq. (11.35) supplies a spectator angular shape, which the thesis
explicitly identifies as inaccurate at low energy. The long-range Born model,
Eq. (11.21b), is an alternative low-energy approximation, not a validated
description of the N2 resonance. A production bundle needs a stated validity
range and a nonnegative inclusive angular decomposition. Failed decompositions
must be reported, not clipped. Available elementary data alone do not establish
rotational coverage through 1 GeV.

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
