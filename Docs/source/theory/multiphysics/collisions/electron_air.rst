.. _theory-iaa-air:
.. _mcc-iaa-sources:

IAA electron--air collision data
================================

The gas-specific IAA manifests in
``warpx-data/MCC_cross_sections/<gas>/IAA/manifest.json``
combine the existing elastic/rotation families with source-resolved
vibrational excitation, electronic excitation, residual dissociation,
RBEQ ionization and oxygen attachment. Neutral internal states and neutral
fragments are not evolved. Each excitation retains its own energy loss.

The default is the completed IAA set, including complementary excitations.
Aggregate scattering, momentum-transfer, excitation and dissociation tables
are not additional collision channels. Nitrogen's residual dissociation
already excludes predissociation; oxygen's dissociative electronic bands
are counted once. Adding complements follows the source's completion
convention and does not change the existing elastic/rotation family.

.. list-table:: Production process inventory
   :header-rows: 1
   :widths: 25 37 38
   :class: rotational-scattering-table

   * - Family
     - N2
     - O2
   * - Elastic and rotation
     - One combined reciprocal family
     - One combined reciprocal family
   * - Bound vibration
     - All 15 supplied ground-state transitions
     - All 41 supplied ground-state transitions
   * - Neutral dissociation
     - Residual dissociation, already excluding predissociation
     - Resonant vibrational dissociation
   * - Electronic excitation
     - 17 states and 7/14 eV complements
     - a, b, Herzberg, two Schumann--Runge losses, longest/second bands,
       and 11 eV complement
   * - Ionization
     - One matched RBEQ total and secondary model
     - One matched RBEQ total and secondary model
   * - Attachment
     - None
     - One dissociative and one three-body source

The default dictionaries have 37 N2 and 54 O2 processes, including each combined
elastic family. The ledger also retains exclusive source alternatives; they
are not simultaneous event rates. Ionization shell selection remains internal
to one process. :ref:`electron-impact-ionization` specifies that sampling.

Integral data and continuations
-------------------------------

Each gas has ``elastic``, ``rotational``, ``vibrational``, ``electronic``,
``dissociation`` and ``ionization`` directories, with ``attachment`` for O2.
Production integral tables are directly inside each family. ``source_bounded``
contains the source-limited counterparts, while ``angular`` contains normalized
profiles and their reader descriptors. Molecular-state and transition names are
explicit in filenames. Table paths in each gas manifest are relative to its ``IAA``
directory. The PICMI loader also accepts the older aggregate manifest layout.

The reader and manifest interfaces are documented in :ref:`usage-collision-data`;
offline reproduction is in :ref:`beam-air-data-preparation`.
``Tools/CrossSections/export_iaa_air.py`` reads the pinned elmolcs archive
without executing package code. NumPy is needed only by the offline exporter.

The source archive SHA256 is
``b864381086120fff37eb8c658dfcb0f150b80969cf9add626a38d89bab0a3c38``.
The exporter requires the existing elastic/rotation and matching RBEQ data.
It writes physical source records, the process ledger, source-bounded tables,
and separately identified production continuations through 1 GeV.
Archive numerical evaluations independently check the analytic excitation
equations. The N2 v=1 source points at 20 and 30 eV are restored. The two
Schumann--Runge loss channels partition the analytic total over
9.955581--10.213826 eV, avoiding the numerical source's artificial zero.

Integral files contain energy in eV and cross section in m2 (m5 for three-body
attachment). Full-line and trailing ``#`` comments are accepted. The adaptive
master grid targets 0.05% interpolation error with a 1e-6-of-peak floor;
individual files retain source knots and omit unnecessary master-grid points.
RBEQ retains its separate 0.02% target and source-consistent secondary model.
Files are bounded above; extrapolation beyond the declared range is an error.

Very narrow threshold regions can need energy knots closer than binary32
spacing. Files marked ``# native_energy_grid = resample`` are interpolated
onto unique representable particle energies during initialization. This
operation uses the original binary64 coordinates and unscaled values, including
m5 attachment data. Source zeros are retained at native coordinates and
bracketed by neighboring representable energies, preserving thresholds and
zero tails. A channel that disappears in particle precision is rejected.
Unmarked legacy files still require a strictly increasing
grid in particle precision. The metadata does not permit unordered or duplicate
source coordinates.

Source-bounded tables stop at their supported endpoints. Production electronic
and ionization tables use the documented analytic models. Positive-ended
resonant tables have endpoint-matched declining power-law continuations fitted
to the final decreasing segment; the N2 14 eV complement uses its dipole-forbidden
asymptote. Documented zero tails remain zero. The continuation is a model
assumption, not an extension of experimental coverage.

Tabulated positive source values interpolate log--log; intervals touching zero
interpolate linearly. Analytic electronic fits are evaluated directly from the
source ``genA``/``genF`` equations and their archived coefficients, not from a
coarse plot. The canonical numerical ``iaa`` inventory supplies bound vibrations;
``iaa*`` alone does not contain every transition. ``sources.json`` retains the
original equation variant, threshold and coefficients, and
``Tools/CrossSections/iaa_air/rates.py`` implements the variants explicitly.

For a positive-ended resonant channel, the production continuation is

.. math::

   \sigma(E)=\sigma(E_h)(E/E_h)^p,\qquad E>E_h,\quad p<0.

The fit uses at least three points from the terminal decreasing segment, at
most its last five points, and matches the endpoint exactly. The fit interval
and exponent are recorded. N2's 14 eV complement fixes :math:`p=-1` from its
dipole-forbidden asymptote. An unsupported or increasing tail is rejected.
Source-bounded files retain the source endpoint instead.

For Schumann--Runge, with :math:`a=9.955581` eV and :math:`b=10.213826` eV,
the partition uses :math:`w=\operatorname{clip}[(E-a)/(b-a),0,1]`:
:math:`\sigma_{6.12}=(1-w)\sigma_{\rm SR}` and
:math:`\sigma_{8.60}=w\sigma_{\rm SR}`. It preserves one analytic total and two
different losses. The longest-band loss is 9.97 eV, correcting the archive's
9.7 eV header. No second fragmentation event is added after a dissociative band.

Angular distributions
---------------------

For channel :math:`j`, the differential cross section is

.. math::

   \frac{d\sigma_j}{d\Omega} = \sigma_j(E) F_j(E,\theta),
   \qquad \int F_j\,d\Omega=1.

The angle uses incoming electron energy in the neutral rest frame and a uniform
azimuth. Borrowing the elastic angular shape does not invoke an additional
elastic or rotational collision. Its normalization is the integral of the DCS,
not the residual elastic integral cross section.

The hierarchy is the explicit thesis prescription, resolved experimental
inelastic DCS, then normalized elastic IAA. For all N2 bound vibrational
channels below 5 eV, thesis section 11.3.2 gives

.. math::

   p(\mu)=\frac{15}{28}\left(1-3\mu^2+\frac{14}{3}\mu^4\right),
   \qquad \mu=\cos\theta.

From 5 to 10 eV this joins the N2 v=0->1 measured distribution, or the elastic
fallback for other vibrational channels. O2 bound vibrations use isotropy
through 20 eV and join elastic IAA over 20--40 eV. These O2 boundaries are an
adopted modeling choice; the thesis does not supply numerical boundaries.
Vibrational dissociation uses elastic fallback.

For N2 v=0->1 above 5 eV, the measured reference uses Linert--Zubek at 5, 10,
15 and 20 eV and Truhlar at 30, 50 and 75 eV. Sparse 35/45 eV and isolated
higher-vibration measurements remain validation references. Higher N2
vibrations join elastic fallback over 5--10 eV where channel data are absent.

N2 electronic profiles use Khakoo 2005/2008 and Malone 2009, with canonical
state mappings. O2 uses Shyn--Sweeney 1993 for a/b, Shyn--Sweeney 2000 for
Herzberg, and Suzuki 2011 for the other three bands. The per-profile provenance
records the actual usable energy rows. Narrow forward-cone measurements are
retained as references, not promoted to complete angular distributions.

.. list-table:: Measured electronic angular intervals before boundary blending
   :header-rows: 1
   :widths: 47 28 25
   :class: rotational-scattering-table

   * - Channels
     - Source
     - Incident energy
   * - N2 A, B, W
     - Khakoo 2005
     - 10--100 eV
   * - N2 B', a', a, w
     - Khakoo 2005
     - 12.5--100 eV
   * - N2 C, E, a''
     - Malone 2009
     - 13--100 eV
   * - N2 b, c3, o3, b', c4', G, F
     - Khakoo 2008
     - 17.5--50 eV
   * - O2 a and b
     - Shyn--Sweeney 1993
     - 5--20 eV
   * - O2 Herzberg
     - Shyn--Sweeney 2000
     - 10--30 eV
   * - O2 Schumann--Runge, longest and second bands
     - Suzuki 2011
     - 15--50 eV

Canonical states, original header labels and corrections are retained together.
In particular, Malone's E and a'' headers map to :math:`E^3\Sigma_g^+` and
:math:`a''{}^1\Sigma_g^+`. Unknown/ambiguous mappings are rejected.
Schumann--Runge losses share one angular distribution. Aggregate angular data
do not add an event rate. The N2 Rydberg 100 eV and O2 100--200 eV rows, and
Green's forward Herzberg measurements, remain partial-angle references.

Two O2 b-state archive entries are corrected against Table I(b), p. 1007 of
`Shyn and Sweeney (1993) <https://doi.org/10.1103/PhysRevA.47.1006>`__:
94 becomes 9.4 at 5 eV and 12 degrees, and 89 becomes 8.9 at 15 eV and
36 degrees, in units of 1e-20 cm2/sr. The rendered paper confirms the decimal
points, which are also missing from its text layer. Provenance retains the
archived values, corrections, citation and verified PDF hash. The corrections
affect the normalized angular profile; the integral rate and loss are unchanged.

Measured angles are interpolated linearly in sin(theta/2). Missing endpoints
use boundary-matched elastic shapes. Accepted rows require eight angles,
75% solid-angle coverage, gaps no larger than 20 degrees and less than 50%
modeled endpoint probability. Each completed row is normalized before linear
incident-energy interpolation. Endpoint probabilities and completion sensitivity
are recorded independently of numerical errors.

Outside a measured interval, boundary profiles join energy-dependent elastic
IAA over one neighboring logarithmic energy interval, limited below by the
configured channel loss :math:`Q`. For usable energies
:math:`E_0<E_1<\cdots<E_n`, the transition ends are
:math:`E_L=\max(Q,E_0^2/E_1)` and :math:`E_H=E_n^2/E_{n-1}`. At least two
usable rows are required. Freeze the nearest measured boundary profile while
evaluating the elastic component at the actual incident energy. Use pure
elastic outside those finite transition intervals. The explicit vibrational
intervals above override this generic rule.
Mixture weights use :math:`3s^2-2s^3`, where :math:`s` is the
fractional logarithmic energy. Shared profiles are tabulated offline against
the existing runtime log-energy interpolation, with CDF and first-four-Legendre
moment errors below 1e-3. Startup prepares inverse CDFs; device sampling uses
bounded lookups and one conditioned uniform draw. It does not interpolate
sampled angles between distributions.

The existing ``differential_cross_section`` parameter also accepts this text
descriptor; paths are relative to the descriptor:

.. code-block:: text

   WARPX_EXCITATION_ANGULAR_V1
   N2
   v0_to_v01_pdf.txt
   ../../elastic/differential_cross_section.txt
   0.00001 112.5

The last line is the inclusive active energy interval. Outside it, the shared
elastic model is sampled. The profile file uses the existing uniform-theta DCS
layout but omits legacy ``SPECIES:`` metadata. Its fallback must be a legacy
elastic table of the specified molecule, not another descriptor. The explicit
fallback includes the existing elastic high-energy continuation.

Attachment
----------

The default dissociative attachment is IAA; MuroranIT is an alternative.
Three-body alternatives are MuroranIT (default), Phelps and Biagi. Their
exclusive groups prevent adding multiple representations of one reaction.
Phelps/Biagi reference-density values are multiplied by 1e-6 to obtain m5.
The supplied MuroranIT m5 values are preserved. Their peak differs by 0.0177724%
from the nominal normalization in Kawaguchi et al. 2025, section 2.6; this source
discrepancy is reported rather than removed by rescaling.

For the default third-body convention,
:math:`\nu_3=n_{O_2}^2 v\sigma_3`. Do not use total air density for the second
density factor. These effective rates approximate stabilization competing with
autodetachment; they are not an exact pressure-dependent resonance model.

The manifest helper and exclusive-source arguments are specified in
:ref:`usage-collision-data`, with a PICMI example in :ref:`usage-beam-air`.
Attachment has no outgoing free-electron angle. Its effective negative ion
does not resolve a dissociation-fragment spectrum.

Verification
------------

Set ``WarpX_IAA_AIR_TEST_DATA`` to the external ``MCC_cross_sections`` directory
to enable production source and angular checks. Analytical angular tests do not
need external data. They check isotropy, exact d-pi moments, mixtures, finite-step
angular-kernel relaxation and an axial very-small-deflection case. These do not
establish complete-MCC relaxation or arbitrary-direction stored-momentum accuracy;
:ref:`beam-air-validation` records the wider audit and remaining gaps.
Synthetic fixtures are generated in the build directory.

``inputs_test_1d_iaa_air_picmi.py`` runs pure-gas and dry-air cooling, a driven
swarm and a homogeneous two-pulse 800 MeV proton collision regression. The
reduced pulse case tests coupled particle production and electron collisions;
it does not reproduce the experiment's transverse geometry or beam propagation.
CPU results do not establish CUDA, HIP or SYCL performance. Use the same portable
tests and benchmark workload on each supported accelerator.

``benchmark_iaa_air.py`` compares compact tables with redundant dense tables
representing identical physical interpolants. It records startup, selector and
angular memory, complete timestep timings, and ensemble variance per cost.
``audit_iaa_air.py`` reports completed-set totals, angular completion and
transition-width sensitivity, and Maxwellian rate sensitivity to uncertain
tails. These diagnostics are written outside the production data repository.

Only production physical data, provenance and usable manifests belong in
``warpx-data``. Exporters, synthetic fixtures, sensitivity experiments, audits,
test outputs and benchmarks belong in WarpX or its build directory.
