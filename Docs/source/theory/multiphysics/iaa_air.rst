.. _theory-iaa-air:

IAA electron--air collision data
===============================

The production IAA manifest in ``warpx-data/MCC_cross_sections/IAA/air``
combines the existing elastic/rotation families with source-resolved
vibrational excitation, electronic excitation, residual dissociation,
RBEQ ionization and oxygen attachment. Neutral internal states and neutral
fragments are not evolved. Each excitation retains its own energy loss.

The default is the completed IAA set, including complementary excitations.
Aggregate scattering, momentum-transfer, excitation and dissociation tables
are not additional collision channels. Nitrogen's residual dissociation
already excludes predissociation; oxygen's dissociative electronic bands
are counted once. Adding complements follows the source's completion
convention and does not change the existing elastic/rotation family.

Data and regeneration
---------------------

``Tools/CrossSections/export_iaa_air.py`` reads the pinned elmolcs archive
without executing package code. NumPy is needed only by the offline exporter;
WarpX does not generate physical tables during startup or time stepping.

.. code-block:: bash

   python Tools/CrossSections/export_iaa_air.py \
       --archive /path/to/elmolcs-main.tar.gz \
       --attachment /path/to/O2_attachment.txt \
       --data-dir /path/to/warpx-data/MCC_cross_sections

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

N2 electronic profiles use Khakoo 2005/2008 and Malone 2009, with canonical
state mappings. O2 uses Shyn--Sweeney 1993 for a/b, Shyn--Sweeney 2000 for
Herzberg, and Suzuki 2011 for the other three bands. The per-profile provenance
records the actual usable energy rows. Narrow forward-cone measurements are
retained as references, not promoted to complete angular distributions.
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
channel loss. Mixture weights use :math:`3s^2-2s^3`, where :math:`s` is the
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
   v01.txt
   ../../elastic_dcs.txt
   0.00001 112.5

The last line is the inclusive active energy interval. Outside it, the shared
elastic model is sampled. The profile file uses the existing uniform-theta DCS
layout but omits legacy ``SPECIES:`` metadata. Its fallback must be a legacy
elastic table of the specified molecule, not another descriptor. The explicit
fallback includes the existing elastic high-energy continuation.

Attachment and PICMI
--------------------

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

.. code-block:: python

   from pywarpx.iaa_air import load_iaa_air_processes

   processes = load_iaa_air_processes(
       data_dir, "O2", ion_species=O2plus,
       attachment_species={"dissociative": Ominus, "three_body": O2minus},
       third_body_density=n_O2,
   )
   collision = picmi.MCCCollisions(
       name="electron_O2", species=electrons,
       background_density=n_O2, background_temperature=300,
       background_mass=31.9988 * 1.66053906660e-27,
       scattering_processes=processes,
   )

The helper reads only the production manifest and uses the supplied 300 K
rotational bundle. Other rotational temperatures require independently prepared
families. The RBEQ alternative selects both its total and secondary model.
No neutral-state species are introduced.

Verification
------------

Set ``WarpX_IAA_AIR_TEST_DATA`` to the external ``MCC_cross_sections`` directory
to enable production source and angular checks. Analytical angular tests do not
need external data. They check isotropy, exact d-pi moments, mixtures, finite-step
Lorentz-gas relaxation and preservation of very small excitation deflections.
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
