.. _usage-collision-data:

Collision data files and model selection
========================================

Physical electron cross sections are distributed separately in
`warpx-data <https://github.com/ssarwar/warpx-data/tree/IAA/MCC_cross_sections>`__.
Use the ``IAA`` branch and record its commit alongside the WarpX commit.
WarpX reads prepared data; it does not import elmolcs or fit physical models
during a simulation. Source and data licensing remains with those sources.

Directory layout
----------------

.. code-block:: text

   MCC_cross_sections/
     N2/IAA/                         O2/IAA/
       manifest.json                   manifest.json
       elastic/                        elastic/
       rotational/                     rotational/
       vibrational/                    vibrational/
       electronic/                     electronic/
       dissociation/                   dissociation/
       ionization/                     ionization/
                                       attachment/

Production integral tables are directly inside the process directory.
``source_bounded/`` contains source-limited counterparts; ``angular/`` contains
normalized profiles and descriptors. ``sources.json`` records citations,
original labels, units, source hashes, corrections and interpolation choices.
The per-gas ``manifest.json`` is the process ledger and runtime entry point:
canonical channel IDs, rates, energy losses, angular models, continuation and
mutually exclusive groups. Its paths are relative to that gas's ``IAA`` directory.
Original archive paths in provenance identify sources, not runtime files.

Production tables extend to 1 GeV with documented models. Source-bounded tables
stop at their own endpoints; combining them narrows a collision object's valid
energy interval to the intersection of its input ranges. A configured energy
loss, physical threshold, fitted onset and source coverage are different
quantities. Preserve each when replacing a source.

Integral cross-section reader
-----------------------------

Each non-comment row has exactly two numbers: electron energy in eV and cross
section in m2, or raw m5 for three-body attachment. Blank lines, full-line ``#``
comments and trailing ``#`` comments are accepted. At least two finite,
nonnegative, strictly increasing energy knots and nonnegative finite cross
sections are required. Nonuniform grids are supported. Ordinary lookup is
piecewise linear in energy and cross section, even when the offline reference
uses log--log interpolation between positive source points.

Recognized metadata use an exact ``# key = value`` line with one value and no
trailing words. Other comment text is provenance only. Important metadata are:

.. list-table:: Reader metadata
   :header-rows: 1
   :widths: 38 62
   :class: rotational-scattering-table

   * - Keys
     - Contract
   * - ``outside_energy_range``, ``energy_min_eV``, ``energy_max_eV``
     - Supply all three, with ``outside_energy_range = error``. Values outside
       the declared supported range are errors. Unmarked legacy tables clamp
       endpoints; that behavior is not a physical continuation model.
   * - ``native_energy_grid = resample``
     - For single-precision particles, prepare unique representable energy
       knots from the original double-precision interpolant, retaining source
       zeros and their neighboring native energies. This does not guarantee
       the offline interpolation budget in particle precision.
   * - ``rbeq_model``, ``rbeq_target``, ``rbeq_normalization``
     - Identify the matching total/secondary model and normalization.
       RBEQ initialization checks these when supplied. Keep the generated
       metadata with the production ionization tables.

The current native-grid path validates cross sections after resampling.
An invalid negative source value can disappear between rounded knots in single
precision; use the independent production-data checks before loading edited
files. This is a known reader-validation gap, not permission for signed rates.
See :ref:`beam-air-validation` for precision limits and follow-up coverage.

Units are configured explicitly for attachment with
``<process>_cross_section_units``. For m5 data, provide
``<process>_third_body_density`` in m-3; WarpX applies it once. Header prose does
not set units. Likewise ``<process>_energy`` sets an excitation's discrete loss;
a threshold comment alone does not configure that loss.

Angular and rotational inputs
-----------------------------

``<process>_differential_cross_section`` selects either an elmolcs-style DCS
table or a ``WARPX_EXCITATION_ANGULAR_V1`` descriptor. A DCS table has one energy
followed by uniformly spaced polar-angle values from zero to pi; all rows have
the same angular count. It is normalized internally and controls angular shape,
not the event rate. Recognized N2/O2 elastic metadata enable the documented
high-energy IAA continuation. Shared path strings reuse prepared angular data.

An excitation descriptor names the molecule, its prepared normalized profile,
an elastic fallback file, and the inclusive active energy interval. Paths are
relative to the descriptor. :ref:`theory-iaa-air` defines its format, interpolation
and transitions; :ref:`multiphysics-collisions-mcc` defines the runtime sampler.

For ``rotation_model = reciprocal_hybrid``, ``rotation_file`` points to
``rotational/elastic_rotation_300K/thermal_rotation.rot`` and its companion
arrays. The bundle supplies a joint angle/rotational-change distribution and
replaces the ordinary elastic rate. Do not add the elementary rotational source
files as additional collision channels. :ref:`rotation-data-format` specifies
the readable format, legacy binary compatibility, units and initialization.

Production manifest helper
--------------------------

``pywarpx.iaa_air.load_iaa_air_processes`` takes the ``MCC_cross_sections`` root,
``"N2"`` or ``"O2"``, and product species, and returns the PICMI process dictionary.
It accepts both the current ``warpx_iaa_gas_v1`` per-gas manifests and older
``warpx_iaa_air_v1`` aggregate manifests. It does not translate literal paths
in a user's existing input deck; update those after moving to the new layout.

.. list-table:: Exclusive choices accepted by the helper
   :header-rows: 1
   :widths: 35 28 37
   :class: rotational-scattering-table

   * - Argument
     - Default
     - Alternative
   * - ``ionization_model``
     - ``iaa_thesis_2023``
     - ``elmolcs_b8643810``; selects both rate and secondary model
   * - ``dissociative_attachment`` (O2)
     - ``iaa``
     - ``muroran``
   * - ``three_body_attachment`` (O2)
     - ``muroran``
     - ``phelps`` or ``biagi``

O2 additionally requires ``attachment_species`` mapping ``dissociative`` to O-
and ``three_body`` to O2-, and a positive ``third_body_density`` equal to the
O2 stabilizer density. Unknown or missing selections and missing files are
errors. The helper selects the completed excitation set and the 300 K
reciprocal family. It introduces no neutral species. See :ref:`usage-beam-air`
for PICMI configuration and :ref:`beam-air-data-preparation` for regeneration.
