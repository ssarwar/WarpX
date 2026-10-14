.. _beam-air-development:

Beam and neutral-gas model development
======================================

These pages cover the beam--air additions to WarpX: prescribed beam and ion
fluids, proton-impact production, electron MCC, IAA physical data, diagnostics
and restart. The user entry point is :ref:`usage-beam-air`. Input parameters
remain in :ref:`running-cpp-parameters`; theory is grouped by physical model.

.. toctree::
   :maxdepth: 1

   runtime
   data_preparation
   validation

Feature and documentation map
-----------------------------

.. list-table:: Maintained specification and independent checks
   :header-rows: 1
   :widths: 25 40 35
   :class: rotational-scattering-table

   * - Feature
     - Theory, algorithm and usage
     - Test location
   * - Relativistic electron MCC, majorants and coupled scheduling
     - :ref:`multiphysics-collisions-mcc`, :ref:`beam-air-runtime`
     - ``Tools/Algorithms/BackgroundMCC``;
       collision selector, timing and coupled-subcycle tests
   * - Elastic/rotational energy exchange and alternatives
     - :ref:`theory-rotational-scattering`, :ref:`rotation-alternatives`
     - ``Tools/CrossSections/reciprocal_rotation``;
       reciprocal/thermal collision tests
   * - Vibrational/electronic excitation, dissociation and attachment
     - :ref:`theory-iaa-air`, :ref:`usage-collision-data`
     - ``Tools/CrossSections/test_iaa_air.py``, angular sampler and
       analytical attachment/cooling tests
   * - Electron ionization and secondary angles
     - :ref:`electron-impact-ionization`
     - Independent RBEQ quadrature, collision energy/momentum and product tests
   * - Calibrated proton/bare-ion spectrum and quiet products
     - :ref:`multiphysics-collisions-proton-impact-ionization`,
       :ref:`proton-ionization-calibration`
     - ``Tools/Algorithms/ProtonImpactIonization``;
       mixed-parent, yield, spectrum and budget tests
   * - Gaussian pulse trains, immobile ion charge and solver coupling
     - :ref:`theory-prescribed-fluids`, :ref:`usage-beam-air`
     - ``Tools/Algorithms/PrescribedFluids``;
       RZ fluid, boundary, source and implicit tests
   * - Fluid-aware diagnostics and persistent restart state
     - :ref:`usage-beam-air`, :ref:`running-cpp-parameters-fluids`,
       the reduced-diagnostic input reference
     - ``Examples/Tests/reduced_diags`` and collision restart tests;
       distributed prescribed-fluid comparisons
   * - Backend and integration fixes
     - :ref:`beam-air-validation`, :ref:`developers-portability`
     - RZ spectral/implicit controls, array-backend and restart checks

The map was checked against all branch changes since upstream ancestor
``58137895feb8f498a6a46c1491c253f4c571e76f``. It groups related commits by
feature rather than creating a documentation page for every historical patch.
Integration fixes include RZ PSATD centering, HIP inverse-FFT workspace sizing,
SYCL queue integration, semi-implicit rollback, collision timing, NumPy 2
buffer conversion and dimensionless Twiss alpha. They retain their tests and
archived evidence; a portability fix or compile result is not an accelerator
physics/performance validation.

Keep model equations, source decisions, reader contracts and user settings in
their owning pages. Link to the commands and evidence beside the tools instead
of copying chronological audit logs into multiple theory pages. Preserve
section/equation labels when moving prose. Physical model changes must update
theory, source provenance, production data and independent tests together.
