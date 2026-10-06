.. _rotation-development:

Continuing development and reproducing results
==============================================

This is the entry point for continuing this work in a new checkout or task.
The physics specification and numerical inputs are maintained in the
repositories; no chat transcript or temporary build directory is required to
understand the model. Local build outputs are useful caches, not authoritative
scientific inputs.

Authoritative material
----------------------

* :ref:`rotation-dcs` is the complete piecewise N2/O2 rotational DCS
  specification: raw source strengths, energy and angular joins, reciprocal
  state-to-state DCS and the channel probabilities used by MCC.
* :ref:`rotation-sources` identifies primary papers, processed measurements,
  source discrepancies and inferred quantities. It distinguishes experimental
  coverage from the energy intervals chosen for interpolation or continuation.
* :ref:`rotation-reciprocity` specifies the normalization equation and its
  implicit reverse term. :ref:`rotation-continuations` specifies cold and
  high-energy limits, bound rotor states and supported energies.
* :ref:`rotation-algorithm` maps equations to host preparation, GPU sampling,
  MCC selection and signed recoil, and lists the responsible source files.
* :ref:`rotation-data-format` specifies readable data and initialization.
  Users select gas, physical model and rotational temperature. WarpX
  identifies the storage layout automatically.
* :ref:`rotation-validation` identifies tests, numerical error budgets,
  measured CPU/GPU behavior and what has not been established physically.

Repository locations and working branches
-----------------------------------------

The implementation is maintained on WarpX's
``codex/rigid-beam-immobile-ions-development-sync`` branch; production numerical
inputs are on warpx-data's ``codex/elmolcs-elastic-rotation`` branch. Before
editing, check each checkout's status and fetch its requested branch. Preserve
local work and other active worktrees. The physical data directories are

.. code-block:: text

   warpx-data/MCC_cross_sections/N2/IAA/
   warpx-data/MCC_cross_sections/O2/IAA/

Each contains ``elastic.txt``, ``elastic_dcs.txt``, short elementary rotational
source curves, ``reciprocal_sources/`` and
``reciprocal_hybrid_300K/thermal_rotation.rot`` with its readable arrays.
Source inputs and production numerical data belong in warpx-data. Exporters,
physics algorithms and tests belong in WarpX. Synthetic fixtures, reference
checkpoints, fit diagnostics, validation records and benchmark outputs belong
in a build directory, not in warpx-data.

The maintained preparation entry points are in
``Tools/CrossSections/reciprocal_rotation/README.rst``. A source environment
needs NumPy, SciPy, Numba and the supplied elmolcs package. Its original
``Data`` directory must remain accessible to elmolcs. The supplied research
resources include the thesis and the Read--Andrick, Jung, Morrison, Gote and
Bhattacharyya papers. Their filenames and source tables are identified in the
source documentation; do not replace them with an arbitrary newer package's
default orbital or rotational data without a fresh consistency audit.

The supplied resource directory uses these filenames (it is outside the source
repository; locate it before running source-level work):

.. list-table:: Source resource inventory
   :header-rows: 1
   :class: rotational-scattering-table
   :widths: 35 65

   * - Resource
     - Supplied filename or directory
   * - IAA thesis
     - ``tesis_iaa_schmalzried_anthony.pdf``
   * - Read--Andrick resonance theory
     - ``read_angular_distributions_rotation.pdf``
   * - Jung experimental branches
     - ``jung_rotational_excitations.pdf``
   * - Morrison low-energy theory
     - ``morrison_rotational_low_energies.pdf``
   * - Gote intermediate-energy measurements
     - ``Gote_rotational_excitation_medium_energies.pdf``
   * - Bhattacharyya oxygen calculation
     - ``bhattacharyya_rotational_oxygen.pdf``
   * - elmolcs implementation and original tables
     - ``elmolcs/src`` and ``elmolcs/Data``

Kutz--Meyer elementary curves are supplied through elmolcs and cited in the
thesis and :ref:`rotation-sources`; consult the paper when revisiting their
interpretation. Numerical paper constraints in warpx-data are sufficient for
normal reproduction with the supplied package. Missing primary resources
should be identified explicitly, rather than replaced by guessed source values.

.. code-block:: bash

   export PYTHONPATH=/path/to/elmolcs/src
   export WARPX_CROSS_SECTION_DATA=/path/to/warpx-data/MCC_cross_sections
   export WARPX_ROTATION_OUTPUT="$PWD/build/reciprocal-rotation"

Preparation and validation sequence
-----------------------------------

1. Reproduce ``prepare_jung.py`` from the readable digitization and run
   ``source_anchor_checks.py``. Source constraints are scientific inputs;
   generated fit/check reports remain in the build directory.
2. For a source-model or temperature change, run ``build_reference.py`` and
   ``test_transition_normalization.py``. Rebuild the reference rather than
   reusing checkpoints generated under a different physical normalization.
3. Run ``adaptive_low_grid.py`` and ``export.py`` for each gas. They construct
   the supported energy domain, threshold rows, angular quantiles, outcome
   palette and high-energy momentum-transfer distribution.
4. Run ``verify.py`` against independently evaluated intermediate source
   energies and ``verify_high.py`` against the angular/rate continuations.
   Use the separate q-bin and phase-average checks when changing those parts.
5. Use ``text_bundle.py`` to encode the validated distributions as readable
   production inputs. It changes their numerical representation, not the
   physical source model. Retain the intermediate binary reference only in
   the build directory for comparison.
6. Build the C++ sampler and run it with ``cell_output=<path>``. Compare every
   native decoded cell using ``check_text_cells.py --reference <index>
   --decoded <path>``. This checks probabilities and signed energy moments
   after float32 alias packing, including rare tails. Run stochastic
   ``sampling_reference.py`` checks as well.
7. Run full MCC thermal, endpoint, precision and integration tests before
   replacing production data. Regenerating data is not itself validation.
   Keep rate, first/second transfer-moment errors below 0.2% and equilibrium
   power imbalance below 0.1% of heating plus cooling.

For a documentation or file-identifier-only change, do not rerun hours of
physical fitting. Verify the documented equations against the implementation,
check numerical payload equality where applicable, test readers and relevant
samplers, and rebuild the documentation. For a numerical encoding change,
compare every reconstructed cell and rerun thermal/sampling tests; for a
physical model change, also redo independent source and equilibrium checks.

Local builds and quick runtime tests
------------------------------------

Use the repository's ``AGENTS.md`` and keep build directories inside its root.
The ``warpx-cpu-mpich-dev`` environment supplies the CPU development tools on
the original workstation. With an appropriate environment and existing build,

.. code-block:: bash

   cmake -S . -B build -DWarpX_DIMS=1 -DWarpX_PYTHON=ON \
       -DWarpX_RECIPROCAL_TEST_DATA="$WARPX_CROSS_SECTION_DATA"
   cmake --build build -j 8 --target test_reciprocal_rotation \
       test_rotation_thresholds benchmark_reciprocal_mcc pyWarpX_1d
   ctest --test-dir build \
       -R 'test_1d_(reciprocal|rotation_thresholds)' --output-on-failure

The CMake configuration must also match the local backend, MPI and optional
library availability; see the installation instructions for a fresh build.
Tests include the real production data, gas-specific and air thermal cases,
malformed/range checks and the direct/cached/fallback selectors. All-single
and double-particle builds have different stored-momentum accuracy even when
the sampling tables use double energy labels. Do not relax physical assertions
to make a failing run pass.

GPU validation and performance
------------------------------

``Tools/CrossSections/perlmutter_reciprocal.sbatch`` gives the maintained
Perlmutter campaign. It requires explicit account and data/environment paths,
verifies an A100 with 40 GB memory, and compiles for ``sm_80``. Keep both gases
resident. Compare alias and cumulative sampling with identical physics, and
compare readable and prepared-array inputs when evaluating startup changes.
Measure initialization, host/device memory, full MCC throughput and complete
PIC timestep cost separately. Use warmed, synchronized timings, alternating
execution order, repeated seeds and mixed electron energies within warps.
A scalar lookup benchmark alone does not establish performance of MCC.

The combined resident rotational-table budget is 1 GiB per GPU, using one MPI
rank per GPU. There must be no per-event neutral-state loop, Boltzmann sum,
Bessel evaluation, rejection loop over rotational transitions or additional
particle pass. Scientific preparation may occur once at initialization only
when its cost is measured and remains reasonable; production source fitting
currently remains offline. Check Compute Sanitizer and meaningful CPU/CUDA
physics cases. Compile/run HIP and SYCL when toolchains and hardware are
available, and distinguish unavailable coverage from a successful check.

Assumptions that require a deliberate physics decision
------------------------------------------------------

* The neutral rotational bath is fixed and Boltzmann distributed, including
  adopted nuclear-spin statistical weights; its population, spin-conversion
  kinetics and energy are not evolved.
* The N2 resonance uses inferred Jung fractions and fitted missing-angle
  corrections. Pure Read shapes are not the full vibrationally elastic DCS.
* The O2 1--20 eV bridge is not a measured resonance law. Its completion uses
  clamped source-energy constraints and the actual-energy elastic marginal.
* The paired low-energy reference is reciprocal in the heavy-target model;
  finite tables approximate it. The high-energy spectator bank neglects
  finite-gap reverse-energy corrections and is not exactly reciprocal at
  finite incident energy.
* The high-J bound rigid rotor omits dissociation, non-rigid spectroscopy,
  excited vibration and resolved O2 electronic-spin structure.
* The nominal-threshold recoil continuation has a documented small energy
  defect. Removing an angle veto is not a claim of exact finite-mass
  kinematics inside that band.
* The sub-eV IAA elastic DCS has sparse source knots and model uncertainty.
  Numerical refinement and equilibrium tests do not supply new transport
  measurements.

If one of these assumptions is changed, update the equations, source rationale,
preparation code, production data, numerical budgets and independent tests
together. Preserve ordinary MCC and other physical models when this feature is
disabled. Follow WarpX commenting/documentation style, make normally sized
validated commits, and push to the intended branches. Do not open a pull
request without explicit instruction.
