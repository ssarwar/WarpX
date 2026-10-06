Reciprocal rotational scattering
===============================

These are offline reference and export tools for a fixed N2 or O2 rotational
bath. They implement the model described in ``IMPLEMENTATION.rst``; source
qualifications are in ``SOURCES.rst``. They are not imported by WarpX.

Install NumPy and SciPy in an environment containing the supplied elmolcs
snapshot. Point ``WARPX_CROSS_SECTION_DATA`` to ``MCC_cross_sections`` in the
warpx-data checkout. Numerical source constraints live there, alongside the
production cross sections, rather than in the WarpX repository. Set
``WARPX_ROTATION_OUTPUT`` to a build directory for intermediate calculations::

    export PYTHONPATH=/path/to/elmolcs/src
    export WARPX_CROSS_SECTION_DATA=/path/to/warpx-data/MCC_cross_sections
    export WARPX_ROTATION_OUTPUT="$PWD/build/reciprocal-rotation"
    python Tools/CrossSections/reciprocal_rotation/source_anchor_checks.py
    python Tools/CrossSections/reciprocal_rotation/build_reference.py --target N2
    python Tools/CrossSections/reciprocal_rotation/build_reference.py --target O2

The reference uses common forward/reverse primitives and explicit cold and
high-energy continuations. Refinement solves additional source rows before
testing interpolation, instead of comparing a coarse table with itself.
Production exports must also validate their packed sampling distributions.

The supported runtime domain is 0--1 GeV. Source curves, modeled continuations,
and the internal reverse-evaluation buffer have distinct meanings. A final
nonzero source entry does not authorize constant extrapolation. An appended
zero does not establish a physical cutoff. Source-only rotational tables end
at 1 keV for N2 and 20 eV for O2; the combined model supplies the continuation.

Generated validation reports, benchmarks, temporary temperature bundles, and
reference checkpoints belong in build directories. Only scientific source
inputs and validated production bundles are exported to warpx-data.

The exporter includes small quantile lookup arrays to shorten exact searches
of the adaptive CDF grids. ``lookup_index.py`` can append these arrays to an
existing V6 bundle without recomputing or changing its physical arrays. Compare
indexed and unindexed alias data with ``benchmark_reciprocal_rotation.py
--reference-dir ... --samplings alias reference``. The ``lookup_check=1``
sampler test requires event-by-event equality for identical random uniforms.
GPU particle sorting can change ordering between separate PIC runs; those
runs need statistical comparisons.

After building the source reference, prepare and verify each gas separately::

    python Tools/CrossSections/reciprocal_rotation/test_transition_normalization.py
    python Tools/CrossSections/reciprocal_rotation/adaptive_low_grid.py --target N2
    python Tools/CrossSections/reciprocal_rotation/export.py --target N2 \
        --reference-dir "$WARPX_ROTATION_OUTPUT" \
        --output build/production/N2/IAA/reciprocal_hybrid_300K
    python Tools/CrossSections/reciprocal_rotation/verify.py \
        --reference-dir "$WARPX_ROTATION_OUTPUT" \
        --bundle build/production/N2/IAA/reciprocal_hybrid_300K/thermal_rotation.rot
    python Tools/CrossSections/reciprocal_rotation/verify_high.py \
        --bundle build/production/N2/IAA/reciprocal_hybrid_300K/thermal_rotation.rot \
        --output build/reciprocal-rotation/high-N2.json

Repeat with O2, then copy only the validated index and binary parts to
warpx-data. ``verify.py`` independently evaluates intermediate source rows
and decodes packed probabilities; it is intentionally more expensive than
quick runtime tests. Never reuse a sampling grid built from an older source
reference version. Use ``cumulative_reference.py`` to prepare an identical
distribution with cumulative sampling outside warpx-data::

    python Tools/CrossSections/reciprocal_rotation/cumulative_reference.py \
        --bundle build/production/N2/IAA/reciprocal_hybrid_300K/thermal_rotation.rot \
        --output build/cumulative/N2/IAA/reciprocal_hybrid_300K

Configure runtime checks with
``-DWarpX_RECIPROCAL_TEST_DATA=/path/to/warpx-data/MCC_cross_sections``.
Build ``test_reciprocal_rotation``, ``test_rotation_thresholds``,
``benchmark_reciprocal_mcc``, ``pyWarpX_1d``, ``pyWarpX_python_sources``, and
``pyAMReX_python_sources``, then run
``ctest --test-dir build -R 'test_1d_(reciprocal|rotation_thresholds)'``.
The sampler and independent-moment tests require NumPy, but do not import
elmolcs or regenerate scientific data. Each quick thermal test uses 8192
electrons and completes within a few seconds on the development CPU.

For the longer thermal campaign, prepare the additional temperatures offline
with ``prepare_temperature_tests.py`` and run
``Examples/Tests/collision/analysis_reciprocal_thermal.py`` with
``--temperature-data build/reciprocal-temperatures`` and ``--integration``.
The bounded 20 eV temporary fixtures exercise real source kernels at
0, 100, 250, 350 and 1000 K. They are test data, not published beam bundles.

``Tools/CrossSections/perlmutter_reciprocal.sbatch`` requires an explicit
GPU allocation and the prepared production/cumulative data paths. It selects
``gpu&hbm40g``, verifies the allocated device capacity, and compiles for
``sm_80``. It runs independent sampler checks, Compute Sanitizer, thermal
tests, and alternating alias/cumulative benchmarks of full MCC operators
and PIC timesteps. Override the task/GPU count for a four-GPU MPI run, using
a separate build directory within the checkout. Single-precision checks use
``WARPX_ROTATION_PARTICLE_PRECISION=SINGLE`` and optionally
``WARPX_ROTATION_REAL_PRECISION=SINGLE``. HIP/SYCL use the same portable
sampler; compilation and execution must be reported only when those
backends are actually available.

The performance report in ``VALIDATION.rst`` includes two particle layouts.
To reproduce the less concentrated layout, pass ``--cells-per-rank 4096
--max-grid-size 1024`` to the benchmark driver. The default 128-cell layout
has much longer per-cell linked lists in the CUDA deposition sorter.
