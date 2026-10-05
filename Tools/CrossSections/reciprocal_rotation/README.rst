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
