#!/usr/bin/env python3
"""Exercise the molecular tail at the lower supported proton energy."""

import runpy
import sys
from pathlib import Path

# The calibration boundary is an incident neutral-rest-frame energy. Use a
# stationary target here; moving-target coverage is exercised inside the range.
sys.argv += ["--energy-keV", "5", "--steps", "2", "--cold"]
runpy.run_path(
    str(
        Path(__file__)
        .resolve()
        .with_name("inputs_test_3d_proton_impact_ionization_picmi.py")
    ),
    run_name="__main__",
)
