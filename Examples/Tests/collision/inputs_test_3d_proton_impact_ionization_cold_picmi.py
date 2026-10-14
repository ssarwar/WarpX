#!/usr/bin/env python3
"""Preserve small constant and parsed neutral temperatures in both gases."""

import runpy
import sys
from pathlib import Path

sys.argv += ["--cold"]
runpy.run_path(
    str(
        Path(__file__)
        .resolve()
        .with_name("inputs_test_3d_proton_impact_ionization_picmi.py")
    ),
    run_name="__main__",
)
