#!/usr/bin/env python3
"""Exercise bounded product creation with the O2 production table."""

import runpy
import sys
from pathlib import Path

sys.argv += ["--target", "O2"]
runpy.run_path(
    str(
        Path(__file__)
        .resolve()
        .with_name("inputs_test_3d_proton_impact_ionization_performance_picmi.py")
    ),
    run_name="__main__",
)
