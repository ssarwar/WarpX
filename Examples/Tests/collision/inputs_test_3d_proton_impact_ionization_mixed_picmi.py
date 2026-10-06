#!/usr/bin/env python3
"""Exercise conditional energy sampling for mixed parents in each gas."""

import runpy
import sys
from pathlib import Path

sys.argv += ["--mixed"]
runpy.run_path(
    str(
        Path(__file__)
        .resolve()
        .with_name("inputs_test_3d_proton_impact_ionization_picmi.py")
    ),
    run_name="__main__",
)
