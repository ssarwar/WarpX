#!/usr/bin/env python3
"""A loose proposal rate must not suppress physical collision probabilities."""

import runpy
import sys
from pathlib import Path

sys.argv += ["--majorant-factor", "10000"]
runpy.run_path(
    str(
        Path(__file__)
        .resolve()
        .with_name("inputs_test_1d_background_mcc_selector_picmi.py")
    ),
    run_name="__main__",
)
