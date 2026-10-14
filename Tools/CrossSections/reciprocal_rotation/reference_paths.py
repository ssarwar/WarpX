# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Keep generated data and validation output in a build directory."""

import os
from pathlib import Path

OUTPUT = Path(
    os.environ.get("WARPX_ROTATION_OUTPUT", "build/rotation-reference")
).resolve()
OUTPUT.mkdir(parents=True, exist_ok=True)


def output_path(name):
    return OUTPUT / name
