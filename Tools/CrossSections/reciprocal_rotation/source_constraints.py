# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Published constraints; see SOURCES.rst for provenance and limitations."""

import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(os.environ["WARPX_CROSS_SECTION_DATA"])
N2 = ROOT / "N2/IAA/reciprocal_sources"
O2 = ROOT / "O2/IAA/reciprocal_sources"
GOTE = {float(k): v for k, v in json.loads((N2 / "n2_gote.json").read_text()).items()}
JUNG = json.loads((N2 / "n2_jung.json").read_text())
B = json.loads((O2 / "o2_bhattacharyya.json").read_text())
ENERGY = np.asarray(B["energy_eV"], float)
TOTAL = np.asarray(B["integrals_a02"], float)
MOMENT = np.asarray(B["momentum_transfer_a02"], float)
