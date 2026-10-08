# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check exported totals against independent differential quadrature.

Usage: python analysis_rbeq_sources.py path/to/export_rbeq
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from rbeq_reference import SNAPSHOT, THESIS, partial

with tempfile.TemporaryDirectory() as directory:
    for name, parameters in [
        ("iaa_thesis_2023", THESIS),
        ("elmolcs_b8643810", SNAPSHOT),
    ]:
        for target, shells in parameters.items():
            file = Path(directory) / "table.txt"
            subprocess.run([sys.argv[1], target, name, str(file)], check=True)
            table = np.loadtxt(file)
            assert np.all(np.diff(table[:, 0].astype(np.float32)) > 0)
            assert np.all(table[:, 1] >= 0)
            assert table[0, 1] == 0 and table[-1, 0] == 1e9
            energies = np.unique(
                np.r_[
                    np.geomspace(shells[-1][0] * 1.0001, 1e9, 51),
                    [s[0] * f for s in shells for f in [0.999, 1.001, 1.01, 1.1]],
                ]
            )
            reference = np.array(
                [
                    sum(
                        max(partial(e, s, name.startswith("elmolcs")), 0)
                        for s in shells
                    )
                    for e in energies
                ]
            )
            actual = np.interp(energies, table[:, 0], table[:, 1])
            np.testing.assert_allclose(actual, reference, rtol=3e-4, atol=1e-28)
            print(target, name, "independent quadrature and float32 knots passed")
