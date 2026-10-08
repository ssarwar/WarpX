# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Measure the detailed-balance error of the compressed high-energy tables.

Unlike the continuous reference, finite bins are not claimed exactly
reciprocal pointwise. Test flux-weighted defects and flag negligible channels
separately. Bin-boundary crossings are included explicitly in the quadrature.
"""

import json

import numpy as np
from hybrid_reference import REST, Elastic, cg_array
from reference_paths import OUTPUT
from rotation_reference import ROTATION

ROOT = OUTPUT
out = []
for target in ["N2", "O2"]:
    data = np.load(ROOT / f"binned-q-moments-{target}.npz")
    edges = data["edges"]
    A = data["rank_coefficients"]
    L = data["ranks"]
    elastic = Elastic(target)
    R = {"N2": 2.068, "O2": 2.281}[target]
    B = ROTATION[target]
    pairs = (
        [(0, 2), (8, 10), (8, 16), (10, 40), (0, 196)]
        if target == "N2"
        else [(1, 3), (9, 11), (9, 17), (11, 41), (1, 167)]
    )
    for E in [1e4, 1e5, 2.5e6, 1e9]:
        p = np.sqrt(E * (E + 2 * REST))
        zmax = R * p / (7.2973525693e-3 * REST)
        for i, f in pairs:
            delta = B * (f * (f + 1) - i * (i + 1))
            ep = E - delta
            pp = np.sqrt(ep * (ep + 2 * REST))
            ratio = pp / p
            coef = cg_array([i], [f], L)[0]
            values = A @ coef
            cuts = np.unique(
                np.r_[
                    0, zmax, edges[edges < zmax], (edges / ratio)[edges / ratio < zmax]
                ]
            )
            x, w = np.polynomial.legendre.leggauss(4)
            z = (cuts[:-1, None] + np.diff(cuts)[:, None] * (x + 1) / 2).ravel()
            dz = (np.diff(cuts)[:, None] * w / 2).ravel()
            y = z / zmax
            measure = dz * 8 * np.pi * z / zmax**2
            k = np.clip(np.searchsorted(edges, z) - 1, 0, len(values) - 1)
            kr = np.clip(np.searchsorted(edges, z * ratio) - 1, 0, len(values) - 1)
            fwd = p * p * elastic.dcs(E, y) * values[k]
            rev = pp * pp * elastic.dcs(ep, y) * values[kr]
            norm = 0.5 * ((fwd + rev) @ measure)
            defect = abs(fwd - rev) @ measure / norm if norm > 0 else 0.0
            strength = (elastic.dcs(E, y) * values[k]) @ measure / elastic.s.elastic(E)
            out.append(
                {
                    "target": target,
                    "E_eV": E,
                    "initial": i,
                    "final": f,
                    "relative_state_specific_rate": float(strength),
                    "flux_weighted_relative_balance_defect": float(defect),
                }
            )
            print(target, E, i, f, strength, defect, flush=True)
(ROOT / "high-table-balance.json").write_text(json.dumps(out, indent=2) + "\n")
