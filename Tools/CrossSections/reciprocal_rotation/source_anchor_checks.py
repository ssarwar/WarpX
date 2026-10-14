# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Checks of the final model against the independent published constraints."""

import json

import numpy as np
from hybrid_reference import REST, C, Hybrid, cg_array
from reference_paths import OUTPUT
from rotation_reference import A0, ROTATION, weights
from source_constraints import ENERGY, MOMENT, TOTAL

ROOT = OUTPUT
out = {}
m = Hybrid("N2", 500, 1000, 2, 48).solve()
phi = m.X @ m.measure
coef = cg_array(m.initial, m.final, m.ranks)
for E in [2.22, 2.47]:
    p = np.sqrt(E * (E + 2 * REST))
    kin = C / (E + REST)
    phase = (
        np.sqrt(np.maximum(E - m.loss, 0) * (np.maximum(E - m.loss, 0) + 2 * REST)) / p
    )
    up = kin * phase * m.pop[m.initial] * (coef @ m.interp(phi, E)[0])
    pf = np.where(m.final <= m.jmax, m.pop[np.minimum(m.final, m.jmax)], 0)
    down = (
        kin
        * pf
        * weights("N2", m.initial)
        / weights("N2", m.final)
        * (coef * m.interp(phi, E + m.loss)).sum(axis=1)
    )
    total = m.evaluate(E)["total_rate"]
    branches = (
        np.array(
            [down[(m.final - m.initial) == j].sum() for j in [6, 4, 2]]
            + [0]
            + [up[(m.final - m.initial) == j].sum() for j in [2, 4, 6]]
        )
        / total
    )
    branches[3] = 1 - branches.sum()
    out[f"N2_{E}eV_500K"] = branches.tolist()
    if E == 2.47:
        published = np.array([1.65, 3.04, 19.1, 3.97, 2.82]) / 30.58
        # Jung's finite-angle reconstruction constrains the model, with the
        # source uncertainty retained separately from numerical interpolation.
        assert np.max(abs(branches[1:6] / published - 1)) < 0.055
    print("N2", E, branches, flush=True)
del m, phi, coef
m = Hybrid("O2", 300, 1000, 2, 48).solve()
phi = m.X @ m.measure
phim = m.X @ (m.measure * 2 * m.y * m.y)
coef = np.r_[cg_array([1], [3], m.ranks)[0], np.zeros(m.ns)]
values = []
for k, E in enumerate(ENERGY):
    p = np.sqrt(E * (E + 2 * REST))
    delta = 10 * ROTATION["O2"]
    phase = np.sqrt((E - delta) * (E - delta + 2 * REST)) / p
    # Remove only the known finite-threshold factor to compare to the paper's ANR values.
    sigma = coef @ m.interp(phi, E)[0] / p / A0**2
    sigma_m = coef @ m.interp(phim, E)[0] / p / A0**2
    row = {
        "E_eV": E,
        "sigma13_a02_without_phase": float(sigma),
        "paper_sigma13_a02": TOTAL[k, 1],
        "sigmaMT13_a02_without_phase": float(sigma_m),
        "paper_sigmaMT13_a02": MOMENT[k, 1],
        "relative_errors": [
            float(sigma / TOTAL[k, 1] - 1),
            float(sigma_m / MOMENT[k, 1] - 1),
        ],
    }
    values.append(row)
    assert max(abs(v) for v in row["relative_errors"]) < 1e-6
    print("O2", row, flush=True)
out["O2_Bhattacharyya"] = values
# Explicit zero-energy endpoint is finite and has no unchanged or excitation rate.
zero = m.evaluate(0.0)
out["O2_zero_energy"] = zero
assert zero["moments"][0] == 0 and zero["moments"][1] == 0 and zero["total_rate"] > 0
(ROOT / "source-anchor-checks.json").write_text(json.dumps(out, indent=2) + "\n")
