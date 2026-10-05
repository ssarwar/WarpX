# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Validate joint-row mixtures, including threshold and angle-weighted balance.

The same selected row must provide the elastic angle and rotational outcome.
Using an unrelated uniform deviate as the global angular CDF coordinate of
an energy-mixture sampler would give the wrong joint distribution.
"""

import argparse
import json
import time

import numpy as np
from hybrid_reference import REST, C, Hybrid
from reference_paths import OUTPUT, output_path
from rotation_reference import KB

p = argparse.ArgumentParser()
p.add_argument("--target", default="N2")
p.add_argument("--tolerance", type=float, default=0.0005)
a = p.parse_args()
ROOT = OUTPUT
start = time.perf_counter()
m = Hybrid(a.target, 300, 1000, 2, 48, angular_resolution=4 if a.target == "N2" else 2)
from build_reference import build_reference

checkpoint = output_path(f"refined-reference-{a.target}.npz")
if checkpoint.exists():
    saved = np.load(checkpoint)
    if "temperature" in saved and float(saved["temperature"]) != 300.0:
        raise ValueError(
            "Reference temperature does not match this sampling-grid check"
        )
    if not np.array_equal(saved["y"], m.y):
        raise ValueError(
            "Reference angular mesh does not match the sampling-grid check"
        )
    m.energy = saved["energy"]
    m.X = saved["X"]
    m.c = saved["c"]
    m.p = np.sqrt(m.energy * (m.energy + 2 * REST))
    m.kin = C / (m.energy + REST)
    m.boundary = np.searchsorted(m.energy, 1001.0)
    anchor = np.flatnonzero(m.energy == 0.001)[0]
    m.cold_anchor = m.X[anchor, 0] / m.p[anchor]
else:
    m, _ = build_reference(a.target)
weight = m.measure[:, None] * np.column_stack(
    (np.ones(len(m.y)), 2 * m.y * m.y, 4 * m.y**4)
)
phi = m.X @ weight
m.X = None
if hasattr(m.prior, "o2_completion"):
    m.prior.o2_completion.cache.clear()
thresholds = np.unique(m.loss[(m.final - m.initial) <= 6])
thresholds = thresholds[thresholds < 1]


def reference(E):
    values = []
    for e in np.atleast_1d(E):
        kin = C / (e + REST)
        if e == 0:
            future = m.interp(phi, m.d_delta)
            down = np.array(
                [
                    kin
                    * np.einsum(
                        "db,dbw->w", m.down * m.d_delta[:, None] ** power, future
                    )
                    for power in [0, 1, 2]
                ]
            )
            row = np.vstack([np.zeros((2, 3)), down[0], np.zeros(3), down[1], down[2]])
        else:
            X = m.interp(phi, e)[0]
            p = np.sqrt(e * (e + 2 * REST))
            phase = (
                np.sqrt(
                    np.maximum(e - m.delta, 0) * (np.maximum(e - m.delta, 0) + 2 * REST)
                )
                / p
            )
            u = []
            d = []
            future = m.interp(phi, e + m.d_delta)
            for power in [0, 1, 2]:
                u.append(kin * (np.asarray((phase * m.delta**power) @ m.up) @ X))
                d.append(
                    kin
                    * np.einsum(
                        "db,dbw->w", m.down * m.d_delta[:, None] ** power, future
                    )
                )
            row = np.vstack([kin * m.diag @ X, u[0], d[0], u[1], d[1], u[2] + d[2]])
        values.append(row.ravel())
    return np.array(values)


def blend(e, lo, hi):
    flags = np.isin(lo, thresholds) | (lo == 0)
    t = (e - lo) / (hi - lo)
    t[flags] = np.sqrt((e[flags] - lo[flags]) / (hi[flags] - lo[flags]))
    return t


anchors = np.array(
    [
        0,
        1e-9,
        0.001,
        1,
        1.25,
        2.22,
        2.47,
        4,
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        100,
        150,
        200,
        220,
        1000,
    ]
)
grid = np.unique(
    np.r_[
        anchors,
        thresholds,
        m.prior.s.knots[m.prior.s.knots <= 1000],
        m.elastic.energy[m.elastic.energy <= 1000],
        np.geomspace(1e-9, 1000, 80),
    ]
)
grid = grid[(grid >= 0) & (grid <= 1000)]
rows = reference(grid)
peak = np.max(abs(rows), axis=0)
pool = m.energy[m.energy <= 1000]
for iteration in range(16):
    f = np.array([0.25, 0.5, 0.75])
    probe = np.unique(
        np.r_[(grid[:-1, None] + np.diff(grid)[:, None] * f).ravel(), pool]
    )
    index = np.clip(np.searchsorted(grid, probe, side="right") - 1, 0, len(grid) - 2)
    exact = reference(probe)
    lo = grid[index]
    hi = grid[index + 1]
    t = blend(probe, lo, hi)
    fit = (1 - t[:, None]) * rows[index] + t[:, None] * rows[index + 1]
    error = abs(exact - fit) / np.maximum(abs(exact), peak * 1e-12)
    per_probe = error.max(axis=1)
    per_interval = np.zeros(len(grid) - 1)
    np.maximum.at(per_interval, index, per_probe)
    bad = per_interval > a.tolerance
    print(
        "refine",
        iteration,
        "rows",
        len(grid),
        "bad intervals",
        int(bad.sum()),
        "max",
        float(error.max()),
        flush=True,
    )
    if not bad.any():
        break
    mid = np.array(
        [
            probe[np.flatnonzero(index == i)[np.argmax(per_probe[index == i])]]
            for i in np.flatnonzero(bad)
        ]
    )
    grid = np.r_[grid, mid]
    rows = np.vstack([rows, reference(mid)])
    order = np.argsort(grid)
    grid = grid[order]
    rows = rows[order]
assert not bad.any()


def tabulated(E):
    index = np.clip(np.searchsorted(grid, E, side="right") - 1, 0, len(grid) - 2)
    t = blend(E, grid[index], grid[index + 1])
    return (1 - t[:, None]) * rows[index] + t[:, None] * rows[index + 1]


# Separate energy integrations of gain and loss for three angular weights.
end = 40 * KB * 300
edges = np.unique(np.r_[0, end, grid[grid <= end]])
x, w = np.polynomial.legendre.leggauss(16)
E = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
dE = (np.diff(edges)[:, None] * w / 2).ravel()
values = tabulated(E).reshape(-1, 6, 3)
density = np.sqrt(E * (E + 2 * REST)) * (E + REST) * np.exp(-E / (KB * 300)) * dE
powers = np.einsum("e,epw->pw", density, values[:, 3:5])
imbalance = (powers[1] - powers[0]) / powers.sum(axis=0)
assert np.max(abs(imbalance)) < 0.001
out = {
    "target": a.target,
    "energy_rows": len(grid),
    "tolerance": a.tolerance,
    "max_probe_relative_error": float(error.max()),
    "angle_weighted_power_imbalance": imbalance.tolist(),
    "seconds": time.perf_counter() - start,
}
(ROOT / f"adaptive-low-table-{a.target}.json").write_text(
    json.dumps(out, indent=2) + "\n"
)
np.savez_compressed(
    ROOT / f"adaptive-low-grid-{a.target}.npz",
    energy=grid,
    moments=rows,
    thresholds=thresholds,
)
print(out, flush=True)
