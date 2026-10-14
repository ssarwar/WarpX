# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Resolve bounded-rotor impulse moments on a momentum-transfer grid.

Uses a narrow energy window to evaluate shifted-energy reverse transitions.
The omitted z>10000 quadrature tail is retained as coarse angular cells; its
rate and maximum possible transfer moments are bounded independently.
"""

import argparse
import json
import time

import numpy as np
from hybrid_reference import REST, Hybrid
from reference_paths import output_path

p = argparse.ArgumentParser()
p.add_argument("--target", default="N2")
p.add_argument("--order", type=int, default=3)
p.add_argument("--step", type=float, default=2.0)
p.add_argument("--energy", type=float, nargs="+", default=[1e4, 1e5, 2.5e6, 1e9])
p.add_argument("--output", default="high-q")
a = p.parse_args()
out = []
for E in a.energy:
    start = time.perf_counter()
    R = {"N2": 2.068, "O2": 2.281}[a.target]
    radius = {"N2": 0.6052, "O2": 0.5677}[a.target]
    p0 = np.sqrt(E * (E + 2 * REST))
    zmax = R * p0 / (7.2973525693e-3 * REST)
    upper = min(10000.0, zmax)
    edges = np.unique(
        np.r_[0, np.geomspace(1e-7, upper, 80), np.arange(0, upper, a.step), upper]
    )
    if zmax > upper:
        edges = np.r_[edges, np.geomspace(upper, zmax, 30)[1:]]
    x, w = np.polynomial.legendre.leggauss(a.order)
    z = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
    weights = (np.diff(edges)[:, None] * w / 2).ravel() * 8 * np.pi * z / zmax**2
    y = z / zmax
    rank = 264 if a.target == "N2" else 232
    m = Hybrid(
        a.target,
        300,
        E + 2,
        1,
        rank,
        angular_override=(y, weights),
        energy_override=E + np.linspace(0, 4, 17),
    ).solve()
    v = m.evaluate(E)
    v["seconds"] = time.perf_counter() - start
    v["angle_nodes"] = len(y)
    # Rutherford tail P(z>cut), exact for its angular marginal.
    eta = 1 / (4 * radius**2 * (p0 / (7.2973525693e-3 * REST)) ** 2)
    cut = upper / zmax
    tail = eta * (1 - cut * cut) / (eta + cut * cut) if zmax > upper else 0.0
    diss = {"N2": 9.759, "O2": 5.116}[a.target]
    v["tail_event_probability_bound"] = tail
    v["tail_moment1_absolute_bound_eV"] = tail * diss
    v["tail_moment2_absolute_bound_eV2"] = tail * diss**2
    out.append(v)
    print(json.dumps(v), flush=True)
    del m
output_path(a.output + "-" + a.target + ".json").write_text(
    json.dumps(out, indent=2) + "\n"
)
