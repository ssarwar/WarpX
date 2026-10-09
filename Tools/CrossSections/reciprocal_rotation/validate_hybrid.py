# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent physics and interpolation checks of the complete source model."""

import argparse
import json
import time

import numpy as np
from hybrid_reference import REST, C, Hybrid, cg_array
from reference_paths import output_path
from rotation_reference import KB
from scipy.special import eval_legendre

p = argparse.ArgumentParser()
p.add_argument("--target", default="N2")
p.add_argument("--temperature", type=float, default=300)
p.add_argument("--resolution", type=int, default=1)
p.add_argument("--output", default="validation")
a = p.parse_args()
start = time.perf_counter()
m = Hybrid(a.target, a.temperature, 1000, a.resolution, 48).solve()
# Test factorial CG coefficients against independent triple-Legendre integrals.
rng = np.random.default_rng(1384)
nodes, w = np.polynomial.legendre.leggauss(256)
errs = []
for _ in range(200):
    i = int(rng.integers(0, 65))
    L = int(2 * rng.integers(0, 25))
    f = int(rng.integers(abs(i - L), i + L + 1))
    f += (i + L + f) % 2
    if f > i + L:
        continue
    direct = (
        (2 * f + 1)
        / 2
        * (w * eval_legendre(i, nodes) * eval_legendre(L, nodes))
        @ eval_legendre(f, nodes)
    )
    stored = cg_array([i], [f], [L])[0, 0]
    errs.append(abs(stored - direct))
assert max(errs) < 1e-10
phi = m.X @ m.measure


def interp(E):
    return m.interp(phi, E)


def moments(E):
    E = np.atleast_1d(E)
    v = []
    for e in E:
        p = np.sqrt(e * (e + 2 * REST))
        kin = C / (e + REST)
        X = interp(e)[0]
        phase = (
            np.sqrt(
                np.maximum(e - m.delta, 0) * (np.maximum(e - m.delta, 0) + 2 * REST)
            )
            / p
        )
        powers = []
        for power in [0, 1, 2]:
            u = kin * (np.asarray((phase * m.delta**power) @ m.up) @ X)
            d = kin * np.einsum(
                "db,db->", m.down * m.d_delta[:, None] ** power, interp(e + m.d_delta)
            )
            powers.append((u, d))
        v.append([kin * m.diag @ X, *powers[0], *powers[1], sum(powers[2])])
    return np.array(v)


# Quadratures are performed separately for excitation and de-excitation;
# forward flux is not algebraically reused as reverse flux.
T = a.temperature
power = []
for order in [8, 16]:
    end = 40 * KB * T
    edges = np.unique(np.r_[0, end, np.geomspace(1e-11, end, 120), m.delta, m.energy])
    edges = edges[(edges >= 0) & (edges <= end)]
    x, w = np.polynomial.legendre.leggauss(order)
    E = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
    dE = (np.diff(edges)[:, None] * w / 2).ravel()
    values = moments(E)[:, 3:5]
    rel = []
    for Te in [0.8 * T, T, 1.2 * T]:
        density = np.sqrt(E * (E + 2 * REST)) * (E + REST) * np.exp(-E / (KB * Te)) * dE
        P = density @ values
        rel.append(float((P[1] - P[0]) / P.sum()))
    power.append({"order": order, "heating_imbalance_for_Te_over_Trot_08_1_12": rel})
assert abs(power[-1]["heating_imbalance_for_Te_over_Trot_08_1_12"][1]) < 1e-3
assert (
    power[-1]["heating_imbalance_for_Te_over_Trot_08_1_12"][0] > 0
    and power[-1]["heating_imbalance_for_Te_over_Trot_08_1_12"][2] < 0
)

# Shape and moment probes deliberately do not coincide with grid nodes.
E = np.unique(
    np.r_[
        np.geomspace(0.001003, 998.0, 180),
        [
            0.025,
            1 - 1e-6,
            1 + 1e-6,
            1.25 - 1e-6,
            1.25 + 1e-6,
            2.22,
            2.47,
            4 - 1e-6,
            4 + 1e-6,
            10 - 1e-6,
            10 + 1e-6,
            20 - 1e-6,
            20 + 1e-6,
            200 - 1e-5,
            200 + 1e-5,
        ],
    ]
)
vals = moments(E)
target_rate = np.sqrt(E * (E + 2 * REST)) * C / (E + REST) * m.prior.s.elastic(E)
rate_error = abs(vals[:, :3].sum(axis=1) / target_rate - 1)
# Signed labels retain their thresholds; a positive interpolated primitive
# cannot open a subthreshold channel because the outgoing momentum is zero.
thresh_err = 0
for delta in rng.choice(m.delta, 200):
    e = max(float(delta) * (1 - 1e-7), 1e-12)
    phase = np.sqrt(max(e - delta, 0) * (max(e - delta, 0) + 2 * REST))
    thresh_err = max(thresh_err, phase)
assert thresh_err == 0

# Round stored primitives, not kinematic labels, to binary32; retain the
# source's double-precision interpolation and compare physical moments.
old = phi
phi = phi.astype(np.float32).astype(float)
single = moments(E)
phi = old
scale = np.maximum(abs(vals), np.max(abs(vals), axis=0) * 1e-12)
sp_error = float(np.max(abs(single - vals) / scale))
zero_X = interp(m.d_delta)
zero_moments = [
    float(C / REST * np.einsum("db,db->", m.down * m.d_delta[:, None] ** k, zero_X))
    for k in [0, 1, 2]
]
report = {
    "target": a.target,
    "temperature_K": T,
    "resolution": a.resolution,
    "energy_nodes": len(m.energy),
    "angle_nodes": len(m.y),
    "seconds": time.perf_counter() - start,
    "maximum_CG_absolute_error": max(errs),
    "minimum_unchanged_background": m.raw_residual_min,
    "zero_energy_gain_rate_power_second": zero_moments,
    "power_checks": power,
    "maximum_offgrid_integral_rate_error": float(rate_error.max()),
    "worst_rate_energy_eV": float(E[rate_error.argmax()]),
    "maximum_float32_primitive_moment_error": sp_error,
    "energy_probes": E.tolist(),
    "moments": vals.tolist(),
}
output_path(a.output + "-" + a.target + f"-{T:g}K.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(
    json.dumps(
        {k: v for k, v in report.items() if k not in ["energy_probes", "moments"]},
        indent=2,
    ),
    flush=True,
)
