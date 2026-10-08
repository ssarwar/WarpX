# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent resolved/phase-averaged checks for the oscillatory high-q tail.

For frozen slowly varying Hankel envelopes the phase integral is analytic:
<vv^T/(v^T M v)> = M^(-1/2)/tr(M^(1/2)). This averages probabilities,
not losses, and retains the discrete spectrum and its second moment.
"""

import argparse
import json

import numpy as np
from reference_paths import OUTPUT
from scipy.special import spherical_jn, spherical_yn

p = argparse.ArgumentParser()
p.add_argument("--target", default="N2")
a = p.parse_args()
ROOT = OUTPUT
from rotor_basis import make_basis

r = make_basis(a.target)
L = r["ranks"]
loss = r["loss"]
coeff = r["coefficients"]
cm = (
    np.array(
        [
            np.ones(len(loss)),
            loss == 0,
            loss > 0,
            loss < 0,
            np.maximum(loss, 0),
            np.maximum(-loss, 0),
            loss**2,
        ]
    )
    @ coeff
)
radius = {"N2": 0.6052, "O2": 0.5677}[a.target]
R = {"N2": 2.068, "O2": 2.281}[a.target]
b = (R / (2 * radius)) ** 2


def averaged(z):
    j = np.sqrt(2 * L + 1) * spherical_jn(L, z)
    y = -np.sqrt(2 * L + 1) * spherical_yn(L, z)
    V = np.column_stack((j, y))
    M = V.T @ (cm[0, :, None] * V)
    ev, U = np.linalg.eigh(M)
    assert np.all(ev > 0)
    root = (U / np.sqrt(ev)) @ U.T / np.sqrt(ev).sum()
    factor = np.linalg.cholesky(root)
    val = np.square(V @ factor).sum(axis=1)
    return val


def integrate(lo, hi, phase):
    if phase:
        x, w = np.polynomial.legendre.leggauss(64)
        z = lo + (hi - lo) * (x + 1) / 2
        W = (hi - lo) * w / 2
        A = np.array([averaged(v) for v in z]).T
    else:
        edge = np.linspace(lo, hi, int(np.ceil((hi - lo) / 2)) + 1)
        x, w = np.polynomial.legendre.leggauss(4)
        z = (edge[:-1, None] + np.diff(edge)[:, None] * (x + 1) / 2).ravel()
        W = (np.diff(edge)[:, None] * w / 2).ravel()
        V = (2 * L[:, None] + 1) * spherical_jn(L[:, None], z[None, :]) ** 2
        A = V / (cm[0] @ V)
    M = cm[1:] @ A
    density = 2 * b * z / (z * z + b) ** 2
    return np.array([M @ (W * density), M @ (W * density * z * z)])


out = []
for lo, hi in [(10000, 12000), (12000, 16000), (16000, 20000), (20000, 30000)]:
    exact = integrate(lo, hi, False)
    avg = integrate(lo, hi, True)
    error = abs(avg / exact - 1)
    out.append(
        {
            "range": [lo, hi],
            "relative_errors_rate_moments_and_z2_weighted": error.tolist(),
        }
    )
    print(lo, hi, error.max(), flush=True)
(ROOT / f"phase-average-{a.target}.json").write_text(json.dumps(out, indent=2) + "\n")
