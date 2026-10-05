# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Offline momentum-transfer spectra, preserving discrete rotational diffusion."""

import numpy as np
from rotation_reference import REST
from rotor_basis import make_basis
from scipy.special import spherical_jn, spherical_yn


def high_tables(target, temperature=300.0, ratio=1.01):
    zdata = make_basis(target, temperature)
    coeff = zdata["coefficients"]
    initial = zdata["initial"]
    final = zdata["final"]
    loss = zdata["loss"]
    L = zdata["ranks"]
    R = {"N2": 2.068, "O2": 2.281}[target]
    radius = {"N2": 0.6052, "O2": 0.5677}[target]
    b = (R / (2 * radius)) ** 2
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
    changing = loss != 0
    palette_loss = np.r_[0.0, loss[changing]]
    palette_initial = (
        np.r_[0.0, 0.0002477204284695341 * initial[changing] * (initial[changing] + 1)]
        if target == "N2"
        else np.r_[
            0.0,
            0.00017828927734694198 * (initial[changing] * (initial[changing] + 1) - 2),
        ]
    )
    count = int(np.ceil(np.log(10000 / 1e-5) / np.log(ratio)))
    edges = np.r_[
        0, np.geomspace(1e-5, 10000, count + 1), np.geomspace(10000, 1e6, 101)[1:]
    ]

    def phase_average(z):
        V = np.column_stack(
            (
                np.sqrt(2 * L + 1) * spherical_jn(L, z),
                -np.sqrt(2 * L + 1) * spherical_yn(L, z),
            )
        )
        M = V.T @ (cm[0, :, None] * V)
        ev, U = np.linalg.eigh(M)
        assert np.all(ev > 0)
        G = (U / np.sqrt(ev)) @ U.T / np.sqrt(ev).sum()
        return np.square(V @ np.linalg.cholesky(G)).sum(axis=1)

    def integrate(lo, hi):
        if lo >= 10000:
            x, w = np.polynomial.legendre.leggauss(8)
            z = lo + (hi - lo) * (x + 1) / 2
            measure = w * (hi - lo) / 2 * 2 * b * z / (z * z + b) ** 2
            shape = np.array([phase_average(v) for v in z]).T
            A = shape @ measure
            return A, cm[1:] @ A, cm[1:] @ (shape @ (measure * z * z))
        order = max(8, int(np.ceil((hi - lo) * 2)) + 4)
        x, w = np.polynomial.legendre.leggauss(order)
        z = lo + (hi - lo) * (x + 1) / 2
        measure = w * (hi - lo) / 2 * 2 * b * z / (z * z + b) ** 2
        shape = (2 * L[:, None] + 1) * spherical_jn(L[:, None], z[None, :]) ** 2
        norm = cm[0] @ shape
        A = (shape * (measure / norm)).sum(axis=1)
        return A, cm[1:] @ A, (cm[1:] @ (shape * (measure * z * z / norm))).sum(axis=1)

    def alias_table(prob, ids):
        n = len(prob)
        q = prob * n
        cut = np.ones(n)
        dest = np.arange(n, dtype=np.uint16)
        small = list(np.flatnonzero(q < 1))
        large = list(np.flatnonzero(q >= 1))
        while small and large:
            lo = small.pop()
            hi = large.pop()
            cut[lo] = q[lo]
            dest[lo] = hi
            q[hi] -= 1 - q[lo]
            (small if q[hi] < 1 else large).append(hi)
        table = np.empty(
            n, dtype=[("cut", "<f4"), ("alias", "<u2"), ("outcome", "<u2")]
        )
        table["cut"] = cut
        table["alias"] = dest
        table["outcome"] = ids
        assert np.all((table["cut"] >= 0) & (table["cut"] <= 1))
        reconstructed = np.asarray(table["cut"], float) / n
        np.add.at(reconstructed, dest, (1 - np.asarray(table["cut"], float)) / n)
        return table, reconstructed

    cells = []
    parts = []
    cdfs = []
    rank_rows = []
    means = []
    mom2angular = []
    masses = []
    momerr = []
    max_pruned = 0.0
    offset = 0
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        A, moment, z2 = integrate(lo, hi)
        rank_rows.append(A / (cm[0] @ A))
        mass = b / (lo * lo + b) - b / (hi * hi + b)
        pfull = coeff @ A / mass
        prob = np.r_[pfull[~changing].sum(), pfull[changing]]
        keep = prob > 1e-15
        keep[0] = True
        removed = prob[~keep].sum()
        max_pruned = max(max_pruned, float(removed))
        prob[0] += removed
        ids = np.flatnonzero(keep).astype(np.uint16)
        prob = prob[keep]
        prob /= prob.sum()
        tab, restored = alias_table(prob, ids)
        features = np.vstack(
            [
                np.ones(len(ids)),
                palette_loss[ids] == 0,
                palette_loss[ids] > 0,
                palette_loss[ids] < 0,
                np.maximum(palette_loss[ids], 0),
                np.maximum(-palette_loss[ids], 0),
                palette_loss[ids] ** 2,
            ]
        )
        exact = features @ prob
        actual = features @ restored
        momerr.append(np.max(abs(actual - exact) / np.maximum(exact, 1e-12)))
        cells.append([offset, len(tab)])
        parts.append(tab)
        cdfs.append(np.cumsum(prob))
        offset += len(tab)
        means.append(exact[1:])
        mom2angular.append(z2)
        masses.append(mass)
    # The final tail uses the last positive distribution, and is tested with an
    # independent probability*D^n bound, not assumed exact.
    means = np.array(means)
    masses = np.array(masses)
    moment_integrals = means * masses[:, None]
    prefix = np.vstack([np.zeros(6), np.cumsum(moment_integrals, axis=0)])
    angle_prefix = np.vstack([np.zeros(6), np.cumsum(mom2angular, axis=0)])
    checks = []
    for E in np.geomspace(1e4, 1e9, 220):
        zmax = R * np.sqrt(E * (E + 2 * REST)) / (7.2973525693e-3 * REST)
        cut = min(zmax, edges[-1])
        j = np.clip(np.searchsorted(edges, cut, side="right") - 1, 0, len(means) - 1)
        mass = b / (edges[j] ** 2 + b) - b / (cut * cut + b)
        tabulated = prefix[j] + mass * means[j]
        _, part, _ = (
            integrate(edges[j], cut) if cut > edges[j] else (None, np.zeros(6), None)
        )
        reference = prefix[j] + part
        if zmax > edges[-1]:
            tail = b / (edges[-1] ** 2 + b) - b / (zmax * zmax + b)
            tabulated += tail * means[-1]
            reference += tail * means[-1]
        norm = 1 + b / zmax**2
        tabulated *= norm
        reference *= norm
        checks.append(
            {
                "E_eV": float(E),
                "relative_moment_error": (
                    abs(tabulated - reference) / np.maximum(reference, 1e-15)
                ).tolist(),
            }
        )
    all_alias = np.concatenate(parts)
    all_cdf = np.concatenate(cdfs)
    return dict(
        edges=edges,
        cells=np.asarray(cells, "<u4"),
        alias=all_alias,
        cdf=all_cdf,
        loss=palette_loss,
        initial_energy=palette_initial,
        initial=np.r_[-1, initial[changing]],
        final=np.r_[-1, final[changing]],
        moments=means,
        angular_moments=angle_prefix,
        partial_error=max(max(c["relative_moment_error"]) for c in checks),
        storage_error=max(momerr),
        pruned_probability=max_pruned,
    )
