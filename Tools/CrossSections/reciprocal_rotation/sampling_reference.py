# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Deterministic packed-distribution references for CPU/GPU sampler tests."""

import argparse
from pathlib import Path

import numpy as np
from bundle import read_bundle
from verify import cell_moments, low_row

REST = 510998.95069
ALPHA = 7.2973525693e-3


def reference(meta, a, moments, maximum, energy, cache, powers=(0, 1, 2)):
    i = np.clip(
        np.searchsorted(a["energies"], energy, side="right") - 1,
        0,
        len(a["energies"]) - 2,
    )
    t = (energy - a["energies"][i]) / (a["energies"][i + 1] - a["energies"][i])
    if a["coordinates"][i]:
        t = np.sqrt(t)
    rate = (1 - t) * a["rates"][i] + t * a["rates"][i + 1]
    if energy < 1000:
        for row in [i, i + 1]:
            if row not in cache:
                cache[row] = low_row(a, row, moments, maximum, powers)
        result = (1 - t) * cache[i] + t * cache[i + 1]
        return rate, result / rate if rate else result
    separation, radius = map(float, meta[1].split()[-2:])
    zmax = separation * np.sqrt(energy * (energy + 2 * REST)) / (ALPHA * REST)
    qedges = np.r_[a["high_edges"][a["high_edges"] < zmax], zmax]
    high = moments[a["high_cells"]]
    values = np.column_stack((1 - high[:, 0] - high[:, 1], high))
    if energy >= 10000:
        x, w = np.polynomial.legendre.leggauss(8)
        q = qedges[:-1, None] + np.diff(qedges)[:, None] * (x + 1) / 2
        b = (separation / (2 * radius)) ** 2
        density = 2 * b * q / (q * q + b) ** 2 * (1 + b / zmax**2)
        measure = np.diff(qedges)[:, None] * w / 2 * density
        d = 2 * (q / zmax) ** 2
        result = np.column_stack(
            [values[: len(q)].T @ (measure * d**j).sum(axis=1) for j in powers]
        )
        return rate, result
    result = np.zeros((6, 3))
    for row, weight in [(i, (1 - t) * a["rates"][i]), (i + 1, t * a["rates"][i + 1])]:
        if weight == 0:
            continue
        lo, hi = a["angular_offsets"][row : row + 2]
        u, d = a["angular_u"][lo:hi], a["deflection"][lo:hi]
        cuts = 2 * (qedges / zmax) ** 2
        grid = np.unique(np.r_[u, np.interp(cuts, d, u)])
        x, w = np.polynomial.legendre.leggauss(3)
        q = (grid[:-1, None] + np.diff(grid)[:, None] * (x + 1) / 2).ravel()
        measure = (np.diff(grid)[:, None] * w / 2).ravel()
        deflection = np.interp(q, u, d)
        z = zmax * np.sqrt(deflection / 2)
        cell = np.clip(
            np.searchsorted(a["high_edges"], z, side="right") - 1, 0, len(high) - 1
        )
        result += weight * np.column_stack(
            [values[cell].T @ (measure * deflection**j) for j in powers]
        )
    return rate, result / rate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--sampler", type=Path, required=True)
    p.add_argument("--samples", type=int, default=32768)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    meta, arrays = read_bundle(args.bundle)
    moments, maximum = cell_moments(arrays)
    squared_moments, _ = cell_moments(arrays, squared=True)
    sampled = np.loadtxt(args.sampler, ndmin=2)
    output, cache, squared_cache = [], {}, {}
    for row in sampled:
        energy, rate = row[:2]
        expected_rate, expected = reference(
            meta, arrays, moments, maximum, energy, cache
        )
        expected = expected.ravel()
        _, expected_square = reference(
            meta,
            arrays,
            squared_moments,
            maximum,
            energy,
            squared_cache,
            powers=(0, 2, 4),
        )
        variance_of_mean = (
            np.maximum(expected_square.ravel() - expected**2, 0) / args.samples
        )
        observed, error = row[2::2], row[3::2]
        assert abs(rate / expected_rate - 1) < 1e-12 if expected_rate else rate == 0
        # The finite-support term covers rare outcomes with fewer than one
        # expected sample. Deterministic quadrature supplies their accuracy test.
        loss = arrays["outcomes"].reshape(-1, 2)[:, 0]
        support = np.repeat([1, 1, 1, max(loss), max(-loss), max(loss * loss)], 3)
        support *= np.tile([1, 2, 4], 6)
        bound = 6 * error + 8 * support / args.samples
        assert np.all(abs(observed - expected) <= bound), (
            energy,
            observed,
            expected,
            bound,
        )
        # When a moment is well sampled, the global rare-tail support bound
        # must not hide a biased energy transfer. Use an independently decoded
        # variance, including fourth energy and angular moments, in this check.
        well_sampled = expected**2 > 400 * variance_of_mean
        assert np.all(
            abs(observed[well_sampled] - expected[well_sampled])
            <= 6 * np.sqrt(variance_of_mean[well_sampled]) + 1e-13
        ), (energy, observed, expected, variance_of_mean)
        output.append(np.r_[energy, expected_rate, expected])
    if args.output:
        np.savetxt(args.output, output, fmt="%.17e")
    print("PASS: independent decoded moments for", len(sampled), "energies", flush=True)


if __name__ == "__main__":
    main()
