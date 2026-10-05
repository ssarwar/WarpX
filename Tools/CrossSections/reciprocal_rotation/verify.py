# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent quadrature of decoded V6 distributions, without random sampling."""

import argparse
import json
from pathlib import Path

import numpy as np
from bundle import read_bundle
from rotation_reference import KB, REST


def cell_moments(a):
    palette = a["outcomes"].reshape(-1, 2)
    loss = palette[:, 0]
    features = np.vstack(
        (loss > 0, loss < 0, np.maximum(loss, 0), np.maximum(-loss, 0), loss**2)
    )
    offsets = a["cell_offsets"]
    result = np.zeros((len(offsets) - 1, 5))
    maximum = np.zeros(len(result))
    for i, (start, end) in enumerate(zip(offsets[:-1], offsets[1:])):
        count = int(end - start)
        if "aliases" in a:
            cell = a["aliases"][start:end]
            assert np.all(cell["alias"] < count) and np.all(cell["outcome"] < len(loss))
            cut = cell["cut"].astype(float)
            assert np.all((cut >= 0) & (cut <= 1))
            probability = (
                cut + np.bincount(cell["alias"], weights=1 - cut, minlength=count)
            ) / count
            ids = cell["outcome"]
        else:
            probability = np.diff(np.r_[0, a["cdf"][start:end]])
            ids = a["outcome_ids"][start:end]
        assert np.all(probability >= 0) and abs(probability.sum() - 1) < 1e-12
        result[i] = features[:, ids] @ probability
        maximum[i] = np.max(loss[ids[probability > 0]])
    return result, maximum


def low_row(a, row, moments, maximum):
    lo, hi = a["angular_offsets"][row : row + 2]
    u = a["angular_u"][lo:hi]
    d = a["deflection"][lo:hi]
    rho = a["changing"][lo:hi]
    lo, hi = a["conditional_offsets"][row : row + 2]
    cu = a["conditional_u"][lo:hi]
    cells = a["conditional_cells"][lo:hi]
    assert u[0] == 0 and u[-1] == 1 and np.all(np.diff(u) > 0)
    assert np.all(np.diff(d) >= 0) and np.all((d >= 0) & (d <= 2))
    assert np.all((rho >= 0) & (rho <= 1))
    if len(cu):
        assert np.max(maximum[cells]) <= a["energies"][row]
    edges = np.unique(np.r_[u, cu])
    x, w = np.polynomial.legendre.leggauss(3)
    q = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
    measure = (np.diff(edges)[:, None] * w / 2).ravel()
    angle = np.interp(q, u, d)
    probability = np.interp(q, u, rho)
    values = np.zeros((6, len(q)))
    values[0] = 1 - probability
    if len(cu):
        for j in range(5):
            values[j + 1] = probability * np.interp(q, cu, moments[cells, j])
    return (
        a["rates"][row] * np.array([values @ (measure * angle**j) for j in range(3)]).T
    )


def main():
    from export import load_reference

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    args = p.parse_args()
    meta, a = read_bundle(args.bundle)
    target = meta[0].split()[0]
    temperature = float(meta[1].split()[0])
    print(
        target, "decoding", len(a.get("aliases", a.get("cdf"))), "entries", flush=True
    )
    moments, maximum = cell_moments(a)
    n = np.searchsorted(a["energies"], 1000, side="right")
    rows = np.array([low_row(a, i, moments, maximum) for i in range(n)])
    reference = load_reference(target, temperature, args.reference_dir)
    exact = np.array(
        [
            reference.angular_moments(e)
            @ (
                reference.measure[:, None]
                * np.column_stack(
                    (np.ones(len(reference.y)), 2 * reference.y**2, 4 * reference.y**4)
                )
            )
            for e in a["energies"][:n]
        ]
    )
    geometry = reference.measure[:, None] * np.column_stack(
        (np.ones(len(reference.y)), 2 * reference.y**2, 4 * reference.y**4)
    )
    primitive = reference.X @ geometry

    def independent_moments(e):
        kin = 299792458.0 / (e + REST)
        current = reference.interp(primitive, e)[0]
        phase = np.zeros(len(reference.delta))
        if e > 0:
            outgoing = np.maximum(e - reference.delta, 0)
            phase = np.sqrt(outgoing * (outgoing + 2 * REST) / (e * (e + 2 * REST)))
        future = reference.interp(primitive, e + reference.d_delta)
        up, down = [], []
        for power in range(3):
            up.append(
                kin
                * (
                    np.asarray((phase * reference.delta**power) @ reference.up)
                    @ current
                )
            )
            down.append(
                kin
                * np.einsum(
                    "db,dbw->w",
                    reference.down * reference.d_delta[:, None] ** power,
                    future,
                )
            )
        return np.vstack(
            (
                kin * reference.diag @ current,
                up[0],
                down[0],
                up[1],
                down[1],
                up[2] + down[2],
            )
        )

    grid = a["energies"][:n]
    probe = (
        grid[:-1, None] + np.diff(grid)[:, None] * np.array([0.25, 0.5, 0.75])
    ).ravel()
    probe_rows = np.repeat(np.arange(n - 1), 3)
    fraction = (probe - grid[probe_rows]) / (grid[probe_rows + 1] - grid[probe_rows])
    fraction = np.where(a["coordinates"][probe_rows] != 0, np.sqrt(fraction), fraction)
    current = (1 - fraction[:, None, None]) * rows[probe_rows] + fraction[
        :, None, None
    ] * rows[probe_rows + 1]
    independent = np.array([independent_moments(e) for e in probe])
    interpolation_error = abs(current - independent) / np.maximum(
        np.maximum(abs(independent), np.max(abs(exact), axis=0)[None, :, :] * 1e-12),
        1e-100,
    )
    scale = np.maximum(abs(exact), np.max(abs(exact), axis=0)[None, :, :] * 1e-12)
    error = abs(rows - exact) / np.maximum(scale, 1e-100)
    worst = np.unravel_index(error.argmax(), error.shape)
    print(
        target,
        "maximum packed-row error",
        error[worst],
        "at",
        a["energies"][worst[0]],
        "component",
        worst[1:],
        flush=True,
    )
    # Separate energy quadratures of heating and cooling, using actual packed rates.
    end = 40 * KB * temperature
    grid = a["energies"][:n]
    edges = np.unique(np.r_[0, end, grid[grid < end]])
    x, w = np.polynomial.legendre.leggauss(16)
    energy = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
    measure = (np.diff(edges)[:, None] * w / 2).ravel()
    i = np.clip(np.searchsorted(grid, energy, side="right") - 1, 0, n - 2)
    t = (energy - grid[i]) / (grid[i + 1] - grid[i])
    t = np.where(a["coordinates"][i] != 0, np.sqrt(t), t)
    power = (1 - t[:, None, None]) * rows[i, 3:5] + t[:, None, None] * rows[i + 1, 3:5]
    density = (
        np.sqrt(energy * (energy + 2 * REST))
        * (energy + REST)
        * np.exp(-energy / (KB * temperature))
        * measure
    )
    integral = np.einsum("e,epw->pw", density, power)
    imbalance = (integral[1] - integral[0]) / integral.sum(axis=0)
    report = dict(
        target=target,
        maximum_row_error=float(error.max()),
        maximum_interpolation_error=float(interpolation_error.max()),
        worst_energy_eV=float(a["energies"][worst[0]]),
        angle_weighted_equilibrium_imbalance=imbalance.tolist(),
    )
    (args.reference_dir / f"verification-{target}.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    np.savez_compressed(
        args.reference_dir / f"packed-moments-{target}.npz",
        energy=grid,
        moments=rows,
        cell_moments=moments,
    )
    print(report, flush=True)
    assert error.max() < 0.002
    assert interpolation_error.max() < 0.002
    assert np.max(abs(imbalance)) < 0.001


if __name__ == "__main__":
    main()
