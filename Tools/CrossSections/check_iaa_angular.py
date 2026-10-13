#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Prepare independent angular quadrature references and run the portable sampler."""

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.special import eval_legendre


def reference(path, energies):
    table = np.loadtxt(path)
    y = np.sin(np.linspace(0, np.pi, table.shape[1] - 1) / 2)
    # Gauss integration is exact for the linear DCS times each Legendre
    # polynomial in mu=1-2*y*y, including the spherical angular measure.
    nodes, weights = np.polynomial.legendre.leggauss(8)
    cdf_knots = np.sin(np.linspace(0, np.pi, 181) / 2)
    edges = np.unique(np.r_[y, 1 / np.sqrt(2), cdf_knots])
    points = (edges[:-1, None] + edges[1:, None]) / 2 + np.diff(edges)[
        :, None
    ] * nodes / 2
    weights = (np.diff(edges)[:, None] * weights / 2).ravel()
    points = points.ravel()
    mu = 1 - 2 * points**2
    rows = []
    for energy in energies:
        i = np.clip(np.searchsorted(table[:, 0], energy) - 1, 0, len(table) - 2)
        f = np.clip(
            np.log(energy / table[i, 0]) / np.log(table[i + 1, 0] / table[i, 0]), 0, 1
        )
        dcs = (1 - f) * table[i, 1:] + f * table[i + 1, 1:]
        probability = weights * points * np.interp(points, y, dcs)
        probability /= probability.sum()
        rows.append(
            [
                energy,
                *(
                    float(probability @ eval_legendre(order, mu))
                    for order in range(1, 5)
                ),
                float(probability[mu < 0].sum()),
                *(float(probability[points <= upper].sum()) for upper in cdf_knots),
            ]
        )
    return rows


def run(program, output, data_dir=None):
    output.mkdir(parents=True, exist_ok=True)
    theta = np.linspace(0, np.pi, 181)
    mu = np.cos(theta)
    dpi = 1 - 3 * mu**2 + 14 * mu**4 / 3
    isotropic_cdf = (1 - mu) / 2
    dpi_cdf = 0.5 - 15 / 28 * (mu - mu**3) - mu**5 / 2
    np.savetxt(
        output / "elastic.txt",
        [[0.01, *np.ones(181)], [1e4, *np.ones(181)]],
        header="SPECIES: e / N2",
    )
    np.savetxt(output / "dpi.txt", [[1, *dpi], [20, *dpi]])
    (output / "dpi.angular").write_text(
        "WARPX_EXCITATION_ANGULAR_V1\nN2\ndpi.txt\nelastic.txt\n1 20\n"
    )
    np.savetxt(
        output / "dpi.reference",
        [
            [0.5, 0, 0, 0, 0, 0.5, *isotropic_cdf],
            [2, 0, 1 / 7, 0, 8 / 63, 0.5, *dpi_cdf],
            [20, 0, 1 / 7, 0, 8 / 63, 0.5, *dpi_cdf],
            [30, 0, 0, 0, 0, 0.5, *isotropic_cdf],
        ],
    )
    commands = [(output / "dpi.angular", output / "dpi.reference", True)]
    np.savetxt(output / "mixture.txt", [[1, *np.ones(181)], [16, *(dpi * 15 / 14)]])
    np.savetxt(
        output / "mixture.reference",
        [[4, 0, 1 / 14, 0, 4 / 63, 0.5, *((isotropic_cdf + dpi_cdf) / 2)]],
    )
    commands.append((output / "mixture.txt", output / "mixture.reference", False))
    np.savetxt(
        output / "native_m5.txt",
        [[9.9999999, 1e-50], [10.0000001, 3e-50], [11, 4e-50]],
        header="native_energy_grid = resample",
    )
    np.savetxt(output / "large.txt", [[0, 1e-20], [20, 1e-20]])
    np.savetxt(output / "small.txt", [[0, 1e-30], [20, 1e-30]])
    np.savetxt(
        output / "unresolved.txt",
        [[9.99999999, 0], [10, 2e-50], [10.00000001, 0]],
        header="native_energy_grid = resample",
    )
    np.savetxt(
        output / "zeros.txt",
        [[9.9999999, 0], [10.0000001, 2e-20], [11.0000001, 0], [1e9, 0]],
        header="native_energy_grid = resample\nenergy_min_eV = 0\n"
        "energy_max_eV = 1000000000\noutside_energy_range = error",
    )
    if data_dir:
        for gas in ["N2", "O2"]:
            root = data_dir / gas / "IAA/air/angular"
            description = json.loads((root / "angular_sources.json").read_text())
            for profile in description["profiles"]:
                path = root / (profile["file_stem"] + ".txt")
                table = np.loadtxt(path)
                indices = np.unique(np.linspace(0, len(table) - 1, 7).astype(int))
                energies = np.unique(
                    np.r_[
                        table[indices, 0],
                        np.sqrt(table[indices[:-1], 0] * table[indices[1:], 0]),
                    ]
                )
                rows = reference(path, energies)
                lo, hi = profile["active_range_eV"]
                rows += reference(
                    data_dir / gas / "IAA/elastic_dcs.txt", [lo / 2, hi * 2]
                )
                ref = output / f"{gas}_{profile['file_stem']}.reference"
                np.savetxt(ref, rows, fmt="%.17e")
                commands.append((path.with_suffix(".angular"), ref, False))
    for file, ref, relaxation in commands:
        result = subprocess.run(
            [
                str(program),
                f"file={file.resolve()}",
                f"reference={ref.resolve()}",
                f"relaxation={int(relaxation)}",
                f"native_grid={(output / 'native_m5.txt').resolve()}",
                f"large_rate={(output / 'large.txt').resolve()}",
                f"small_rate={(output / 'small.txt').resolve()}",
                f"zero_grid={(output / 'zeros.txt').resolve()}",
            ],
            text=True,
            capture_output=True,
        )
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        print(
            file.name,
            next(x for x in result.stdout.splitlines() if x.startswith("IAA_ANGULAR")),
            flush=True,
        )
    for name, contents in {
        "bad_species": "O2\ndpi.txt\nelastic.txt\n1 20",
        "bad_bounds": "N2\ndpi.txt\nelastic.txt\n1 21",
        "recursive": "N2\ndpi.txt\ndpi.angular\n1 20",
    }.items():
        path = output / (name + ".angular")
        path.write_text("WARPX_EXCITATION_ANGULAR_V1\n" + contents + "\n")
        result = subprocess.run(
            [
                str(program),
                f"file={path.resolve()}",
                f"reference={(output / 'dpi.reference').resolve()}",
                "amrex.throw_exception=1",
                "amrex.signal_handling=0",
            ],
            capture_output=True,
        )
        assert result.returncode != 0, name
    command = [
        str(program),
        f"file={(output / 'dpi.angular').resolve()}",
        f"reference={(output / 'dpi.reference').resolve()}",
        f"native_grid={(output / 'unresolved.txt').resolve()}",
        f"large_rate={(output / 'large.txt').resolve()}",
        f"small_rate={(output / 'small.txt').resolve()}",
        f"zero_grid={(output / 'zeros.txt').resolve()}",
        "amrex.throw_exception=1",
        "amrex.signal_handling=0",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        assert "unresolved in particle precision" in result.stdout + result.stderr
    else:
        assert "particle_bits=64" in result.stdout


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path)
    args = parser.parse_args()
    run(args.program, args.output, args.data_dir)
