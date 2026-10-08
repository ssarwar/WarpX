# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check high-energy rate and angular interpolation against the source DCS."""

import argparse
import json
from pathlib import Path

import numpy as np
from bundle import read_bundle
from hybrid_reference import Elastic
from sampling_reference import ALPHA, REST, reference
from verify import cell_moments


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    metadata, a = read_bundle(args.bundle)
    target = metadata[0].split()[0]
    source = Elastic(target)
    moments, maximum = cell_moments(a)
    separation = float(metadata[1].split()[-2])
    high = moments[a["high_cells"]]
    probabilities = np.column_stack((1 - high[:, 0] - high[:, 1], high))
    energies = np.unique(
        np.r_[
            np.geomspace(1000, 10000, 161),
            np.linspace(7900, 10100, 100),
            np.geomspace(1e4, 1e9, 41),
        ]
    )
    x, w = np.polynomial.legendre.leggauss(6)
    out = []
    for energy in energies:
        rate, packed = reference(metadata, a, moments, maximum, energy, {})
        zmax = separation * np.sqrt(energy * (energy + 2 * REST)) / (ALPHA * REST)
        edges = np.unique(
            np.r_[0, 1, source.y, a["high_edges"][a["high_edges"] < zmax] / zmax]
        )
        y = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
        measure = (np.diff(edges)[:, None] * w / 2).ravel() * 8 * np.pi * y
        shape = source.shape(energy, y)
        bins = np.clip(
            np.searchsorted(a["high_edges"], zmax * y, side="right") - 1,
            0,
            len(high) - 1,
        )
        exact = np.column_stack(
            [
                probabilities[bins].T @ (shape * measure * (2 * y * y) ** j)
                for j in range(3)
            ]
        )
        expected_rate = (
            299792458
            * np.sqrt(energy * (energy + 2 * REST))
            / (energy + REST)
            * source.s.elastic(energy)
        )
        error = abs(rate * packed - expected_rate * exact) / np.maximum(
            abs(expected_rate * exact), 1e-100
        )
        out.append(
            dict(
                energy_eV=float(energy),
                maximum_error=float(error.max()),
                rate_error=float(abs(rate / expected_rate - 1)),
            )
        )
    result = dict(
        target=target, maximum_error=max(r["maximum_error"] for r in out), checks=out
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        target,
        result["maximum_error"],
        max(out, key=lambda r: r["maximum_error"]),
        flush=True,
    )
    assert result["maximum_error"] < 0.0015


if __name__ == "__main__":
    main()
