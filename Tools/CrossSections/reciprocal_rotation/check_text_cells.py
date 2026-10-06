# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare every native decoded cell with the original production probabilities.

Use test_reciprocal_rotation cell_output=<path> with a readable V7 input.
The reference is its validated V6 input retained outside warpx-data during
conversion. This deterministic check includes float32 alias packing.
"""

import argparse
from pathlib import Path

import numpy as np
from bundle import read_bundle


def check(reference, decoded):
    _, arrays = read_bundle(reference)
    loss = arrays["outcomes"].reshape(-1, 2)[:, 0]
    features = np.array(
        [
            np.ones(len(loss)),
            loss == 0,
            loss > 0,
            loss < 0,
            np.maximum(loss, 0),
            np.maximum(-loss, 0),
            loss**2,
            loss**4,
        ]
    )
    observed = np.loadtxt(decoded)
    offsets = arrays["cell_offsets"]
    assert observed.shape == (len(offsets) - 1, len(features))
    worst = 0.0
    for cell, (lo, hi) in enumerate(zip(offsets[:-1], offsets[1:])):
        values = arrays["aliases"][lo:hi]
        probability = values["cut"].astype(float) / len(values)
        np.add.at(
            probability,
            values["alias"],
            (1 - values["cut"].astype(float)) / len(values),
        )
        expected = features[:, values["outcome"]] @ probability
        # Keep positive/negative energy powers separate. The mean signed loss
        # can vanish at equilibrium and is not a useful relative-error scale.
        error = np.max(abs(observed[cell] - expected) / np.maximum(expected, 1e-100))
        worst = max(worst, float(error))
    # Encoding permits 2e-6; allow 1e-7 additionally for float32 alias cutoffs.
    assert worst < 2.1e-6, worst
    print("PASS: every cell, maximum relative moment error", worst)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--decoded", type=Path, required=True)
    args = parser.parse_args()
    check(args.reference, args.decoded)
