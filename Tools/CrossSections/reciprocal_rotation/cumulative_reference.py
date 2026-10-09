# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Prepare an offline cumulative reference from independently decoded aliases."""

import argparse
import math
from pathlib import Path

import numpy as np
from bundle import read_bundle, write_bundle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    metadata, data = read_bundle(args.bundle)
    aliases = data.pop("aliases")
    cdf = np.empty(len(aliases))
    for start, end in zip(data["cell_offsets"][:-1], data["cell_offsets"][1:]):
        values = aliases[start:end]
        count = int(end - start)
        cut = values["cut"].astype(float)
        probability = (
            cut + np.bincount(values["alias"], weights=1 - cut, minlength=count)
        ) / count
        # Normalize the accumulated roundoff before conversion to binary64.
        # Forcing only the last element to one can create a descending final
        # interval when the preceding partial sum rounds just above one.
        partial = np.cumsum(probability, dtype=np.longdouble)
        assert abs(math.fsum(probability) - 1) < 16 * np.finfo(float).eps
        cdf[start:end] = partial / partial[-1]
    arrays = []
    for name, value in data.items():
        kind = {"f": {4: "f32", 8: "f64"}, "u": {2: "u16", 4: "u32"}}[value.dtype.kind][
            value.dtype.itemsize
        ]
        arrays.append((name, kind, value))
    arrays.extend([("cdf", "f64", cdf), ("outcome_ids", "u16", aliases["outcome"])])
    metadata[0] = metadata[0].replace(" alias", " cumulative")
    print(write_bundle(args.output, metadata, arrays), "bytes", flush=True)


if __name__ == "__main__":
    main()
