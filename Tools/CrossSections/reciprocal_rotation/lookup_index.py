# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Add exact search bounds to existing, unchanged rotational distributions."""

import argparse
from pathlib import Path

import numpy as np
from bundle import read_bundle, write_bundle

BINS = 256


def lookup_arrays(data):
    """Uniform quantile bins bound a search; they never approximate a CDF."""
    result = []
    for prefix in ["angular", "conditional"]:
        offsets = np.asarray(data[prefix + "_offsets"], dtype=np.int64)
        coordinates = np.asarray(data[prefix + "_u"])
        lookup = np.zeros((len(offsets) - 1, BINS + 1), dtype=np.uint32)
        for row, (start, end) in enumerate(zip(offsets[:-1], offsets[1:])):
            if end == start:
                continue
            values = coordinates[start:end]
            assert len(values) >= 2 and values[0] == 0 and values[-1] == 1
            indices = (
                np.searchsorted(values, np.arange(BINS + 1) / BINS, side="right") - 1
            )
            lookup[row] = start + np.clip(indices, 0, len(values) - 2)
        result.append((prefix + "_lookup", "u32", lookup.ravel()))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    metadata, data = read_bundle(args.bundle)
    arrays = []
    for name, value in data.items():
        if name in ["angular_lookup", "conditional_lookup"]:
            continue
        kind = (
            "alias"
            if name == "aliases"
            else {"f": {4: "f32", 8: "f64"}, "u": {2: "u16", 4: "u32"}}[
                value.dtype.kind
            ][value.dtype.itemsize]
        )
        arrays.append((name, kind, value))
    arrays.extend(lookup_arrays(data))
    print(write_bundle(args.output, metadata, arrays), "bytes", flush=True)


if __name__ == "__main__":
    main()
