#!/usr/bin/env python3
"""Report native-state differences without changing any acceptance criterion."""

import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("actual", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = {}
    with np.load(args.reference) as reference, np.load(args.actual) as actual:
        assert set(reference) == set(actual)
        for name in reference:
            expected, value = reference[name], actual[name]
            assert expected.shape == value.shape, name
            difference = np.abs(value - expected)
            nonzero = expected != 0
            peak = float(np.max(np.abs(expected), initial=0))
            error = float(np.max(difference, initial=0))
            rows[name] = dict(
                shape=list(expected.shape),
                changed_elements=int(np.count_nonzero(difference)),
                maximum_absolute_difference=error,
                reference_peak=peak,
                difference_over_peak=error / peak if peak else None,
                maximum_pointwise_relative_difference=float(
                    np.max(difference[nonzero] / np.abs(expected[nonzero]), initial=0)
                ),
            )
    args.output.write_text(json.dumps(rows, indent=2) + "\n")
    print(
        json.dumps(
            {key: row for key, row in rows.items() if row["changed_elements"]}, indent=2
        )
    )


if __name__ == "__main__":
    main()
