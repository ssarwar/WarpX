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
    parser.add_argument(
        "--rigid-source-criteria",
        action="store_true",
        help="Also evaluate the unchanged rigid-source fixture assertions.",
    )
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
            if args.rigid_source_criteria:
                if name.startswith(("Efield_fp", "Bfield_fp")):
                    kind = name.rsplit("_", 1)[0]
                    scale = max(
                        np.max(np.abs(reference[kind + "_" + direction]))
                        for direction in ("r", "theta", "z")
                    )
                    atol = float(64 * np.finfo(value.dtype).eps * scale)
                    rtol = 0.0
                else:
                    rtol, atol = 2e-14, 1e-300
                matches = np.isclose(value, expected, rtol=rtol, atol=atol)
                rows[name].update(
                    original_rtol=rtol,
                    original_atol=atol,
                    failed_elements=int(np.count_nonzero(~matches)),
                    passes_original_criterion=bool(np.all(matches)),
                )
    args.output.write_text(json.dumps(rows, indent=2) + "\n")
    print(
        json.dumps(
            {key: row for key, row in rows.items() if row["changed_elements"]}, indent=2
        )
    )
    if args.rigid_source_criteria:
        raise SystemExit(
            int(any(not row["passes_original_criterion"] for row in rows.values()))
        )


if __name__ == "__main__":
    main()
