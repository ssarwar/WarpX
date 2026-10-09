# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Prepare real, bounded low-energy temperature fixtures outside warpx-data."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", choices=["N2", "O2"], required=True)
    parser.add_argument(
        "--temperatures", type=float, nargs="+", default=[0, 100, 250, 350, 1000]
    )
    args = parser.parse_args()
    tools = Path(__file__).resolve().parent
    results = []
    for temperature in args.temperatures:
        root = args.output.resolve() / f"{temperature:g}K"
        reference = root / "reference"
        reference.mkdir(parents=True, exist_ok=True)
        env = os.environ | {"WARPX_ROTATION_OUTPUT": str(reference)}
        shared = [
            "--target",
            args.target,
            "--temperature",
            str(temperature),
            "--maximum-energy",
            "20",
        ]
        commands = [
            [sys.executable, str(tools / "build_reference.py"), *shared],
            [sys.executable, str(tools / "adaptive_low_grid.py"), *shared],
            [
                sys.executable,
                str(tools / "export.py"),
                *shared,
                "--reference-dir",
                str(reference),
                "--output",
                str(
                    root
                    / "data"
                    / args.target
                    / "IAA"
                    / f"reciprocal_hybrid_{temperature:g}K"
                ),
            ],
        ]
        passed = True
        with (reference / f"prepare-{args.target}.log").open("w") as log:
            for command in commands:
                result = subprocess.run(
                    command, env=env, stdout=log, stderr=subprocess.STDOUT
                )
                if result.returncode:
                    passed = False
                    break
        results.append(dict(temperature=temperature, passed=passed))
        print(args.target, results[-1], flush=True)
    (args.output / f"preparation-{args.target}.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    if not all(r["passed"] for r in results):
        raise SystemExit(
            "Some temperature fixtures failed; inspect their source checks."
        )


if __name__ == "__main__":
    main()
