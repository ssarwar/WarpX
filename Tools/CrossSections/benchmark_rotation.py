# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare complete MCC timestep cost and estimator noise with identical physics.

Run in a WarpX Python environment. Results include all setup and timestep
measurements and variance times cost; no hardware-dependent speed assertion is
used in CI. The scalar particle reduction in the input fences GPU execution.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--data-dir",
    type=Path,
    default=Path(__file__).resolve().parents[2]
    / "build/Tools/CrossSections/rotation_reference",
)
parser.add_argument("--particles", type=int, default=65536)
parser.add_argument("--steps", type=int, default=8)
parser.add_argument("--repeats", type=int, default=5)
parser.add_argument("--warmup-runs", type=int, default=1)
args = parser.parse_args()
if args.repeats < 2:
    raise ValueError("At least two independent seeds are required to measure variance")
if args.warmup_runs < 0:
    raise ValueError("The number of warmup runs must be nonnegative")
script = (
    Path(__file__).resolve().parents[2]
    / "Examples/Tests/collision/inputs_test_1d_background_mcc_rotation_picmi.py"
)
report = {
    "particles_per_case": args.particles,
    "steps": args.steps,
    "repeats": args.repeats,
    "warmup_runs_per_mode": args.warmup_runs,
}
records = {"alias": [], "cumulative": []}
for repeat in range(-args.warmup_runs, args.repeats):
    # Alternate order to reduce thermal/throttling bias.
    for mode in ["alias", "cumulative"] if repeat % 2 == 0 else ["cumulative", "alias"]:
        name = f"{mode}_{repeat}" if repeat >= 0 else f"warmup_{mode}_{-repeat}"
        directory = (args.output / name).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable,
            str(script),
            "--data-dir",
            str(args.data_dir.resolve()),
            "--particles",
            str(args.particles),
            "--steps",
            str(args.steps),
            "--seed",
            str(2026 + repeat),
        ]
        if mode == "cumulative":
            command.append("--cumulative")
        with (directory / "run.log").open("w") as log:
            subprocess.run(
                command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True
            )
        # Warm both paths in separate processes, then measure fresh ensembles.
        # This avoids mixing one-time driver/cache costs into one sampler's mean.
        if repeat < 0:
            continue
        with np.load(directory / "background_mcc_rotation_results.npz") as data:
            records[mode].append(
                {
                    "startup_seconds": float(data["startup_seconds"]),
                    "step_seconds": float(data["step_seconds"]),
                    "means": {
                        key: float(data[key][0]) for key in data if key.startswith("e_")
                    },
                    "cases": int(data["cases"]),
                }
            )
for mode, runs in records.items():
    seconds = np.array([run["step_seconds"] for run in runs])
    average = float(seconds.mean())
    report[mode] = {
        "runs": runs,
        "mean_step_seconds": average,
        "median_step_seconds": float(np.median(seconds)),
        "standard_deviation_step_seconds": float(seconds.std(ddof=1)),
        "particle_steps_per_second": args.particles
        * args.steps
        * runs[0]["cases"]
        / average,
        "variance_times_seconds": {
            key: float(np.var([r["means"][key] for r in runs], ddof=1) * average)
            for key in runs[0]["means"]
        },
    }
report["cumulative_over_alias_time"] = (
    report["cumulative"]["mean_step_seconds"] / report["alias"]["mean_step_seconds"]
)
(args.output / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
print(
    json.dumps(
        {key: value for key, value in report.items() if key not in records}, indent=2
    )
)
