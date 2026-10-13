#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare compact and redundant dense representations of identical IAA physics."""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


def densify(source, destination, logarithmic=False):
    table = np.loadtxt(source)
    x, values = table[:, 0], table[:, 1:]
    f = np.arange(4) / 4
    if logarithmic:
        energy = np.exp(np.log(x[:-1, None]) + f * np.diff(np.log(x))[:, None])
        fraction = (np.log(energy) - np.log(x[:-1, None])) / np.diff(np.log(x))[:, None]
    else:
        energy = x[:-1, None] + f * np.diff(x)[:, None]
        fraction = (energy - x[:-1, None]) / np.diff(x)[:, None]
    rows = (1 - fraction[:, :, None]) * values[:-1, None, :] + fraction[
        :, :, None
    ] * values[1:, None, :]
    result = np.column_stack(
        (
            np.r_[energy.ravel(), x[-1]],
            np.vstack((rows.reshape(-1, values.shape[1]), values[-1])),
        )
    )
    # Adjacent binary64 threshold knots have no representable midpoint.
    result = result[np.r_[True, np.diff(result[:, 0]) > 0]]
    header = "\n".join(
        line[1:].strip()
        for line in source.read_text().splitlines()
        if line.startswith("#")
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(destination, result, fmt="%.17e", header=header)


def prepare_dense(source, destination):
    if destination == source or source in destination.parents:
        raise ValueError(
            "Dense benchmark data must be outside the production repository"
        )
    manifest_path = source / "IAA/air/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    (destination / "IAA/air").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(manifest_path, destination / "IAA/air/manifest.json")
    angular = set()
    for target in manifest["targets"].values():
        for key in ["elastic", "elastic_dcs", "rotation"]:
            relative = Path(target[key])
            if key == "rotation":
                relative = relative.parent
            path = destination / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            if not path.exists():
                path.symlink_to(
                    source / relative, target_is_directory=key == "rotation"
                )
        for channel in target["channels"]:
            relative = Path(channel["cross_section"])
            densify(source / relative, destination / relative)
            if "differential_cross_section" in channel:
                angular.add(Path(channel["differential_cross_section"]))
    for relative in angular:
        if relative.suffix != ".angular":
            continue
        descriptor = source / relative
        words = descriptor.read_text().split()
        destination_file = destination / relative
        destination_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(descriptor, destination_file)
        densify(
            descriptor.parent / words[2],
            destination_file.parent / words[2],
            logarithmic=True,
        )


def benchmark(args):
    args.output.mkdir(parents=True, exist_ok=True)
    dense = args.output / "dense-data"
    prepare_dense(args.data_dir, dense)
    script = (
        Path(__file__).resolve().parents[2]
        / "Examples/Tests/collision/inputs_test_1d_iaa_air_picmi.py"
    )
    results = {}
    for label, data in [("compact", args.data_dir), ("dense", dense)]:
        trials = []
        for repeat in range(-1, args.repeats):
            directory = args.output / f"{label}-{repeat + 1}"
            directory.mkdir(exist_ok=True)
            command = [
                sys.executable,
                str(script),
                "--data-dir",
                str(data),
                "--target",
                args.target,
                "--mode",
                args.mode,
                "--particles",
                str(args.particles),
                "--steps",
                str(args.steps),
                "--warmup",
                "4",
                "--seed",
                str(1234 + max(0, repeat)),
                "--output",
                str(directory / "results.json"),
            ]
            run = subprocess.run(command, cwd=directory, capture_output=True, text=True)
            (directory / "run.log").write_text(run.stdout + run.stderr)
            if run.returncode:
                raise RuntimeError(f"Benchmark failed; see {directory / 'run.log'}")
            result = json.loads((directory / "results.json").read_text())
            result["selector_bytes_per_copy"] = sum(
                int(x)
                for x in re.findall(
                    r"cumulative prefix table \(\d+ energy points,\s*(\d+) bytes",
                    run.stdout,
                )
            )
            result["angular_bytes_per_copy"] = sum(
                int(x)
                for x in re.findall(
                    r"Background MCC angular arrays: (\d+) bytes", run.stdout
                )
            )
            assert "exact per-process fallback" not in run.stdout
            if repeat >= 0:
                trials.append(result)
            print(
                label,
                "warmup" if repeat < 0 else repeat,
                result["step_seconds"],
                flush=True,
            )
        elapsed = np.array([r["step_seconds"] for r in trials])
        energy = np.array([r["final"]["mean_energy_eV"] for r in trials])
        results[label] = {
            "trials": trials,
            "mean_step_seconds": float(elapsed.mean()),
            "std_step_seconds": float(elapsed.std(ddof=1)),
            "mean_energy_eV": float(energy.mean()),
            "variance_energy_eV2": float(energy.var(ddof=1)),
            "energy_variance_times_seconds": float(energy.var(ddof=1) * elapsed.mean()),
        }
    compact, dense_result = results["compact"], results["dense"]
    error = abs(compact["mean_energy_eV"] - dense_result["mean_energy_eV"])
    standard_error = np.sqrt(
        (compact["variance_energy_eV2"] + dense_result["variance_energy_eV2"])
        / args.repeats
    )
    results["physics_difference_standard_errors"] = float(
        error / max(standard_error, 1e-30)
    )
    timing_error = np.sqrt(
        (compact["std_step_seconds"] ** 2 + dense_result["std_step_seconds"] ** 2)
        / args.repeats
    )
    regression = compact["mean_step_seconds"] - dense_result["mean_step_seconds"]
    results["compact_regression_beyond_uncertainty"] = bool(
        regression > 0.05 * dense_result["mean_step_seconds"]
        and regression > 2 * timing_error
    )
    (args.output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    assert error < 6 * standard_error + 1e-3
    assert not results["compact_regression_beyond_uncertainty"], (
        "Investigate the recorded timing regression"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", choices=["N2", "O2", "air"], default="air")
    parser.add_argument(
        "--mode", choices=["cooling", "swarm", "pulse"], default="cooling"
    )
    parser.add_argument("--particles", type=int, default=8192)
    parser.add_argument("--steps", type=int, default=64)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    assert args.repeats >= 3
    args.data_dir = args.data_dir.resolve()
    args.output = args.output.resolve()
    benchmark(args)
