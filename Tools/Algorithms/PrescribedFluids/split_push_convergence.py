#!/usr/bin/env python3
"""Resolve a marginal Coulomb split-push energy failure using seeds and timestep.

The original analysis and its tolerances remain unchanged. Refinement halves
the timestep and doubles the step count, preserving the physical duration.
The Slurm return code describes execution only; every analysis return code and
energy error is retained, including failures.
"""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from discharge_convergence import run_simulation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--build-name", default="build_pm_gpu_sync")
    parser.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    parser.add_argument(
        "--versions",
        nargs="+",
        choices=["branch", "stock"],
        default=["branch", "stock"],
    )
    args = parser.parse_args()
    args.root, args.output = args.root.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(
        interpretation=__doc__,
        seeds=args.seeds,
        original_input_sha256={},
        original_analysis_sha256={},
        cases=[],
    )
    base = "test_2d_collisions_split_momentum_push_electromagnetic"
    for version in args.versions:
        root = args.root if version == "branch" else args.root / "build/stock-warpx"
        build = root / args.build_name
        example = root / "Examples/Tests/collision"
        report["original_input_sha256"][version] = {
            name: hashlib.sha256((example / name).read_bytes()).hexdigest()
            for name in [
                "inputs_" + base,
                "inputs_base_2d_collisions_split_momentum_push",
            ]
        }
        report["original_analysis_sha256"][version] = hashlib.sha256(
            (
                example / "analysis_test_2d_collisions_split_momentum_push.py"
            ).read_bytes()
        ).hexdigest()
        inventory = json.loads(
            subprocess.check_output(
                ["ctest", "--test-dir", str(build), "--show-only=json-v1"], text=True
            )
        )
        tests = {test["name"]: test for test in inventory["tests"]}
        configurations = [("original-repeat", None, 1)] + [
            (f"seed{seed}-dt{divisor}", seed, divisor)
            for divisor in [1, 2]
            for seed in args.seeds
        ]
        for label, seed, divisor in configurations:
            directory = args.output / version / label
            directory.mkdir(parents=True)
            record = dict(
                version=version, case=label, seed=seed, timestep_divisor=divisor
            )
            record["stages"] = []
            for stage in ["run", "analysis"]:
                test = tests[base + "." + stage]
                properties = {
                    item["name"]: item["value"] for item in test["properties"]
                }
                overrides = dict(
                    item.split("=", 1) for item in properties.get("ENVIRONMENT", [])
                )
                environment = os.environ.copy()
                environment.update(overrides)
                environment["MPLBACKEND"] = "Agg"
                command = list(test["command"])
                if stage == "run" and seed is not None:
                    command.extend(
                        [
                            f"warpx.random_seed={seed}",
                            f"my_constants.dt={0.1 / divisor}/wpe",
                            f"max_step={2000 * divisor}",
                            f"diag1.intervals={2000 * divisor}",
                            "warpx.verbose=0",
                        ]
                    )
                returncode = run_simulation(
                    command, directory, environment, 600, log_name=stage + ".log"
                )
                record["stages"].append(
                    dict(
                        stage=stage,
                        command=command,
                        environment=overrides,
                        returncode=returncode,
                    )
                )
                if stage == "run" and returncode:
                    break
            if record["stages"][0]["returncode"] == 0:
                field = np.loadtxt(directory / "diags/reducedfiles/field_energy.txt")
                kinetic = np.loadtxt(
                    directory / "diags/reducedfiles/particle_energy.txt"
                )
                assert np.array_equal(field[:, :2], kinetic[:, :2])
                total = field[:, 2] + kinetic[:, 2]
                error = (total - total[0]) / total[0]
                record.update(
                    time_start=float(field[0, 1]),
                    time_end=float(field[-1, 1]),
                    initial_field_energy=float(field[0, 2]),
                    initial_kinetic_energy=float(kinetic[0, 2]),
                    final_field_energy=float(field[-1, 2]),
                    relative_energy_max=float(np.max(np.abs(error))),
                    relative_energy_final=float(error[-1]),
                )
            report["cases"].append(record)
            (args.output / "summary.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
            print(version, label, record.get("relative_energy_max"), flush=True)
    raise SystemExit(
        int(any(row["stages"][0]["returncode"] for row in report["cases"]))
    )


if __name__ == "__main__":
    main()
