#!/usr/bin/env python3
"""Reproduce acceleration restart failures with repeat and native-state controls.

Run inside a Perlmutter GPU allocation. Original CTest launch commands and
parameters are retained, while all new output uses an isolated directory.
Reported component errors retain the original 1e-12 acceptance criterion.
"""

import argparse
import json
import os
import subprocess
import sys
import traceback
from pathlib import Path

from report_plotfile_difference import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    inventory = json.loads(
        subprocess.check_output(
            ["ctest", "--test-dir", str(args.build), "--show-only=json-v1"], text=True
        )
    )
    tests = {test["name"]: test for test in inventory["tests"]}
    scripts = Path(__file__).resolve().parent
    records = []
    for suffix in ("", "_psatd", "_psatd_time_avg"):
        name = "test_3d_acceleration" + suffix
        test = tests[name + ".run"]
        properties = {item["name"]: item["value"] for item in test["properties"]}
        env = os.environ.copy()
        env.update(item.split("=", 1) for item in properties.get("ENVIRONMENT", []))
        command = test["command"]
        application = next(
            i
            for i, item in enumerate(command)
            if Path(item).name.startswith("warpx.3d.")
        )
        parameters = command[application + 2 :]
        # AMReX prepends this prefix even to absolute input names. Preserve
        # CTest's relative name and its environment for both executables.
        input_path = command[application + 1]
        record = dict(test=name, runs=[], comparisons={})
        try:
            for mode in ("cpp", "native"):
                variants = (
                    ("full", "restart")
                    if mode == "cpp"
                    else ("full", "repeat", "restart")
                )
                for variant in variants:
                    directory = args.output / name / (mode + "_" + variant)
                    directory.mkdir(parents=True)
                    extra = list(parameters)
                    if variant == "restart":
                        checkpoint = directory.parent / (mode + "_full/diags/chk000005")
                        extra.append("amr.restart=" + str(checkpoint))
                    if mode == "cpp":
                        invocation = [
                            *command[: application + 1],
                            str(input_path),
                            *extra,
                        ]
                    else:
                        invocation = [
                            *command[:application],
                            sys.executable,
                            str(scripts / "probe_restart_fields.py"),
                            str(input_path),
                            *("--parameter=" + value for value in extra),
                        ]
                    with (directory / "run.log").open("w") as stream:
                        result = subprocess.run(
                            invocation,
                            cwd=directory,
                            env=env,
                            stdout=stream,
                            stderr=subprocess.STDOUT,
                            timeout=300,
                        )
                    record["runs"].append(
                        dict(
                            directory=str(directory),
                            command=invocation,
                            returncode=result.returncode,
                            environment={
                                key: env[key]
                                for key in (
                                    "PYTHONPATH",
                                    "OMP_NUM_THREADS",
                                    "AMREX_INPUTS_FILE_PREFIX",
                                )
                                if key in env
                            },
                        )
                    )
                    if result.returncode:
                        raise RuntimeError(
                            f"Simulation failed: {directory / 'run.log'}"
                        )
                base = args.output / name
                reference = base / (mode + "_full/diags/diag1000010")
                for variant in variants[1:]:
                    actual = base / (mode + "_" + variant + "/diags/diag1000010")
                    record["comparisons"][mode + "_" + variant] = compare(
                        reference, actual
                    )
            base = args.output / name
            record["comparisons"]["cpp_vs_native_full"] = compare(
                base / "cpp_full/diags/diag1000010",
                base / "native_full/diags/diag1000010",
            )
            for variant in ("repeat", "restart"):
                for step in (5, 6, 10):
                    output = base / f"{variant}-native-{step}.json"
                    command = [
                        sys.executable,
                        str(scripts / "report_state_difference.py"),
                        str(base / f"native_full/native_{step}.npz"),
                        str(base / f"native_{variant}/native_{step}.npz"),
                        "--output",
                        str(output),
                    ]
                    with output.with_suffix(".log").open("w") as stream:
                        subprocess.run(
                            command, stdout=stream, stderr=subprocess.STDOUT, check=True
                        )
                    record["comparisons"][f"{variant}_native_{step}"] = json.loads(
                        output.read_text()
                    )
        except Exception:
            record["execution_error"] = traceback.format_exc()
        records.append(record)
        (args.output / "comparison.json").write_text(
            json.dumps(records, indent=2) + "\n"
        )
        print(name, record.get("execution_error", "recorded"), flush=True)
    failed = any(
        record.get("execution_error")
        or any(
            comparison.get("failed_fields")
            for comparison in record["comparisons"].values()
        )
        for record in records
    )
    raise SystemExit(int(failed))


if __name__ == "__main__":
    main()
