# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare complete MCC steps and estimator variance for identical physics."""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--data-dir", type=Path, required=True)
p.add_argument("--cumulative-dir", type=Path, required=True)
p.add_argument("--source-dir", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
p.add_argument("--particles", type=int, nargs="+", default=[262144, 1048576])
p.add_argument("--targets", nargs="+", choices=["N2", "O2", "air"], default=["air"])
p.add_argument(
    "--cases",
    nargs="+",
    choices=["thermal", "resonance", "intermediate", "broad"],
    default=["thermal", "resonance", "intermediate", "broad"],
)
p.add_argument("--steps", type=int, default=64)
p.add_argument("--repeats", type=int, default=5)
p.add_argument("--mpi-ranks", type=int, default=0)
p.add_argument(
    "--mcc-program", type=Path, help="Time the complete MCC operator separately"
)
p.add_argument(
    "--disabled", action="store_true", help="Also measure the legacy elastic baseline"
)
args = p.parse_args()
if args.repeats < 2:
    raise ValueError("Variance comparisons require independent repeated seeds")
script = (
    Path(__file__).resolve().parents[2]
    / "Examples/Tests/collision/inputs_test_1d_reciprocal_rotation_picmi.py"
)
args.output.mkdir(parents=True, exist_ok=True)
records = []
for count in args.particles:
    for target in args.targets:
        for case in args.cases:
            modes = ["alias", "cumulative"] + (["disabled"] if args.disabled else [])
            for repeat in range(-1, args.repeats):
                order = modes if repeat % 2 == 0 else modes[::-1]
                for mode in order:
                    directory = args.output / f"{count}-{target}-{case}-{mode}-{repeat}"
                    directory.mkdir(exist_ok=True)
                    destination = directory / "results.json"
                    command = [
                        sys.executable,
                        str(script),
                        "--data-dir",
                        str(
                            (
                                args.cumulative_dir
                                if mode == "cumulative"
                                else args.data_dir
                            ).resolve()
                        ),
                        "--source-dir",
                        str(args.source_dir.resolve()),
                        "--target",
                        target,
                        "--particles",
                        str(count),
                        "--cells",
                        str(128 * max(1, args.mpi_ranks)),
                        "--steps",
                        str(args.steps),
                        "--warmup",
                        "4",
                        "--pic",
                        "--seed",
                        str(2026 + repeat),
                        "--output",
                        str(destination.resolve()),
                    ]
                    if mode == "cumulative":
                        command.append("--cumulative")
                    elif mode == "disabled":
                        command.append("--disabled")
                    if case in ["thermal", "broad"]:
                        command.extend(["--mode", case])
                    else:
                        command.extend(
                            [
                                "--mode",
                                "mono",
                                "--energy",
                                "2.47" if case == "resonance" else "50",
                            ]
                        )
                    launcher = []
                    if args.mpi_ranks:
                        launcher = [
                            "srun",
                            "-n",
                            str(args.mpi_ranks),
                            "--gpus-per-task=1",
                        ]
                    native = directory / "inputs"
                    if args.mcc_program:
                        with (directory / "input.log").open("w") as log:
                            subprocess.run(
                                [*command, "--write-input", str(native.resolve())],
                                stdout=log,
                                stderr=subprocess.STDOUT,
                                check=True,
                            )
                    with (directory / "run.log").open("w") as log:
                        subprocess.run(
                            [*launcher, *command],
                            stdout=log,
                            stderr=subprocess.STDOUT,
                            check=True,
                            cwd=directory,
                        )
                    mcc_seconds = None
                    device_bytes = None
                    if args.mcc_program:
                        mcc = subprocess.run(
                            [
                                *launcher,
                                str(args.mcc_program.resolve()),
                                str(native.resolve()),
                                f"benchmark.steps={args.steps}",
                                "benchmark.warmup=4",
                                f"benchmark.broad={int(case == 'broad')}",
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            text=True,
                            check=True,
                            cwd=directory,
                        )
                        (directory / "mcc.log").write_text(mcc.stdout)
                        mcc_seconds = float(
                            re.search(r"MCC_OPERATOR_SECONDS ([\d.eE+-]+)", mcc.stdout)[
                                1
                            ]
                        )
                        device_bytes = int(
                            re.search(r"DEVICE_USED_BYTES (\d+)", mcc.stdout)[1]
                        )
                    if repeat < 0:
                        continue
                    result = json.loads(destination.read_text())
                    result.update(
                        sampling=mode,
                        case=case,
                        repeat=repeat,
                        mcc_seconds=mcc_seconds,
                        device_used_bytes=device_bytes,
                    )
                    records.append(result)
                    print(
                        count,
                        target,
                        case,
                        mode,
                        repeat,
                        result["step_seconds"],
                        flush=True,
                    )
summary = []
for count in args.particles:
    for target in args.targets:
        for case in args.cases:
            for mode in modes:
                rows = [
                    r
                    for r in records
                    if (r["particles"], r["target"], r["case"], r["sampling"])
                    == (count, target, case, mode)
                ]
                seconds = np.array([r["step_seconds"] for r in rows])
                changes = np.array([r["mean_change_eV"] for r in rows])
                summary.append(
                    dict(
                        particles=count,
                        target=target,
                        case=case,
                        sampling=mode,
                        median_step_seconds=float(np.median(seconds)),
                        mean_step_seconds=float(seconds.mean()),
                        step_standard_deviation=float(seconds.std(ddof=1)),
                        particle_steps_per_second=float(
                            count * args.steps / np.median(seconds)
                        ),
                        mean_change_eV=float(changes.mean()),
                        variance_times_step_seconds=float(
                            changes.var(ddof=1) * seconds.mean()
                        ),
                        median_startup_seconds=float(
                            np.median([r["startup_seconds"] for r in rows])
                        ),
                        median_mcc_seconds=(
                            float(np.median([r["mcc_seconds"] for r in rows]))
                            if args.mcc_program
                            else None
                        ),
                        peak_device_bytes=(
                            max(r["device_used_bytes"] for r in rows)
                            if args.mcc_program
                            else None
                        ),
                        peak_host_bytes=max(r["peak_host_bytes"] for r in rows),
                        transport_variance_times_seconds=float(
                            np.var([r["final"][5] / count for r in rows], ddof=1)
                            * seconds.mean()
                        ),
                    )
                )
(args.output / "benchmark.json").write_text(
    json.dumps(dict(summary=summary, runs=records), indent=2) + "\n"
)
print(json.dumps(summary, indent=2), flush=True)
