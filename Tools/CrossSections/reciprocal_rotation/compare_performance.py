# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare repeated, same-physics MCC/PIC timings and transfer variance costs."""

import argparse
import json
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--benchmark", type=Path, required=True)
p.add_argument(
    "--reference", choices=["cumulative", "reference", "disabled"], default="cumulative"
)
p.add_argument("--output", type=Path, required=True)
args = p.parse_args()
data = json.loads(args.benchmark.read_text())
runs = data["runs"] if isinstance(data, dict) else data
keys = sorted({(r["particles"], r["target"], r["case"]) for r in runs})
reports = []
rng = np.random.default_rng(314159)
for key in keys:
    paired = []
    rows = [r for r in runs if (r["particles"], r["target"], r["case"]) == key]
    for repeat in sorted({r["repeat"] for r in rows}):
        by_mode = {r["sampling"]: r for r in rows if r["repeat"] == repeat}
        if "alias" in by_mode and args.reference in by_mode:
            paired.append((by_mode["alias"], by_mode[args.reference]))
    if len(paired) < 3:
        continue
    report = dict(
        particles=key[0],
        target=key[1],
        case=key[2],
        reference=args.reference,
        repeats=len(paired),
    )
    indices = rng.integers(len(paired), size=(10000, len(paired)))
    for field in ["mcc_seconds", "step_seconds"]:
        values = np.array([[a[field], b[field]] for a, b in paired])
        ratio = np.median(values[:, 0]) / np.median(values[:, 1])
        boot = np.median(values[indices, 0], axis=1) / np.median(
            values[indices, 1], axis=1
        )
        report[field + "_ratio"] = float(ratio)
        report[field + "_ratio_interval"] = np.quantile(boot, [0.025, 0.975]).tolist()
        report[field + "_above_five_percent"] = bool(ratio > 1.05)
    times = np.array([[a["step_seconds"], b["step_seconds"]] for a, b in paired]).mean(
        axis=0
    )
    for name, measure in [
        ("energy", lambda r: r["mean_change_eV"]),
        ("momentum", lambda r: (r["final"][4] - r["initial"][4]) / r["particles"]),
    ]:
        values = np.array([[measure(a), measure(b)] for a, b in paired])
        cost = values.var(axis=0, ddof=1) * times
        report[name + "_variance_cost"] = cost.tolist()
        report[name + "_variance_cost_ratio"] = (
            float(cost[0] / cost[1]) if cost[1] else None
        )
    if args.reference == "reference":
        report["identical_particle_moments"] = all(
            a["initial"] == b["initial"] and a["final"] == b["final"] for a, b in paired
        )
    reports.append(report)
args.output.write_text(json.dumps(reports, indent=2) + "\n")
for report in reports:
    print(
        report["particles"],
        report["case"],
        "MCC ratio",
        report["mcc_seconds_ratio"],
        "PIC ratio",
        report["step_seconds_ratio"],
        flush=True,
    )
