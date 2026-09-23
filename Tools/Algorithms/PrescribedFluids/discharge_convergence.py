#!/usr/bin/env python3
"""Repeat the unchanged Turner analysis across seeds and independent refinements.

The stock input is copied into each run with explicit, recorded parameter
changes. Physical duration and diagnostic averaging duration stay fixed.
No collision kernel, reference profile, or assertion tolerance is changed.
Use separate invocations for stock and branch libraries and retain both.
"""

import argparse
import ast
import concurrent.futures
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Expected one occurrence of {old!r}")
    return source.replace(old, new)


def run_case(args, model, seed, refinement, source, reference):
    # (mesh, particles per cell, time steps): vary errors independently before
    # the joint refinement. More cells also increase total particle count.
    cells, particles, timesteps = {
        "base": (1, 1, 1),
        "timestep": (1, 1, 2),
        "particles": (1, 4, 1),
        "joint": (2, 2, 2),
    }[refinement]
    directory = args.output / f"{model}-{refinement}-seed{seed}"
    directory.mkdir()
    adapted = replace_once(
        source,
        "        self.setup_run()",
        f"""        self.nz *= {cells}
        self.seed_nppc *= {particles}
        self.dt /= {timesteps}
        self.max_steps *= {timesteps}
        self.diag_steps *= {timesteps}
        self.setup_run()
        self.sim.random_seed = {seed}
        self.sim.verbose = 0""",
    )
    adapted = replace_once(
        adapted,
        '"../../../../warpx-data/MCC_cross_sections/He/"',
        repr(str(args.data / "MCC_cross_sections/He") + "/"),
    )
    script = directory / "inputs.py"
    script.write_text(adapted)
    command = [*shlex.split(args.launcher), sys.executable, str(script), "--test"]
    if model == "dsmc":
        command.append("--dsmc")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(args.build / "lib/site-packages")
    env["MPLBACKEND"] = "Agg"
    env["OMP_NUM_THREADS"] = "1"
    record = dict(
        model=model,
        seed=seed,
        refinement=refinement,
        cells=32 * cells,
        particles_per_cell=256 * particles,
        timestep_divisor=timesteps,
        command=command,
        library_path=env["PYTHONPATH"],
        adapted_input_sha256=hashlib.sha256(adapted.encode()).hexdigest(),
    )
    started = time.monotonic()
    with (directory / "run.log").open("w") as log:
        result = subprocess.run(
            command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT
        )
    record["run_returncode"] = result.returncode
    record["seconds"] = time.monotonic() - started
    if result.returncode == 0:
        with (directory / "analysis.log").open("w") as log:
            analysis = subprocess.run(
                [sys.executable, str(args.example / "analysis_1d.py")],
                cwd=directory,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        record["analysis_returncode"] = analysis.returncode
        density = np.load(directory / "ion_density_case_1.npy")
        interpolated = np.interp(
            np.linspace(0.0, 1.0, len(density)),
            np.linspace(0.0, 1.0, len(reference)),
            reference,
        )
        error = (density[1:-1] - interpolated[1:-1]) / interpolated[1:-1]
        record["rms_relative_error"] = float(np.sqrt(np.mean(error**2)))
        record["mean_signed_relative_error"] = float(np.mean(error))
        record["density"] = density.tolist()
    (directory / "summary.json").write_text(json.dumps(record, indent=2) + "\n")
    print(
        directory.name,
        record.get("rms_relative_error", "execution failed"),
        flush=True,
    )
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--example", type=Path)
    parser.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4])
    parser.add_argument(
        "--refinements",
        nargs="+",
        choices=["base", "timestep", "particles", "joint"],
        default=["base"],
    )
    parser.add_argument(
        "--models", nargs="+", choices=["mcc", "dsmc"], default=["mcc", "dsmc"]
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--launcher", default="mpiexec -n 1")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    args.example = (
        args.example or root / "Examples/Physics_applications/capacitive_discharge"
    ).resolve()
    args.build, args.output, args.data = (
        path.resolve() for path in (args.build, args.output, args.data)
    )
    args.output.mkdir(parents=True, exist_ok=False)
    source = (args.example / "inputs_base_1d_picmi.py").read_text()
    analysis = (args.example / "analysis_1d.py").read_text()
    # Read the literal published profile, without executing the assertion or
    # needing a simulation output before this independent comparison exists.
    tree = ast.parse(analysis)
    reference = next(
        np.array(ast.literal_eval(node.value.args[0]))
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "ref_density"
            for target in node.targets
        )
    )
    tolerance = float(re.search(r"^tolerance = ([0-9.]+)$", analysis, re.MULTILINE)[1])
    report = dict(
        driver_revision=subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip(),
        input_sha256=hashlib.sha256(source.encode()).hexdigest(),
        analysis_sha256=hashlib.sha256(analysis.encode()).hexdigest(),
        library_sha256={
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((args.build / "lib/site-packages/pywarpx").glob("*1d*.so"))
        },
        original_tolerance=tolerance,
        interpretation="Individual original assertions retained; inspect convergence and seed distributions separately.",
        cases=[],
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [
            pool.submit(run_case, args, model, seed, refinement, source, reference)
            for refinement in args.refinements
            for seed in args.seeds
            for model in args.models
        ]
        for future in concurrent.futures.as_completed(futures):
            report["cases"].append(future.result())
            report["cases"].sort(
                key=lambda row: (row["refinement"], row["model"], row["seed"])
            )
            (args.output / "summary.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
    raise SystemExit(int(any(case["run_returncode"] for case in report["cases"])))


if __name__ == "__main__":
    main()
