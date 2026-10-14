#!/usr/bin/env python3
"""Run reproducible particle/fluid comparisons serially on one CPU/GPU allocation."""

import argparse
import itertools
import json
import shlex
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite", choices=["fields", "source", "coupled", "push"], required=True
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells", nargs=2, type=int, default=[32, 128])
    parser.add_argument("--ppc", nargs="+", type=int, default=[4, 16, 64, 256])
    parser.add_argument(
        "--beams",
        nargs="+",
        choices=["fluid", "quiet", "random"],
        default=["fluid", "quiet", "random"],
    )
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--solver", default="Yee")
    parser.add_argument(
        "--launcher", default="", help="For example: srun -n 1 -c 32 --gpus-per-task=1"
    )
    parser.add_argument("--checkpoint", action="store_true")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args, extra = parser.parse_known_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    program = Path(__file__).with_name("benchmark.py")
    cases = []
    common = [
        "--mode",
        args.suite,
        "--cells",
        *map(str, args.cells),
        "--steps",
        str(args.steps),
        "--solver",
        args.solver,
    ]
    if args.suite == "coupled":
        common += ["--mcc"]
    if args.suite == "push":
        common += ["--no-self-fields"]
    if args.checkpoint:
        common += ["--checkpoint"]
    if args.profile:
        common += ["--profile"]
    common += extra
    ions = ["fluid", "frozen"] if args.suite in ["source", "coupled"] else ["fluid"]
    for beam, ion, seed in itertools.product(args.beams, ions, range(args.seeds)):
        for ppc in [1] if beam == "fluid" else args.ppc:
            name = f"{beam}-{ion}-ppc{ppc}-seed{41 + seed}"
            command = (
                shlex.split(args.launcher) + [sys.executable, str(program)] + common
            )
            command += [
                "--beam",
                beam,
                "--ions",
                ion,
                "--ppc",
                str(ppc),
                "--seed",
                str(41 + seed),
            ]
            cases.append(dict(name=name, command=command))
    manifest = output / "manifest.json"
    if manifest.exists() and json.loads(manifest.read_text()) != cases:
        parser.error("The existing manifest differs; use a new output directory")
    manifest.write_text(json.dumps(cases, indent=2) + "\n")
    if args.dry_run:
        print(f"Wrote {len(cases)} cases to {manifest}")
        return
    for index, case in enumerate(cases):
        directory = output / case["name"]
        directory.mkdir(exist_ok=True)
        if args.resume and (directory / "result.json").exists():
            continue
        print(f"[{index + 1}/{len(cases)}] {case['name']}", flush=True)
        with (directory / "run.log").open("w") as log:
            subprocess.run(
                case["command"],
                cwd=directory,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=1800,
            )


if __name__ == "__main__":
    main()
