#!/usr/bin/env python3
"""Write and run independent one-parameter convergence studies.

Each directory retains the exact command, JSON measurements and native arrays.
Timings from this study are not performance measurements: use ensemble.py on an
otherwise idle allocation for those. Source cases share a 10 ps physical duration.
"""

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
        "study", choices=["deposition", "source", "joint", "solvers", "continuum"]
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=int, default=6)
    parser.add_argument("--launcher", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    program = Path(__file__).with_name("benchmark.py")
    cases = []

    def add(name, *options):
        command = shlex.split(args.launcher) + [sys.executable, str(program)]
        cases.append(dict(name=name, command=command + list(map(str, options))))

    if args.study == "deposition":
        # Compare direct and charge-conserving particle currents separately.
        # Particle quadrature and the direct-current shape both converge with h.
        for nr, solver, beam, deposition in itertools.product(
            [16, 32, 64],
            ["semi_implicit_em", "semi_implicit_mm"],
            ["fluid", "quiet"],
            ["direct", "villasenor"],
        ):
            if beam == "fluid" and deposition == "direct":
                continue
            add(
                f"{solver}-n{nr}-{beam}-{deposition}",
                "--mode",
                "fields",
                "--solver",
                solver,
                "--beam",
                beam,
                "--implicit-deposition",
                deposition,
                "--ppc",
                256,
                "--cells",
                nr,
                4 * nr,
                "--steps",
                6,
            )
    elif args.study == "continuum":
        # Hold mesh spacing fixed while enlarging the domain, then refine the
        # enlarged domain. This separates finite-wall bias from mesh error.
        for name, nr, nz, radial, longitudinal in [
            ("mesh16", 16, 64, 8, 12),
            ("mesh32", 32, 128, 8, 12),
            ("mesh64", 64, 256, 8, 12),
            ("domain16", 64, 256, 16, 24),
            ("domain32", 128, 512, 32, 48),
            ("expanded_fine", 256, 1024, 32, 48),
        ]:
            add(
                name,
                "--mode",
                "fields",
                "--beam",
                "fluid",
                "--cells",
                nr,
                nz,
                "--radial-sigmas",
                radial,
                "--longitudinal-sigmas",
                longitudinal,
                "--steps",
                4,
                "--dt",
                1e-13,
            )
    elif args.study == "solvers":
        for solver, dt, seed in itertools.product(
            ["Yee", "PSATD", "semi_implicit_em", "semi_implicit_mm"],
            [5e-13, 2.5e-13],
            range(41, 41 + args.seeds),
        ):
            add(
                f"{solver}-dt{dt}-seed{seed}",
                "--mode",
                "coupled",
                "--mcc",
                "--beam",
                "fluid",
                "--solver",
                solver,
                "--cells",
                32,
                128,
                "--dt",
                dt,
                "--steps",
                round(1e-11 / dt),
                "--seed",
                seed,
            )
    else:
        variants = {
            "base": [],
            "mesh16": ["--cells", 16, 64],
            "mesh64": ["--cells", 64, 256],
            "dt_large": ["--dt", 1e-12, "--steps", 10],
            "dt_small": ["--dt", 2.5e-13, "--steps", 40],
            "weight_large": ["--weight", 400],
            "weight_small": ["--weight", 25],
            "cap2": ["--cap", 2],
            "cap32": ["--cap", 32],
            "resolution2": ["--source-resolution", 2],
            "resolution32": ["--source-resolution", 32],
            "subcycles4": ["--subcycles", 4],
            "thermal": ["--ions", "thermal"],
        }
        if args.study == "joint":
            # Weight decreases faster than cell area: the global pending
            # population tends to zero while the spatial mesh is refined.
            variants = {
                "mesh16_weight100": ["--cells", 16, 64, "--weight", 100],
                "mesh32_weight12.5": ["--cells", 32, 128, "--weight", 12.5],
                "mesh64_weight1.5625": ["--cells", 64, 256, "--weight", 1.5625],
            }
        for name, options in variants.items():
            for seed in range(41, 41 + args.seeds):
                add(
                    f"{name}-seed{seed}",
                    "--mode",
                    "coupled",
                    "--mcc",
                    "--beam",
                    "fluid",
                    "--cells",
                    32,
                    128,
                    "--steps",
                    20,
                    "--seed",
                    seed,
                    *options,
                )

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
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
