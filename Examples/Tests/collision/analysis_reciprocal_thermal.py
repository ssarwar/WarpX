# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Thermal, timestep and restart checks with offline molecular distributions."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--data-dir", type=Path, required=True)
p.add_argument("--source-dir", type=Path, required=True)
p.add_argument("--temperature-data", type=Path)
p.add_argument("--temperatures", type=float, nargs="+", default=[300])
p.add_argument("--targets", nargs="+", default=["N2", "O2", "air"])
p.add_argument("--output", type=Path, required=True)
p.add_argument("--particles", type=int, default=32768)
p.add_argument("--steps", type=int, default=512)
p.add_argument("--mpi-ranks", type=int, default=0)
p.add_argument("--integration", action="store_true")
args = p.parse_args()
args.output = args.output.resolve()
args.output.mkdir(parents=True, exist_ok=True)
script = Path(__file__).with_name("inputs_test_1d_reciprocal_rotation_picmi.py")
launcher = (
    ["srun", "-n", str(args.mpi_ranks), "--gpus-per-task=1"] if args.mpi_ranks else []
)
records = []


def run(name, gas="air", temperature=300, electron_temperature=300, **options):
    directory = args.output / name
    directory.mkdir(exist_ok=True)
    if temperature == 300:
        data = args.data_dir
    else:
        if args.temperature_data is None:
            raise ValueError("Other temperatures require prepared offline fixtures")
        data = args.temperature_data / f"{temperature:g}K" / "data"
    settings = dict(
        data_dir=data.resolve(),
        source_dir=args.source_dir.resolve(),
        target=gas,
        temperature=temperature,
        electron_temperature=electron_temperature,
        particles=args.particles,
        steps=args.steps,
        density=1e28,
        dt=1e-16,
        check_thermal=True,
        output=directory / "results.json",
    )
    if gas == "O2":
        settings.update(density=1e30, dt=1e-17)
    settings.update(options)
    command = [*launcher, sys.executable, str(script.resolve())]
    for key, value in settings.items():
        if value is False:
            continue
        command.append("--" + key.replace("_", "-"))
        if value is not True:
            command.append(str(value))
    with (directory / "run.log").open("w") as log:
        subprocess.run(
            command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True
        )
    result = json.loads((directory / "results.json").read_text())
    records.append(dict(name=name, **result))
    print(name, result["mean_change_eV"], result["standard_error_bound_eV"], flush=True)
    return result


for temperature in args.temperatures:
    for gas in args.targets:
        electrons = (
            [0, 300]
            if temperature == 0
            else [temperature / 4, temperature, 4 * temperature]
        )
        for te in electrons:
            run(f"{gas}-{temperature:g}K-{te:g}K", gas, temperature, te)

if args.integration:
    # Very cold electrons heat even when the neutral translational bath is cold.
    run("unequal-cold-translation", electron_temperature=0, translational_temperature=0)
    # Hot electrons cool with both rotational and translational baths colder.
    run(
        "unequal-warm-translation",
        electron_temperature=1200,
        translational_temperature=100,
    )
    steps = args.steps
    convergence = [
        run(
            f"timestep-{factor}",
            electron_temperature=75,
            dt=1e-16 / factor,
            steps=steps * factor,
            seed=2026 + factor,
        )
        for factor in (1, 2, 4)
    ]
    coarse, fine = convergence[1:]
    sigma = np.hypot(coarse["standard_error_bound_eV"], fine["standard_error_bound_eV"])
    assert abs(coarse["mean_change_eV"] - fine["mean_change_eV"]) < 6 * sigma
    subcycled = run("subcycles", electron_temperature=75, subcycles=2, seed=27182)
    sigma = np.hypot(
        subcycled["standard_error_bound_eV"], coarse["standard_error_bound_eV"]
    )
    assert abs(subcycled["mean_change_eV"] - coarse["mean_change_eV"]) < 6 * sigma
    checkpoint_dir = args.output / "checkpoint" / "diags"
    first = run("checkpoint", steps=8, checkpoint=checkpoint_dir, check_thermal=False)
    restarted = run(
        "restart", steps=8, restart=checkpoint_dir / "chk000008", check_thermal=False
    )
    assert np.allclose(first["final"], restarted["initial"], rtol=2e-12, atol=0)
    assert restarted["mean_initial_eV"] != restarted["mean_final_eV"]

(args.output / "summary.json").write_text(json.dumps(records, indent=2) + "\n")
