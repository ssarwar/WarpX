#!/usr/bin/env python3
"""Check angular-model compatibility for fusion and DSMC collisions."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

CASES = {
    "forward": None,
    "backward": None,
    "isotropic": None,
    "IAA": "IAA scattering is only supported for electron Background MCC.",
    "legendre": "legendre requires <collision_name>.legendre_angular_distribution_coefficients",
}


def initialize(model, collision):
    from pywarpx import geometry, libwarpx

    geometry.dims = 3
    deck = Path("collision_angles.inputs").resolve()
    deck.write_text(
        """
max_step = 0
amrex.throw_exception = 1
amrex.signal_handling = 0
amr.n_cell = 8 8 8
amr.max_grid_size = 8
amr.blocking_factor = 8
amr.max_level = 0
geometry.dims = 3
geometry.prob_lo = 0 0 0
geometry.prob_hi = 1 1 1
boundary.field_lo = periodic periodic periodic
boundary.field_hi = periodic periodic periodic
algo.maxwell_solver = Yee
algo.particle_shape = 1
particles.species_names = d t neutron alpha
d.species_type = deuterium
t.species_type = tritium
neutron.species_type = neutron
alpha.species_type = helium4
d.injection_style = none
t.injection_style = none
neutron.injection_style = none
alpha.injection_style = none
collisions.collision_names = test
test.species = d t
"""
        + (
            "test.type = nuclearfusion\n"
            "test.product_species = neutron alpha\n"
            f"test.scattering_angle_model = {model}\n"
            if collision == "fusion"
            else "test.type = dsmc\n"
            "test.scattering_processes = elastic\n"
            "test.elastic_cross_section = elastic.txt\n"
            f"test.elastic_scattering_angle_model = {model}\n"
        )
    )
    Path("elastic.txt").write_text("0 1e-20\n1e6 1e-20\n")
    libwarpx.initialize(["collision_angles", str(deck)])
    libwarpx.finalize()


def check(collision):
    for model, expected in CASES.items():
        environment = os.environ.copy()
        # The child initializes its own MPI singleton and reads a generated deck.
        for name in list(environment):
            if (
                name.startswith(("PMI_", "PMIX_", "OMPI_"))
                or name == "AMREX_INPUTS_FILE_PREFIX"
            ):
                environment.pop(name)
        result = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--model",
                model,
                "--collision",
                collision,
            ],
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=60,
        )
        Path(f"{collision}_{model}.log").write_text(result.stdout)
        if collision == "dsmc" and model == "legendre":
            expected = "legendre is not supported for DSMC/MCC collisions."
        if expected is None:
            assert result.returncode == 0, result.stdout
        else:
            message = " ".join(
                " ".join(line.lstrip("# ").split())
                for line in result.stdout.splitlines()
            )
            assert result.returncode != 0 and expected in message, result.stdout
        print(f"PASS: {collision} model {model}")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", choices=CASES)
parser.add_argument("--collision", choices=("fusion", "dsmc"), default="fusion")
args = parser.parse_args()
if args.model:
    initialize(args.model, args.collision)
else:
    for collision in ("fusion", "dsmc"):
        check(collision)
