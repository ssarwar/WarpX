#!/usr/bin/env python3
"""Reject missing or incompatible prescribed state before advancing a restart."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

source = Path(__file__).with_name("inputs_test_rz_rigid_source_picmi.py").resolve()
checkpoint = Path("../test_rz_rigid_source_Yee_picmi/diags/chk000002").resolve()
cases = {
    "density": (
        "Level_0/fluid_density_i_N2_immobile[level=0]_H",
        "Missing required prescribed-fluid checkpoint density",
    ),
    "remainder": (
        "Level_0/N2_immobile_product_weight_remainder[level=0]_H",
        "Missing required prescribed-source checkpoint state",
    ),
    "counter": (
        "Level_0/N2_immobile_sampling_counter[level=0]_H",
        "Missing required prescribed-source checkpoint state",
    ),
    "budget": (
        "Level_0/N2_immobile_source_budget[level=0]_H",
        "Missing required prescribed-source checkpoint state",
    ),
    "beam": (None, "prescribed-fluid species, models"),
    "source": (None, "immutable physics/sampling configuration changed"),
    "particle_diagnostic": (None, "requires kinetic species"),
}
for name, (missing, expected) in cases.items():
    directory = Path(name).resolve()
    directory.mkdir(exist_ok=True)
    restart = checkpoint
    if missing:
        restart = directory / "checkpoint"
        shutil.copytree(checkpoint, restart, dirs_exist_ok=True)
        (restart / missing).unlink()
    environment = os.environ.copy()
    for key in list(environment):
        if key.startswith(("PMI_", "PMIX_", "OMPI_")):
            environment.pop(key)
    environment["WARPX_TEST_RESTART_MUTATION"] = name
    result = subprocess.run(
        [
            sys.executable,
            str(source),
            "--solver",
            "Yee",
            "--subcycles",
            "2",
            "--restart",
            str(restart),
        ],
        cwd=directory,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=30,
        check=False,
    )
    assert result.returncode != 0, f"Invalid restart {name} unexpectedly succeeded"
    assert expected in " ".join(result.stdout.replace("#", " ").split()), result.stdout[
        -6000:
    ]
    print(f"PASS: {name} rejected before advancing")
