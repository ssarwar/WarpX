#!/usr/bin/env python3
"""Reject source/product aliasing and invalid parser values, including Release GPUs."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--rz", action="store_true")
args = parser.parse_args()
source = (
    Path(__file__)
    .resolve()
    .with_name(
        "inputs_test_rz_rigid_source_picmi.py"
        if args.rz
        else "inputs_test_3d_proton_impact_ionization_picmi.py"
    )
)
cases = (
    []
    if args.rz
    else [
        (
            "--alias-product",
            "Proton-impact product species must differ from the projectile",
        )
    ]
)
for quantity in ["density", "temperature"]:
    cases.append(
        (
            f"--bad-{quantity}",
            f"Rigid-source neutral {quantity} must be finite and non-negative."
            if args.rz
            else "Proton-impact ionization requires finite, non-negative neutral",
        )
    )
environment = {
    key: value
    for key, value in os.environ.items()
    if not key.startswith(("PMI_", "PMIX_", "OMPI_"))
}
for argument, message in cases:
    result = subprocess.run(
        [sys.executable, str(source), argument],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=60,
        env=environment,
    )
    assert result.returncode != 0, f"Invalid input {argument} was accepted"
    assert message in result.stdout, result.stdout
    print(f"PASS: {argument} rejected")
