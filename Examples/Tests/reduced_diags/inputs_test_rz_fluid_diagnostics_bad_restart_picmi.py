#!/usr/bin/env python3
"""Reject missing, incompatible and corrupt cumulative probe checkpoints."""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("checkpoint", type=Path)
args = parser.parse_args()
source = Path(__file__).with_name("inputs_test_rz_fluid_diagnostics_picmi.py").resolve()
environment = {
    key: value
    for key, value in os.environ.items()
    if not key.startswith(("PMI_", "PMIX_", "OMPI_"))
}
for case in [
    "missing_header",
    "bad_version",
    "changed_configuration",
    "future_time",
    "future_step",
]:
    with tempfile.TemporaryDirectory(prefix=case + "-", dir=Path.cwd()) as temporary:
        directory = Path(temporary)
        checkpoint = directory / "checkpoint"
        shutil.copytree(args.checkpoint, checkpoint)
        header = checkpoint / "ReducedDiags/integral/FieldProbeHeader"
        lines = header.read_text().splitlines()
        if case == "missing_header":
            header.unlink()
        else:
            if case == "bad_version":
                lines[0] = "FieldProbe_unknown"
            elif case == "changed_configuration":
                lines[1] += " incompatible"
            elif case == "future_time":
                lines[2] = "1 1.0"
            else:
                lines[2] = "100 2e-13"
            header.write_text("\n".join(lines) + "\n")
        result = subprocess.run(
            [
                sys.executable,
                str(source),
                "--restart",
                str(checkpoint),
                "--record-only",
            ],
            cwd=directory,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=60,
            check=False,
        )
        Path(case + ".log").write_text(result.stdout)
        assert result.returncode != 0, f"Invalid checkpoint {case} was accepted"
        assert (
            "Missing, invalid or incompatible FieldProbe checkpoint state"
            in result.stdout
        ), result.stdout
        print(f"PASS: {case} rejected")
