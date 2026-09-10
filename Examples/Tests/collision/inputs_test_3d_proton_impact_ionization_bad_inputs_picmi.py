#!/usr/bin/env python3
"""Reject source/product aliasing before product-tile growth can invalidate pointers."""

import os
import subprocess
import sys
from pathlib import Path

source = (
    Path(__file__)
    .resolve()
    .with_name("inputs_test_3d_proton_impact_ionization_picmi.py")
)
environment = {
    key: value
    for key, value in os.environ.items()
    if not key.startswith(("PMI_", "PMIX_", "OMPI_"))
}
result = subprocess.run(
    [sys.executable, str(source), "--alias-product"],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    text=True,
    check=False,
    timeout=60,
    env=environment,
)
assert result.returncode != 0, "Projectile/product aliasing was accepted"
assert (
    "Proton-impact product species must differ from the projectile" in result.stdout
), result.stdout
print("PASS: aliased projectile/product rejected")
