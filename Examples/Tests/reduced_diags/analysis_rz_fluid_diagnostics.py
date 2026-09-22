#!/usr/bin/env python3
"""Compare diagnostics with a stationary Maxwell solution and a restarted run."""

import argparse
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("full", type=Path)
parser.add_argument("restart", type=Path)
parser.add_argument("--nonzero-flux", action="store_true")
args = parser.parse_args()
reference = np.load(args.full / "diagnostics.npz")
resumed = np.load(args.restart / "diagnostics.npz")
errors = {}
prefix_file = args.restart / "continuation_prefix.json"
prefix_rows = json.loads(prefix_file.read_text()) if prefix_file.exists() else {}
for label, data in [("full", reference), ("restart", resumed)]:
    for name in data.files:
        assert np.isfinite(data[name]).all(), (label, name)
        np.testing.assert_allclose(
            data[name][:, 1], data[name][:, 0] * 1e-13, rtol=2e-13
        )
    if args.nonzero_flux:
        continue
    # Independent volume integral B^2*pi*R^2*L/(2*mu0).
    energy = 0.01**2 * np.pi * 0.016**2 * 0.064 / (2 * 1.2566370612685e-6)
    errors[label + "_energy_relative"] = float(
        np.max(np.abs(data["FieldEnergy"][:, 4] / energy - 1))
    )
    for name in ["instant", "integral"]:
        values = data[name]
        expected = 0.01 * (values[:, 1] if name == "integral" else 1)
        scale = 0.01 * (4e-13 if name == "integral" else 1)
        errors[label + "_" + name + "_scaled"] = float(
            np.max(np.abs(values[:, 10] - expected)) / scale
        )
    errors[label + "_poynting_absolute"] = float(
        np.max(np.abs(data["FieldPoyntingFlux"][:, 2:]))
    )
if args.nonzero_flux:
    values = reference["FieldPoyntingFlux"]
    expected = 2 * np.pi * 0.016 * 0.064 * 1000 * 0.01 / 1.2566370612685e-6
    errors["initial_power_relative"] = float(abs(values[0, 4] / expected - 1))
    errors["initial_integral_scaled"] = float(
        np.max(np.abs(values[0, 6:])) / (expected * 1e-13)
    )
for name in reference.files:
    full, restart = reference[name], resumed[name]
    assert np.isfinite(full).all() and np.isfinite(restart).all(), name
    if name in prefix_rows:
        count = prefix_rows[name]
        np.testing.assert_array_equal(restart[:count], full[:count], err_msg=name)
        restart = restart[count:]
    full = full[full[:, 0] >= restart[0, 0]]
    assert full.shape == restart.shape, (name, full.shape, restart.shape)
    scale = np.maximum(np.max(np.abs(full), axis=0), 1e-30)
    errors["continuation_" + name] = float(np.max(np.abs(full - restart) / scale))
Path("comparison.json").write_text(json.dumps(errors, indent=2) + "\n")
print(json.dumps(errors, indent=2))
assert max(errors.values()) < 2e-13, errors
