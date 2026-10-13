#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent source anchors, interpolation and accounting checks for IAA air."""

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from iaa_air.rates import PEAK_FLOOR, SR_JOIN, TOLERANCE, Rate
from iaa_air.sources import file_stem
from scipy.integrate import quad


def check(data_dir):
    manifest = json.loads((data_dir / "IAA/air/manifest.json").read_text())
    result = {}
    for gas in ["N2", "O2"]:
        root = data_dir / gas / "IAA/air"
        sources = json.loads((root / "sources.json").read_text())["records"]
        assert (
            sum(r["family"] == "vibrational" for r in sources)
            == {"N2": 15, "O2": 41}[gas]
        )
        assert len({r["id"] for r in sources}) == len(sources)
        maximum_error = 0.0
        point_count = 0
        for record in sources:
            rate = Rate(record)
            table = np.loadtxt(root / "production" / (file_stem(record["id"]) + ".txt"))
            energy, sigma = table.T
            assert np.all(np.diff(energy) > 0) and energy[-1] == 1e9
            assert np.all(np.isfinite(table)) and np.all(sigma >= 0)
            fractions = np.array([0.113, 0.327, 0.613, 0.887])
            probes = (energy[:-1, None] + np.diff(energy)[:, None] * fractions).ravel()
            expected = rate(probes)
            actual = np.interp(probes, energy, sigma)
            errors = np.abs(actual - expected) / np.maximum(
                expected, rate.peak * PEAK_FLOOR
            )
            maximum_error = max(maximum_error, float(errors.max()))
            assert errors.max() <= TOLERANCE * 1.01, (
                gas,
                record["id"],
                errors.max(),
                probes[errors.argmax()],
            )
            assert np.all(rate(np.array([record["loss_eV"] * 0.99])) == 0)
            if record["model"] == "analytic" and "partition" not in record:
                # Archived numerical evaluations are independent of this exporter.
                x, y = np.asarray(record["table"]).T
                assert np.max(np.abs(rate(x) - y)) < 2e-5 * rate.peak, (
                    gas,
                    record["id"],
                )
            point_count += len(energy)
        for profile in json.loads((root / "angular/angular_sources.json").read_text())[
            "profiles"
        ]:
            table = np.loadtxt(root / "angular" / (profile["file_stem"] + ".txt"))
            y = np.sin(np.linspace(0, np.pi, table.shape[1] - 1) / 2)
            # Adaptive quadrature over each original linear segment, independent
            # of the exporter's closed-form integral and the C++ CDF builder.
            for row in table[:: max(1, len(table) // 5), 1:]:
                integral = sum(
                    quad(lambda z: z * np.interp(z, y, row), lo, hi, epsabs=1e-13)[0]
                    for lo, hi in zip(y[:-1], y[1:])
                )
                assert abs(8 * np.pi * integral - 1) < 2e-10
            for row in profile["rows"]:
                if row["accepted"]:
                    assert row["completion"]["modeled_probability"] < 0.5
        result[gas] = {
            "maximum_relative_interpolation_error": maximum_error,
            "rate_points": point_count,
        }
    n2 = np.loadtxt(data_dir / "N2/IAA/air/production/v01.txt")
    assert np.interp(20, n2[:, 0], n2[:, 1]) == 1.225e-21
    assert np.interp(30, n2[:, 0], n2[:, 1]) == 2.66e-22
    # Verify the archived m5 values without renormalizing the physical source.
    # Its peak differs by 0.0177724% from the nominal paper normalization
    # (1e-20 cm2 at 3.2188e16 cm^-3); report that source discrepancy separately.
    o2 = np.loadtxt(data_dir / "O2/IAA/air/production/three_body_muroran.txt")
    assert o2[:, 1].max() == 3.1073e-47
    result["O2"]["muroran_nominal_peak_relative_offset"] = float(
        o2[:, 1].max() * 3.2188e22 / 1e-24 - 1
    )
    sr0 = np.loadtxt(data_dir / "O2/IAA/air/production/schumann_runge_onset.txt")
    sr1 = np.loadtxt(data_dir / "O2/IAA/air/production/schumann_runge_bulk.txt")
    e = np.linspace(*SR_JOIN, 301)
    total = np.interp(e, sr0[:, 0], sr0[:, 1]) + np.interp(e, sr1[:, 0], sr1[:, 1])
    assert np.all(total > 0) and np.all(np.diff(total) > 0)
    for gas, description in manifest["targets"].items():
        selected = [c for c in description["channels"] if c["default"]]
        assert sum(c["type"] == "ionization" for c in selected) == 1
        assert sum(c["type"] == "attachment" for c in selected) == (
            2 if gas == "O2" else 0
        )
    # Read the public loader without importing pywarpx's compiled extension.
    path = Path(__file__).resolve().parents[2] / "Python/pywarpx/iaa_air.py"
    spec = importlib.util.spec_from_file_location("air_loader", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for gas in ["N2", "O2"]:
        processes = module.load_iaa_air_processes(
            data_dir,
            gas,
            ion_species="ion",
            attachment_species={"dissociative": "Ominus", "three_body": "O2minus"},
            third_body_density=1e24,
        )
        assert len(processes) == {"N2": 37, "O2": 54}[gas]
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(check(args.data_dir), indent=2))
