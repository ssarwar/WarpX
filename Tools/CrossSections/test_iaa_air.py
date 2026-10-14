#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent source anchors, interpolation and accounting checks for IAA air."""

import argparse
import copy
import importlib.util
import itertools
import json
import os
import tempfile
from pathlib import Path

import numpy as np
from iaa_air.layout import load_angular_sources, load_manifest, load_sources
from iaa_air.rates import PEAK_FLOOR, SR_JOIN, TOLERANCE, Rate
from scipy.integrate import quad


def check_shyn_sweeney(root, profiles):
    # Independent transcription of Table I(a,b), p. 1007, Phys. Rev. A 47,
    # 1006 (1993), checked against the rendered paper rather than its text layer.
    # The parenthesized 168-degree points are extrapolations, not measurements.
    angles = np.arange(12, 157, 12)
    energies = [5, 7, 10, 15, 20]
    measured = {
        "a": [
            [9.9, 9.5, 7.8, 5.8, 5.3, 5.7, 5.7, 6.2, 6.1, 4.9, 5.3, 6.3, 6.1],
            [12.3, 11, 9.3, 6.9, 5.8, 6.4, 6.5, 7.4, 7.9, 8.6, 9.9, 11, 12],
            [14, 11, 8.2, 7.4, 6.3, 6.2, 6.7, 7.0, 6.3, 6.9, 7.2, 7.6, 7.6],
            [4.4, 2.9, 2.7, 2.4, 2.7, 2.9, 2.8, 3.2, 3.5, 3.7, 4.2, 4.3, 4.8],
            [2.2, 1.7, 1.5, 1.6, 1.7, 1.7, 1.7, 1.7, 1.9, 1.9, 2.1, 2.1, 2.4],
        ],
        "b": [
            [9.4, 15, 18, 17, 14, 16, 14, 15, 17, 19, 16, 19, 16],
            [8.8, 23, 29, 31, 28, 27, 29, 26, 27, 30, 27, 25, 20],
            [9.3, 16, 21, 24, 18, 19, 14, 14, 13, 11, 11, 11, 9.1],
            [3.9, 8.5, 8.9, 10, 7.9, 7.4, 6.2, 5.5, 5.3, 4.5, 4.7, 3.5, 3.1],
            [2.6, 4.7, 6.0, 6.1, 5.8, 5.7, 5.5, 5.0, 4.2, 3.3, 2.5, 1.8, 0.9],
        ],
    }
    for key, expected in measured.items():
        profile = next(p for p in profiles if p["profile"] == key)
        rows = profile["rows"]
        np.testing.assert_array_equal([r["energy_eV"] for r in rows], energies)
        np.testing.assert_array_equal([r["values"] for r in rows], expected)
        table = np.loadtxt(root / profile["profile_file"])
        assert table.shape[1] == 182
        for row, values in zip(rows, expected):
            np.testing.assert_array_equal(row["angles_deg"], angles)
            assert row["accepted"]
            # Normalization must preserve the measured shape at every source
            # angle, including both repaired values, in the exported profile.
            actual = table[table[:, 0] == row["energy_eV"], 1 + angles]
            assert actual.shape == angles.shape
            np.testing.assert_allclose(
                actual / actual[3], np.array(values) / values[3], rtol=2e-12, atol=0
            )
            repairs = row.get("source_repairs", [])
            if key == "b" and row["energy_eV"] in (5, 15):
                assert len(repairs) == 1
                repair = repairs[0]
                expected_repair = (
                    (12, 94, 9.4) if row["energy_eV"] == 5 else (36, 89, 8.9)
                )
                assert (
                    tuple(
                        repair[k]
                        for k in ("angle_deg", "original_value", "corrected_value")
                    )
                    == expected_repair
                )
                assert repair["units"] == "1e-20 cm2/sr"
                assert repair["doi"] == "10.1103/PhysRevA.47.1006"
                assert repair["location"] == "Table I(b), p. 1007"
                assert repair["source_pdf_sha256"] == (
                    "81ae455b82ba89b7a39d34364d10a7feb539f8bb3ac27070d699b226f17fc3c9"
                )
            else:
                assert not repairs


def check(data_dir):
    manifest = load_manifest(data_dir)
    result = {}
    for gas in ["N2", "O2"]:
        root = data_dir / gas / "IAA"
        sources = load_sources(data_dir, gas)
        channels = {r["id"]: r for r in manifest["targets"][gas]["channels"]}
        assert (
            sum(r["family"] == "vibrational" for r in sources)
            == {"N2": 15, "O2": 41}[gas]
        )
        assert len({r["id"] for r in sources}) == len(sources)
        maximum_error = 0.0
        point_count = 0
        for record in sources:
            rate = Rate(record)
            table = np.loadtxt(data_dir / channels[record["id"]]["cross_section"])
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
        profiles = load_angular_sources(data_dir, gas)
        if gas == "O2":
            check_shyn_sweeney(root, profiles)
        for profile in profiles:
            table = np.loadtxt(root / profile["profile_file"])
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
    n2 = np.loadtxt(data_dir / "N2/IAA/vibrational/v0_to_v01.txt")
    assert np.interp(20, n2[:, 0], n2[:, 1]) == 1.225e-21
    assert np.interp(30, n2[:, 0], n2[:, 1]) == 2.66e-22
    # Verify the archived m5 values without renormalizing the physical source.
    # Its peak differs by 0.0177724% from the nominal paper normalization
    # (1e-20 cm2 at 3.2188e16 cm^-3); report that source discrepancy separately.
    o2 = np.loadtxt(data_dir / "O2/IAA/attachment/three_body_muroran.txt")
    assert o2[:, 1].max() == 3.1073e-47
    result["O2"]["muroran_nominal_peak_relative_offset"] = float(
        o2[:, 1].max() * 3.2188e22 / 1e-24 - 1
    )
    sr0 = np.loadtxt(data_dir / "O2/IAA/electronic/schumann_runge_6.12_eV.txt")
    sr1 = np.loadtxt(data_dir / "O2/IAA/electronic/schumann_runge_8.60_eV.txt")
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
    # Exercise both public manifest formats with the same real physical files.
    # Temporary metadata stays outside warpx-data; no duplicate tables are needed.
    with tempfile.TemporaryDirectory(prefix="iaa-manifest-") as temporary:
        legacy_root = Path(temporary).resolve()
        targets = copy.deepcopy(manifest["targets"])
        for description in targets.values():
            for key in ("elastic", "elastic_dcs", "rotation"):
                description[key] = os.path.relpath(
                    data_dir / description[key], legacy_root
                )
            for channel in description["channels"]:
                for key in (
                    "cross_section",
                    "source_cross_section",
                    "differential_cross_section",
                ):
                    if key in channel:
                        channel[key] = os.path.relpath(
                            data_dir / channel[key], legacy_root
                        )
        legacy = legacy_root / "IAA/air/manifest.json"
        legacy.parent.mkdir(parents=True)
        legacy.write_text(
            json.dumps({"format": "warpx_iaa_air_v1", "targets": targets})
        )
        for gas in ("N2", "O2"):
            choices = itertools.product(
                ("iaa_thesis_2023", "elmolcs_b8643810"),
                ("iaa", "muroran") if gas == "O2" else ("iaa",),
                ("muroran", "phelps", "biagi") if gas == "O2" else ("muroran",),
            )
            for ionization, dissociative, three_body in choices:
                configurations = []
                for root in (data_dir, legacy_root):
                    processes = module.load_iaa_air_processes(
                        root,
                        gas,
                        ion_species="ion",
                        attachment_species={
                            "dissociative": "Ominus",
                            "three_body": "O2minus",
                        },
                        third_body_density=1e24,
                        ionization_model=ionization,
                        dissociative_attachment=dissociative,
                        three_body_attachment=three_body,
                    )
                    for process in processes.values():
                        for key in (
                            "cross_section",
                            "differential_cross_section",
                            "rotation_file",
                        ):
                            if key in process:
                                process[key] = str(Path(process[key]).resolve())
                    configurations.append(processes)
                assert configurations[0] == configurations[1], (
                    gas,
                    ionization,
                    dissociative,
                    three_body,
                )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(check(args.data_dir), indent=2))
