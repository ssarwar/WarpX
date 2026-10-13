#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Report completed-set totals and model sensitivities outside warpx-data."""

import argparse
import json
from pathlib import Path

import numpy as np
from iaa_air.angular import THETA, Elastic, Profile
from iaa_air.rates import REST, Rate
from iaa_air.sources import Archive, blocks, file_stem
from scipy.integrate import trapezoid


def vector(path):
    text = "\n".join(x for x in path.read_text().splitlines() if not x.startswith("#"))
    words = text.split()
    assert words[0] == "values"
    values = np.asarray(words[2:], dtype=float)
    assert len(values) == int(words[1])
    return values


def audit(data, archive):
    manifest = json.loads((data / "IAA/air/manifest.json").read_text())
    totals = [
        r
        for r in blocks(archive.read("Data/cs/Cross section.txt"))
        if "COMMENT: TOTAL  " in r["header"]
    ]
    report = {
        "interpretation": "Physical model differences are reported, not removed by retuning rates.",
        "third_body_density_m3": 5.25e24,
        "targets": {},
    }
    for gas in ["N2", "O2"]:
        description = manifest["targets"][gas]
        energies = np.array([5, 10, 20, 50, 100, 300, 1000.0])
        speed = (
            299792458 * np.sqrt(energies * (energies + 2 * REST)) / (energies + REST)
        )
        root = data / gas / "IAA/air"
        rotation = (data / description["rotation"]).parent
        knots = vector(rotation / "energies.txt")
        rates = vector(rotation / "rates.txt")
        coordinates = vector(rotation / "coordinates.txt")
        indices = np.clip(np.searchsorted(knots, energies) - 1, 0, len(knots) - 2)
        fraction = (energies - knots[indices]) / (knots[indices + 1] - knots[indices])
        fraction = np.where(coordinates[indices] != 0, np.sqrt(fraction), fraction)
        combined = (
            (1 - fraction) * rates[indices] + fraction * rates[indices + 1]
        ) / speed
        for channel in description["channels"]:
            if not channel["default"]:
                continue
            table = np.loadtxt(data / channel["cross_section"])
            multiplier = (
                report["third_body_density_m3"] if channel["units"] == "m5" else 1
            )
            combined += multiplier * np.interp(energies, table[:, 0], table[:, 1])
        reference = np.asarray(next(r for r in totals if r["target"] == gas)["table"])
        reference_values = np.interp(energies, reference[:, 0], reference[:, 1])
        records = json.loads((root / "sources.json").read_text())["records"]
        tails = []
        for record in records:
            rate = Rate(record)
            if rate.tail is None:
                continue
            energy = min(10 * rate.maximum, 1e9)
            base = float(rate(energy))
            ratio = energy / rate.maximum
            thermal = []
            grid = np.unique(np.r_[rate.knots, np.geomspace(1e-5, 1e6, 5000)])
            speed_grid = 299792458 * np.sqrt(grid * (grid + 2 * REST)) / (grid + REST)
            beyond = np.maximum(grid / rate.maximum, 1)
            reference_rate = rate(grid) * speed_grid
            for temperature in [0.1, 1, 10, 100]:
                distribution = (
                    2
                    / np.sqrt(np.pi)
                    * np.sqrt(grid)
                    * np.exp(-grid / temperature)
                    / temperature**1.5
                )
                baseline = trapezoid(reference_rate * distribution, grid)
                variants = [
                    trapezoid(reference_rate * distribution * factor, grid)
                    for factor in [1 / beyond, np.sqrt(beyond), grid <= rate.maximum]
                ]
                assert variants[2] <= baseline * (1 + 1e-12)
                assert (
                    variants[0] <= baseline * (1 + 1e-12) <= variants[1] * (1 + 2e-12)
                )
                thermal.append(
                    {
                        "electron_temperature_eV": temperature,
                        "baseline_coefficient": float(baseline),
                        "relative_changes_steeper_shallower_zero": [
                            float(v / baseline - 1) for v in variants
                        ],
                    }
                )
            tails.append(
                {
                    "channel": record["id"],
                    "energy_eV": energy,
                    "baseline": base,
                    "steeper_exponent_by_one": base / ratio,
                    "shallower_exponent_by_half": base * np.sqrt(ratio),
                    "zero_tail": 0,
                    "maxwellian_rate_sensitivity": thermal,
                }
            )
        angles = json.loads((root / "angular/angular_sources.json").read_text())
        transition_sensitivity = []
        elastic = Elastic(data / gas / "IAA/elastic_dcs.txt")
        measure = 2 * np.pi * np.sin(THETA)
        cosine = np.cos(THETA)
        for entry in angles["profiles"]:
            key = entry["profile"]
            matching = next(
                r
                for r in records
                if r["id"] == key
                or (key == "vibrational" and r["family"] == "vibrational")
                or (key == "schumann_runge" and r["id"].startswith(key))
            )
            profile = Profile(gas, key, elastic, entry["rows"], matching["loss_eV"])
            if key in ("v01", "vibrational"):
                lo, hi = (5, 10) if gas == "N2" else (20, 40)
                energy = np.sqrt(lo * hi)
                low = profile(lo)
                high = profile.measured(energy) if profile.rows else elastic(energy)
            else:
                lo, hi = profile.energy[-1], profile.hi
                energy = np.sqrt(lo * hi)
                low, high = profile.pdfs[-1], elastic(energy)
            delta = []
            for width_factor in [0.5, 2]:
                coordinate = np.clip(0.5 / width_factor, 0, 1)
                weight = coordinate**2 * (3 - 2 * coordinate)
                difference = (weight - 0.5) * (high - low) * measure
                delta.append(
                    [
                        float(trapezoid(difference * p, THETA))
                        for p in [cosine, (3 * cosine * cosine - 1) / 2]
                    ]
                )
            transition_sensitivity.append(
                {
                    "profile": key,
                    "energy_eV": float(energy),
                    "log_width_factors": [0.5, 2],
                    "P1_P2_changes": delta,
                }
            )
        cap_mass = [
            r["completion"]["modeled_probability"]
            for p in angles["profiles"]
            for r in p["rows"]
            if r["accepted"]
        ]
        sensitivity = [
            r["completion"]["moment_sensitivity_to_constant_caps"]
            for p in angles["profiles"]
            for r in p["rows"]
            if r["accepted"]
        ]
        # Check loss moments at independent energies across the entire production range.
        probes = np.geomspace(0.01, 1e9, 5001)
        exact = np.zeros((3, len(probes)))
        actual = np.zeros_like(exact)
        for record in records:
            if record["type"] != "excitation":
                continue
            table = np.loadtxt(root / "production" / (file_stem(record["id"]) + ".txt"))
            for order in range(3):
                weight = record["loss_eV"] ** order
                exact[order] += weight * Rate(record)(probes)
                actual[order] += weight * np.interp(probes, table[:, 0], table[:, 1])
        errors = np.max(
            np.abs(actual - exact)
            / np.maximum(exact, exact.max(axis=1)[:, None] * 1e-6),
            axis=1,
        )
        assert np.all(errors < 5e-4), (gas, errors)
        report["targets"][gas] = {
            "energy_eV": energies.tolist(),
            "completed_total_m2": combined.tolist(),
            "archive_total_m2": reference_values.tolist(),
            "relative_total_difference": (combined / reference_values - 1).tolist(),
            "excitation_rate_and_loss_moment_interpolation_errors": errors.tolist(),
            "maximum_modeled_angular_probability": max(cap_mass),
            "maximum_cap_completion_moment_sensitivity": np.max(
                sensitivity, axis=0
            ).tolist(),
            "tail_sensitivity": tails,
            "transition_width_sensitivity": transition_sensitivity,
        }
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.data_dir.resolve() in args.output.resolve().parents:
        raise ValueError("Audit output must be outside the production repository")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(audit(args.data_dir, Archive(args.archive)), indent=2) + "\n"
    )
