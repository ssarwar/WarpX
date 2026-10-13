#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Export only physical production data; diagnostics go to stdout or a build directory."""

import argparse
import json
from pathlib import Path

import numpy as np
from iaa_air.angular import Elastic, select_profiles, write_profile
from iaa_air.rates import MAX_ENERGY, Rate, master_grid, write_rate
from iaa_air.sources import ARCHIVE_SHA256, Archive, digest, file_stem, read_sources


def write_json(path, content):
    path.write_text(
        json.dumps(content, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )


def export(archive_path, attachment, data_dir, integral_only=False):
    archive = Archive(archive_path)
    sources = read_sources(archive, attachment)
    manifest = {
        "format": "warpx_iaa_air_v1",
        "archive_sha256": ARCHIVE_SHA256,
        "attachment_sha256": digest(attachment.read_bytes()),
        "maximum_energy_eV": MAX_ENERGY,
        "defaults": {
            "ionization": "iaa_thesis_2023",
            "dissociative_attachment": "iaa",
            "three_body_attachment": "muroran",
            "third_body": "O2",
        },
        "assumptions": [
            "completed IAA; ground-state targets; neutral states are not evolved",
            "existing elastic/rotation family is preserved",
            "three-body stabilization uses O2 density, not total air density",
        ],
        "targets": {},
    }
    for gas, records in sources.items():
        root = data_dir / gas / "IAA"
        out = root / "air"
        for sub in ["source", "production", "angular"]:
            (out / sub).mkdir(parents=True, exist_ok=True)
        rates = [Rate(r) for r in records]
        elastic_file = root / "elastic.txt"
        grid = master_grid(rates, np.loadtxt(elastic_file)[:, 0])
        profiles = (
            {}
            if integral_only
            else select_profiles(
                archive, gas, Elastic(root / "elastic_dcs.txt"), records
            )
        )
        angular_sources = [write_profile(out / "angular", p) for p in profiles.values()]
        write_json(
            out / "angular" / "angular_sources.json",
            {
                "format": "warpx_iaa_angular_sources_v1",
                "archive_sha256": ARCHIVE_SHA256,
                "thesis": "Schmalzried 2023, section 11.3.2, equation 11.40",
                "vibrational_transition_eV": [5, 10] if gas == "N2" else [20, 40],
                "minimum_angles": 8,
                "minimum_solid_angle_coverage": 0.75,
                "maximum_angular_gap_degrees": 20,
                "maximum_completion_probability": 0.5,
                "profiles": angular_sources,
            },
        )
        channels = []
        for rate in rates:
            r = rate.record
            source_path = out / "source" / (file_stem(r["id"]) + ".txt")
            production_path = out / "production" / (file_stem(r["id"]) + ".txt")
            write_rate(source_path, rate, grid, source=True)
            points = write_rate(production_path, rate, grid)
            channel = {
                k: v for k, v in r.items() if k not in ("table", "header", "fit")
            }
            channel.update(
                cross_section=str(production_path.relative_to(data_dir)),
                source_cross_section=str(source_path.relative_to(data_dir)),
                continuation=rate.continuation,
                points=points,
                default=r["id"]
                not in ["attachment_muroran", "three_body_phelps", "three_body_biagi"],
            )
            if r["type"] == "excitation":
                key = r["id"]
                if r["family"] == "vibrational":
                    key = "v01" if gas == "N2" and key == "v01" else "vibrational"
                if key.startswith("schumann_runge"):
                    key = "schumann_runge"
                if key in profiles:
                    angle_file = out / "angular" / (file_stem(key) + ".angular")
                    channel["angular_model"] = key
                else:
                    angle_file = root / "elastic_dcs.txt"
                    channel["angular_model"] = "normalized_elastic_IAA"
                channel["differential_cross_section"] = str(
                    angle_file.relative_to(data_dir)
                )
            channels.append(channel)
        write_json(
            out / "sources.json",
            {
                "format": "warpx_iaa_rate_sources_v1",
                "archive_sha256": ARCHIVE_SHA256,
                "records": records,
            },
        )
        for model in ["iaa_thesis_2023", "elmolcs_b8643810"]:
            path = data_dir / "IAA" / "RBEQ" / f"{gas}_{model}.txt"
            if not path.is_file():
                raise ValueError(f"Generate the matching RBEQ table first: {path}")
            original = path.read_text()
            if "outside_energy_range" not in original:
                path.write_text(
                    "# energy_min_eV = 0\n# energy_max_eV = 1000000000\n"
                    "# outside_energy_range = error\n" + original
                )
            channels.append(
                {
                    "id": "ionization_" + model,
                    "type": "ionization",
                    "family": "ionization",
                    "cross_section": str(path.relative_to(data_dir)),
                    "rbeq_model": model,
                    "rbeq_target": gas,
                    "alternative_group": "ionization",
                    "loss_eV": float(np.loadtxt(path)[0, 0]),
                    "units": "m2",
                    "default": model == "iaa_thesis_2023",
                }
            )
        manifest["targets"][gas] = {
            "elastic": str(elastic_file.relative_to(data_dir)),
            "elastic_dcs": str((root / "elastic_dcs.txt").relative_to(data_dir)),
            "rotation": str(
                (root / "reciprocal_hybrid_300K" / "thermal_rotation.rot").relative_to(
                    data_dir
                )
            ),
            "master_energy_points": len(grid),
            "native_float32_energy_points": len(np.unique(grid.astype(np.float32))),
            "channels": channels,
        }
        print(
            f"{gas}: {len(channels)} channels including alternatives; {len(grid)} union knots; {len(profiles)} shared angular profiles",
            flush=True,
        )
    destination = data_dir / "IAA" / "air"
    destination.mkdir(parents=True, exist_ok=True)
    write_json(destination / "manifest.json", manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--attachment", type=Path, required=True)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="warpx-data MCC_cross_sections directory",
    )
    parser.add_argument(
        "--integral-only",
        action="store_true",
        help="Development export; do not publish without angular profiles",
    )
    args = parser.parse_args()
    export(args.archive, args.attachment, args.data_dir, args.integral_only)


if __name__ == "__main__":
    main()
