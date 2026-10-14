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
from iaa_air.layout import (
    ELASTIC_DCS,
    ELASTIC_RATE,
    ROTATION_BUNDLE,
    angular_family,
    angular_paths,
    family_name,
    rate_path,
)
from iaa_air.rates import MAX_ENERGY, Rate, master_grid, write_rate
from iaa_air.sources import ARCHIVE_SHA256, Archive, digest, read_sources


def write_json(path, content):
    path.write_text(
        json.dumps(content, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )


def export(archive_path, attachment, data_dir):
    archive = Archive(archive_path)
    sources = read_sources(archive, attachment)
    manifest = {
        "format": "warpx_iaa_gas_v1",
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
    }
    for gas, records in sources.items():
        root = data_dir / gas / "IAA"
        for family in sorted({family_name(r) for r in records}):
            (root / family / "source_bounded").mkdir(parents=True, exist_ok=True)
        for family in ("vibrational", "electronic"):
            (root / family / "angular").mkdir(parents=True, exist_ok=True)
        rates = [Rate(r) for r in records]
        elastic_file = root / ELASTIC_RATE
        grid = master_grid(rates, np.loadtxt(elastic_file)[:, 0])
        profiles = select_profiles(archive, gas, Elastic(root / ELASTIC_DCS), records)
        angular_sources = [write_profile(root, p) for p in profiles.values()]
        for family in ("vibrational", "electronic"):
            write_json(
                root / family / "angular/sources.json",
                {
                    "format": "warpx_iaa_angular_sources_v1",
                    "archive_sha256": ARCHIVE_SHA256,
                    "thesis": "Schmalzried 2023, section 11.3.2, equation 11.40",
                    "vibrational_transition_eV": [5, 10] if gas == "N2" else [20, 40],
                    "minimum_angles": 8,
                    "minimum_solid_angle_coverage": 0.75,
                    "maximum_angular_gap_degrees": 20,
                    "maximum_completion_probability": 0.5,
                    "target": gas,
                    "family": family,
                    "profiles": [
                        p
                        for p in angular_sources
                        if angular_family(p["profile"]) == family
                    ],
                },
            )
        channels = []
        for rate in rates:
            r = rate.record
            source_path = root / rate_path(gas, r, source=True)
            production_path = root / rate_path(gas, r)
            write_rate(source_path, rate, grid, source=True)
            points = write_rate(production_path, rate, grid)
            channel = {
                k: v for k, v in r.items() if k not in ("table", "header", "fit")
            }
            channel.update(
                cross_section=str(production_path.relative_to(root)),
                source_cross_section=str(source_path.relative_to(root)),
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
                    angle_file = root / angular_paths(gas, key)[0]
                    channel["angular_model"] = key
                else:
                    angle_file = root / ELASTIC_DCS
                    channel["angular_model"] = "normalized_elastic_IAA"
                channel["differential_cross_section"] = str(
                    angle_file.relative_to(root)
                )
            channels.append(channel)
        for family in sorted({family_name(r) for r in records}):
            write_json(
                root / family / "sources.json",
                {
                    "format": "warpx_iaa_rate_sources_v1",
                    "archive_sha256": ARCHIVE_SHA256,
                    "target": gas,
                    "family": family,
                    "records": [r for r in records if family_name(r) == family],
                },
            )
        for model in ["iaa_thesis_2023", "elmolcs_b8643810"]:
            path = root / "ionization" / f"rbeq_{model}.txt"
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
                    "cross_section": str(path.relative_to(root)),
                    "rbeq_model": model,
                    "rbeq_target": gas,
                    "alternative_group": "ionization",
                    "loss_eV": float(np.loadtxt(path)[0, 0]),
                    "units": "m2",
                    "default": model == "iaa_thesis_2023",
                }
            )
        description = {
            "elastic": str(elastic_file.relative_to(root)),
            "elastic_dcs": str((root / ELASTIC_DCS).relative_to(root)),
            "rotation": str(ROTATION_BUNDLE),
            "master_energy_points": len(grid),
            "rounded_master_energy_points_float32": len(
                np.unique(grid.astype(np.float32))
            ),
            "channels": channels,
        }
        gas_manifest = dict(manifest, target=gas, **description)
        if gas == "N2":
            gas_manifest["defaults"] = {
                "ionization": manifest["defaults"]["ionization"]
            }
            gas_manifest["assumptions"] = manifest["assumptions"][:2]
            gas_manifest.pop("attachment_sha256")
        write_json(root / "manifest.json", gas_manifest)
        print(
            f"{gas}: {len(channels)} channels including alternatives; {len(grid)} union knots; {len(profiles)} shared angular profiles",
            flush=True,
        )


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
    args = parser.parse_args()
    export(args.archive, args.attachment, args.data_dir)


if __name__ == "__main__":
    main()
