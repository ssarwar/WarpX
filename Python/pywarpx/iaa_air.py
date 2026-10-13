# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Load the production electron-air process manifest from warpx-data.

This module reads prepared physical data. It does not generate cross sections
and has no dependency on elmolcs, NumPy, or the offline exporters.
"""

import json
from pathlib import Path


def load_iaa_air_processes(
    data_dir,
    target,
    *,
    ion_species,
    attachment_species=None,
    third_body_density=None,
    ionization_model="iaa_thesis_2023",
    dissociative_attachment="iaa",
    three_body_attachment="muroran",
):
    """Return ``scattering_processes`` for a PICMI ``MCCCollisions`` object.

    ``data_dir`` is the warpx-data MCC_cross_sections directory. For O2,
    ``attachment_species`` maps ``dissociative`` and ``three_body`` to the O-
    and O2- product species. ``third_body_density`` is the O2 density in m^-3,
    independently of the total gas density. The supplied rotation bundle is
    for a 300 K rotational bath. Neutral populations are not evolved.
    """
    root = Path(data_dir).resolve()
    manifest = json.loads((root / "IAA/air/manifest.json").read_text())
    if manifest.get("format") != "warpx_iaa_air_v1":
        raise ValueError("Unsupported IAA air manifest format")
    if target not in ("N2", "O2"):
        raise ValueError("IAA air target must be N2 or O2")
    selected = {
        "ionization": "ionization_" + ionization_model,
        "dissociative_attachment": "attachment_" + dissociative_attachment,
        "three_body_attachment": "three_body_" + three_body_attachment,
    }
    description = manifest["targets"][target]

    def path(name):
        value = root / name
        if not value.is_file():
            raise FileNotFoundError(value)
        return str(value)

    processes = {
        "elastic": {
            "cross_section": path(description["elastic"]),
            "scattering_angle_model": "IAA",
            "rotation_model": "reciprocal_hybrid",
            "rotation_file": path(description["rotation"]),
            "rotational_temperature": 300.0,
        }
    }
    groups = set()
    for channel in description["channels"]:
        group = channel.get("alternative_group")
        if group and channel["id"] != selected[group]:
            continue
        if group:
            if group in groups:
                raise ValueError(f"Duplicate IAA alternative group: {group}")
            groups.add(group)
        process = {"cross_section": path(channel["cross_section"])}
        kind = channel["type"]
        if kind == "excitation":
            process.update(
                energy=channel["loss_eV"],
                scattering_angle_model="IAA",
                differential_cross_section=path(channel["differential_cross_section"]),
            )
        elif kind == "ionization":
            process.update(
                energy=channel["loss_eV"],
                species=ion_species,
                energy_sharing_model="RBEQ",
                scattering_angle_model="IAA",
                rbeq_model=channel["rbeq_model"],
                rbeq_target=target,
            )
        elif kind == "attachment":
            product = "three_body" if channel["units"] == "m5" else "dissociative"
            if not attachment_species or product not in attachment_species:
                raise ValueError(
                    f"O2 requires the {product} attachment product species"
                )
            process.update(
                species=attachment_species[product],
                cross_section_units=channel["units"],
            )
            if product == "three_body":
                if third_body_density is None:
                    raise ValueError("Supply the O2 third-body density explicitly")
                process["third_body_density"] = third_body_density
        else:
            raise ValueError(f"Unsupported IAA process: {kind}")
        processes[kind + "_" + channel["id"]] = process
    required = {"ionization"} if target == "N2" else set(selected)
    if groups != required:
        raise ValueError(f"Unknown or missing IAA model selection: {selected}")
    return processes
