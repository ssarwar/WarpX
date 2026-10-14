# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Gas-local paths for physical IAA data, independent of the export models."""

import json
from pathlib import Path

ELASTIC_RATE = Path("elastic/cross_section.txt")
ELASTIC_DCS = Path("elastic/differential_cross_section.txt")
ROTATION_BUNDLE = Path("rotational/elastic_rotation_300K/thermal_rotation.rot")

N2_ELECTRONIC = {
    "A": "A_3Sigma_u_plus",
    "B": "B_3Pi_g",
    "W": "W_3Delta_u",
    "B_prime": "B_prime_3Sigma_u_minus",
    "a_prime": "a_prime_1Sigma_u_minus",
    "a": "a_1Pi_g",
    "w": "w_1Delta_u",
    "C": "C_3Pi_u",
    "E": "E_3Sigma_g_plus",
    "a_double_prime": "a_double_prime_1Sigma_g_plus",
    "b": "b_1Pi_u",
    "c3": "c3_1Pi_u",
    "o3": "o3_1Pi_u",
    "b_prime": "b_prime_1Sigma_u_plus",
    "c4_prime": "c4_prime_1Sigma_u_plus",
    "G": "G_3Pi_u",
    "F": "F_3Pi_u",
}
O2_ELECTRONIC = {
    "a": "a_1Delta_g",
    "b": "b_1Sigma_g_plus",
    "herzberg": "herzberg_continuum",
    "schumann_runge_onset": "schumann_runge_6.12_eV",
    "schumann_runge_bulk": "schumann_runge_8.60_eV",
}


def family_name(record):
    """Complementary excitations belong with the electronic bands."""
    family = record["family"]
    return "electronic" if family == "complementary" else family


def channel_stem(gas, channel):
    if gas not in ("N2", "O2"):
        raise ValueError(f"Unsupported IAA target: {gas}")
    if channel.startswith("v") and channel[1:].isdigit():
        return "v0_to_v" + channel[1:]
    if channel.startswith("complement_"):
        return channel + "_eV"
    if channel.startswith("attachment_"):
        return "dissociative_" + channel.removeprefix("attachment_")
    if channel.startswith("ionization_"):
        return "rbeq_" + channel.removeprefix("ionization_")
    states = N2_ELECTRONIC if gas == "N2" else O2_ELECTRONIC
    return states.get(channel, channel)


def rate_path(gas, record, *, source=False):
    root = Path(family_name(record))
    if source:
        root /= "source_bounded"
    return root / (channel_stem(gas, record["id"]) + ".txt")


def angular_family(profile):
    return "vibrational" if profile in ("v01", "vibrational") else "electronic"


def angular_paths(gas, profile):
    name = (
        "shared_bound_states"
        if profile == "vibrational"
        else channel_stem(gas, profile)
    )
    directory = Path(angular_family(profile)) / "angular"
    return directory / (name + ".angular"), directory / (name + "_pdf.txt")


def load_manifest(data_dir, gases=("N2", "O2")):
    """Read gas-local manifests and expose MCC-root-relative paths to tools.

    The legacy aggregate manifest remains readable for pinned older datasets.
    Runtime input loading is independently implemented in pywarpx.iaa_air.
    """
    data_dir = Path(data_dir)
    targets = {}
    for gas in gases:
        path = data_dir / gas / "IAA/manifest.json"
        if not path.is_file():
            legacy = json.loads((data_dir / "IAA/air/manifest.json").read_text())
            if legacy.get("format") != "warpx_iaa_air_v1":
                raise ValueError("Unsupported IAA manifest format")
            targets[gas] = legacy["targets"][gas]
            continue
        entry = json.loads(path.read_text())
        if entry.get("format") != "warpx_iaa_gas_v1" or entry.get("target") != gas:
            raise ValueError(f"Invalid IAA gas manifest: {path}")
        prefix = Path(gas) / "IAA"
        for key in ("elastic", "elastic_dcs", "rotation"):
            entry[key] = str(prefix / entry[key])
        for channel in entry["channels"]:
            for key in (
                "cross_section",
                "source_cross_section",
                "differential_cross_section",
            ):
                if key in channel:
                    channel[key] = str(prefix / channel[key])
        targets[gas] = entry
    return {"targets": targets}


def load_sources(data_dir, gas):
    """Return source records in manifest order, preserving process ordering."""
    root = Path(data_dir) / gas / "IAA"
    legacy = root / "air/sources.json"
    if legacy.is_file():
        return json.loads(legacy.read_text())["records"]
    records = {}
    for path in sorted(root.glob("*/sources.json")):
        for record in json.loads(path.read_text())["records"]:
            if record["id"] in records:
                raise ValueError(f"Duplicate IAA source record: {record['id']}")
            records[record["id"]] = record
    manifest = load_manifest(data_dir, (gas,))["targets"][gas]
    expected = [c["id"] for c in manifest["channels"] if c["type"] != "ionization"]
    if set(expected) != set(records):
        raise ValueError(f"IAA source records do not match the {gas} manifest")
    return [records[channel] for channel in expected]


def load_angular_sources(data_dir, gas):
    root = Path(data_dir) / gas / "IAA"
    legacy = root / "air/angular/angular_sources.json"
    if legacy.is_file():
        profiles = json.loads(legacy.read_text())["profiles"]
        for profile in profiles:
            stem = profile["file_stem"]
            profile["profile_file"] = f"air/angular/{stem}.txt"
            profile["descriptor_file"] = f"air/angular/{stem}.angular"
        return profiles
    profiles = []
    for family in ("vibrational", "electronic"):
        source = root / family / "angular/sources.json"
        profiles.extend(json.loads(source.read_text())["profiles"])
    return profiles
