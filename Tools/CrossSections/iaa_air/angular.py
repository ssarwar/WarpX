# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Positive, source-resolved excitation angular profiles prepared offline."""

import re

import numpy as np

from .layout import angular_paths
from .sources import N2_STATES, O2_STATES, digest, file_stem

THETA = np.linspace(0, np.pi, 181)
HALF_SINE = np.sin(THETA / 2)


def central(token):
    token = token.replace("−", "-")
    match = re.fullmatch(r"\(([^±)]+)±[^)]+\)(e[+-]?\d+)", token)
    if match:
        return float(match[1] + match[2])
    try:
        return float(token.split("±")[0].split("(")[0])
    except ValueError:
        return np.nan


def read_dcs(text):
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if not line.startswith("PROCESS:") or not re.search(
            "Inelastic|Vibrational", line
        ):
            continue
        j = i
        while j < len(lines) and not re.match(r"-{5,}", lines[j]):
            j += 1
        k = j + 1
        while k < len(lines) and not re.match(r"-{5,}", lines[k]):
            k += 1
        header = "\n".join(lines[i:j])
        columns = next((x for x in lines[i:j] if x.startswith("COLUMNS:")), "").split(
            "|"
        )[1:]
        energies = [re.search(r"^\s*(?:DCS\s+)?([0-9]+\.?[0-9]*)", c) for c in columns]
        rows = []
        for row in lines[j + 1 : k]:
            fields = row.partition("#")[0].split()
            if len(fields) != len(columns) + 1:
                continue
            try:
                angle = float(fields[0])
            except ValueError:
                continue
            rows.append([angle, *(central(t) for t in fields[1:])])
        if not rows:
            continue
        array = np.asarray(rows)
        for c, energy in enumerate(energies, 1):
            if energy is None:
                continue
            valid = np.isfinite(array[:, c]) & (array[:, c] >= 0)
            # The aggregate inherits extrapolated component endpoints, without
            # the parentheses used to mark them in the component source rows.
            if "Herzberg Pseudo-continuum" in header:
                valid &= array[:, 0] <= 156
            yield {
                "original_process": line[9:],
                "header": header,
                "line": i + 1,
                "energy_eV": float(energy[1]),
                "angles_deg": array[valid, 0].tolist(),
                "values": array[valid, c].tolist(),
            }


def repair_shyn_sweeney_b(rows):
    """Restore two decimal points verified visually in Table I(b), p. 1007."""
    for energy, angle, original, corrected in [(5, 12, 94, 9.4), (15, 36, 89, 8.9)]:
        matches = [row for row in rows if row["energy_eV"] == energy]
        if len(matches) != 1 or matches[0]["angles_deg"].count(angle) != 1:
            raise ValueError("Missing or ambiguous Shyn-Sweeney O2 b-state repair")
        row = matches[0]
        index = row["angles_deg"].index(angle)
        if row["values"][index] != original:
            raise ValueError("Unexpected archived Shyn-Sweeney O2 b-state value")
        row["values"][index] = corrected
        row["source_repairs"] = [
            {
                "kind": "decimal_transcription",
                "angle_deg": angle,
                "original_value": original,
                "corrected_value": corrected,
                "units": "1e-20 cm2/sr",
                "reference": "Shyn and Sweeney, Phys. Rev. A 47, 1006 (1993)",
                "doi": "10.1103/PhysRevA.47.1006",
                "location": "Table I(b), p. 1007",
                "verification": "visual inspection of the published table; "
                "the PDF text layer also omits these decimal points",
                "source_pdf_sha256": (
                    "81ae455b82ba89b7a39d34364d10a7feb539f8bb3ac27070d699b226f17fc3c9"
                ),
            }
        ]


def normalize(values):
    # Exact integral of a DCS linear in y=sin(theta/2), dOmega=8*pi*y*dy.
    y = HALF_SINE
    h = np.diff(y)
    slope = np.diff(values, axis=-1) / h
    integrals = (
        y[:-1] * values[..., :-1] * h
        + (y[:-1] * slope + values[..., :-1]) * h * h / 2
        + slope * h * h * h / 3
    )
    total = 8 * np.pi * integrals.sum(axis=-1)
    return values / total[..., None]


class Elastic:
    def __init__(self, path):
        table = np.loadtxt(path)
        self.energy = table[:, 0]
        self.values = table[:, 1:]
        self.y = np.sin(np.linspace(0, np.pi, self.values.shape[1]) / 2)

    def __call__(self, energy):
        i = np.clip(np.searchsorted(self.energy, energy) - 1, 0, len(self.energy) - 2)
        f = np.clip(
            np.log(energy / self.energy[i])
            / np.log(self.energy[i + 1] / self.energy[i]),
            0,
            1,
        )
        values = (1 - f) * self.values[i] + f * self.values[i + 1]
        return normalize(np.interp(HALF_SINE, self.y, values))


def complete(row, elastic):
    angles = np.asarray(row["angles_deg"])
    values = np.asarray(row["values"])
    if len(angles) < 8:
        return None
    coverage = (np.cos(np.deg2rad(angles[0])) - np.cos(np.deg2rad(angles[-1]))) / 2
    if coverage < 0.75 or np.diff(angles).max() > 20 or not values.max() > 0:
        return None
    y = np.sin(np.deg2rad(angles) / 2)
    base = elastic(row["energy_eV"])
    pdf = np.interp(HALF_SINE, y, values)
    constant = pdf.copy()
    for low, endpoint in [(True, 0), (False, -1)]:
        mask = HALF_SINE < y[endpoint] if low else HALF_SINE > y[endpoint]
        boundary = np.interp(y[endpoint], HALF_SINE, base)
        if boundary <= 0:
            return None
        pdf[mask] = values[endpoint] * base[mask] / boundary
    pdf = normalize(pdf)
    constant = normalize(constant)
    measure = 2 * np.pi * np.sin(THETA)
    cap = 0.0
    for mask in [HALF_SINE <= y[0], HALF_SINE >= y[-1]]:
        cap += np.trapezoid((pdf * measure)[mask], THETA[mask])
    mu = np.cos(THETA)
    sensitivity = [
        float(abs(np.trapezoid((pdf - constant) * measure * p, THETA)))
        for p in [mu, (3 * mu**2 - 1) / 2]
    ]
    row["completion"] = {
        "model": "boundary_matched_elastic",
        "solid_angle_coverage": float(coverage),
        "modeled_probability": float(cap),
        "moment_sensitivity_to_constant_caps": sensitivity,
    }
    if cap >= 0.5:
        return None
    return pdf


def smooth(energy, lo, hi):
    x = np.clip(np.log(energy / lo) / np.log(hi / lo), 0, 1)
    return x * x * (3 - 2 * x)


class Profile:
    def __init__(self, gas, key, elastic, records, loss):
        self.gas, self.key, self.elastic = gas, key, elastic
        self.rows, pdfs = [], []
        self.references = []
        for row in sorted(records, key=lambda r: r["energy_eV"]):
            pdf = complete(row, elastic)
            row["accepted"] = pdf is not None
            self.references.append(row)
            if pdf is not None:
                self.rows.append(row)
                pdfs.append(pdf)
        if self.rows:
            self.energy = np.array([r["energy_eV"] for r in self.rows])
            if len(self.energy) < 2 or np.any(np.diff(self.energy) <= 0):
                raise ValueError(
                    f"Insufficient or ambiguous angular series for {gas}/{key}"
                )
            self.pdfs = np.asarray(pdfs)
            self.lo = max(loss, self.energy[0] ** 2 / self.energy[1])
            self.hi = self.energy[-1] ** 2 / self.energy[-2]
            self.knots = np.r_[self.lo, self.energy, self.hi]
        elif key not in ("vibrational", "v01"):
            raise ValueError(f"Missing angular mapping for {gas}/{key}")
        if key in ("vibrational", "v01"):
            self.lo = 1e-5
            if gas == "N2":
                self.hi = self.hi if self.rows else 10.0
                self.knots = np.unique(
                    np.r_[self.lo, 5.0, 10.0, self.knots if self.rows else [self.hi]]
                )
            else:
                self.hi = 40.0
                self.knots = np.array([self.lo, 20.0, 40.0])

    def measured(self, energy):
        i = np.clip(np.searchsorted(self.energy, energy) - 1, 0, len(self.energy) - 2)
        f = np.clip(
            (energy - self.energy[i]) / (self.energy[i + 1] - self.energy[i]), 0, 1
        )
        return (1 - f) * self.pdfs[i] + f * self.pdfs[i + 1]

    def __call__(self, energy):
        base = self.elastic(energy)
        if self.key in ("vibrational", "v01"):
            if self.gas == "N2":
                mu = np.cos(THETA)
                low = normalize(1 - 3 * mu**2 + 14 * mu**4 / 3)
                if energy <= 5:
                    return low
                if energy < 10:
                    w = smooth(energy, 5, 10)
                    high = self.measured(energy) if self.rows else base
                    return (1 - w) * low + w * high
                if not self.rows:
                    return base
            else:
                low = np.full(len(THETA), 1 / (4 * np.pi))
                w = smooth(energy, 20, 40)
                return (1 - w) * low + w * base
        if energy < self.energy[0]:
            w = smooth(energy, self.lo, self.energy[0])
            return (1 - w) * base + w * self.pdfs[0]
        if energy > self.energy[-1]:
            w = smooth(energy, self.energy[-1], self.hi)
            return (1 - w) * self.pdfs[-1] + w * base
        return self.measured(energy)


def select_profiles(archive, gas, elastic, records):
    text = archive.read(f"Data/dcs/{gas}/Diff cross section.e-{gas}")
    data = list(read_dcs(text))
    for row in data:
        row["source_sha256"] = digest(text.encode())
    selection = {}
    for record in records:
        key = record["id"]
        if record["family"] == "vibrational":
            key = "v01" if gas == "N2" and key == "v01" else "vibrational"
        elif key.startswith("schumann_runge"):
            key = "schumann_runge"
        elif record["family"] != "electronic":
            continue
        if key in selection:
            continue
        rows = []
        if key == "vibrational":
            pass
        elif key == "v01":
            rows = [
                r
                for r in data
                if "v1)" in r["original_process"]
                and ("Linert" in r["header"] or "Truhlar et al.(1977)" in r["header"])
            ]
        else:
            states = N2_STATES if gas == "N2" else O2_STATES
            state = next(s for s, name in states.items() if name == key)
            raw_state = state.replace("Σ", "SIG").replace("Π", "PI").replace("Δ", "DEL")
            source = (
                (
                    "Malone et al (2009)"
                    if key in ("C", "E", "a_double_prime")
                    else "Khakoo et al (2008)"
                    if key in ("b", "c3", "o3", "b_prime", "c4_prime", "G", "F")
                    else "Khakoo et al (2005)"
                )
                if gas == "N2"
                else (
                    "Shyn&Sweeney (1993)"
                    if key in ("a", "b")
                    else "Herzberg Pseudo-continuum"
                    if key == "herzberg"
                    else "Suzuki et al (2010)"
                )
            )
            if gas == "O2" and key == "b":
                raw_state = "b1SIGg"  # Source omits the reflection-symmetry suffix.
            # These two source headers contain a documented g/u transcription error.
            if source == "Malone et al (2009)" and key in ("E", "a_double_prime"):
                raw_state = raw_state.replace("SIGg", "SIGu")
            rows = [
                r
                for r in data
                if source in r["header"] and f"({raw_state})" in r["original_process"]
            ]
            if gas == "O2" and key == "b":
                repair_shyn_sweeney_b(rows)
        selection[key] = Profile(gas, key, elastic, rows, record["loss_eV"])
    return selection


def tabulate(profile):
    # Runtime interpolation is in log(E); refine its error against the linear-E
    # reference and explicitly modeled bridges, retaining measured energy rows.
    knots = np.unique(
        np.r_[
            profile.knots,
            profile.elastic.energy[
                (profile.elastic.energy > profile.lo)
                & (profile.elastic.energy < profile.hi)
            ],
        ]
    )
    measure = 2 * np.pi * np.sin(THETA)
    for _ in range(20):
        pdfs = np.array([profile(e) for e in knots])
        split = []
        for i, (lo, hi) in enumerate(zip(knots[:-1], knots[1:])):
            for f in (0.2113248654, 0.5, 0.7886751346):
                e = np.exp((1 - f) * np.log(lo) + f * np.log(hi))
                delta = ((1 - f) * pdfs[i] + f * pdfs[i + 1] - profile(e)) * measure
                # Total variation also bounds CDF error and any |P_l|<=1 moment.
                if np.trapezoid(np.abs(delta), THETA) > 3e-4:
                    split.append(np.sqrt(lo * hi))
                    break
        if not split:
            return knots, pdfs
        knots = np.sort(np.r_[knots, split])
    raise ValueError(f"Angular refinement did not converge: {profile.key}")


def write_profile(directory, profile):
    energy, pdf = tabulate(profile)
    descriptor, pdf_file = angular_paths(profile.gas, profile.key)
    name = descriptor.stem
    np.savetxt(
        directory / pdf_file,
        np.column_stack((energy, pdf)),
        fmt=["%.17e"] + ["%.12e"] * len(THETA),
        header=f"IAA {profile.gas} {file_stem(profile.key)}: normalized excitation angular PDF (sr^-1)\n"
        "Uniform theta grid 0..pi; runtime interpolation linear in sin(theta/2) and log(E).\n"
        "Measured coverage, endpoint completion and transition assumptions: sources.json",
    )
    (directory / descriptor).write_text(
        f"WARPX_EXCITATION_ANGULAR_V1\n{profile.gas}\n{pdf_file.name}\n../../elastic/differential_cross_section.txt\n"
        f"{profile.lo:.17g} {profile.hi:.17g}\n"
    )
    return {
        "profile": profile.key,
        "file_stem": name,
        "profile_file": str(pdf_file),
        "descriptor_file": str(descriptor),
        "active_range_eV": [profile.lo, profile.hi],
        "reference_energy_interpolation": "linear",
        "runtime_energy_interpolation": "logarithmic",
        "transition": "smoothstep_log_energy",
        "rows": profile.references,
        "energy_rows": len(energy),
        "angle_points": len(THETA),
    }
