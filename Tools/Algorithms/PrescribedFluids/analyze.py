#!/usr/bin/env python3
"""Summarize ensemble error, noise, populations, spectra and measured cost.

The field baseline is the analytic fluid on the same mesh. Independent continuum
and PJG quadratures are reported separately; finite-domain and mesh errors are
not absorbed into a fitted tolerance. Confidence intervals describe seed means,
not uncertainty in the underlying air cross sections.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from gaussian_reference import gaussian_fields
from scipy.constants import c, e, proton_mass
from scipy.stats import t as student_t


def weighted_norm(array, cells):
    """Cylindrical volume norm, excluding a duplicated periodic endpoint."""
    nr, nz = cells
    r = np.arange(array.shape[0], dtype=float)
    if array.shape[0] == nr + 1:
        r[0] = 1 / 6  # Verboncoeur volume divided by 2*pi*dr*dr*dz.
    else:
        r += 0.5
    return np.sqrt(np.sum(array[:, :nz] ** 2 * r[:, None]))


def interval(values):
    values = np.asarray(values, dtype=float)
    n = values.size
    half = (
        student_t.ppf(0.995, n - 1) * values.std(ddof=1) / np.sqrt(n) if n > 1 else None
    )
    return {
        "mean": float(values.mean()),
        "ci99_half_width": None if half is None else float(half),
        "samples": n,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    cases = []
    for path in sorted(args.directory.glob("*/result.json")):
        data = json.loads(path.read_text())
        data["arrays"] = dict(np.load(path.with_suffix(".npz")))
        data["path"] = str(path)
        cases.append(data)
    if not cases:
        parser.error("No result.json files found")
    reference = next(
        case
        for case in cases
        if case["parameters"]["beam"] == "fluid"
        and case["parameters"]["ions"] == "fluid"
    )
    settings = reference["parameters"]
    sigma_r, sigma_t = 0.002, 25e-12
    gamma = 1 + 800e6 * e / (proton_mass * c**2)
    velocity = c * np.sqrt(1 - gamma**-2)
    sigma_z = velocity * sigma_t
    charge = 0.6 * np.sqrt(2 * np.pi) * sigma_t
    rmax = settings.get("radial_sigmas", 8) * sigma_r
    zmax = settings.get("longitudinal_sigmas", 12) * sigma_z
    continuum = {}
    nr, nz = settings["cells"]
    for component, index in [("Er", 0), ("Ez", 1), ("Btheta", 2)]:
        shape = reference["arrays"]["0_" + component].shape
        r = (np.arange(shape[0]) + (0 if shape[0] == nr + 1 else 0.5)) * rmax / nr
        z = (
            np.arange(shape[1]) + (0 if shape[1] == nz + 1 else 0.5)
        ) * 2 * zmax / nz - zmax
        mask = (r[:, None] < 3 * sigma_r) & (np.abs(z[None, :]) < 3 * sigma_z)
        field = gaussian_fields(
            r[:, None], z[None, :], sigma_r, sigma_z, charge, velocity
        )[index]
        continuum[component] = (field, mask)
    # Reject mixing grids, physical inputs or solver methods in one ensemble.
    for case in cases:
        for key in ["ranks", "backend", "revision"]:
            if case[key] != reference[key]:
                raise ValueError(f"Incompatible {key}: {case['path']}")
        for key in [
            "radial_sigmas",
            "longitudinal_sigmas",
            "implicit_deposition",
            "mcc_all",
            "max_grid_size",
        ]:
            if case["parameters"].get(key) != settings.get(key):
                raise ValueError(f"Incompatible {key}: {case['path']}")
        for key in [
            "cells",
            "steps",
            "dt",
            "mode",
            "solver",
            "gas_density",
            "weight",
            "cap",
            "shape",
            "mcc",
            "no_self_fields",
            "subcycles",
            "source_resolution",
            "electron_ppc",
            "temperature",
        ]:
            if case["parameters"][key] != settings[key]:
                raise ValueError(f"Incompatible {key}: {case['path']}")
    groups = defaultdict(list)
    for case in cases:
        p = case["parameters"]
        groups[(p["beam"], p["ions"], p["ppc"])].append(case)
    rows = []
    step = settings["steps"]
    for (beam, ions, ppc), members in groups.items():
        row = dict(beam=beam, ions=ions, ppc=ppc)
        row["timestep_s"] = interval(
            [case["ordinary_step_median_s"] for case in members]
        )
        row["initialization_s"] = interval(
            [case["initialization_s"] for case in members]
        )
        row["electron_macroparticles"] = interval(
            [case["history"][-1]["macroparticles"]["electrons"] for case in members]
        )
        row["particle_payload_bytes"] = interval(
            [sum(case["particle_payload_bytes"].values()) for case in members]
        )
        row["fluid_density_bytes"] = interval(
            [sum(case["fluid_density_bytes"].values()) for case in members]
        )
        for key in ["particle_payload_bytes", "fluid_density_bytes"]:
            row[key + "_by_species"] = {
                species: interval([case[key][species] for case in members])
                for species in members[0][key]
            }
        for key in ["checkpoint_step_s", "checkpoint_bytes", "device_reserved_bytes"]:
            values = [case[key] for case in members]
            if all(value is not None for value in values):
                row[key] = interval(values)
        if not settings["no_self_fields"]:
            for component, (exact, mask) in continuum.items():
                scale = weighted_norm(exact * mask, settings["cells"])
                row["continuum_" + component + "_relative_rms_error"] = interval(
                    [
                        weighted_norm(
                            (case["arrays"]["0_" + component] - exact) * mask,
                            settings["cells"],
                        )
                        / scale
                        for case in members
                    ]
                )
        for field in ["beam", "electrons", "Er", "Ez", "Btheta"]:
            key = f"{step}_{field}"
            baseline = reference["arrays"][key]
            scale = weighted_norm(baseline, settings["cells"])
            samples = np.array([case["arrays"][key] for case in members])
            mean = samples.mean(axis=0)
            if scale == 0:
                continue
            row[field + "_relative_rms_error"] = interval(
                [
                    weighted_norm(array - baseline, settings["cells"]) / scale
                    for array in samples
                ]
            )
            row[field + "_relative_seed_noise"] = (
                float(
                    np.sqrt(
                        sum(
                            weighted_norm(array - mean, settings["cells"]) ** 2
                            for array in samples
                        )
                        / (len(samples) - 1)
                    )
                    / scale
                )
                if len(samples) > 1
                else None
            )
        for species in reference["history"][-1]["physical"]:
            row[species + "_number"] = interval(
                [case["history"][-1]["physical"][species] for case in members]
            )
        row["electron_energy_J"] = interval(
            [case["history"][-1]["electron_energy_J"] for case in members]
        )
        row["pair_charge_residual"] = interval(
            [
                (
                    case["history"][-1]["physical"]["N2plus"]
                    + case["history"][-1]["physical"]["O2plus"]
                    - case["history"][-1]["physical"]["Ominus"]
                    - case["history"][-1]["physical"].get("O2minus", 0.0)
                    - case["history"][-1]["physical"]["electrons"]
                )
                / max(1.0, case["history"][-1]["physical"]["electrons"])
                for case in members
            ]
        )
        histograms = [case["arrays"]["electron_histogram"] for case in members]
        total = sum(histograms)
        if total.sum() > 0:
            row["electron_cdf"] = (np.cumsum(total) / total.sum()).tolist()
        if settings["mode"] == "source":
            sys.path.insert(
                0, str(Path(__file__).resolve().parents[1] / "ProtonImpactIonization")
            )
            from calibrated_pjg import total_cross_section

            qe, mp = e, proton_mass
            gamma = 1 + 800e6 * qe / (mp * c**2)
            speed = c * np.sqrt(1 - gamma**-2)
            number = 0.6 * 25e-12 * np.sqrt(2 * np.pi) / qe
            for target, fraction in [("N2", 0.79), ("O2", 0.21)]:
                expected = (
                    number
                    * settings["gas_density"]
                    * fraction
                    * speed
                    * step
                    * settings["dt"]
                )
                expected *= float(total_cross_section(target, 800e6)) * 1e-4
                row[target + "_primary_yield_relative_error"] = interval(
                    [
                        (
                            case["history"][-1]["physical"][target + "plus"]
                            + case["history"][-1]["pending"]["p_" + target]
                        )
                        / expected
                        - 1
                        for case in members
                    ]
                )
        rows.append(row)
    output = {
        "settings": settings,
        "revision": sorted({case["revision"] for case in cases}),
        "backend": sorted({case["backend"] for case in cases}),
        "ranks": sorted({case["ranks"] for case in cases}),
        "groups": rows,
    }
    (args.directory / "summary.json").write_text(json.dumps(output, indent=2) + "\n")
    # Scientific plots are standalone artifacts, not screenshots of an interactive dashboard.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for group_index, (beam, ions) in enumerate(
        dict.fromkeys((row["beam"], row["ions"]) for row in rows)
    ):
        selected = sorted(
            [row for row in rows if row["beam"] == beam and row["ions"] == ions],
            key=lambda row: row["ppc"],
        )
        label = f"{beam} beam, {ions} ions"
        color = f"C{group_index}"
        x = [row["timestep_s"]["mean"] for row in selected]
        xerr = [row["timestep_s"]["ci99_half_width"] or 0 for row in selected]
        y = [
            row.get("beam_relative_rms_error", {"mean": 0})["mean"] for row in selected
        ]
        yerr = [
            row.get("beam_relative_rms_error", {}).get("ci99_half_width", 0) or 0
            for row in selected
        ]
        axes[0].errorbar(x, y, xerr=xerr, yerr=yerr, fmt="o-", color=color, label=label)
        if beam == "fluid":
            axes[1].axhline(x[0], linestyle="--", color=color, label=label)
        else:
            axes[1].errorbar(
                [row["ppc"] for row in selected],
                x,
                yerr=xerr,
                fmt="o-",
                color=color,
                label=label,
            )
    axes[0].set(
        xlabel="Median timestep wall time (s)",
        ylabel="Beam density relative RMS error",
        xscale="log",
    )
    axes[0].set_yscale("symlog", linthresh=1e-5, linscale=0.4)
    axes[0].set_ylim(bottom=0)
    axes[1].set(
        xlabel="Beam macroparticles per cell (input)",
        ylabel="Median timestep wall time (s)",
        xscale="log",
        yscale="log",
    )
    for axis in axes:
        axis.grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, title="99% confidence intervals", title_fontsize=8)
    fig.savefig(args.directory / "noise_cost.png", dpi=180)
    fig.savefig(args.directory / "noise_cost.pdf")
    plt.close(fig)
    print(f"Summarized {len(cases)} runs in {args.directory / 'summary.json'}")


if __name__ == "__main__":
    main()
