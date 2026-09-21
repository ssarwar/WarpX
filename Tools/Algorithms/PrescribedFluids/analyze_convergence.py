#!/usr/bin/env python3
"""Analyze a convergence.py study without mixing discretizations or ensembles."""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from analyze import interval, weighted_norm
from gaussian_reference import gaussian_fields
from scipy.constants import c, e, proton_mass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "study", choices=["deposition", "source", "joint", "solvers", "continuum"]
    )
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    cases = {}
    for path in args.directory.glob("*/result.json"):
        result = json.loads(path.read_text())
        cases[path.parent.name] = (result, dict(np.load(path.with_suffix(".npz"))))
    rows = []
    if args.study == "deposition":
        for name, (result, arrays) in sorted(cases.items()):
            p = result["parameters"]
            if p["beam"] == "fluid":
                continue
            base = f"{p['solver']}-n{p['cells'][0]}-fluid-villasenor"
            reference = cases[base][1]
            errors = {}
            for field in ["beam", "Er", "Ez", "Btheta"]:
                key = f"{p['steps']}_{field}"
                errors[field] = weighted_norm(
                    arrays[key] - reference[key], p["cells"]
                ) / weighted_norm(reference[key], p["cells"])
            rows.append(dict(case=name, errors=errors))
    elif args.study == "continuum":
        sigma_r, sigma_t = 0.002, 25e-12
        gamma = 1 + 800e6 * e / (proton_mass * c**2)
        velocity = c * np.sqrt(1 - gamma**-2)
        sigma_z = velocity * sigma_t
        charge = 0.6 * np.sqrt(2 * np.pi) * sigma_t
        for name, (result, arrays) in sorted(cases.items()):
            p = result["parameters"]
            nr, nz = p["cells"]
            rmax, zmax = (
                p["radial_sigmas"] * sigma_r,
                p["longitudinal_sigmas"] * sigma_z,
            )
            errors = {}
            for index, field in enumerate(["Er", "Ez", "Btheta"]):
                actual = arrays["0_" + field]
                r = (
                    (np.arange(actual.shape[0]) + (actual.shape[0] == nr) / 2)
                    * rmax
                    / nr
                )
                z = (
                    np.arange(actual.shape[1]) + (actual.shape[1] == nz) / 2
                ) * 2 * zmax / nz - zmax
                # Restrict expensive continuum quadrature to the core being measured.
                ir = r < 3 * sigma_r
                iz = np.abs(z) < 3 * sigma_z
                exact = np.zeros_like(actual)
                exact[np.ix_(ir, iz)] = gaussian_fields(
                    r[ir, None], z[None, iz], sigma_r, sigma_z, charge, velocity
                )[index]
                mask = ir[:, None] & iz[None, :]
                errors[field] = weighted_norm(
                    (actual - exact) * mask, p["cells"]
                ) / weighted_norm(exact, p["cells"])
            rows.append(dict(case=name, errors=errors))
    else:
        sys.path.insert(
            0, str(Path(__file__).resolve().parents[1] / "ProtonImpactIonization")
        )
        from calibrated_pjg import total_cross_section

        groups = defaultdict(list)
        for name, case in cases.items():
            groups[name.rsplit("-seed", 1)[0]].append(case)
        gamma = 1 + 800e6 * e / (proton_mass * c**2)
        velocity = c * np.sqrt(1 - gamma**-2)
        number = 0.6 * 25e-12 * np.sqrt(2 * np.pi) / e
        sigma = {
            target: float(total_cross_section(target, 800e6)) * 1e-4
            for target in ["N2", "O2"]
        }
        for name, members in sorted(groups.items()):
            p = members[0][0]["parameters"]
            final = [result["history"][-1] for result, _ in members]
            row = dict(case=name, parameters=p)
            for species in ["electrons", "N2plus", "O2plus", "Ominus"]:
                row[species] = interval([state["physical"][species] for state in final])
            row["electron_energy_J"] = interval(
                [state["electron_energy_J"] for state in final]
            )
            for target, fraction in [("N2", 0.79), ("O2", 0.21)]:
                expected = (
                    number
                    * p["gas_density"]
                    * fraction
                    * velocity
                    * sigma[target]
                    * p["dt"]
                    * p["steps"]
                )
                errors = [
                    (
                        state["primary_source_budgets"]["p_" + target][0]
                        + state["pending"]["p_" + target]
                    )
                    / expected
                    - 1
                    for state in final
                ]
                # Existing PJG interpolation bound, set independently of this sweep.
                assert max(map(abs, errors)) < 1e-3, (name, target, errors)
                row[target + "_primary_yield_relative_error"] = interval(errors)
                row[target + "_pending_fraction"] = interval(
                    [state["pending"]["p_" + target] / expected for state in final]
                )
            hist = sum(arrays["electron_histogram"] for _, arrays in members)
            row["electron_cdf"] = (np.cumsum(hist) / hist.sum()).tolist()
            rows.append(row)
    output = args.directory / "comparison.json"
    output.write_text(json.dumps(rows, indent=2) + "\n")
    print(f"Analyzed {len(cases)} runs in {output}")


if __name__ == "__main__":
    main()
