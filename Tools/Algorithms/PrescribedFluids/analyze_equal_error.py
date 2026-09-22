#!/usr/bin/env python3
"""Compare beam cost at predefined continuum-field accuracy targets.

The targets follow the independent expanded-domain convergence study, before
running this ensemble: 1% core RMS error for Er/Btheta and 2% for Ez. They do
not change any regression assertion. Both initialization and the final output
are compared with the translating free-space Gaussian solution. This is a
vacuum beam-field benchmark, not an equal-error claim for nonlinear chemistry.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from analyze import interval, weighted_norm
from gaussian_reference import gaussian_fields
from scipy.constants import c, e, proton_mass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    targets = {"Er": 0.01, "Ez": 0.02, "Btheta": 0.01}
    sigma_r, sigma_t = 0.002, 25e-12
    gamma = 1 + 800e6 * e / (proton_mass * c**2)
    velocity = c * np.sqrt(1 - gamma**-2)
    sigma_z = velocity * sigma_t
    charge = 0.6 * np.sqrt(2 * np.pi) * sigma_t
    exact_fields, groups = {}, defaultdict(list)
    configurations = set()
    parameters = None
    for path in sorted(args.directory.glob("*/result.json")):
        result = json.loads(path.read_text())
        configurations.add(
            (
                result["revision"],
                result["backend"],
                result["ranks"],
                result["benchmark_sha256"],
            )
        )
        p = result["parameters"]
        assert p["mode"] == "fields" and p["solver"] == "Yee"
        assert p["cells"] == [256, 1024]
        assert p["radial_sigmas"] == 32 and p["longitudinal_sigmas"] == 48
        assert not p["no_self_fields"] and not p["mcc"]
        settings = {
            key: value for key, value in p.items() if key not in ["beam", "ppc", "seed"]
        }
        if parameters is None:
            parameters = settings
        assert settings == parameters, path
        errors = defaultdict(list)
        nr, nz = p["cells"]
        with np.load(path.with_suffix(".npz")) as arrays:
            for state in [result["history"][0], result["history"][-1]]:
                time, step = state["time"], state["step"]
                for index, component in enumerate(targets):
                    actual = arrays[f"{step}_{component}"]
                    assert np.isfinite(actual).all(), (path, component)
                    key = (time, component, actual.shape)
                    if key not in exact_fields:
                        r = np.arange(actual.shape[0]) + (actual.shape[0] == nr) / 2
                        r *= 32 * sigma_r / nr
                        z = np.arange(actual.shape[1]) + (actual.shape[1] == nz) / 2
                        z = z * 96 * sigma_z / nz - 48 * sigma_z - velocity * time
                        ir, iz = r < 3 * sigma_r, np.abs(z) < 3 * sigma_z
                        exact = np.zeros_like(actual)
                        exact[np.ix_(ir, iz)] = gaussian_fields(
                            r[ir, None], z[None, iz], sigma_r, sigma_z, charge, velocity
                        )[index]
                        exact_fields[key] = (exact, ir[:, None] & iz[None, :])
                    exact, mask = exact_fields[key]
                    errors[component].append(
                        weighted_norm((actual - exact) * mask, p["cells"])
                        / weighted_norm(exact, p["cells"])
                    )
        groups[p["beam"], p["ppc"]].append(
            dict(
                seed=p["seed"],
                timestep_s=result["ordinary_step_median_s"],
                initialization_s=result["initialization_s"],
                worst_output_errors={
                    key: max(values) for key, values in errors.items()
                },
            )
        )
    assert len(configurations) == 1, (
        "Cannot mix builds, benchmark versions or MPI sizes"
    )
    assert set(groups) == {("fluid", 1)} | {
        (beam, ppc) for beam in ["quiet", "random"] for ppc in [4, 16, 64, 256]
    }, "The prescribed ensemble is incomplete"
    rows = []
    for (beam, ppc), samples in sorted(groups.items()):
        assert len(samples) >= 6 and len({sample["seed"] for sample in samples}) == len(
            samples
        )
        errors = {
            key: interval([sample["worst_output_errors"][key] for sample in samples])
            for key in targets
        }
        rows.append(
            dict(
                beam=beam,
                ppc=ppc,
                errors=errors,
                qualified=all(
                    errors[key]["mean"] + errors[key]["ci99_half_width"] <= limit
                    for key, limit in targets.items()
                ),
                timestep_s=interval([sample["timestep_s"] for sample in samples]),
                initialization_s=interval(
                    [sample["initialization_s"] for sample in samples]
                ),
                samples=samples,
            )
        )
    fluid = next(row for row in rows if row["beam"] == "fluid")
    comparisons = {}
    for beam in ["quiet", "random"]:
        eligible = [row for row in rows if row["beam"] == beam and row["qualified"]]
        best = (
            min(eligible, key=lambda row: row["timestep_s"]["mean"])
            if eligible
            else None
        )
        comparisons[beam] = (
            None
            if best is None
            else dict(
                ppc=best["ppc"],
                particle_over_fluid_timestep_ratio=best["timestep_s"]["mean"]
                / fluid["timestep_s"]["mean"],
            )
        )
    report = dict(
        targets=targets, parameters=parameters, groups=rows, comparisons=comparisons
    )
    (args.directory / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    assert fluid["qualified"], (
        "Fluid reference did not meet the predefined continuum targets"
    )
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
