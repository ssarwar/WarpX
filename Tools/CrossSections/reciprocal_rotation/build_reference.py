# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Build and independently refine the low-energy source reference."""

import argparse
import json

import numpy as np
from hybrid_reference import REFERENCE_VERSION, REST, C, Hybrid
from reference_paths import output_path
from reference_refinement import inserted_moments


def build_reference(
    target, temperature=300.0, maximum_energy=1000.0, tolerance=2e-4, seed=None
):
    """Return the positive reference and source-grid convergence diagnostics."""
    model = Hybrid(
        target,
        temperature,
        maximum_energy,
        resolution=2,
        rank_max=48,
        angular_resolution=4,
        energy_override=None if seed is None else np.load(seed)["energy"],
    ).solve()
    history = []
    for iteration in range(12):
        ids = np.flatnonzero(
            (model.energy[:-1] >= 0.001) & (model.energy[1:] <= maximum_energy)
        )
        probes = (
            model.energy[ids, None]
            + (model.energy[ids + 1] - model.energy[ids])[:, None]
            * np.array([0.25, 0.5, 0.75])
        ).ravel()
        added = []
        worst = 0.0
        worst_energy = 0.0
        for energy in probes:
            fresh, _ = inserted_moments(model, energy)
            interpolated = model.angular_moments(energy)
            measure = model.measure[:, None] * np.column_stack(
                (np.ones(len(model.y)), 2 * model.y**2, 4 * model.y**4)
            )
            exact = fresh @ measure
            current = interpolated @ measure
            error = float(np.max(abs(exact - current) / np.maximum(abs(exact), 1e-40)))
            if error > worst:
                worst, worst_energy = error, float(energy)
            if error > tolerance:
                added.append(energy)
        history.append(
            dict(
                iteration=iteration,
                rows=len(model.energy),
                added=len(added),
                maximum_error=worst,
                worst_energy_eV=worst_energy,
            )
        )
        print(history[-1], flush=True)
        if not added:
            return model, history
        model.energy = np.unique(np.r_[model.energy, added])
        model.p = np.sqrt(model.energy * (model.energy + 2 * REST))
        model.kin = C / (model.energy + REST)
        model.boundary = np.searchsorted(model.energy, maximum_energy + 1)
        model.X = np.zeros((len(model.energy), model.nb, len(model.y)))
        model.c = np.ones((len(model.energy), len(model.y)))
        model.cold_anchor = None
        model.raw_residual_min = 1.0
        model.fail = []
        model.solve()
    raise RuntimeError("Source-grid refinement did not converge")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", choices=["N2", "O2"], required=True)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--maximum-energy", type=float, default=1000.0)
    parser.add_argument(
        "--seed", type=str, help="Refined energy mesh; solve all angular rows anew"
    )
    args = parser.parse_args()
    model, history = build_reference(
        args.target, args.temperature, args.maximum_energy, seed=args.seed
    )
    np.savez_compressed(
        output_path(f"refined-reference-{args.target}.npz"),
        energy=model.energy,
        y=model.y,
        measure=model.measure,
        X=model.X,
        c=model.c,
        temperature=model.T,
        target=model.target,
        maximum_energy=args.maximum_energy,
        reference_version=REFERENCE_VERSION,
        initial_jmax=model.jmax,
        ranks=model.ranks,
    )
    output_path(f"reference-refinement-{args.target}.json").write_text(
        json.dumps(history, indent=2) + "\n"
    )
