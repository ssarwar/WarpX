#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""End-to-end combined-family rates, signed losses and MCC timing."""

import argparse
import time
from pathlib import Path

import numpy as np
from analysis_rotation_reference import analytic_bundle

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--particles", type=int, default=65536)
parser.add_argument("--steps", type=int, default=1)
parser.add_argument("--cumulative", action="store_true")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

grid = picmi.Cartesian1DGrid(
    number_of_cells=[1],
    lower_bound=[0],
    upper_bound=[100],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)
density = 1e29
dt = 1e-16
c = picmi.constants.c
me = picmi.constants.m_e
qe = picmi.constants.q_e
rest = me * c * c / qe
species, collisions, references = [], [], []
for target in ["N2", "O2"]:
    bundle = analytic_bundle(target)
    file = Path(f"{target}_analytic.rot").resolve()
    bundle.write(file)
    inclusive = Path(f"{target}_inclusive.txt").resolve()
    np.savetxt(inclusive, [[0, 2e-20], [10, 2e-20]])
    # An ordinary channel exercises the retained selector alongside the family.
    ordinary = Path(f"{target}_ordinary.txt").resolve()
    np.savetxt(ordinary, [[0, 1e-21], [10, 1e-21]])
    for temperature in [0, 100, 300, 1000]:
        rates = bundle.at_temperature(temperature).sum(axis=1)
        for energy in [0, 0.5]:
            if energy == 0 and temperature == 0:
                continue
            name = f"e_{target}_{temperature}_{energy}".replace(".", "_")
            u = c * np.sqrt(energy * (energy + 2 * rest)) / rest
            electron = picmi.Species(
                particle_type="electron",
                name=name,
                initial_distribution=picmi.UniformDistribution(
                    density=1, directed_velocity=[0, 0, u]
                ),
                warpx_do_not_deposit=True,
                warpx_do_not_gather=True,
            )
            rotation = {
                "cross_section": str(inclusive),
                "rotation_file": str(file),
                "rotation_model": "analytic_test",
                "rotational_temperature": temperature,
                "rotation_sampling": "cumulative" if args.cumulative else "alias",
            }
            collisions.append(
                picmi.MCCCollisions(
                    name=f"mcc_{name}",
                    species=electron,
                    background_density=density,
                    background_temperature=0,
                    scattering_processes={
                        "elastic": rotation,
                        "elastic_ordinary": {"cross_section": str(ordinary)},
                    },
                )
            )
            species.append(electron)
            interpolated = np.array(
                [np.interp(energy, bundle.energies, column) for column in rates.T]
            )
            v = c * np.sqrt(energy * (energy + 2 * rest)) / (energy + rest)
            total = interpolated.sum() + v * 1e-21
            probability = -np.expm1(-density * total * dt)
            mean = probability * np.dot(interpolated, bundle.losses) / total
            second = probability * np.dot(interpolated, bundle.losses**2) / total
            references.append((name, target, energy, mean, second, probability))

sim = picmi.Simulation(
    solver=solver,
    time_step_size=dt,
    max_steps=args.steps,
    warpx_collisions=collisions,
    warpx_random_seed=args.seed,
    verbose=0,
)
for electron in species:
    sim.add_species(
        electron,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[args.particles], grid=grid
        ),
    )
sim.initialize_inputs()
start = time.perf_counter()
sim.initialize_warpx()
startup = time.perf_counter() - start
start = time.perf_counter()
sim.step(args.steps)
# A scalar reduction fences GPU work without copying the particle arrays.
sim.particles.get(species[0].name).sum_particle_weight(local=True)
elapsed = time.perf_counter() - start


def component(container, key):
    arrays = [tile[key] for tile in container.iterator(level=0)]
    return np.concatenate(
        [a.get() if hasattr(a, "get") else np.asarray(a) for a in arrays]
    ).astype(float)


output = {
    "startup_seconds": startup,
    "step_seconds": elapsed,
    "particles": args.particles,
    "steps": args.steps,
    "cases": len(references),
}
for name, target, energy, mean, second, probability in references:
    container = sim.particles.get(name)
    ux, uy, uz = (component(container, key) for key in ["ux", "uy", "uz"])
    u2 = ux**2 + uy**2 + uz**2
    final_energy = me * u2 / (qe * (1 + np.sqrt(1 + u2 / c**2)))
    # Recover internal loss using the independently reconstructed neutral recoil.
    incoming = c * np.sqrt(energy * (energy + 2 * rest)) / rest
    mass = (28.0134 if target == "N2" else 31.9988) * 1.66053906660e-27
    recoil_u2 = (me / mass) ** 2 * (ux**2 + uy**2 + (uz - incoming) ** 2)
    recoil = mass * recoil_u2 / (qe * (1 + np.sqrt(1 + recoil_u2 / c**2)))
    # One event permits a recoil reconstruction. In multi-step benchmarks use
    # the physical electron energy change, without combining distinct neutrals.
    losses = energy - final_energy - (recoil if args.steps == 1 else 0)
    changed = (ux != 0) | (uy != 0) | (np.abs(uz - incoming) > 1e-6 * max(incoming, 1))
    output[name] = [
        losses.mean(),
        (losses**2).mean(),
        mean,
        second,
        changed.mean(),
        probability,
    ]
if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez("background_mcc_rotation_results.npz", **output)
sim.finalize()
