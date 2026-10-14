#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Analytical finite-step cooling and attachment limits of the actual MCC loop."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import binom, poisson

from pywarpx import amrex, picmi

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mode", choices=["cooling", "attachment"], default="cooling")
parser.add_argument("--steps", type=int, default=32)
parser.add_argument("--third-body-factor", type=float, default=1)
args, _ = parser.parse_known_args()
count, cells = 32768, 64
nu, density, third_body = 1e11, 1e28, 3e24
c, me, qe = picmi.constants.c, picmi.constants.m_e, picmi.constants.q_e
rest = me * c * c / qe
initial_energy = 5.5 if args.mode == "cooling" else 0.1
energy = np.unique(
    np.r_[np.linspace(0.01, 8, 4000), 0.1, np.arange(0.5, 7), 1, 1.00001]
)
speed = c * np.sqrt(energy * (energy + 2 * rest)) / (energy + rest)
sigma = nu / (density * speed)
grid = picmi.Cartesian1DGrid(
    number_of_cells=[cells],
    lower_bound=[0],
    upper_bound=[1],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=cells,
    warpx_blocking_factor=1,
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(
        density=1,
        directed_velocity=[
            0,
            0,
            c * np.sqrt(initial_energy * (initial_energy + 2 * rest)) / rest,
        ],
    ),
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
    warpx_do_not_deposit=True,
)
products = []
if args.mode == "cooling":
    sigma[energy <= 1] = 0
    np.savetxt(
        "cooling.txt",
        np.column_stack((energy, sigma)),
        header="Synthetic constant-frequency fixed-loss ladder",
    )
    np.savetxt("isotropic.txt", [[0.01, 1, 1, 1], [8, 1, 1, 1]])
    processes = {
        "excitation": {
            "cross_section": str(Path("cooling.txt").resolve()),
            "energy": 1,
            "scattering_angle_model": "IAA",
            "differential_cross_section": str(Path("isotropic.txt").resolve()),
        }
    }
    total_nu = nu
else:
    processes = {}
    for name, fraction, units in [("binary", 0.7, "m2"), ("ternary", 0.3, "m5")]:
        product = picmi.Species(
            name=name,
            charge=-qe,
            mass=1.0,
            warpx_do_not_push=True,
            warpx_do_not_gather=True,
            warpx_do_not_deposit=True,
        )
        products.append(product)
        values = fraction * sigma / (third_body if units == "m5" else 1)
        np.savetxt(
            name + ".txt",
            np.column_stack((energy, values)),
            header="Synthetic attachment; raw m5 values below binary32 range",
        )
        processes["attachment_" + name] = {
            "cross_section": str(Path(name + ".txt").resolve()),
            "cross_section_units": units,
            "species": product,
        }
        if units == "m5":
            processes["attachment_" + name]["third_body_density"] = (
                third_body * args.third_body_factor
            )
    total_nu = nu * (0.7 + 0.3 * args.third_body_factor)
dt = 4 / (nu * args.steps)
collision = picmi.MCCCollisions(
    name="analytical",
    species=electrons,
    background_density=density,
    background_temperature=0,
    background_mass=1.0,
    scattering_processes=processes,
    nu_max=2 * total_nu,
)
simulation = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9),
    time_step_size=dt,
    warpx_collisions=[collision],
    warpx_random_seed=23451,
    verbose=0,
)
simulation.add_species(
    electrons,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[count // cells], grid=grid),
)
for product in products:
    simulation.add_species(product, layout=None)
simulation.initialize_inputs()
amrex.the_arena_init_size = 8 * 1024 * 1024
simulation.initialize_warpx()
simulation.step(args.steps)


def array(values):
    return np.asarray(values.get() if hasattr(values, "get") else values)


if args.mode == "cooling":
    energies = []
    for tile in simulation.particles.get("electrons").iterator(level=0):
        u2 = sum(array(tile[k]).astype(float) ** 2 for k in ["ux", "uy", "uz"])
        energies.extend(me * u2 / (qe * (1 + np.sqrt(1 + u2 / c**2))))
    coordinate = initial_energy - np.asarray(energies)
    levels = np.rint(coordinate).astype(int)
    assert np.max(np.abs(coordinate - levels)) < 3e-5
    assert np.all((levels >= 0) & (levels <= 5))
    measured = np.bincount(levels, minlength=6) / count
    probability = -math.expm1(-nu * dt)
    expected = np.r_[
        binom.pmf(np.arange(5), args.steps, probability),
        binom.sf(4, args.steps, probability),
    ]
    continuous = np.r_[poisson.pmf(np.arange(5), 4), poisson.sf(4, 4)]
    half_step = -math.expm1(-nu * dt / 2)
    refined = np.r_[
        binom.pmf(np.arange(5), 2 * args.steps, half_step),
        binom.sf(4, 2 * args.steps, half_step),
    ]
    assert np.linalg.norm(refined - continuous) < np.linalg.norm(expected - continuous)
else:
    names = ["electrons", "binary", "ternary"]
    measured = (
        np.array(
            [
                sum(
                    len(tile["w"])
                    for tile in simulation.particles.get(name).iterator(level=0)
                )
                for name in names
            ]
        )
        / count
    )
    survival = np.exp(-total_nu * args.steps * dt)
    expected = np.array(
        [
            survival,
            (1 - survival) * 0.7 * nu / total_nu,
            (1 - survival) * 0.3 * args.third_body_factor * nu / total_nu,
        ]
    )
    assert measured.sum() == 1
assert np.all(
    np.abs(measured - expected)
    < 7 * np.sqrt(expected * (1 - expected) / count) + 1 / count
)
Path("analytical_results.json").write_text(
    json.dumps(
        {
            "mode": args.mode,
            "steps": args.steps,
            "observed": measured.tolist(),
            "expected": expected.tolist(),
        },
        indent=2,
    )
    + "\n"
)
print("IAA analytical", args.mode, "PASS", measured, expected)
