#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""End-to-end combined-family rates, signed losses and MCC timing."""

import argparse
import shlex
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "Tools/CrossSections"))
from rotation_reference import Bundle

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--data-dir", type=Path, required=True)
parser.add_argument("--particles", type=int, default=65536)
parser.add_argument("--cells", type=int, default=1)
parser.add_argument("--steps", type=int, default=1)
parser.add_argument("--cumulative", action="store_true")
parser.add_argument("--anisotropic", action="store_true")
parser.add_argument(
    "--elmolcs",
    action="store_true",
    help="read prepared N2/IAA and O2/IAA source families",
)
parser.add_argument("--seed", type=int, default=42)
parser.add_argument(
    "--spectator", action="store_true", help="kinetic IAA spectator closure"
)
parser.add_argument("--spectator-reference", type=Path)
parser.add_argument(
    "--broad-spectrum",
    action="store_true",
    help="benchmark spectator scattering with log-uniform 0.0021--900 eV electrons",
)
args = parser.parse_args()
if args.cells <= 0 or args.particles % args.cells:
    raise ValueError("The particle count must be divisible by the positive cell count")
if args.broad_spectrum and (not args.spectator or args.steps <= 1):
    raise ValueError(
        "The broad-spectrum benchmark requires --spectator and --steps > 1"
    )
spectator_reference = {}
if args.spectator:
    if args.spectator_reference is None:
        raise ValueError(
            "Spectator tests require the separately prepared quadrature reference"
        )
    for line in args.spectator_reference.read_text().splitlines():
        fields = shlex.split(line)
        temperature, energy = map(float, fields[1:3])
        spectator_reference[temperature, energy] = np.array(fields[4:], float)

grid = picmi.Cartesian1DGrid(
    number_of_cells=[args.cells],
    lower_bound=[0],
    upper_bound=[100],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=args.cells,
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
for target in ["N2"] if args.spectator else ["N2", "O2"]:
    directory = (
        args.data_dir / target / "IAA"
        if args.elmolcs or args.spectator
        else args.data_dir
    )
    file = (
        directory
        / (
            "thermal_spectator.rot"
            if args.spectator
            else "thermal_rotation.rot"
            if args.elmolcs
            else f"{target}_analytic.rot"
        )
    ).resolve()
    bundle = None if args.spectator else Bundle.read(file)
    inclusive = (
        directory
        / (
            "elastic.txt"
            if args.spectator
            else "thermal_rotation_elastic.txt"
            if args.elmolcs or args.spectator
            else f"{target}_inclusive.txt"
        )
    ).resolve()
    angular_file = "elastic_dcs_anisotropic" if args.anisotropic else "elastic_dcs"
    elastic_dcs = (
        directory
        / (
            "elastic_dcs.txt"
            if args.elmolcs or args.spectator
            else f"{target}_{angular_file}.txt"
        )
    ).resolve()
    ordinary = (args.data_dir / f"{target}_ordinary.txt").resolve()
    for temperature in [0, 100, 300, 1000]:
        rates = (
            None if args.spectator else bundle.at_temperature(temperature).sum(axis=1)
        )
        energies = (
            [0.0021, 0.01, 0.1, 0.5, 2.3, 2.6, 10, 100, 900]
            if args.spectator
            else [0, 0.0021, 0.01, 0.5, 2.3, 10]
            if args.elmolcs
            else [0, 0.5]
        )
        if args.spectator and args.steps == 1:
            energies.append(1000)
        if args.broad_spectrum:
            # One ensemble per temperature; overwrite its nominal momentum
            # below with a spectrum spanning the loaded source table.
            energies = [0.1]
        for energy in energies:
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
                "rotation_model": "iaa_spectator"
                if args.spectator
                else "elastic_dcs"
                if args.elmolcs
                else "analytic_test",
                "scattering_angle_model": "IAA",
                "differential_cross_section": str(elastic_dcs),
                "rotational_temperature": temperature,
                "rotation_sampling": "cumulative" if args.cumulative else "alias",
            }
            collisions.append(
                picmi.MCCCollisions(
                    name=f"mcc_{name}",
                    species=electron,
                    background_density=density,
                    background_temperature=0,
                    scattering_processes={"elastic": rotation}
                    if args.elmolcs or args.spectator
                    else {
                        "elastic": rotation,
                        "elastic_ordinary": {"cross_section": str(ordinary)},
                    },
                )
            )
            species.append(electron)
            v = c * np.sqrt(energy * (energy + 2 * rest)) / (energy + rest)
            if args.spectator:
                table = np.loadtxt(inclusive)
                total = v * np.interp(energy, table[:, 0], table[:, 1])
                probability = -np.expm1(-density * total * dt)
                expected = spectator_reference[temperature, energy]
                mean, second, fourth = probability * expected[[0, 1, 9]]
            else:
                interpolated = np.array(
                    [np.interp(energy, bundle.energies, column) for column in rates.T]
                )
                total = interpolated.sum() + (0 if args.elmolcs else v * 1e-21)
                probability = -np.expm1(-density * total * dt)
                mean = probability * np.dot(interpolated, bundle.losses) / total
                second = probability * np.dot(interpolated, bundle.losses**2) / total
                fourth = probability * np.dot(interpolated, bundle.losses**4) / total
            references.append((name, target, energy, mean, second, fourth, probability))

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
            n_macroparticle_per_cell=[args.particles // args.cells], grid=grid
        ),
    )
sim.initialize_inputs()
start = time.perf_counter()
sim.initialize_warpx()
startup = time.perf_counter() - start
initial_momentum = {}
initial_mean_energy = {}
spectrum_rng = np.random.default_rng(args.seed)
for electron in species:
    energy_sum, particle_count = 0.0, 0
    # Exhaust each AMReX iterator before opening the next species' iterator.
    for tile in sim.particles.get(electron.name).iterator(level=0):
        if args.broad_spectrum:
            kinetic = np.exp(
                spectrum_rng.uniform(np.log(0.0021), np.log(900), len(tile["uz"]))
            )
            stored = (c * np.sqrt(kinetic * (kinetic + 2 * rest)) / rest).astype(
                tile["uz"].dtype
            )
            if hasattr(tile["uz"], "set"):
                tile["uz"].set(stored)
                tile["uz"].device.synchronize()
            else:
                tile["uz"][:] = stored
            u2 = stored.astype(float) ** 2
            energy_sum += (me * u2 / (qe * (1 + np.sqrt(1 + u2 / c**2)))).sum()
            particle_count += len(stored)
        else:
            value = tile["uz"][0]
            initial_momentum[electron.name] = float(
                value.get() if hasattr(value, "get") else value
            )
    if args.broad_spectrum:
        initial_mean_energy[electron.name] = energy_sum / particle_count
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
    "cells": args.cells,
    "broad_spectrum": args.broad_spectrum,
}
for name, target, energy, mean, second, fourth, probability in references:
    container = sim.particles.get(name)
    ux, uy, uz = (component(container, key) for key in ["ux", "uy", "uz"])
    assert len(ux) == args.particles, (name, len(ux), args.particles)
    u2 = ux**2 + uy**2 + uz**2
    final_energy = me * u2 / (qe * (1 + np.sqrt(1 + u2 / c**2)))
    if args.broad_spectrum:
        # A difference of ensemble means remains valid if particles reorder.
        # There is no fixed-energy moment reference for this timing-only input.
        output[name] = [initial_mean_energy[name] - final_energy.mean()]
        assert np.isfinite(output[name]).all(), name
        continue
    # Recover internal loss using the independently reconstructed neutral recoil.
    # Reconstruct small rotational changes from the momentum actually stored.
    # At keV energies, float32 initialization roundoff can bias a meV signal.
    incoming = initial_momentum[name]
    initial_energy = me * incoming**2 / (qe * (1 + np.sqrt(1 + incoming**2 / c**2)))
    mass = (28.0134 if target == "N2" else 31.9988) * 1.66053906660e-27
    recoil_u2 = (me / mass) ** 2 * (ux**2 + uy**2 + (uz - incoming) ** 2)
    recoil = mass * recoil_u2 / (qe * (1 + np.sqrt(1 + recoil_u2 / c**2)))
    # One event permits a recoil reconstruction. In multi-step benchmarks use
    # the physical electron energy change, without combining distinct neutrals.
    losses = initial_energy - final_energy - (recoil if args.steps == 1 else 0)
    changed = (ux != 0) | (uy != 0) | (np.abs(uz - incoming) > 1e-6 * max(incoming, 1))
    if args.anisotropic:
        assert args.steps == 1
        rotated = np.abs(losses) > 1e-4
        mu = uz[rotated] / np.sqrt(u2[rotated])
        assert len(mu) > 0
        expected_mu = 0 if energy == 0 else -2 / 15
        expected_mu2 = 1 / 3 if energy == 0 else (2 / 3 + 44 * 3 / 105) / 6
        for values, expected in [(mu, expected_mu), (mu**2, expected_mu2)]:
            assert (
                abs(values.mean() - expected)
                < 7 * np.sqrt(values.var() / len(values)) + 2e-4
            )
        correlation = (mu - expected_mu) * losses[rotated]
        assert abs(correlation.mean()) < 7 * np.sqrt(correlation.var() / len(mu)) + 2e-7
    output[name] = [
        losses.mean(),
        (losses**2).mean(),
        mean,
        second,
        changed.mean(),
        probability,
    ]
    if (args.elmolcs or args.spectator) and args.steps == 1:
        # Test with independently reconstructed neutral recoil. The source
        # bundle is loaded only; neither a source fit nor a test table is built.
        first_error = 7 * np.sqrt(
            max(losses.var(), second - mean**2, 0) / args.particles
        )
        second_error = 7 * np.sqrt(
            max((losses**2).var(), fourth - second**2, 0) / args.particles
        )
        assert abs(losses.mean() - mean) <= first_error + 0.002 * max(
            abs(mean), 1e-10
        ), (name, losses.mean(), mean)
        assert abs((losses**2).mean() - second) <= second_error + 0.002 * max(
            second, 1e-10
        ), (name, (losses**2).mean(), second)
        assert (
            abs(changed.mean() - probability)
            <= 7 * np.sqrt(probability * (1 - probability) / args.particles) + 1e-5
        ), (name, changed.mean(), probability)
        if args.spectator:
            temperature = int(name.split("_")[2])
            expected = spectator_reference[temperature, energy]
            cosine = np.divide(uz, np.sqrt(u2), out=np.ones_like(uz), where=u2 > 0)
            mixed = losses * cosine
            predicted = probability * expected[3]
            error = 7 * np.sqrt(
                max(mixed.var(), probability * expected[11] - predicted**2, 0)
                / args.particles
            )
            assert abs(mixed.mean() - predicted) <= error + 0.002 * max(
                abs(predicted), 1e-10
            ), (name, mixed.mean(), predicted)
if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez("background_mcc_rotation_results.npz", **output)
sim.finalize()
