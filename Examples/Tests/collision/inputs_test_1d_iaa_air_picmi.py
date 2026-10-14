#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Production IAA cooling, driven swarm and reduced proton-pulse regressions.

The pulse case isolates particle creation and electron-air collisions in a
homogeneous cell ensemble. It is not a model of the experiment's beam geometry.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from pywarpx import amrex, picmi
from pywarpx.iaa_air import load_iaa_air_processes

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data-dir", type=Path, required=True)
parser.add_argument("--target", choices=["N2", "O2", "air"], default="air")
parser.add_argument("--mode", choices=["cooling", "swarm", "pulse"], default="cooling")
parser.add_argument("--particles", type=int, default=4096)
parser.add_argument("--cells", type=int, default=16)
parser.add_argument("--steps", type=int, default=32)
parser.add_argument("--warmup", type=int, default=0)
parser.add_argument("--dt", type=float, default=1e-13)
parser.add_argument("--density", type=float, default=2.5e25)
parser.add_argument("--seed", type=int, default=2026)
parser.add_argument("--output", type=Path, default=Path("iaa_air_results.json"))
args, _ = parser.parse_known_args()
assert args.particles >= args.cells and args.particles % args.cells == 0

c, me, qe = picmi.constants.c, picmi.constants.m_e, picmi.constants.q_e
amu = 1.66053906660e-27
rest = me * c * c / qe
proper = c * np.sqrt(20 * (20 + 2 * rest)) / rest
grid = picmi.Cartesian1DGrid(
    number_of_cells=[args.cells],
    lower_bound=[0],
    upper_bound=[1],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=args.cells,
    warpx_blocking_factor=1,
)
distribution = picmi.UniformDistribution(
    density=1e10,
    directed_velocity=[0, 0, 0 if args.mode == "swarm" else proper],
    rms_velocity=[proper / np.sqrt(3)] * 3 if args.mode == "swarm" else [0] * 3,
)
electron = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=distribution,
    warpx_do_not_deposit=True,
    warpx_do_not_gather=args.mode != "swarm",
)
species = []
positive, negative = [], []


def product(name, charge, mass):
    result = picmi.Species(
        name=name,
        charge=charge,
        mass=mass,
        warpx_do_not_push=True,
        warpx_do_not_gather=True,
        warpx_do_not_deposit=True,
    )
    species.append(result)
    (positive if charge > 0 else negative).append(name)
    return result


fractions = {"N2": 0.79, "O2": 0.21} if args.target == "air" else {args.target: 1.0}
collisions = []
beam = None
if args.mode == "pulse":
    proton_rest = picmi.constants.m_p * c * c / qe
    beam = picmi.Species(
        particle_type="proton",
        name="protons",
        initial_distribution=picmi.UniformDistribution(
            density=1e12,
            directed_velocity=[
                0,
                0,
                c * np.sqrt(8e8 * (8e8 + 2 * proton_rest)) / proton_rest,
            ],
        ),
        warpx_do_not_push=True,
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )
for gas, fraction in fractions.items():
    mass = {"N2": 28.0134, "O2": 31.9988}[gas] * amu
    ion = product(gas + "plus", qe, mass - me)
    attachments = None
    if gas == "O2":
        attachments = {
            "dissociative": product("Ominus", -qe, mass / 2 + me),
            "three_body": product("O2minus", -qe, mass + me),
        }
    processes = load_iaa_air_processes(
        args.data_dir,
        gas,
        ion_species=ion,
        attachment_species=attachments,
        third_body_density=args.density * fraction,
    )
    if beam:
        collisions.append(
            picmi.ProtonImpactIonizationCollisions(
                name="proton_" + gas,
                species=beam,
                product_species=[electron, ion],
                ionization_target=gas,
                background_density=args.density * fraction,
                background_temperature=300,
                fixed_product_weight=1e10 / args.particles,
                max_products_per_cell=128,
            )
        )
    collisions.append(
        picmi.MCCCollisions(
            name="air_" + gas,
            species=electron,
            background_density=args.density * fraction,
            background_temperature=300,
            background_mass=mass,
            scattering_processes=processes,
            ndt_subcycle=4,
        )
    )
simulation = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9),
    time_step_size=args.dt,
    warpx_collisions=collisions,
    warpx_random_seed=args.seed,
    warpx_serialize_initial_conditions=True,
    verbose=0,
)
simulation.add_species(
    electron,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[args.particles // args.cells], grid=grid
    ),
)
for item in species:
    simulation.add_species(item, layout=None)
if beam:
    simulation.add_species(
        beam, layout=picmi.GriddedLayout(n_macroparticle_per_cell=[4], grid=grid)
    )
if args.mode == "swarm":
    simulation.add_applied_field(
        picmi.AnalyticAppliedField(
            Ex_expression="0",
            Ey_expression="0",
            Ez_expression=str(50e-21 * args.density),
        )
    )
simulation.initialize_inputs()
amrex.the_arena_init_size = 8 * 1024 * 1024
begin = time.perf_counter()
simulation.initialize_warpx()
initialization = time.perf_counter() - begin


def array(value):
    return np.asarray(value.get() if hasattr(value, "get") else value)


def weight(name):
    return sum(
        float(array(tile["w"]).sum(dtype=np.float64))
        for tile in simulation.particles.get(name).iterator(level=0)
    )


def statistics():
    sums = np.zeros(4)
    count = 0
    bits = 64
    for tile in simulation.particles.get("electrons").iterator(level=0):
        bits = 8 * tile["ux"].dtype.itemsize
        u = [array(tile[k]).astype(float) for k in ["ux", "uy", "uz"]]
        w = array(tile["w"]).astype(float)
        gamma = np.sqrt(1 + sum(v * v for v in u) / c**2)
        energy = me * sum(v * v for v in u) / (qe * (1 + gamma))
        assert np.all(np.isfinite(energy)) and np.all(energy >= 0)
        sums += [w.sum(), w @ energy, w @ (energy * energy), w @ (u[2] / gamma)]
        count += len(w)
    return {
        "electrons": sums[0],
        "mean_energy_eV": sums[1] / sums[0],
        "second_energy_moment_eV2": sums[2] / sums[0],
        "drift_velocity_m_s": sums[3] / sums[0],
        "particles": count,
        "particle_precision_bits": bits,
        "positive_ions": sum(weight(name) for name in positive),
        "negative_ions": sum(weight(name) for name in negative),
    }


initial = statistics()
base_weights = []
if beam:
    base_weights = [
        array(tile["w"]).copy()
        for tile in simulation.particles.get("protons").iterator(level=0)
    ]
history = []
elapsed = 0.0
for step in range(args.warmup + args.steps):
    if beam:
        phase = step % args.steps
        on = phase < max(
            1, args.steps // 8
        ) or args.steps // 2 <= phase < args.steps // 2 + max(1, args.steps // 8)
        for tile, original in zip(
            simulation.particles.get("protons").iterator(level=0), base_weights
        ):
            values = original * int(on)
            if hasattr(tile["w"], "set"):
                tile["w"].set(values)
            else:
                tile["w"][:] = values
    begin = time.perf_counter()
    simulation.step(1)
    duration = time.perf_counter() - begin
    if step >= args.warmup:
        elapsed += duration
    if beam and (step + 1) % max(1, args.steps // 4) == 0:
        history.append(statistics())
final = statistics()
charge_error = (
    final["positive_ions"]
    - final["negative_ions"]
    - final["electrons"]
    + initial["electrons"]
)
assert abs(charge_error) < 2e-6 * max(initial["electrons"], final["positive_ions"])
if args.mode == "cooling":
    assert final["mean_energy_eV"] < initial["mean_energy_eV"]
elif args.mode == "swarm":
    assert final["drift_velocity_m_s"] < 0
else:
    assert history[0]["positive_ions"] > 0
    assert history[2]["positive_ions"] > history[1]["positive_ions"]
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(
    json.dumps(
        {
            "target": args.target,
            "mode": args.mode,
            "seed": args.seed,
            "steps": args.steps,
            "dt_s": args.dt,
            "density_m3": args.density,
            "initialization_seconds": initialization,
            "step_seconds": elapsed,
            "initial": initial,
            "final": final,
            "history": history,
            "charge_error_electrons": charge_error,
        },
        indent=2,
    )
    + "\n"
)
print(f"IAA_AIR {args.mode} {args.target}: {final['mean_energy_eV']:.8g} eV; PASS")
