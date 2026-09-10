#!/usr/bin/env python3
"""Generate N2 and O2 products from rigid monoenergetic or mixed beams."""

import argparse
import time

import numpy as np

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--energy-keV", type=float, default=50.0)
parser.add_argument("--mixed", action="store_true")
parser.add_argument("--steps", type=int, default=1)
args = parser.parse_args()

PARTICLES_PER_DIRECTION = 16
PARTICLES_PER_BEAM = PARTICLES_PER_DIRECTION**3
PROJECTILE_DENSITY = 1.0e8
PROJECTILE_ENERGY = args.energy_keV * 1e3
BACKGROUND_DENSITY = 2.0e20
TIME_STEP = 1.0e-9
FIXED_PRODUCT_WEIGHT = 15.0 if args.energy_keV < 10 else 200.0
MAX_PRODUCTS_PER_CELL = 30000

C = picmi.constants.c
M_E = picmi.constants.m_e
M_P = picmi.constants.m_p
Q_E = picmi.constants.q_e
M_U = 1.660_539_066_60e-27
PROJECTILE_REST_ENERGY = M_P * C**2 / Q_E
PROJECTILE_GAMMA = 1.0 + PROJECTILE_ENERGY / PROJECTILE_REST_ENERGY
PROJECTILE_PROPER_SPEED = C * np.sqrt(
    PROJECTILE_ENERGY
    / PROJECTILE_REST_ENERGY
    * (PROJECTILE_ENERGY / PROJECTILE_REST_ENERGY + 2)
)

grid = picmi.Cartesian3DGrid(
    number_of_cells=[1, 1, 1],
    lower_bound=[0.0, 0.0, 0.0],
    upper_bound=[1.0, 1.0, 1.0],
    lower_boundary_conditions=["periodic"] * 3,
    upper_boundary_conditions=["periodic"] * 3,
    lower_boundary_conditions_particles=["periodic"] * 3,
    upper_boundary_conditions_particles=["periodic"] * 3,
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)


def make_beam(name, direction):
    return picmi.Species(
        particle_type="proton",
        name=f"beam_{name}",
        initial_distribution=picmi.UniformDistribution(
            density=PROJECTILE_DENSITY,
            directed_velocity=np.asarray(direction) * PROJECTILE_PROPER_SPEED,
        ),
        warpx_do_not_push=True,
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


def make_electrons(name):
    return picmi.Species(
        particle_type="electron",
        name=f"electrons_{name}",
        warpx_do_not_push=True,
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


def make_ions(name, neutral_mass):
    return picmi.Species(
        name=f"ions_{name}",
        charge="q_e",
        mass=neutral_mass - M_E,
        warpx_do_not_push=True,
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


cases = {
    "N2": {
        "direction": [0.0, 0.0, 1.0],
        "neutral_mass": 28.0134 * M_U,
        "temperature": 300.0,
        "density": "2.0e20*(1.0+0.0*x)",
    },
    "O2": {
        "direction": [1.0, 0.0, 0.0],
        "neutral_mass": 31.9988 * M_U,
        "temperature": "600.0+0.0*t",
        "density": BACKGROUND_DENSITY,
    },
}

collisions = []
for name, case in cases.items():
    case["beam"] = make_beam(name, case["direction"])
    case["electrons"] = make_electrons(name)
    case["ions"] = make_ions(name, case["neutral_mass"])
    collisions.append(
        picmi.ProtonImpactIonizationCollisions(
            name=f"pjg_{name}",
            species=case["beam"],
            product_species=[case["electrons"], case["ions"]],
            ionization_target=name,
            background_density=case["density"],
            background_temperature=case["temperature"],
            fixed_product_weight=FIXED_PRODUCT_WEIGHT,
            max_products_per_cell=MAX_PRODUCTS_PER_CELL,
        )
    )

sim = picmi.Simulation(
    solver=solver,
    time_step_size=TIME_STEP,
    max_steps=args.steps,
    warpx_collisions=collisions,
    warpx_random_seed=42,
    warpx_serialize_initial_conditions=True,
    verbose=1,
)
for case in cases.values():
    sim.add_species(
        case["beam"],
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[PARTICLES_PER_DIRECTION] * 3, grid=grid
        ),
    )
    for product in [case["electrons"], case["ions"]]:
        sim.add_species(
            product,
            layout=picmi.GriddedLayout(n_macroparticle_per_cell=[0, 0, 0], grid=grid),
        )


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def component(container, name):
    return np.concatenate(
        [to_numpy(tile[name]) for tile in container.iterator(level=0)]
    ).astype(np.float64)


def ids(container):
    return np.concatenate(
        [
            to_numpy(libwarpx.amr.unpack_ids(tile["idcpu"]))
            for tile in container.iterator(level=0)
        ]
    )


def particle_data(container):
    return {
        name: component(container, name)
        for name in ["x", "y", "z", "ux", "uy", "uz", "w"]
    } | {"id": ids(container)}


sim.initialize_inputs()
initialization_start = time.perf_counter()
sim.initialize_warpx()
initialization_elapsed = time.perf_counter() - initialization_start

results = {
    "projectile_energy": PROJECTILE_ENERGY,
    "projectile_rest_energy": PROJECTILE_REST_ENERGY,
    "projectile_density": PROJECTILE_DENSITY,
    "background_density": BACKGROUND_DENSITY,
    "time_step": TIME_STEP,
    "steps": args.steps,
    "max_products_per_cell": MAX_PRODUCTS_PER_CELL,
    "fixed_product_weight": FIXED_PRODUCT_WEIGHT,
    "particles_per_beam": PARTICLES_PER_BEAM,
    "initialization_elapsed": initialization_elapsed,
}
for name in cases:
    beam = sim.particles.get(f"beam_{name}")
    if args.mixed:
        # Keep distinct parent energies in the same cell, in opposite orders
        # for the two gases. This detects parent/energy-quantile correlation.
        for tile in beam.iterator(level=0):
            subset = (
                slice(None, len(tile["w"]) // 2)
                if name == "N2"
                else slice(len(tile["w"]) // 2, None)
            )
            other_energy = 500e3
            other_speed = C * np.sqrt(
                other_energy
                / PROJECTILE_REST_ENERGY
                * (other_energy / PROJECTILE_REST_ENERGY + 2)
            )
            for component_name, direction in zip(
                ("ux", "uy", "uz"), cases[name]["direction"], strict=True
            ):
                tile[component_name][subset] = direction * other_speed
    for component_name, values in particle_data(beam).items():
        results[f"{name}_beam_initial_{component_name}"] = values.copy()

step_start = time.perf_counter()
sim.step(args.steps)
results["step_elapsed"] = time.perf_counter() - step_start

for name, case in cases.items():
    for species_name in ["beam", "electrons", "ions"]:
        container = sim.particles.get(f"{species_name}_{name}")
        for component_name, values in particle_data(container).items():
            results[f"{name}_{species_name}_{component_name}"] = values
    results[f"{name}_neutral_mass"] = case["neutral_mass"]
    results[f"{name}_temperature"] = float(300.0 if name == "N2" else 600.0)
    results[f"{name}_direction"] = np.asarray(case["direction"])

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez("proton_impact_ionization_results.npz", **results)
    total_products = sum(results[f"{name}_electrons_w"].size for name in cases)
    print(
        "Proton-impact ionization physics test: "
        f"products={total_products}, initialization={initialization_elapsed:.6f} s, "
        f"step={results['step_elapsed']:.6f} s"
    )

sim.finalize()
