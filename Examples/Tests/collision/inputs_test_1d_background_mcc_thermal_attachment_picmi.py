#!/usr/bin/env python3
"""Exercise near-rest electrons colliding with thermal neutral molecules."""

import math
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

PARTICLE_COUNT = 131072
INITIAL_PRODUCT_COUNT = 1
BACKGROUND_DENSITY = 1.0e20
BACKGROUND_TEMPERATURE = 300.0
BACKGROUND_MASS = 32.0 * picmi.constants.m_p
CROSS_SECTION = 1.0e-22
MAJORANT_THERMAL_STD = 8.0
MAJORANT_OPTICAL_DEPTH = 0.5
COLLISION_SUBCYCLES = 16

NEUTRAL_VELOCITY_STD = math.sqrt(
    picmi.constants.kb * BACKGROUND_TEMPERATURE / BACKGROUND_MASS
)
NU_MAX = (
    BACKGROUND_DENSITY * CROSS_SECTION * MAJORANT_THERMAL_STD * NEUTRAL_VELOCITY_STD
)
DT = MAJORANT_OPTICAL_DEPTH / NU_MAX

source_dir = Path(__file__).resolve().parent

grid = picmi.Cartesian1DGrid(
    number_of_cells=[1],
    lower_bound=[0.0],
    upper_bound=[100000.0],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)

electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, 0.0],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)


def make_negative_ions(name):
    return picmi.Species(
        name=name,
        charge="-q_e",
        mass=BACKGROUND_MASS + picmi.constants.m_e,
        initial_distribution=picmi.UniformDistribution(density=1.0),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


negative_a = make_negative_ions("negative_a")
negative_b = make_negative_ions("negative_b")


def make_collision(name, product_species):
    return picmi.MCCCollisions(
        name=name,
        species=electrons,
        background_density=BACKGROUND_DENSITY,
        background_temperature=BACKGROUND_TEMPERATURE,
        background_mass=BACKGROUND_MASS,
        scattering_processes={
            "attachment": {
                "cross_section": str(source_dir / "background_mcc_attachment_m2.txt"),
                "cross_section_units": "m2",
                "species": product_species,
            }
        },
        nu_max=NU_MAX,
        ndt_subcycle=COLLISION_SUBCYCLES,
    )


collision_a = make_collision("mcc_a", negative_a)
collision_b = make_collision("mcc_b", negative_b)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=DT,
    max_steps=1,
    warpx_collisions=[collision_a, collision_b],
    verbose=1,
)
sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid),
)
for species in [negative_a, negative_b]:
    sim.add_species(
        species,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[INITIAL_PRODUCT_COUNT], grid=grid
        ),
    )

sim.initialize_inputs()
sim.initialize_warpx()
sim.step(1)


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def product_data(name):
    container = sim.particles.get(name)
    count = container.number_of_particles(only_local=True)
    proper_velocity = {
        component: np.concatenate(
            [to_numpy(pti[component]) for pti in container.iterator(level=0)]
        )[INITIAL_PRODUCT_COUNT:]
        for component in ["ux", "uy", "uz"]
    }
    return count - INITIAL_PRODUCT_COUNT, proper_velocity


negative_a_events, negative_a_u = product_data("negative_a")
negative_b_events, negative_b_u = product_data("negative_b")
electron_count = sim.particles.get("electrons").number_of_particles(only_local=True)

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "background_mcc_thermal_attachment_results.npz",
        electron_count=electron_count,
        negative_a_events=negative_a_events,
        negative_b_events=negative_b_events,
        negative_a_ux=negative_a_u["ux"],
        negative_a_uy=negative_a_u["uy"],
        negative_a_uz=negative_a_u["uz"],
        negative_b_ux=negative_b_u["ux"],
        negative_b_uy=negative_b_u["uy"],
        negative_b_uz=negative_b_u["uz"],
        neutral_velocity_std=NEUTRAL_VELOCITY_STD,
    )

sim.finalize()
