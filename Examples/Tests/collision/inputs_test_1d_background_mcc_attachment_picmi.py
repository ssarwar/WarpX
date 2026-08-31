#!/usr/bin/env python3
"""Exercise multiple MCC product channels across several particle tiles."""

import math
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

PARTICLE_COUNT = 65536
GRID_CELL_COUNT = 8
PARTICLES_PER_CELL = PARTICLE_COUNT // GRID_CELL_COUNT
INITIAL_PRODUCT_COUNT_PER_CELL = 1
BACKGROUND_DENSITY = 1.0e20
BACKGROUND_TEMPERATURE = 300.0
BACKGROUND_MASS = 1000.0 * picmi.constants.m_e
CROSS_SECTION = 1.0e-22
THIRD_BODY_DENSITY = 1.0e25
THREE_BODY_CROSS_SECTION = CROSS_SECTION / THIRD_BODY_DENSITY
ELECTRON_ENERGY_EV = 1000.0
COLLISION_OPTICAL_DEPTH = 0.4
MAJORANT_SAFETY = 1.01

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e
KB = picmi.constants.kb

GAMMA = 1.0 + ELECTRON_ENERGY_EV * Q_E / (M_E * C**2)
ELECTRON_SPEED = C * math.sqrt(1.0 - 1.0 / GAMMA**2)
NU_MAX = MAJORANT_SAFETY * BACKGROUND_DENSITY * 4.0 * CROSS_SECTION * ELECTRON_SPEED
DT = COLLISION_OPTICAL_DEPTH / NU_MAX
NEUTRAL_VELOCITY_STD = math.sqrt(KB * BACKGROUND_TEMPERATURE / BACKGROUND_MASS)

source_dir = Path(__file__).resolve().parent

grid = picmi.Cartesian1DGrid(
    number_of_cells=[GRID_CELL_COUNT],
    lower_bound=[0.0],
    upper_bound=[100.0 * GRID_CELL_COUNT],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=2,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)


def make_product_species(name, charge, mass):
    return picmi.Species(
        name=name,
        charge=charge,
        mass=mass,
        initial_distribution=picmi.UniformDistribution(
            density=1.0,
            directed_velocity=[0.0, 0.0, 0.0],
        ),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, ELECTRON_SPEED],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
ion_a = make_product_species("ion_a", "q_e", BACKGROUND_MASS - M_E)
ion_b = make_product_species("ion_b", "q_e", BACKGROUND_MASS - M_E)
negative_a = make_product_species("negative_a", "-q_e", BACKGROUND_MASS + M_E)
negative_b = make_product_species("negative_b", "-q_e", BACKGROUND_MASS + M_E)

collision = picmi.MCCCollisions(
    name="mcc",
    species=electrons,
    background_density=BACKGROUND_DENSITY,
    background_temperature=BACKGROUND_TEMPERATURE,
    background_mass=BACKGROUND_MASS,
    scattering_processes={
        "ionization_a": {
            "cross_section": str(
                source_dir / "background_mcc_attachment_ionization.txt"
            ),
            "energy": 15.0,
            "species": ion_a,
        },
        "ionization_b": {
            "cross_section": str(
                source_dir / "background_mcc_attachment_ionization.txt"
            ),
            "energy": 25.0,
            "species": ion_b,
        },
        "attachment_dissociative": {
            "cross_section": str(source_dir / "background_mcc_attachment_m2.txt"),
            "cross_section_units": "m2",
            "species": negative_a,
        },
        "attachment_three_body": {
            "cross_section": str(source_dir / "background_mcc_attachment_m5.txt"),
            "cross_section_units": "m5",
            "third_body_density": THIRD_BODY_DENSITY,
            "species": negative_b,
        },
    },
    nu_max=NU_MAX,
)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=DT,
    max_steps=1,
    warpx_collisions=[collision],
    verbose=1,
)
sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[PARTICLES_PER_CELL], grid=grid
    ),
)
for species in [ion_a, ion_b, negative_a, negative_b]:
    sim.add_species(
        species,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[INITIAL_PRODUCT_COUNT_PER_CELL], grid=grid
        ),
    )


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def get_component(container, component):
    return np.concatenate(
        [to_numpy(pti[component]) for pti in container.iterator(level=0)]
    )


def get_ids(container):
    return np.concatenate(
        [
            to_numpy(libwarpx.amr.unpack_ids(pti["idcpu"]))
            for pti in container.iterator(level=0)
        ]
    )


sim.initialize_inputs()
sim.initialize_warpx()
electron_container = sim.particles.get("electrons")
initial_electron_w = get_component(electron_container, "w")
initial_electron_ids = get_ids(electron_container)
initial_product_ids = {
    name: get_ids(sim.particles.get(name))
    for name in ["ion_a", "ion_b", "negative_a", "negative_b"]
}
sim.step(1)


def product_data(name):
    container = sim.particles.get(name)
    ids = get_ids(container)
    created = ~np.isin(ids, initial_product_ids[name])
    components = {
        component: get_component(container, component)[created]
        for component in ["ux", "uy", "uz", "w"]
    }
    return np.count_nonzero(created), components, ids[created]


electron_count = electron_container.number_of_particles(only_local=True)
electron_w = get_component(electron_container, "w")
electron_ids = get_ids(electron_container)
ion_a_count, ion_a_data, ion_a_ids = product_data("ion_a")
ion_b_count, ion_b_data, ion_b_ids = product_data("ion_b")
negative_a_count, negative_a_data, negative_a_ids = product_data("negative_a")
negative_b_count, negative_b_data, negative_b_ids = product_data("negative_b")

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "background_mcc_attachment_results.npz",
        electron_count=electron_count,
        ion_a_events=ion_a_count,
        ion_b_events=ion_b_count,
        negative_a_events=negative_a_count,
        negative_b_events=negative_b_count,
        initial_electron_w=initial_electron_w,
        initial_electron_ids=initial_electron_ids,
        initial_product_ids=np.concatenate(list(initial_product_ids.values())),
        electron_w=electron_w,
        electron_ids=electron_ids,
        ion_a_w=ion_a_data["w"],
        ion_a_ids=ion_a_ids,
        ion_a_ux=ion_a_data["ux"],
        ion_a_uy=ion_a_data["uy"],
        ion_a_uz=ion_a_data["uz"],
        ion_b_w=ion_b_data["w"],
        ion_b_ids=ion_b_ids,
        ion_b_ux=ion_b_data["ux"],
        ion_b_uy=ion_b_data["uy"],
        ion_b_uz=ion_b_data["uz"],
        negative_a_w=negative_a_data["w"],
        negative_a_ids=negative_a_ids,
        negative_a_ux=negative_a_data["ux"],
        negative_a_uy=negative_a_data["uy"],
        negative_a_uz=negative_a_data["uz"],
        negative_b_w=negative_b_data["w"],
        negative_b_ids=negative_b_ids,
        negative_b_ux=negative_b_data["ux"],
        negative_b_uy=negative_b_data["uy"],
        negative_b_uz=negative_b_data["uz"],
        neutral_velocity_std=NEUTRAL_VELOCITY_STD,
    )

sim.finalize()
