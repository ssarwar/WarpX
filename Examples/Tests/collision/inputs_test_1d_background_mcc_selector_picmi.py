#!/usr/bin/env python3
"""Exercise one-draw Background MCC selection and the PICMI nu_max input."""

import math
from pathlib import Path

import numpy as np

from pywarpx import Collisions, libwarpx, picmi

PARTICLE_COUNT = 131072
INITIAL_ION_COUNT = 1
BACKGROUND_DENSITY = 1.0e20
CROSS_SECTION = 1.0e-20
ELECTRON_ENERGY_EV = 1000.0
COLLISION_OPTICAL_DEPTH = 0.5

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e

GAMMA = 1.0 + ELECTRON_ENERGY_EV*Q_E/(M_E*C**2)
ELECTRON_SPEED = C*math.sqrt(1.0 - 1.0/GAMMA**2)
ELECTRON_PROPER_SPEED = GAMMA*ELECTRON_SPEED
NU_MAX = BACKGROUND_DENSITY*2.0*CROSS_SECTION*ELECTRON_SPEED
DT = COLLISION_OPTICAL_DEPTH/NU_MAX

source_dir = Path(__file__).resolve().parent

grid = picmi.Cartesian1DGrid(
    number_of_cells=[1],
    lower_bound=[0.0],
    upper_bound=[100.0],
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
        directed_velocity=[0.0, 0.0, ELECTRON_SPEED],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
ion_product = picmi.Species(
    name="ion_product",
    charge="q_e",
    mass=1000.0*M_E,
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, 0.0],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)

processes = {
    "elastic": {
        "cross_section": str(source_dir / "background_mcc_selector_elastic.txt"),
        "scattering_angle_model": "backward",
    },
    "ionization": {
        "cross_section": str(source_dir / "background_mcc_selector_ionization.txt"),
        "energy": 15.0,
        "species": ion_product,
    },
}
collision = picmi.MCCCollisions(
    name="mcc",
    species=electrons,
    background_density=BACKGROUND_DENSITY,
    background_temperature=0.0,
    background_mass=100.0*M_E,
    scattering_processes=processes,
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
        n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid
    ),
)
sim.add_species(
    ion_product,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[INITIAL_ION_COUNT], grid=grid
    ),
)

sim.initialize_inputs()
collision_bucket = next(
    bucket for bucket in Collisions.collisions_list if bucket.instancename == "mcc"
)
assert math.isclose(collision_bucket.nu_max, NU_MAX, rel_tol=1.0e-15)

sim.initialize_warpx()
sim.step(1)

electron_container = sim.particles.get("electrons")
ion_container = sim.particles.get("ion_product")

electron_count = electron_container.number_of_particles(only_local=True)
ion_count = ion_container.number_of_particles(only_local=True)


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


uz = np.concatenate(
    [to_numpy(pti["uz"]) for pti in electron_container.iterator(level=0)]
)
elastic_events = int(np.count_nonzero(uz < -0.9*ELECTRON_PROPER_SPEED))
ionization_events = int(ion_count - INITIAL_ION_COUNT)

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "background_mcc_selector_results.npz",
        electron_count=electron_count,
        elastic_events=elastic_events,
        ionization_events=ionization_events,
    )

sim.finalize()
