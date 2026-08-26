#!/usr/bin/env python3
"""Performance smoke test for many attachment channels sharing one product."""

import math
import time
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

PARTICLE_COUNT = 32768
INITIAL_NEGATIVE_COUNT = 1
PROCESS_COUNT = 64
BACKGROUND_DENSITY = 1.0e20
BACKGROUND_MASS = 1000.0 * picmi.constants.m_e
CROSS_SECTION_PER_PROCESS = 1.0e-24
ELECTRON_ENERGY_EV = 1000.0
MAJORANT_OPTICAL_DEPTH = 0.2

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e

GAMMA = 1.0 + ELECTRON_ENERGY_EV * Q_E / (M_E * C**2)
ELECTRON_SPEED = C * math.sqrt(1.0 - 1.0 / GAMMA**2)
RATE_UPPER_BOUND = (
    BACKGROUND_DENSITY
    * PROCESS_COUNT
    * CROSS_SECTION_PER_PROCESS
    * C
)
DT = MAJORANT_OPTICAL_DEPTH / RATE_UPPER_BOUND

cross_section_path = Path("background_mcc_many_attachment.txt").resolve()
energies = np.geomspace(1.0e-6, 2.5e6, 1025)
sigmas = np.full_like(energies, CROSS_SECTION_PER_PROCESS)
sigmas[-1] = 0.0
np.savetxt(cross_section_path, np.column_stack((energies, sigmas)))
cross_section_file = str(cross_section_path)

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
negative_ions = picmi.Species(
    name="negative_ions",
    charge="-q_e",
    mass=BACKGROUND_MASS + M_E,
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, 0.0],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)

processes = {
    f"attachment{index}": {
        "cross_section": cross_section_file,
        "species": negative_ions,
    }
    for index in range(PROCESS_COUNT)
}
collision = picmi.MCCCollisions(
    name="mcc",
    species=electrons,
    background_density=BACKGROUND_DENSITY,
    background_temperature=0.0,
    background_mass=BACKGROUND_MASS,
    scattering_processes=processes,
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
    negative_ions,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[INITIAL_NEGATIVE_COUNT], grid=grid
    ),
)

sim.initialize_warpx()
start = time.perf_counter()
sim.step(1)
elapsed = time.perf_counter() - start

electron_count = sim.particles.get("electrons").number_of_particles(
    only_local=True
)
negative_count = sim.particles.get("negative_ions").number_of_particles(
    only_local=True
)

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "background_mcc_many_attachment_results.npz",
        electron_count=electron_count,
        attachment_events=negative_count - INITIAL_NEGATIVE_COUNT,
        elapsed=elapsed,
    )

sim.finalize()
