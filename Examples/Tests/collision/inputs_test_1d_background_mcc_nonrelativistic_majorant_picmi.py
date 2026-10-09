#!/usr/bin/env python3
"""Exercise the automatic MCC majorant for non-electron projectiles."""

import math
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

PARTICLE_COUNT = 131072
BACKGROUND_DENSITY = 1.0e20
PARTICLE_MASS = 4.0 * picmi.constants.m_p
PEAK_ENERGY_EV = 10.0
PEAK_CROSS_SECTION = 1.0e-20
MAJORANT_OPTICAL_DEPTH = 0.1
TAIL_RELATIVE_PROPER_SPEED = 100.0 * picmi.constants.c

C = picmi.constants.c
Q_E = picmi.constants.q_e

# Invert ParticleUtils::getCollisionEnergy for a stationary equal-mass target.
energy_mass = PEAK_ENERGY_EV * Q_E / C**2
tau = energy_mass * (4.0 * PARTICLE_MASS + energy_mass) / (
    2.0 * PARTICLE_MASS**2
)
relative_proper_speed = C * math.sqrt(tau * (tau + 2.0))
collision_speed = relative_proper_speed / (tau + 1.0)
collision_frequency = (
    BACKGROUND_DENSITY * PEAK_CROSS_SECTION * collision_speed
)
time_step = MAJORANT_OPTICAL_DEPTH / collision_frequency
tail_cross_section = PEAK_CROSS_SECTION * collision_speed / C

cross_section_path = Path("background_mcc_nonrelativistic_majorant.txt").resolve()
np.savetxt(
    cross_section_path,
    np.array(
        [
            [0.0, 0.0],
            [PEAK_ENERGY_EV, PEAK_CROSS_SECTION],
            [2.0 * PEAK_ENERGY_EV, 0.0],
        ]
    ),
)
tail_cross_section_path = Path(
    "background_mcc_nonrelativistic_majorant_tail.txt"
).resolve()
np.savetxt(
    tail_cross_section_path,
    np.array(
        [
            [0.0, tail_cross_section],
            [PEAK_ENERGY_EV, tail_cross_section],
        ]
    ),
)

grid = picmi.Cartesian1DGrid(
    number_of_cells=[1],
    lower_bound=[0.0],
    upper_bound=[1.0],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)

peak_ions = picmi.Species(
    name="peak_ions",
    charge="q_e",
    mass=PARTICLE_MASS,
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, relative_proper_speed],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
tail_ions = picmi.Species(
    name="tail_ions",
    charge="q_e",
    mass=PARTICLE_MASS,
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, TAIL_RELATIVE_PROPER_SPEED],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
peak_collision = picmi.MCCCollisions(
    name="peak_mcc",
    species=peak_ions,
    background_density=BACKGROUND_DENSITY,
    background_temperature=0.0,
    background_mass=PARTICLE_MASS,
    scattering_processes={
        "elastic": {
            "cross_section": str(cross_section_path),
            "scattering_angle_model": "backward",
        }
    },
)
tail_collision = picmi.MCCCollisions(
    name="tail_mcc",
    species=tail_ions,
    background_density=BACKGROUND_DENSITY,
    background_temperature=0.0,
    background_mass=PARTICLE_MASS,
    scattering_processes={
        "elastic": {
            "cross_section": str(tail_cross_section_path),
            "scattering_angle_model": "backward",
        }
    },
)
simulation = picmi.Simulation(
    solver=solver,
    time_step_size=time_step,
    max_steps=1,
    warpx_collisions=[peak_collision, tail_collision],
    verbose=1,
)
simulation.add_species(
    peak_ions,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid
    ),
)
simulation.add_species(
    tail_ions,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid
    ),
)

simulation.step(1)


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


peak_container = simulation.particles.get("peak_ions")
peak_uz = np.concatenate(
    [to_numpy(pti["uz"]) for pti in peak_container.iterator(level=0)]
).astype(np.float64)
peak_collided = peak_uz < 0.5 * relative_proper_speed

tail_container = simulation.particles.get("tail_ions")
tail_uz = np.concatenate(
    [to_numpy(pti["uz"]) for pti in tail_container.iterator(level=0)]
).astype(np.float64)
tail_collided = tail_uz < 0.5 * TAIL_RELATIVE_PROPER_SPEED
tail_gamma = math.sqrt(1.0 + (TAIL_RELATIVE_PROPER_SPEED / C) ** 2)
tail_collision_speed = TAIL_RELATIVE_PROPER_SPEED / tail_gamma

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "background_mcc_nonrelativistic_majorant_results.npz",
        particle_count=PARTICLE_COUNT,
        optical_depth=MAJORANT_OPTICAL_DEPTH,
        peak_event_count=np.count_nonzero(peak_collided),
        peak_maximum_collided_speed=np.max(np.abs(peak_uz[peak_collided])),
        peak_relative_proper_speed=relative_proper_speed,
        tail_event_count=np.count_nonzero(tail_collided),
        tail_maximum_collided_speed=np.max(np.abs(tail_uz[tail_collided])),
        tail_relative_proper_speed=TAIL_RELATIVE_PROPER_SPEED,
        tail_collision_speed=tail_collision_speed,
        tail_collision_speed_over_c=tail_collision_speed / C,
    )

simulation.finalize()
