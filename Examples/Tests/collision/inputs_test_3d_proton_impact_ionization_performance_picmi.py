#!/usr/bin/env python3
"""Bounded-creation performance smoke test for proton-impact ionization."""

import argparse
import time

import numpy as np

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--cells-per-direction", type=int, default=8)
parser.add_argument("--particles-per-cell-direction", type=int, default=4)
parser.add_argument("--steps", type=int, default=5)
parser.add_argument("--max-products-per-cell", type=int, default=8)
args = parser.parse_args()

assert args.cells_per_direction > 0
assert args.particles_per_cell_direction > 0
assert args.steps > 0
assert args.max_products_per_cell > 0

cell_count = args.cells_per_direction**3
particles_per_cell = args.particles_per_cell_direction**3
particle_count = cell_count * particles_per_cell

C = picmi.constants.c
M_P = picmi.constants.m_p
Q_E = picmi.constants.q_e
PROJECTILE_ENERGY = 800.0e6
PROJECTILE_REST_ENERGY = M_P * C**2 / Q_E
PROJECTILE_GAMMA = 1.0 + PROJECTILE_ENERGY / PROJECTILE_REST_ENERGY
PROJECTILE_PROPER_SPEED = C * np.sqrt(PROJECTILE_GAMMA**2 - 1.0)

grid = picmi.Cartesian3DGrid(
    number_of_cells=[args.cells_per_direction] * 3,
    lower_bound=[0.0, 0.0, 0.0],
    upper_bound=[float(args.cells_per_direction)] * 3,
    lower_boundary_conditions=["periodic"] * 3,
    upper_boundary_conditions=["periodic"] * 3,
    lower_boundary_conditions_particles=["periodic"] * 3,
    upper_boundary_conditions_particles=["periodic"] * 3,
    warpx_max_grid_size=args.cells_per_direction,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)

beam = picmi.Species(
    particle_type="proton",
    name="beam",
    initial_distribution=picmi.UniformDistribution(
        density=1.0e8,
        directed_velocity=[0.0, 0.0, PROJECTILE_PROPER_SPEED],
    ),
    warpx_do_not_push=True,
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    warpx_do_not_push=True,
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
ions = picmi.Species(
    name="ions",
    charge="q_e",
    mass=28.0134 * 1.660_539_066_60e-27 - picmi.constants.m_e,
    warpx_do_not_push=True,
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
collision = picmi.ProtonImpactIonizationCollisions(
    name="pjg",
    species=beam,
    product_species=[electrons, ions],
    ionization_target="N2",
    background_density=1.0e21,
    background_temperature=300.0,
    fixed_product_weight=3.0e5,
    max_products_per_cell=args.max_products_per_cell,
)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=1.0e-9,
    max_steps=args.steps,
    warpx_collisions=[collision],
    warpx_random_seed=314159,
    verbose=1,
)
sim.add_species(
    beam,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[args.particles_per_cell_direction] * 3,
        grid=grid,
    ),
)
for product in [electrons, ions]:
    sim.add_species(
        product,
        layout=picmi.GriddedLayout(n_macroparticle_per_cell=[0, 0, 0], grid=grid),
    )

sim.initialize_inputs()
initialization_start = time.perf_counter()
sim.initialize_warpx()
initialization_elapsed = time.perf_counter() - initialization_start
step_start = time.perf_counter()
sim.step(args.steps)
step_elapsed = time.perf_counter() - step_start

electron_count = sim.particles.get("electrons").number_of_particles(only_local=False)
ion_count = sim.particles.get("ions").number_of_particles(only_local=False)
particle_calls = particle_count * args.steps

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "proton_impact_ionization_performance_results.npz",
        cells_per_direction=args.cells_per_direction,
        cell_count=cell_count,
        particles_per_cell=particles_per_cell,
        particle_count=particle_count,
        steps=args.steps,
        max_products_per_cell=args.max_products_per_cell,
        electron_count=electron_count,
        ion_count=ion_count,
        initialization_elapsed=initialization_elapsed,
        step_elapsed=step_elapsed,
        particle_calls_per_second=particle_calls / step_elapsed,
        product_pairs_per_second=electron_count / step_elapsed,
    )
    print(
        "Proton-impact ionization performance: "
        f"cells={cell_count}, particles={particle_count}, pairs={electron_count}, "
        f"initialization={initialization_elapsed:.6f} s, step={step_elapsed:.6f} s, "
        f"throughput={particle_calls / step_elapsed:.3e} particle-calls/s"
    )

sim.finalize()
