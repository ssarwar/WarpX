#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check the physical time interval sampled by collision subcycling."""

import argparse
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import picmi

parser = argparse.ArgumentParser()
parser.add_argument(
    "--solver", choices=["Yee", "semi_implicit_em", "semi_implicit_mm"], default="Yee"
)
parser.add_argument("--subcycles", type=int, default=4)
parser.add_argument("--early", action="store_true")
args = parser.parse_args()
dt = 1e-12
grid = picmi.Cartesian1DGrid(
    number_of_cells=[8],
    lower_bound=[0],
    upper_bound=[1],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    warpx_max_grid_size=4,
    warpx_blocking_factor=4,
)
electrons = picmi.Species(
    name="electrons",
    particle_type="electron",
    initial_distribution=picmi.UniformDistribution(
        density=1e8, directed_velocity=[0, 0, 1e6]
    ),
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
    warpx_do_not_deposit=True,
)
if MPI.COMM_WORLD.rank == 0:
    Path("elastic.txt").write_text("0 1e-20\n1000000 1e-20\n")
MPI.COMM_WORLD.Barrier()
# Future gas is identically absent over the entire simulated interval [0,dt].
# Early gas instead checks that implicit subcycling samples the elapsed step,
# rather than evaluating every substep at the end of the field push.
condition = f"t < {0.625 * dt:.17g}" if args.early else f"t > {1.125 * dt:.17g}"
collision = picmi.MCCCollisions(
    name="scatter",
    species=electrons,
    background_density=f"if({condition},1e30,0)",
    max_background_density=1e30,
    background_temperature=0,
    background_mass=28 * 1.66053906660e-27,
    ndt_subcycle=args.subcycles,
    scattering_processes={"elastic": {"cross_section": "elastic.txt"}},
)
implicit = args.solver != "Yee"
mass_matrices = args.solver == "semi_implicit_mm"
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    max_steps=1,
    time_step_size=dt,
    particle_shape=1,
    verbose=0,
    warpx_collisions=[collision],
    warpx_current_deposition_algo="direct" if implicit else None,
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.NewtonNonlinearSolver(
            relative_tolerance=1e-12,
            use_mass_matrices_jacobian=mass_matrices,
            use_mass_matrices_pc=mass_matrices,
            pc_type=picmi.JacobiPreconditioner() if mass_matrices else None,
            linear_solver=picmi.GMRESLinearSolver(relative_tolerance=1e-12),
        )
    )
    if implicit
    else None,
)
sim.add_species(
    electrons, picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[64])
)
sim.step(1)
scattered = 0
for tile in sim.particles.get("electrons").iterator(level=0):
    values = tile["ux"]
    values = values.get() if hasattr(values, "get") else np.asarray(values)
    scattered += np.count_nonzero(values)
scattered = MPI.COMM_WORLD.allreduce(scattered)
sim.finalize()
if args.early:
    assert scattered > 0, "Substeps did not sample gas present during the elapsed step"
else:
    assert scattered == 0, "Collisions sampled gas after the simulated interval"
if MPI.COMM_WORLD.rank == 0:
    print("PASS: collision substep placement stays inside the physical step")
