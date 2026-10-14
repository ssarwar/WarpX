#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare empty cylindrical and annular Poisson solves with exact solutions."""

import argparse

import numpy as np

from pywarpx import picmi

parser = argparse.ArgumentParser()
parser.add_argument("--annular", action="store_true")
args = parser.parse_args()
rmin, rmax, cells = (1.0 if args.annular else 0.0), 2.0, 64
grid = picmi.CylindricalGrid(
    number_of_cells=[cells, 8],
    lower_bound=[rmin, 0],
    upper_bound=[rmax, 1],
    lower_boundary_conditions=["dirichlet" if args.annular else "none", "periodic"],
    upper_boundary_conditions=["dirichlet", "periodic"],
    lower_boundary_conditions_particles=[
        "absorbing" if args.annular else "none",
        "periodic",
    ],
    upper_boundary_conditions_particles=["absorbing", "periodic"],
    n_azimuthal_modes=1,
    warpx_max_grid_size=16,
    warpx_blocking_factor=8,
    warpx_potential_lo_r=0,
    warpx_potential_hi_r=1,
)
sim = picmi.Simulation(
    solver=picmi.ElectrostaticSolver(
        grid=grid, method="Multigrid", required_precision=1e-13
    ),
    time_step_size=1e-12,
    max_steps=1,
    particle_shape=1,
    verbose=0,
)
# Boundary voltages trigger the initial Poisson solve. Inspect that static
# solution without an unnecessary second solve of an already-converged zero RHS.
sim.initialize_inputs()
sim.initialize_warpx()
phi = sim.fields.get("phi_fp", level=0)
values = phi[...]
dr = (rmax - rmin) / cells
radius = rmin + dr * phi.imesh(0)
if args.annular:
    # Laplace's equation: r*d(phi)/dr is constant. Dirichlet data give
    # phi(r) = log(r/rmin)/log(rmax/rmin), independent of z.
    exact = np.log(radius / rmin) / np.log(rmax / rmin)
    # The numerical radial operator is second order, not an exact logarithm.
    tolerance = (dr / rmin) ** 2
    np.testing.assert_allclose(values[0], 0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(values[-1], 1, rtol=0, atol=1e-12)
else:
    # Regularity at the axis and phi(rmax)=1 give the constant solution.
    exact = np.ones_like(radius)
    tolerance = 1e-10
np.testing.assert_allclose(
    values, np.broadcast_to(exact[:, None], values.shape), rtol=0, atol=tolerance
)
sim.finalize()
print("PASS: radial Dirichlet potential agrees with the analytical solution")
