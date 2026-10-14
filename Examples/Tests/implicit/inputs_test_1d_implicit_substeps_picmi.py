#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Force implicit retries; check free streaming and the integrated Lorentz law."""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument("--solver", choices=["theta", "semi"], required=True)
args = parser.parse_args()
comm = MPI.COMM_WORLD
if comm.rank == 0:
    # The nonlinear diagnostic appends to existing files, including on reruns.
    Path("picard.txt").unlink(missing_ok=True)
    Path("substep_result.json").unlink(missing_ok=True)
comm.Barrier()
c, me, qe, eps0, mu0 = (
    picmi.constants.c,
    picmi.constants.m_e,
    picmi.constants.q_e,
    picmi.constants.ep0,
    picmi.constants.mu0,
)
dt, electric, magnetic = 1e-10, 1.0, 1e-7
# The full step is too large for the prescribed Picard iteration budget;
# the implicit solver must retry with smaller steps. The eventual solution
# stays nonrelativistic, although the failed iterates need not do so.
density = (8 / dt) ** 2 * eps0 * me / qe**2
grid = picmi.Cartesian1DGrid(
    number_of_cells=[16],
    lower_bound=[0],
    upper_bound=[1],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=8,
    warpx_blocking_factor=8,
)
nonlinear = picmi.PicardNonlinearSolver(
    max_iterations=20,
    relative_tolerance=1e-9,
    absolute_tolerance=0,
    require_convergence=True,
    diagnostic_file="picard.txt",
    verbose=False,
)
scheme = (
    picmi.ThetaImplicitEMEvolveScheme(nonlinear_solver=nonlinear, theta=0.5)
    if args.solver == "theta"
    else picmi.SemiImplicitEMEvolveScheme(nonlinear_solver=nonlinear)
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    warpx_evolve_scheme=scheme,
    warpx_current_deposition_algo="direct",
    warpx_use_filter=False,
    time_step_size=dt,
    max_steps=1,
    particle_shape=1,
    verbose=0,
)
electrons = picmi.Species(
    name="electrons",
    particle_type="electron",
    initial_distribution=picmi.UniformDistribution(density=density),
)
ions = picmi.Species(
    name="ions",
    particle_type="proton",
    initial_distribution=picmi.UniformDistribution(density=density),
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
)
tracers = picmi.Species(
    name="tracers",
    mass=me,
    charge=0,
    initial_distribution=picmi.UniformDistribution(
        density=1, directed_velocity=[0, 0, 1e6]
    ),
    warpx_do_not_deposit=True,
)
for species in [electrons, ions, tracers]:
    sim.add_species(
        species, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[4])
    )
sim.add_diagnostic(
    picmi.ReducedDiagnostic(diag_type="FieldPoyntingFlux", name="flux", period=1)
)
sim.initialize_inputs()
warpx.E_ext_grid_init_style = "constant"
warpx.E_external_grid = [electric, 0, 0]
warpx.B_ext_grid_init_style = "constant"
warpx.B_external_grid = [0, magnetic, 0]
sim.initialize_warpx()


def host(values):
    return values.get() if hasattr(values, "get") else np.asarray(values)


def particles(name):
    local = []
    for tile in sim.particles.get(name).iterator(level=0):
        local.append(
            {
                key: host(tile[key]).copy()
                for key in ["idcpu", "z", "ux", "uy", "uz", "w"]
            }
        )
    tiles = [tile for rank in comm.allgather(local) for tile in rank]
    arrays = {key: np.concatenate([tile[key] for tile in tiles]) for key in tiles[0]}
    order = np.argsort(arrays["idcpu"])
    return {key: values[order] for key, values in arrays.items()}


initial = {name: particles(name) for name in ["electrons", "tracers"]}
sim.step(1)
final = {name: particles(name) for name in initial}
for name in initial:
    np.testing.assert_array_equal(initial[name]["idcpu"], final[name]["idcpu"])
u0 = initial["tracers"]
speed = u0["uz"] / np.sqrt(1 + (u0["ux"] ** 2 + u0["uy"] ** 2 + u0["uz"] ** 2) / c**2)
np.testing.assert_allclose(
    final["tracers"]["z"], u0["z"] + speed * dt, rtol=0, atol=2e-13
)
for axis in ["ux", "uy", "uz"]:
    np.testing.assert_allclose(final["tracers"][axis], u0[axis], rtol=2e-13, atol=1e-30)

# Uniform fields obey du_x/dt = (q/m)(E_x - v_z B_y), so the exact
# integrated Lorentz law gives int(E_x dt) = (m/q) Delta u_x + B_y Delta z.
# It also holds for the implicit midpoint update. This independently checks
# the time integral of the face Poynting flux, including every accepted
# substep and excluding the initial/end-of-step duplicate samples.
before, after = initial["electrons"], final["electrons"]
integral_e = np.average(
    -me / qe * (after["ux"] - before["ux"]) + magnetic * (after["z"] - before["z"]),
    weights=before["w"],
)
expected_flux_integral = magnetic / mu0 * integral_e
sim.finalize()
if comm.rank == 0:
    history = np.loadtxt("picard.txt", skiprows=1, ndmin=2)
    accepted = history[(history[:, 2] < 20) & (history[:, 5] < 1e-9)]
    assert len(accepted) > 1, "The regression must exercise actual implicit substeps"
    accepted_dt = np.diff(np.r_[0, accepted[:, 1]])
    assert np.all(accepted_dt > 0)
    np.testing.assert_allclose(accepted_dt.sum(), dt, rtol=2e-13)
    flux = np.loadtxt("diags/reducedfiles/flux.txt", skiprows=1, ndmin=2)
    np.testing.assert_array_equal(flux[0, 4:], 0)
    np.testing.assert_allclose(
        flux[-1, 4:],
        [-expected_flux_integral, expected_flux_integral],
        rtol=3e-8,
        atol=1e-25,
    )
    Path("substep_result.json").write_text(
        json.dumps(
            dict(
                solver=args.solver,
                accepted_substeps=len(accepted),
                accepted_dt=accepted_dt.tolist(),
                expected_flux_integral=expected_flux_integral,
                measured_flux_integral=flux[-1, 5],
            ),
            indent=2,
        )
        + "\n"
    )
    print(
        "PASS: implicit substep duration, free streaming and integrated Lorentz/Poynting law"
    )
