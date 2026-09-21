#!/usr/bin/env python3
"""Independent normalization, propagation and continuity checks for a rigid beam."""

import argparse
import math

import numpy as np

from pywarpx import algo, picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument(
    "--solver",
    choices=["Yee", "PSATD", "semi_implicit_em", "semi_implicit_mm"],
    default="Yee",
)
parser.add_argument("--shape", type=int, default=3)
parser.add_argument("--self-fields", action="store_true")
parser.add_argument("--negative", action="store_true")
args = parser.parse_args()

nr, nz = 32, 128
rmax, zmin, zmax = 0.016, -0.08, 0.08
dr, dz = rmax / nr, (zmax - zmin) / nz
sigma_r, sigma_t, current = 0.002, 25e-12, 0.6
times, amplitudes = np.array([0.0, 1e-11]), np.array([0.75, 0.25])
qe, mp, c = picmi.constants.q_e, picmi.constants.m_p, picmi.constants.c
gamma = 1 + 800e6 * qe / (mp * c**2)
velocity = c * np.sqrt(1 - gamma**-2) * (-1 if args.negative else 1)
dt = 5e-13
grid = picmi.CylindricalGrid(
    number_of_cells=[nr, nz],
    lower_bound=[0.0, zmin],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    lower_boundary_conditions_particles=["none", "periodic"],
    upper_boundary_conditions_particles=["absorbing", "periodic"],
    warpx_max_grid_size=32,
    warpx_blocking_factor=8,
    n_azimuthal_modes=1,
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(
        grid=grid,
        method="PSATD" if args.solver == "PSATD" else "Yee",
        stencil_order=[16, 16],
        **({"warpx_current_correction": True} if args.solver == "PSATD" else {}),
    ),
    time_step_size=dt,
    max_steps=3,
    particle_shape=args.shape,
    warpx_use_filter=False,
    warpx_current_deposition_algo="direct" if args.solver.startswith("semi_implicit") else None,
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.NewtonNonlinearSolver(
            relative_tolerance=1e-12,
            use_mass_matrices_jacobian=args.solver == "semi_implicit_mm",
            use_mass_matrices_pc=args.solver == "semi_implicit_mm",
            pc_type=picmi.JacobiPreconditioner()
            if args.solver == "semi_implicit_mm"
            else None,
            linear_solver=picmi.GMRESLinearSolver(relative_tolerance=1e-12),
        )
    )
    if args.solver.startswith("semi_implicit")
    else None,
    verbose=0,
)
sim.initialize_inputs()
algo.particle_shape = args.shape
warpx.get_bucket("fluids").species_names = ["beam"]
beam = warpx.get_bucket("beam")
beam.model = "rigid_beam"
beam.species_type = "proton"
beam.velocity_z = velocity
beam.sigma_r = sigma_r
beam.sigma_t = sigma_t
beam.peak_current = current
beam.pulse_times = times
beam.pulse_amplitudes = amplitudes
beam.initialize_self_fields = args.self_fields
sim.initialize_warpx()


def host(value):
    return value.get() if hasattr(value, "get") else np.asarray(value)


def density():
    mf = sim.fields.get("fluid_density_beam", level=0)
    return np.squeeze(host(mf[...])).copy()


mf = sim.fields.get("fluid_density_beam", level=0)
r, z = host(mf.mesh("r")), host(mf.mesh("z"))
volume = 2 * np.pi * r * dr * dz
if args.solver != "PSATD":
    volume[0] = np.pi * dr**2 * dz / 3
    # The periodic high node duplicates the low node.
    z = z[:-1]
expected_charge = (
    current * math.sqrt(2 * math.pi) * sigma_t * math.erf(8 / math.sqrt(2))
)


def check(time):
    number = density()[:, :nz] * volume[:, None]
    charge = qe * np.sum(number)
    # The outer boundary cuts only a Gaussian tail below 2e-12 of the population.
    np.testing.assert_allclose(charge, expected_charge, rtol=2e-11)
    mean = np.sum(number * z[None, :]) / np.sum(number)
    expected_mean = velocity * (time - np.dot(times, amplitudes))
    np.testing.assert_allclose(mean, expected_mean, rtol=0, atol=2e-13)
    assert not any("fluid_momentum_density_beam" in key for key in sim.fields.list())


check(0.0)
old = density()
for step in range(1, 4):
    sim.step(1)
    check(step * dt)
    new = density()
    if args.solver != "PSATD":
        jz = np.squeeze(host(sim.fields.get("current_fp", dir="z", level=0)[...]))
        change = qe * (new[:, 1:-1] - old[:, 1:-1])
        divergence = dt / dz * (jz[:, 1:] - jz[:, :-1])
        np.testing.assert_allclose(
            change + divergence, 0, rtol=0, atol=1e-12 * qe * np.max(new)
        )
    old = new
assert np.max(np.abs(host(sim.fields.get("Efield_fp", dir="r", level=0)[...]))) > 0
print("PASS: Gaussian charge, centroid, propagation and discrete continuity")
sim.finalize()
