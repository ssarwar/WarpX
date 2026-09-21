#!/usr/bin/env python3
"""Independent normalization, propagation and continuity checks for a rigid beam."""

import argparse
import math

import numpy as np

from pywarpx import picmi

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
    warpx_current_deposition_algo="direct"
    if args.solver.startswith("semi_implicit")
    else None,
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
normalization = {"peak_current": current}
if args.shape == 1:
    normalization = {
        "peak_density": current
        / (qe * abs(velocity) * 2 * np.pi * sigma_r**2 * (-np.expm1(-32)))
    }
elif args.shape == 2:
    normalization = {
        "bunch_charge": current
        * np.sqrt(2 * np.pi)
        * sigma_t
        * math.erf(8 / np.sqrt(2))
    }
beam = picmi.FluidSpecies(
    name="beam",
    model="rigid_beam",
    particle_type="proton",
    velocity_z=velocity,
    r_rms=math.sqrt(2) * sigma_r,
    sigma_t=sigma_t,
    **normalization,
    pulse_times=times,
    pulse_amplitudes=amplitudes,
    initialize_self_fields=args.self_fields,
)
sim.add_fluid_species(beam)
sim.initialize_inputs()
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
if args.self_fields:
    er = np.squeeze(host(sim.fields.get("Efield_fp", dir="r", level=0)[...]))
    ez = np.squeeze(host(sim.fields.get("Efield_fp", dir="z", level=0)[...]))
    bt = np.squeeze(host(sim.fields.get("Bfield_fp", dir="theta", level=0)[...]))
    if args.solver != "PSATD":
        # Independent cylindrical Gauss law on the Yee mesh, including the axis.
        # The outer Dirichlet potential node is excluded from the Poisson equation.
        radial_divergence = np.empty((nr, nz))
        radial_divergence[0] = 4 * er[0, :nz] / dr
        radius = np.arange(1, nr) * dr
        radial_divergence[1:] = (
            (radius[:, None] + dr / 2) * er[1:, :nz]
            - (radius[:, None] - dr / 2) * er[:-1, :nz]
        ) / (radius[:, None] * dr)
        divergence = radial_divergence + (ez[:nr] - np.roll(ez[:nr], 1, axis=1)) / dz
        charge = qe * old[:nr, :nz]
        np.testing.assert_allclose(
            picmi.constants.ep0 * divergence, charge, rtol=0, atol=1e-8 * charge.max()
        )
    # A rigid axial source satisfies B_theta = v_z E_r/c^2. Yee needs an axial
    # average to B_theta's staggering; PSATD already stores collocated fields.
    centered_er = er if args.solver == "PSATD" else (er[:, :-1] + er[:, 1:]) / 2
    expected_bt = velocity / c**2 * centered_er
    np.testing.assert_allclose(bt, expected_bt, rtol=0, atol=2e-13 * abs(bt).max())
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
