#!/usr/bin/env python3
"""Check continuity at the domain faces while a prescribed beam enters and exits."""

import argparse

import numpy as np
from scipy.integrate import quad
from scipy.special import erf

from pywarpx import picmi

parser = argparse.ArgumentParser()
parser.add_argument("--solver", default="Yee")
parser.add_argument("--negative", action="store_true")
parser.add_argument("--filter", action="store_true")
args = parser.parse_args()
implicit = args.solver.startswith("semi_implicit")
velocity = (
    (-1 if args.negative else 1)
    * picmi.constants.c
    * np.sqrt(
        1
        - (
            1
            + 800e6 * picmi.constants.q_e / (picmi.constants.m_p * picmi.constants.c**2)
        )
        ** -2
    )
)
grid = picmi.CylindricalGrid(
    number_of_cells=[32, 64],
    lower_bound=[0, 0],
    upper_bound=[0.016, 0.032],
    lower_boundary_conditions=["none", "none"],
    upper_boundary_conditions=["none", "none"],
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    n_azimuthal_modes=1,
    warpx_max_grid_size=32,
    warpx_blocking_factor=8,
)
dt = 5e-13
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(
        grid=grid,
        method="Yee",
        source_smoother=picmi.BinomialSmoother(n_pass=[0, 3]) if args.filter else None,
    ),
    time_step_size=dt,
    max_steps=3,
    particle_shape=3,
    warpx_use_filter=args.filter,
    warpx_current_deposition_algo="direct" if implicit else None,
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.NewtonNonlinearSolver(
            relative_tolerance=1e-12,
            use_mass_matrices_jacobian=args.solver == "semi_implicit_mm",
        )
    )
    if implicit
    else None,
    verbose=0,
)
sim.add_fluid_species(
    picmi.FluidSpecies(
        name="beam",
        model="rigid_beam",
        particle_type="proton",
        velocity_z=velocity,
        sigma_r=0.002,
        sigma_t=25e-12,
        peak_current=0.6,
        pulse_times=[0],
    )
)
sim.initialize_inputs()
sim.initialize_warpx()


def host(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def density():
    mf = sim.fields.get("fluid_density_beam", level=0)
    if args.filter:
        raw = host(mf[...])
        # WarpX's RZ volume filter closes smoothing flux at charge boundaries.
        for _ in range(3):
            result = raw.copy()
            result[:, 1:-1] = (raw[:, :-2] + 2 * raw[:, 1:-1] + raw[:, 2:]) / 4
            result[:, 0] = (3 * raw[:, 0] + raw[:, 1]) / 4
            result[:, -1] = (3 * raw[:, -1] + raw[:, -2]) / 4
            raw = result
        return raw
    return host(mf[...]).copy()


def check_charge(number, time):
    sigma = abs(velocity) * 25e-12
    dz = dr = 0.0005

    def projected_interval(u):
        shape = 0.75 - u * u if abs(u) < 0.5 else 0.5 * (1.5 - abs(u)) ** 2
        a = (-dz / 2 - velocity * time - dz * u) / (np.sqrt(2) * sigma)
        b = (0.032 + dz / 2 - velocity * time - dz * u) / (np.sqrt(2) * sigma)
        return shape * sigma * np.sqrt(np.pi / 2) * (erf(b) - erf(a))

    # Boundary nodes own a full dual cell. Integrate their projected Gaussian
    # support independently; the clipped beam is never renormalized.
    expected = (
        0.6
        / abs(velocity)
        * quad(projected_interval, -1.5, 1.5, points=[-0.5, 0.5], epsabs=1e-14)[0]
    )
    volume = 2 * np.pi * np.arange(33) * dr**2 * dz
    volume[0] = np.pi * dr**2 * dz / 3
    actual = picmi.constants.q_e * np.sum(number * volume[:, None])
    np.testing.assert_allclose(actual, expected, rtol=2e-11)


old = density()
check_charge(old, 0)
for step in range(3):
    sim.step(1)
    new = density()
    check_charge(new, (step + 1) * dt)
    # Include the two current faces bounding the physical boundary nodes.
    current = host(sim.fields.get("current_fp", "z", level=0)[:, -1j:2j])
    change = picmi.constants.q_e * (new - old)
    divergence = dt / 0.0005 * np.diff(current, axis=1)
    np.testing.assert_allclose(
        change + divergence,
        0,
        rtol=0,
        atol=1e-12 * picmi.constants.q_e * np.max(new),
    )
    old = new
print("PASS: local continuity includes both physical domain faces")
sim.finalize()
