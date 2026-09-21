#!/usr/bin/env python3
"""Independent clipped pulse-train yield, quiet gaps and fractional-state restart."""

import argparse
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI
from scipy.integrate import quad
from scipy.special import erf

from pywarpx import amr, picmi

sys.path.insert(
    0,
    str(
        Path(__file__).resolve().parents[3] / "Tools/Algorithms/ProtonImpactIonization"
    ),
)
from calibrated_pjg import total_cross_section

parser = argparse.ArgumentParser()
parser.add_argument("--solver", default="Yee")
parser.add_argument("--restart")
args = parser.parse_args()
qe, mp, c = picmi.constants.q_e, picmi.constants.m_p, picmi.constants.c
gamma = 1 + 8e8 * qe / (mp * c**2)
speed = c * np.sqrt(1 - gamma**-2)
sigma_t, sigma_r = 1e-12, 0.002
sigma_z = speed * sigma_t
dt, steps, checkpoint = 1e-13, 180, 75
ngas, cutoff = 1e21, 2.0
rmax, zmax = 4 * sigma_r, 4 * sigma_z
times, amplitudes = [0.0, 15e-12], [1.0, 0.7]
grid = picmi.CylindricalGrid(
    number_of_cells=[16, 32],
    lower_bound=[0, -zmax],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "none"],
    upper_boundary_conditions=["none", "none"],
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    n_azimuthal_modes=1,
    warpx_max_grid_size=16,
    warpx_blocking_factor=8,
)
beam = picmi.FluidSpecies(
    name="beam",
    particle_type="proton",
    model="rigid_beam",
    kinetic_energy=8e8,
    sigma_r=sigma_r,
    sigma_t=sigma_t,
    peak_current=0.6,
    cutoff_z=cutoff,
    pulse_times=times,
    pulse_amplitudes=amplitudes,
)
electrons, ions, collisions = [], [], []
for name, density, quadrature in [
    ("uniform", ngas, 2),
    ("varying", None, 8),
    ("zero", 0.0, 2),
]:
    electron = picmi.Species(
        name="e_" + name,
        particle_type="electron",
        warpx_do_not_push=True,
        warpx_do_not_gather=True,
    )
    ion = picmi.FluidSpecies(
        name="i_" + name,
        model="immobile",
        charge=qe,
        mass=28.0134 * 1.66053906660e-27 - picmi.constants.m_e,
    )
    # Linear gas variation in z and t exercises space-time source quadrature.
    gas = (
        density
        if density is not None
        else f"{ngas}*(1+0.1*z/{zmax}+0.1*t/{steps * dt})"
    )
    collisions.append(
        picmi.ProtonImpactIonizationCollisions(
            name=name,
            species=beam,
            product_species=[electron, ion],
            ionization_target="N2",
            background_density=gas,
            background_temperature=0.0,
            fixed_product_weight=1.0,
            max_products_per_cell=2,
            gas_quadrature_points=quadrature,
            ndt_subcycle=2,
        )
    )
    electrons.append(electron)
    ions.append(ion)
implicit = args.solver.startswith("semi_implicit")
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(
        grid=grid,
        method="PSATD" if args.solver == "PSATD" else "Yee",
        stencil_order=[8, 8],
    ),
    time_step_size=dt,
    max_steps=steps,
    particle_shape=2,
    verbose=0,
    warpx_use_filter=False,
    warpx_collisions=collisions,
    warpx_current_deposition_algo="direct" if implicit else None,
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.NewtonNonlinearSolver(
            relative_tolerance=1e-10,
            use_mass_matrices_jacobian=args.solver == "semi_implicit_mm",
        )
    )
    if implicit
    else None,
)
sim.add_fluid_species(beam)
for electron, ion in zip(electrons, ions):
    sim.add_species(electron, layout=None)
    sim.add_fluid_species(ion)
sim.add_diagnostic(picmi.Checkpoint(name="chk", period=checkpoint))
sim.initialize_inputs()
if args.restart:
    amr.restart = args.restart
sim.initialize_warpx()


def host(value):
    return value.get() if hasattr(value, "get") else np.asarray(value)


def snapshot():
    state = {}
    for collision, ion in zip(collisions, ions):
        for suffix in ["product_weight_remainder", "sampling_counter", "source_budget"]:
            name = collision.name + "_" + suffix
            state[name] = host(sim.fields.get(name, level=0)[...]).copy()
        name = "fluid_density_" + ion.name
        state[name] = host(sim.fields.get(name, level=0)[...]).copy()
    for electron in electrons:
        state[electron.name] = np.array(
            sim.particles.get(electron.name).sum_particle_weight(False)
        )
    for kind in ["Efield_fp", "Bfield_fp", "current_fp"]:
        for direction in ["r", "theta", "z"]:
            state[kind + direction] = host(
                sim.fields.get(kind, direction, level=0)[...]
            ).copy()
    return state


def compare(state, saved, resumed_steps=0):
    for name, array in state.items():
        scale = np.max(np.abs(saved[name]), initial=0)
        operations = 64
        if name.startswith(("Efield", "Bfield", "current")):
            # MPI redistribution changes reduction order. Bound accumulated
            # roundoff in the two curl updates (four operations each) per step;
            # immediate restoration and all chemistry retain the strict bound.
            operations += 8 * resumed_steps
        error = operations * np.finfo(float).eps
        np.testing.assert_allclose(
            array,
            saved[name],
            rtol=0,
            atol=error / (1 - error) * scale,
            err_msg=name,
        )


start = sim.extension.warpx.getistep(lev=0)
if args.restart:
    assert start == checkpoint
    with np.load(Path(args.restart).parents[1] / "gap.npz") as saved:
        compare(snapshot(), saved)
gap = None
for step in range(start + 1, steps + 1):
    sim.step(1)
    if step in [65, 85]:
        state = snapshot()
        if gap is not None:
            # No production, remainder loss or sequence advance while both pulses are absent.
            for name in state:
                if (
                    "source_budget" in name
                    or "sampling_counter" in name
                    or "remainder" in name
                ):
                    np.testing.assert_array_equal(state[name], gap[name])
        gap = state
    if step == checkpoint:
        saved_state = snapshot()
        if MPI.COMM_WORLD.rank == 0:
            np.savez_compressed("gap.npz", **saved_state)

state = snapshot()
# Integrate the physical clipped Gaussian directly; no mesh or production CDF is reused.
radial = 2 * np.pi * sigma_r**2 * (-np.expm1(-(rmax**2) / (2 * sigma_r**2)))
n0 = 0.6 / (qe * speed * 2 * np.pi * sigma_r**2 * (-np.expm1(-32)))


def spatial(time, varying):
    result = 0.0
    for arrival, amplitude in zip(times, amplitudes):
        center = speed * (time - arrival)
        a, b = (
            max(-zmax, center - cutoff * sigma_z),
            min(zmax, center + cutoff * sigma_z),
        )
        if b <= a:
            continue
        za, zb = (a - center) / sigma_z, (b - center) / sigma_z
        integral = (
            sigma_z * np.sqrt(np.pi / 2) * (erf(zb / np.sqrt(2)) - erf(za / np.sqrt(2)))
        )
        first = center * integral + sigma_z**2 * (
            np.exp(-(za**2) / 2) - np.exp(-(zb**2) / 2)
        )
        result += amplitude * (
            integral * (1 + 0.1 * time / (steps * dt)) + 0.1 * first / zmax
            if varying
            else integral
        )
    return result


cross_section = float(total_cross_section("N2", 8e8)) * 1e-4
for name in ["uniform", "varying"]:
    # Integrate in ps to keep the quadrature's absolute tolerance dimensionless.
    integrated = (
        quad(
            lambda ps: spatial(ps * 1e-12, name == "varying"),
            0,
            steps * dt / 1e-12,
            points=[2, 6, 9, 13, 17],
            epsabs=1e-14,
            epsrel=1e-10,
        )[0]
        * 1e-12
    )
    expected = n0 * radial * ngas * cross_section * speed * integrated
    actual = (
        state[name + "_source_budget"][..., 0].sum()
        + state[name + "_product_weight_remainder"].sum()
    )
    # Existing PJG interpolation error is <1e-3; gas midpoint error at q=8 is <1e-5 here.
    np.testing.assert_allclose(actual, expected, rtol=1.01e-3)
assert state["e_zero"] == 0
assert np.max(state["zero_source_budget"]) == 0
assert np.max(state["zero_product_weight_remainder"]) == 0
if args.restart:
    with np.load(Path(args.restart).parents[1] / "final.npz") as saved:
        compare(state, saved, steps - checkpoint)
elif MPI.COMM_WORLD.rank == 0:
    np.savez_compressed("final.npz", **state)
sim.finalize()
print("PASS: clipped pulse train, varying gas, zero source and gap restart")
