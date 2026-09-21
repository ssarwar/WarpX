#!/usr/bin/env python3
"""Check integrated rigid-source budgets and immediate checkpoint restoration."""

import argparse
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI
from scipy.special import erf

from pywarpx import picmi, warpx

sys.path.insert(
    0,
    str(
        Path(__file__).resolve().parents[3] / "Tools/Algorithms/ProtonImpactIonization"
    ),
)
from calibrated_pjg import total_cross_section

parser = argparse.ArgumentParser()
parser.add_argument(
    "--solver",
    default="Yee",
    choices=["Yee", "PSATD", "semi_implicit_em", "semi_implicit_mm"],
)
parser.add_argument("--restart")
parser.add_argument("--subcycles", type=int, default=1)
args = parser.parse_args()
qe, mp, me, c = (
    picmi.constants.q_e,
    picmi.constants.m_p,
    picmi.constants.m_e,
    picmi.constants.c,
)
energy, sigma_r, sigma_t, peak_current = 8e8, 0.002, 25e-12, 0.6
gamma = 1 + energy * qe / (mp * c * c)
speed = c * np.sqrt(1 - gamma**-2)
sigma_z = speed * sigma_t
cutoff, ngas, dt, steps = 2.0, 1e21, 1e-12, 6
rmax, zmin, zmax = 8 * sigma_r, -4 * sigma_z, 4 * sigma_z
grid = picmi.CylindricalGrid(
    number_of_cells=[16, 32],
    lower_bound=[0, zmin],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    lower_boundary_conditions_particles=["none", "periodic"],
    upper_boundary_conditions_particles=["absorbing", "periodic"],
    warpx_max_grid_size=16,
    warpx_blocking_factor=8,
    n_azimuthal_modes=1,
)
solver = picmi.ElectromagneticSolver(
    grid=grid, method="PSATD" if args.solver == "PSATD" else "Yee", stencil_order=[8, 8]
)
frozen = dict(
    warpx_do_not_push=True, warpx_do_not_gather=True, warpx_do_not_deposit=True
)
beam = picmi.Species(name="beam", particle_type="proton")
cases, collisions = {}, []
for target, mass in [("N2", 28.0134), ("O2", 31.9988)]:
    for mode in ["kinetic", "immobile", "fine"]:
        name = f"{target}_{mode}"
        electron = picmi.Species(name="e_" + name, particle_type="electron", **frozen)
        ion = picmi.Species(
            name="i_" + name, charge=qe, mass=mass * 1.66053906660e-27 - me, **frozen
        )
        weight = 0.125 if mode == "fine" else 1.0
        cases[name] = (electron, ion, target, mode, weight)
        collisions.append(
            picmi.ProtonImpactIonizationCollisions(
                name=name,
                species=beam,
                product_species=[electron, ion],
                ionization_target=target,
                background_density=ngas,
                background_temperature=0,
                fixed_product_weight=weight,
                max_products_per_cell=2,
                ndt_subcycle=args.subcycles,
            )
        )
sim = picmi.Simulation(
    solver=solver,
    time_step_size=dt,
    max_steps=steps,
    particle_shape=3,
    warpx_collisions=collisions,
    warpx_amr_restart=args.restart,
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
for electron, ion, _, mode, _ in cases.values():
    sim.add_species(electron, layout=None)
    if mode == "kinetic":
        sim.add_species(ion, layout=None)
sim.add_diagnostic(picmi.Checkpoint(name="chk", period=2, write_dir="diags"))
sim.initialize_inputs()
fluid_names = ["beam"]
for _, ion, _, mode, _ in cases.values():
    if mode != "kinetic":
        fluid_names.append(ion.name)
        bucket = warpx.get_bucket(ion.name)
        bucket.model, bucket.mass, bucket.charge = "immobile", ion.mass, qe
warpx.get_bucket("fluids").species_names = fluid_names
bucket = warpx.get_bucket("beam")
bucket.model, bucket.species_type = "rigid_beam", "proton"
bucket.kinetic_energy, bucket.sigma_r, bucket.sigma_t = energy, sigma_r, sigma_t
bucket.peak_current, bucket.cutoff_z = peak_current, cutoff
bucket.pulse_times, bucket.pulse_amplitudes = [0.0, 2e-12], [0.75, 0.25]
sim.initialize_warpx()


def host(value):
    return value.get() if hasattr(value, "get") else np.asarray(value)


def field(name):
    return host(sim.fields.get(name, level=0)[...]).copy()


def population(name):
    local = sum(
        float(host(tile["w"]).sum())
        for tile in sim.particles.get(name).iterator(level=0)
    )
    return MPI.COMM_WORLD.allreduce(local)


def state():
    result = {}
    for name, (electron, ion, _, mode, _) in cases.items():
        for suffix in ["product_weight_remainder", "source_budget", "sampling_counter"]:
            result[name + "_" + suffix] = field(name + "_" + suffix)
        if mode != "kinetic":
            result[ion.name] = field("fluid_density_" + ion.name)
        result[electron.name] = np.array(population(electron.name))
    return result


start = sim.extension.warpx.getistep(lev=0)
if args.restart:
    assert start == 2
    reference = Path(args.restart).parents[1] / "state_2.npz"
    with np.load(reference) as saved:
        for name, value in state().items():
            np.testing.assert_allclose(value, saved[name], rtol=2e-14, atol=0)

number = peak_current * np.sqrt(2 * np.pi) * sigma_t * erf(cutoff / np.sqrt(2)) / qe
previous = {name: population(values[0].name) for name, values in cases.items()}
for step in range(start + 1, steps + 1):
    sim.step(1)
    for name, (electron, ion, target, mode, weight) in cases.items():
        emitted = population(electron.name)
        remaining = field(name + "_product_weight_remainder")
        expected = (
            number
            * ngas
            * float(total_cross_section(target, energy))
            * 1e-4
            * speed
            * step
            * dt
        )
        # The independent SDCS quadrature differs from the existing table by
        # <1e-3 (tested over its full energy range in test_pjg_model.cpp).
        np.testing.assert_allclose(emitted + remaining.sum(), expected, rtol=1e-3)
        assert np.all(remaining >= 0) and np.max(remaining) < weight
        totals = field(name + "_source_budget")
        np.testing.assert_allclose(totals[..., 0].sum(), emitted, rtol=2e-14)
        assert np.all(totals[..., 1:] >= 0)
        erho = host(
            sim.particles.get(electron.name).get_charge_density(lev=0, local=False)[...]
        )
        if mode == "kinetic":
            irho = host(
                sim.particles.get(ion.name).get_charge_density(lev=0, local=False)[...]
            )
            np.testing.assert_allclose(population(ion.name), emitted, rtol=2e-14)
        else:
            irho = qe * field("fluid_density_" + ion.name)
        np.testing.assert_allclose(
            irho + erho, 0, atol=3e-14 * np.max(np.abs(erho)), rtol=0
        )
        assert emitted > previous[name]
        previous[name] = emitted
    for target in ["N2", "O2"]:
        budgets = [
            previous[target + "_" + mode]
            + field(target + "_" + mode + "_product_weight_remainder").sum()
            for mode in ["kinetic", "immobile", "fine"]
        ]
        np.testing.assert_allclose(budgets, budgets[0], rtol=2e-14)
    saved_state = state()
    if MPI.COMM_WORLD.rank == 0:
        np.savez(f"state_{step}.npz", **saved_state)
if args.restart:
    with np.load(Path(args.restart).parents[1] / f"state_{steps}.npz") as saved:
        for name, value in state().items():
            np.testing.assert_allclose(value, saved[name], rtol=2e-14, atol=1e-300)
print(
    "PASS: N2/O2 rigid-source budgets, caps, ion charge footprints and deterministic restart"
)
sim.finalize()
