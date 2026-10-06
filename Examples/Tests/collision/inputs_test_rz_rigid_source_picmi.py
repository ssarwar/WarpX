#!/usr/bin/env python3
"""Check integrated rigid-source budgets and immediate checkpoint restoration."""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy.special import erf

from pywarpx import amrex, libwarpx, picmi, warpx

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
parser.add_argument("--openpmd", action="store_true")
parser.add_argument("--subcycles", type=int, default=1)
parser.add_argument("--alpha", action="store_true")
parser.add_argument("--load-balance", action="store_true")
parser.add_argument("--max-grid-size", type=int, default=16)
parser.add_argument("--cells-z", type=int, default=64)
parser.add_argument("--bad-density", action="store_true")
parser.add_argument("--bad-temperature", action="store_true")
args = parser.parse_args()
if os.environ.get("WARPX_TEST_RESTART_MUTATION") or args.bad_density or args.bad_temperature:
    amrex.throw_exception = 1
    amrex.signal_handling = 0
qe, mp, me, c = (
    picmi.constants.q_e,
    picmi.constants.m_p,
    picmi.constants.m_e,
    picmi.constants.c,
)
energy, sigma_r, sigma_t, peak_current = 8e8, 0.002, 25e-12, 0.6
charge_state = 2 if args.alpha else 1
beam_mass = (4 if args.alpha else 1) * mp
beam_energy = energy * beam_mass / mp
peak_current *= charge_state
gamma = 1 + beam_energy * qe / (beam_mass * c * c)
speed = c * np.sqrt(1 - gamma**-2)
sigma_z = speed * sigma_t
cutoff, ngas, dt, steps = 2.0, 1e21, 1e-12, 6
rmax, zmin, zmax = 8 * sigma_r, -4 * sigma_z, 4 * sigma_z
grid = picmi.CylindricalGrid(
    number_of_cells=[16, args.cells_z],
    lower_bound=[0, zmin],
    upper_bound=[rmax, zmax],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    lower_boundary_conditions_particles=["none", "periodic"],
    upper_boundary_conditions_particles=["absorbing", "periodic"],
    warpx_max_grid_size=args.max_grid_size,
    warpx_blocking_factor=8,
    n_azimuthal_modes=1,
)
solver = picmi.ElectromagneticSolver(
    grid=grid, method="PSATD" if args.solver == "PSATD" else "Yee", stencil_order=[8, 8]
)
frozen = dict(
    warpx_do_not_push=True, warpx_do_not_gather=True, warpx_do_not_deposit=True
)
beam = picmi.FluidSpecies(
    name="beam",
    model="rigid_beam",
    mass=beam_mass,
    charge=charge_state * qe,
    kinetic_energy=beam_energy,
    sigma_r=sigma_r,
    sigma_t=sigma_t,
    peak_current=peak_current,
    cutoff_z=cutoff,
    pulse_times=[0.0, 2e-12],
    pulse_amplitudes=[0.75, 0.25],
)
cases, collisions = {}, []
for target, mass in [("N2", 28.0134), ("O2", 31.9988)]:
    for mode in ["kinetic", "immobile", "fine"]:
        name = f"{target}_{mode}"
        electron = picmi.Species(name="e_" + name, particle_type="electron", **frozen)
        if mode == "kinetic":
            ion = picmi.Species(
                name="i_" + name,
                charge=qe,
                mass=mass * 1.66053906660e-27 - me,
                **frozen,
            )
        else:
            ion = picmi.FluidSpecies(
                name="i_" + name,
                model="immobile",
                charge=qe,
                mass=mass * 1.66053906660e-27 - me,
            )
        weight = 0.125 if mode == "fine" else 1.0
        cases[name] = (electron, ion, target, mode, weight)
        collisions.append(
            picmi.ProtonImpactIonizationCollisions(
                name=name,
                species=beam,
                product_species=[electron, ion],
                ionization_target=target,
                background_density="-(1+t)" if args.bad_density else ngas,
                background_temperature="-(1+t)" if args.bad_temperature else 0,
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
    warpx_load_balance_intervals=1 if args.load_balance else None,
    warpx_load_balance_efficiency_ratio_threshold=0.1 if args.load_balance else None,
    warpx_load_balance_costs_update="heuristic" if args.load_balance else None,
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
sim.add_fluid_species(beam)
for electron, ion, _, mode, _ in cases.values():
    sim.add_species(electron, layout=None)
    if mode == "kinetic":
        sim.add_species(ion, layout=None)
    else:
        sim.add_fluid_species(ion)
sim.add_diagnostic(picmi.Checkpoint(name="chk", period=2, write_dir="diags"))
diagnostic_types = [
    "ParticleNumber",
    "ParticleCharge",
    "ParticleEnergy",
    "ParticleMomentum",
    "PrescribedSourceBudget",
]
for kind in diagnostic_types:
    sim.add_diagnostic(picmi.ReducedDiagnostic(diag_type=kind, name=kind, period=1))
plot_fields = [
    "rho",
    "rho_beam",
    "rho_i_N2_immobile",
    "rho_e_N2_immobile",
    "fluid_density_beam",
    "fluid_current_beamz",
    "fluid_density_i_N2_immobile",
    "N2_immobile_product_weight_remainder",
    "N2_immobile_emitted_number",
    "N2_immobile_electron_energy",
    "N2_immobile_binding_energy",
    "part_per_cell_beam",
    "part_per_cell_i_N2_immobile",
    "part_per_cell_e_N2_immobile",
]
if args.bad_density or args.bad_temperature:
    # Isolate parser validation from unrelated output diagnostics.
    plot_fields = [name for name in plot_fields if not name.startswith("part_per_cell")]
sim.add_diagnostic(
    picmi.FieldDiagnostic(
        name="fields", grid=grid, period=2, data_list=plot_fields, write_dir="diags"
    )
)
if args.openpmd:
    sim.add_diagnostic(
        picmi.FieldDiagnostic(
            name="pmd",
            grid=grid,
            period=2,
            data_list=plot_fields,
            write_dir="diags/pmd",
            warpx_format="openpmd",
            warpx_openpmd_backend="h5",
            warpx_dump_rz_modes=True,
        )
    )
sim.initialize_inputs()
mutation = os.environ.get("WARPX_TEST_RESTART_MUTATION")
if mutation == "beam":
    beam.fluid.peak_current *= 2
elif mutation == "source":
    warpx.get_bucket("N2_immobile").fixed_product_weight = 0.75
elif mutation == "remove":
    warpx.get_bucket("fluids").species_names = ["beam"]
elif mutation == "particle_diagnostic":
    warpx.get_bucket("ParticleEnergy").species = ["beam"]
elif mutation == "load_balance":
    warpx.get_bucket("algo").load_balance_intervals = 1
sim.initialize_warpx()
# The source/recoil checks also run in a serial WarpX build. Initialize mpi4py
# only when AMReX uses MPI, preserving the same collective checks in MPI runs.
mpi_enabled = libwarpx.amr.Config.have_mpi
if mpi_enabled:
    from mpi4py import MPI


def host(value):
    return value.get() if hasattr(value, "get") else np.asarray(value)


def field(name):
    return host(sim.fields.get(name, level=0)[...]).copy()


def population(name):
    local = sum(
        float(host(tile["w"]).sum())
        for tile in sim.particles.get(name).iterator(level=0)
    )
    return MPI.COMM_WORLD.allreduce(local) if mpi_enabled else local


def particle_state(name):
    arrays = []
    for tile in sim.particles.get(name).iterator(level=0):
        arrays.append(
            np.column_stack(
                [host(tile[key]) for key in ["r", "theta", "z", "w", "ux", "uy", "uz"]]
            )
        )
    local = np.concatenate(arrays) if arrays else np.empty((0, 7))
    values = np.concatenate(MPI.COMM_WORLD.allgather(local)) if mpi_enabled else local
    return values[np.lexsort(values.T[::-1])]


def state():
    result = {}
    for name, (electron, ion, _, mode, _) in cases.items():
        for suffix in ["product_weight_remainder", "source_budget", "sampling_counter"]:
            result[name + "_" + suffix] = field(name + "_" + suffix)
        if mode != "kinetic":
            result[ion.name] = field("fluid_density_" + ion.name)
        result[electron.name] = np.array(population(electron.name))
        result[electron.name + "_phase_space"] = particle_state(electron.name)
        if electron.name == "e_N2_immobile":
            positions = result[electron.name + "_phase_space"]
            result[electron.name + "_cell_counts"] = np.histogram2d(
                positions[:, 0],
                positions[:, 2],
                bins=[np.linspace(0, rmax, 17), np.linspace(zmin, zmax, args.cells_z + 1)],
            )[0]
    for kind in ["Efield_fp", "Bfield_fp", "current_fp"]:
        for direction in ["r", "theta", "z"]:
            result[kind + "_" + direction] = host(
                sim.fields.get(kind, direction, level=0)[...]
            ).copy()
    result["beam"] = field("fluid_density_beam")
    result["beam_current"] = host(
        sim.fields.get("fluid_current_beam", "z", level=0)[...]
    ).copy()
    return result


start = sim.extension.warpx.getistep(lev=0)
initial_distribution = list(sim.extension.warpx.DistributionMap(0).ProcessorMap())
redistributed = False
if args.restart:
    assert start == 2
    reference = Path(args.restart).parents[1] / "state_2.npz"
    with np.load(reference) as saved:
        for name, value in state().items():
            np.testing.assert_allclose(value, saved[name], rtol=2e-14, atol=0)

number = (
    peak_current
    * np.sqrt(2 * np.pi)
    * sigma_t
    * erf(cutoff / np.sqrt(2))
    / (charge_state * qe)
)
previous = {name: population(values[0].name) for name, values in cases.items()}
for step in range(start + 1, steps + 1):
    sim.step(1)
    redistributed |= (
        list(sim.extension.warpx.DistributionMap(0).ProcessorMap())
        != initial_distribution
    )
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
            * charge_state**2
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
    beam_population = number
    expected_populations = {
        electron.name: population(electron.name)
        for electron, _, _, _, _ in cases.values()
    }
    for electron, ion, _, _, _ in cases.values():
        expected_populations[ion.name] = expected_populations[electron.name]
    expected_populations["beam"] = beam_population
    source_totals = {
        name: (
            field(name + "_product_weight_remainder").sum(),
            field(name + "_source_budget").reshape(-1, 4).sum(axis=0),
        )
        for name in cases
    }
    import re

    def reduced(kind):
        path = Path("diags/reducedfiles") / (kind + ".txt")
        columns = [
            re.sub(r"^\[\d+\]", "", key)
            for key in path.read_text().splitlines()[0].lstrip("#").split()
        ]
        row = np.loadtxt(path, ndmin=2)[-1]
        np.testing.assert_allclose(row[1], step * dt, rtol=2e-14)
        return dict(zip(columns, row))

    counts = reduced("ParticleNumber")
    charges = reduced("ParticleCharge")
    energies = reduced("ParticleEnergy")
    momenta = reduced("ParticleMomentum")
    for name, population_value in expected_populations.items():
        np.testing.assert_allclose(
            counts[name + "_weight()"], population_value, rtol=3e-12
        )
        charge = -qe if name.startswith("e_") else qe
        if name == "beam":
            charge *= charge_state
        np.testing.assert_allclose(
            charges[name + "(C)"], charge * population_value, rtol=3e-12
        )
    for name in [species.name for species in sim.fluid_species]:
        assert counts[name + "_macroparticles()"] == 0
        if name != "beam":
            assert energies[name + "(J)"] == 0
            assert momenta[name + "_z(kg*m/s)"] == 0
    np.testing.assert_allclose(
        energies["beam(J)"], number * beam_energy * qe, rtol=3e-12
    )
    np.testing.assert_allclose(
        momenta["beam_z(kg*m/s)"], number * gamma * beam_mass * speed, rtol=3e-12
    )
    np.testing.assert_allclose(
        energies["total_mean(J)"],
        energies["total(J)"] / sum(expected_populations.values()),
        rtol=3e-12,
    )
    budgets = reduced("PrescribedSourceBudget")
    for name, (pending, totals) in source_totals.items():
        np.testing.assert_allclose(budgets[name + "_pending()"], pending, rtol=3e-12)
        np.testing.assert_allclose(budgets[name + "_emitted()"], totals[0], rtol=3e-12)
        for comp, suffix in enumerate(
            ["electron_energy", "binding_energy", "discarded_ion_energy"], 1
        ):
            np.testing.assert_allclose(
                budgets[name + "_" + suffix + "(J)"], totals[comp], rtol=3e-12
            )
    saved_state = state()
    if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
        np.savez(f"state_{step}.npz", **saved_state)
if args.load_balance:
    assert redistributed, (
        "The test must redistribute persistent fields across MPI ranks"
    )
if args.restart or args.load_balance:
    reference_directory = (
        Path(args.restart).parents[1]
        if args.restart
        else Path(f"../test_rz_rigid_source_{args.solver}_picmi")
    )
    with np.load(reference_directory / f"state_{steps}.npz") as saved:
        for name, value in state().items():
            if name.startswith(("Efield_fp", "Bfield_fp")):
                # MPI redistribution changes reduction and transform ordering.
                # Bound roundoff in the vector field norm: a whole component
                # can vanish by symmetry, so its own norm is not an error scale.
                # Chemistry and persistent budgets retain stricter checks below.
                kind = name.rsplit("_", 1)[0]
                scale = max(
                    np.max(np.abs(saved[kind + "_" + direction]))
                    for direction in ["r", "theta", "z"]
                )
                rounding = 64 * np.finfo(value.dtype).eps * scale
                np.testing.assert_allclose(
                    value, saved[name], rtol=0, atol=rounding, err_msg=name
                )
            else:
                np.testing.assert_allclose(
                    value, saved[name], rtol=2e-14, atol=1e-300, err_msg=name
                )
print(
    "PASS: N2/O2 rigid-source budgets, caps, ion charge footprints and deterministic restart"
)
sim.finalize()
