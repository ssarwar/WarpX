#!/usr/bin/env python3
"""Manufactured constant-field diagnostics, including cumulative probe restart."""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument("--restart")
parser.add_argument("--continue-output", action="store_true")
parser.add_argument("--stock", action="store_true")
parser.add_argument("--record-only", action="store_true")
parser.add_argument("--nonzero-flux", action="store_true")
parser.add_argument(
    "--solver",
    default="Yee",
    choices=["Yee", "PSATD", "semi_implicit_em", "semi_implicit_mm"],
)
args = parser.parse_args()
if args.continue_output:
    assert args.restart, "Output continuation needs a checkpoint"
    if MPI.COMM_WORLD.rank == 0:
        checkpoint = Path(args.restart).resolve()
        checkpoint_step = int(checkpoint.name[-6:])
        output = Path("diags/reducedfiles")
        output.mkdir(parents=True, exist_ok=True)
        prefix = {}
        for path in (checkpoint.parents[1] / "diags/reducedfiles").glob("*.txt"):
            lines = path.read_text().splitlines(keepends=True)
            kept = [
                line for line in lines[1:] if int(line.split()[0]) <= checkpoint_step
            ]
            (output / path.name).write_text(lines[0] + "".join(kept))
            prefix[path.stem] = len(kept)
        assert prefix
        Path("continuation_prefix.json").write_text(json.dumps(prefix) + "\n")
    MPI.COMM_WORLD.Barrier()
dt, magnetic_field = 1e-13, 0.01
rmax, length = 0.016, 0.064
grid = picmi.CylindricalGrid(
    number_of_cells=[16, 128],
    lower_bound=[0, -length / 2],
    upper_bound=[rmax, length / 2],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    warpx_max_grid_size=16,
    warpx_blocking_factor=8,
    n_azimuthal_modes=1,
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(
        grid=grid,
        method="PSATD" if args.solver == "PSATD" else "Yee",
        stencil_order=[8, 8],
    ),
    time_step_size=dt,
    max_steps=4,
    particle_shape=2,
    warpx_amr_restart=args.restart,
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
electron = picmi.Species(
    name="electrons",
    particle_type="electron",
    initial_distribution=picmi.UniformDistribution(
        density=1e10,
        directed_velocity=[0, 0, 1e6],
        rms_velocity=[1e5, 1e5, 1e5],
    ),
    warpx_do_not_push=True,
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
sim.add_species(
    electron,
    layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[1, 1, 1]),
)
if not args.stock:
    sim.add_fluid_species(
        picmi.FluidSpecies(
            name="beam",
            model="rigid_beam",
            particle_type="proton",
            kinetic_energy=800e6,
            sigma_r=0.002,
            sigma_t=25e-12,
            peak_current=0.6,
            pulse_times=[0],
            warpx_do_not_deposit=True,
        )
    )
    for name, charge in [("positive", 1), ("negative", -1)]:
        sim.add_fluid_species(
            picmi.FluidSpecies(
                name=name,
                model="immobile",
                mass=28 * 1.66053906660e-27,
                charge=charge * picmi.constants.q_e,
                initial_density=1e10,
            )
        )
for kind in [
    "FieldEnergy",
    "FieldPoyntingFlux",
    "ParticleEnergy",
    "ParticleMomentum",
    "ParticleNumber",
    "ParticleCharge",
    "Timestep",
]:
    if args.stock and kind == "ParticleCharge":
        # This reduction is part of the prescribed-fluid branch.
        continue
    sim.add_diagnostic(picmi.ReducedDiagnostic(diag_type=kind, name=kind, period=1))
for kind in ["BeamRelevant", "ParticleExtrema"]:
    sim.add_diagnostic(
        picmi.ReducedDiagnostic(diag_type=kind, name=kind, period=1, species=electron)
    )
sim.add_diagnostic(
    picmi.ReducedDiagnostic(
        diag_type="ParticleHistogram",
        name="histogram",
        period=1,
        species=electron,
        bin_number=16,
        bin_min=0,
        bin_max=rmax,
        histogram_function="sqrt(x*x+y*y)",
    )
)
for integrate in [False, True]:
    sim.add_diagnostic(
        picmi.ReducedDiagnostic(
            diag_type="FieldProbe",
            name="integral" if integrate else "instant",
            period=1,
            probe_geometry="Line",
            x_probe=0.003,
            z_probe=-0.01,
            x1_probe=0.013,
            z1_probe=0.01,
            resolution=9,
            interp_order=2,
            integrate=integrate,
        )
    )
sim.add_diagnostic(picmi.Checkpoint(name="chk", period=2, write_dir="diags"))
sim.initialize_inputs()
warpx.write_diagnostics_on_restart = True
warpx.B_ext_grid_init_style = "constant"
warpx.B_external_grid = [0, 0, magnetic_field]
if args.nonzero_flux:
    warpx.E_ext_grid_init_style = "constant"
    warpx.E_external_grid = [0, 1000, 0]
sim.initialize_warpx()
start = sim.extension.warpx.getistep(lev=0)
sim.step(4 - start)
sim.finalize()
if MPI.COMM_WORLD.rank == 0:
    directory = Path("diags/reducedfiles")
    data = {
        path.stem: np.loadtxt(path, skiprows=1, ndmin=2)
        for path in directory.glob("*.txt")
    }
    np.savez("diagnostics.npz", **data)
    # A line probe writes one row per probe and step. Retain every probe,
    # including the initialization/restart sample, when checking exact fields.
    # The separate analysis also verifies the pre-checkpoint output prefix.
    data = {name: values[values[:, 0] >= start] for name, values in data.items()}
    if not args.record_only and not args.nonzero_flux:
        # A uniform axial magnetic field is a stationary exact Maxwell solution.
        energy = (
            magnetic_field**2 * np.pi * rmax**2 * length / (2 * picmi.constants.mu0)
        )
        np.testing.assert_allclose(data["FieldEnergy"][:, 4], energy, rtol=2e-13)
        np.testing.assert_allclose(data["FieldEnergy"][:, 3], 0, atol=1e-30)
        np.testing.assert_allclose(data["FieldPoyntingFlux"][:, 2:], 0, atol=1e-30)
        for name in ["instant", "integral"]:
            # Each probe record contains x,y,z,Ex,Ey,Ez,Bx,By,Bz,S (no id).
            values = data[name]
            expected = magnetic_field * (values[:, 1] if name == "integral" else 1)
            for column in range(10, values.shape[1], 10):
                np.testing.assert_allclose(
                    values[:, column], expected, rtol=2e-13, atol=1e-30
                )
        for values in data.values():
            assert np.isfinite(values).all()
            np.testing.assert_allclose(values[:, 1], values[:, 0] * dt, rtol=2e-13)
        print("PASS: constant-field diagnostics and cumulative probe timestamps")
    if args.nonzero_flux and not args.record_only:
        values = data["FieldPoyntingFlux"]
        if not args.restart:
            # Exact initial surface integral; subsequent fields need not be stationary.
            expected = (
                2 * np.pi * rmax * length * 1000 * magnetic_field / picmi.constants.mu0
            )
            np.testing.assert_allclose(values[0, 4], expected, rtol=2e-13)
            np.testing.assert_array_equal(values[0, 6:], 0)
        assert np.isfinite(values).all()
