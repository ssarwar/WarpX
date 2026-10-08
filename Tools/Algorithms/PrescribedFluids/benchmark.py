#!/usr/bin/env python3
"""One reproducible RZ particle/fluid comparison; see README.md for ensembles.

Timing excludes Python measurements and synchronizes the device and MPI ranks.
The optional MCC channels are synthetic rate fixtures for representation tests,
not a quantitative air chemistry dataset.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import amrex, picmi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beam", choices=["fluid", "quiet", "random"], default="fluid")
    parser.add_argument(
        "--ions", choices=["fluid", "frozen", "thermal"], default="fluid"
    )
    parser.add_argument(
        "--mode", choices=["fields", "source", "coupled", "push"], default="source"
    )
    parser.add_argument(
        "--solver",
        choices=["Yee", "PSATD", "semi_implicit_em", "semi_implicit_mm"],
        default="Yee",
    )
    parser.add_argument("--cells", type=int, nargs=2, default=[32, 128])
    parser.add_argument("--max-grid-size", type=int, default=64)
    parser.add_argument("--radial-sigmas", type=float, default=8.0)
    parser.add_argument("--longitudinal-sigmas", type=float, default=12.0)
    parser.add_argument("--ppc", type=int, default=16)
    parser.add_argument("--electron-ppc", type=int, default=4)
    parser.add_argument("--seed", type=int, default=41)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--dt", type=float, default=5e-13)
    parser.add_argument("--weight", type=float, default=100.0)
    parser.add_argument("--cap", type=int, default=8)
    parser.add_argument("--gas-density", type=float, default=1e23)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--shape", type=int, default=3)
    parser.add_argument(
        "--implicit-deposition", choices=["direct", "villasenor"], default="villasenor"
    )
    parser.add_argument("--source-resolution", type=int, default=8)
    parser.add_argument("--subcycles", type=int, default=1)
    parser.add_argument("--mcc", action="store_true")
    parser.add_argument(
        "--mcc-all",
        action="store_true",
        help="Include IAA elastic/excitation and three-body attachment fixtures",
    )
    parser.add_argument("--no-self-fields", action="store_true")
    parser.add_argument("--checkpoint", action="store_true")
    parser.add_argument("--checkpoint-period", type=int)
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--check-restored", type=Path)
    parser.add_argument("--snapshot-particles", action="store_true")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("result.json"))
    args = parser.parse_args()
    args.mcc = args.mcc or args.mcc_all
    if args.steps < 3 or args.ppc < 1 or args.electron_ppc < 1:
        parser.error("Use at least three steps and positive particle counts")
    nquiet = int(np.sqrt(args.ppc))
    if args.beam == "quiet" and nquiet * nquiet != args.ppc:
        parser.error("Quiet particle ppc must be a perfect square")
    if args.mcc and args.mode != "coupled":
        parser.error("MCC is a coupled-physics option")
    comm = MPI.COMM_WORLD
    root = Path(__file__).resolve().parents[3]
    qe, me, mp, c = (
        picmi.constants.q_e,
        picmi.constants.m_e,
        picmi.constants.m_p,
        picmi.constants.c,
    )
    amu = 1.66053906660e-27
    energy, sigma_r, sigma_t, current = 800e6, 0.002, 25e-12, 0.6
    gamma = 1 + energy * qe / (mp * c * c)
    speed = c * np.sqrt(1 - gamma**-2)
    sigma_z = speed * sigma_t
    rmax, zmax = args.radial_sigmas * sigma_r, args.longitudinal_sigmas * sigma_z
    n0 = current / (qe * speed * 2 * np.pi * sigma_r**2 * (-np.expm1(-32)))
    if args.steps * args.dt * speed + 8 * sigma_z >= zmax:
        parser.error(
            "Increase longitudinal extent before running through a beam boundary"
        )
    grid = picmi.CylindricalGrid(
        number_of_cells=args.cells,
        lower_bound=[0, -zmax],
        upper_bound=[rmax, zmax],
        lower_boundary_conditions=["none", "periodic"],
        upper_boundary_conditions=[
            "none" if args.solver == "PSATD" else "dirichlet",
            "periodic",
        ],
        lower_boundary_conditions_particles=["none", "periodic"],
        upper_boundary_conditions_particles=["absorbing", "periodic"],
        n_azimuthal_modes=1,
        warpx_max_grid_size=args.max_grid_size,
        warpx_blocking_factor=8,
    )
    if args.beam == "fluid":
        beam = picmi.FluidSpecies(
            name="beam",
            model="rigid_beam",
            particle_type="proton",
            kinetic_energy=energy,
            sigma_r=sigma_r,
            sigma_t=sigma_t,
            peak_current=current,
            pulse_times=[0],
            initialize_self_fields=not args.no_self_fields,
        )
    else:
        beam = picmi.Species(
            name="beam",
            particle_type="proton",
            warpx_do_not_gather=True,
            warpx_random_theta=args.beam == "random",
            initial_distribution=picmi.AnalyticDistribution(
                density_expression=(
                    f"{n0:.17g}*exp(-(x*x+y*y)/(2*{sigma_r:.17g}^2)-z*z/(2*{sigma_z:.17g}^2))"
                    f"*(x*x+y*y < {8 * sigma_r:.17g}^2)*(abs(z) < {8 * sigma_z:.17g})"
                ),
                directed_velocity=[0, 0, gamma * speed],
            ),
        )
    electrons = picmi.Species(
        name="electrons",
        particle_type="electron",
        warpx_do_not_push=args.mode == "source",
        warpx_do_not_gather=args.mode == "source",
        initial_distribution=picmi.UniformDistribution(density=1e10)
        if args.mode == "push"
        else None,
    )
    ions = []
    ion_definitions = [
        ("N2plus", 28.0134 * amu - me, 1),
        ("O2plus", 31.9988 * amu - me, 1),
        ("Ominus", 15.9994 * amu + me, -1),
    ]
    if args.mcc_all:
        ion_definitions.append(("O2minus", 31.9988 * amu + me, -1))
    for name, mass, sign in ion_definitions:
        if args.ions == "fluid":
            ion = picmi.FluidSpecies(
                name=name, model="immobile", mass=mass, charge=sign * qe
            )
        else:
            ion = picmi.Species(
                name=name,
                mass=mass,
                charge=sign * qe,
                warpx_do_not_push=args.ions == "frozen",
                warpx_do_not_gather=args.ions == "frozen",
            )
        ions.append(ion)
    collisions = []
    if args.mode in ["source", "coupled"]:
        for target, fraction, ion in [("N2", 0.79, ions[0]), ("O2", 0.21, ions[1])]:
            options = (
                dict(
                    source_sampling_points=args.source_resolution,
                    sampling_seed=args.seed,
                )
                if args.beam == "fluid"
                else {}
            )
            collisions.append(
                picmi.ProtonImpactIonizationCollisions(
                    name="p_" + target,
                    species=beam,
                    product_species=[electrons, ion],
                    ionization_target=target,
                    background_density=args.gas_density * fraction,
                    background_temperature=args.temperature,
                    fixed_product_weight=args.weight,
                    max_products_per_cell=args.cap,
                    ndt_subcycle=args.subcycles,
                    **options,
                )
            )
    if args.mcc:
        extra_processes = {"N2": {}, "O2": {}}
        if args.mcc_all:
            # Manufactured tables exercise every collision path together. Their
            # values are numerical fixtures, not an evaluated air dataset.
            fixture = Path("mcc-fixture").resolve()
            if comm.rank == 0:
                fixture.mkdir(exist_ok=True)
                np.savetxt(fixture / "elastic.txt", [[0, 2e-20], [1e9, 2e-20]])
                np.savetxt(
                    fixture / "excitation.txt",
                    [[0, 0], [6, 0], [8, 1e-20], [1e9, 1e-20]],
                )
                np.savetxt(fixture / "attachment-m5.txt", [[0, 1e-43], [1e9, 1e-43]])
                theta = np.linspace(0, np.pi, 361)
                for target in extra_processes:
                    np.savetxt(
                        fixture / f"DCS.e-{target}",
                        np.column_stack(
                            (
                                [1, 10, 100, 1000, 10000, 1e9],
                                np.tile(1 + 0.5 * np.cos(theta), (6, 1)),
                            )
                        ),
                        header=(
                            "Synthetic DCS with IAA/elmolcs row layout\n"
                            f"SPECIES: e / {target}\n"
                            "COLUMNS: theta = linspace(0, 180, 361) (deg)"
                        ),
                        comments="",
                    )
            comm.Barrier()
            for target, processes in extra_processes.items():
                for kind in ["elastic", "excitation"]:
                    processes[kind] = dict(
                        cross_section=str(fixture / f"{kind}.txt"),
                        scattering_angle_model="IAA",
                        differential_cross_section=str(fixture / f"DCS.e-{target}"),
                    )
                processes["excitation"]["energy"] = 6.0
            extra_processes["O2"]["attachment_three_body"] = dict(
                cross_section=str(fixture / "attachment-m5.txt"),
                cross_section_units="m5",
                third_body_density=args.gas_density,
                species=ions[3],
            )
        # Intentionally synthetic: isolate destination handling and finite-mass
        # kinematics with the same fixed rate functions in every representation.
        collisions.append(
            picmi.MCCCollisions(
                name="e_N2",
                species=electrons,
                background_density=args.gas_density * 0.79,
                background_temperature=args.temperature,
                background_mass=28.0134 * amu,
                scattering_processes={
                    **extra_processes["N2"],
                    "ionization": dict(
                        cross_section=str(
                            root
                            / "Examples/Tests/collision/background_mcc_immobile_ionization.txt"
                        ),
                        energy=15.58,
                        species=ions[0],
                        energy_sharing_model="RBEQ",
                        scattering_angle_model="IAA",
                        rbeq_target="N2",
                    ),
                },
                ndt_subcycle=args.subcycles,
            )
        )
        collisions.append(
            picmi.MCCCollisions(
                name="e_O2",
                species=electrons,
                background_density=args.gas_density * 0.21,
                background_temperature=args.temperature,
                background_mass=31.9988 * amu,
                scattering_processes={
                    **extra_processes["O2"],
                    "ionization": dict(
                        cross_section=str(
                            root / "Examples/Tests/collision/background_mcc_rbeq_o2.txt"
                        ),
                        energy=12.07,
                        species=ions[1],
                        energy_sharing_model="RBEQ",
                        scattering_angle_model="IAA",
                        rbeq_target="O2",
                    ),
                    "attachment": dict(
                        cross_section=str(
                            root
                            / "Examples/Tests/collision/background_mcc_immobile_attachment.txt"
                        ),
                        cross_section_units="m2",
                        species=ions[2],
                    ),
                },
                ndt_subcycle=args.subcycles,
            )
        )
    implicit = args.solver.startswith("semi_implicit")
    mm = args.solver == "semi_implicit_mm"
    sim = picmi.Simulation(
        solver=picmi.ElectromagneticSolver(
            grid=grid,
            method="PSATD" if args.solver == "PSATD" else "Yee",
            stencil_order=[16, 16],
        ),
        time_step_size=args.dt,
        max_steps=args.steps,
        particle_shape=args.shape,
        warpx_collisions=collisions,
        warpx_amr_restart=str(args.restart) if args.restart else None,
        warpx_random_seed=args.seed,
        warpx_use_filter=False,
        warpx_current_deposition_algo=args.implicit_deposition if implicit else None,
        warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
            nonlinear_solver=picmi.NewtonNonlinearSolver(
                relative_tolerance=1e-10,
                linear_solver=picmi.GMRESLinearSolver(relative_tolerance=1e-10),
                use_mass_matrices_jacobian=mm,
                use_mass_matrices_pc=mm,
                pc_type=picmi.JacobiPreconditioner() if mm else None,
            )
        )
        if implicit
        else None,
        warpx_do_device_synchronize=args.profile,
        warpx_self_fields_required_precision=1e-11,
        verbose=0,
    )
    if args.beam == "fluid":
        sim.add_fluid_species(beam)
    else:
        layout = (
            picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[nquiet, 1, nquiet])
            if args.beam == "quiet"
            else picmi.PseudoRandomLayout(grid=grid, n_macroparticles_per_cell=args.ppc)
        )
        sim.add_species(
            beam, layout=layout, initialize_self_field=not args.no_self_fields
        )
    sim.add_species(
        electrons,
        layout=picmi.GriddedLayout(
            grid=grid, n_macroparticle_per_cell=[args.electron_ppc, 1, 1]
        )
        if args.mode == "push"
        else None,
    )
    for ion in ions:
        if args.ions == "fluid":
            sim.add_fluid_species(ion)
        else:
            sim.add_species(ion, layout=None)
    if args.checkpoint:
        sim.add_diagnostic(
            picmi.Checkpoint(name="chk", period=args.checkpoint_period or args.steps)
        )
    sim.initialize_inputs()
    amrex.the_arena_init_size = 0
    sim.extension.load_library()
    backend = sim.extension.Config.gpu_backend
    if backend in ["CUDA", "HIP"]:
        import cupy

        synchronize = cupy.cuda.runtime.deviceSynchronize
    else:
        # SYCL builds synchronize their work at the end of a step when requested.
        def synchronize():
            pass

        if backend == "SYCL":
            from pywarpx import warpx

            warpx.do_device_synchronize = True

    def timed(function):
        synchronize()
        comm.Barrier()
        start = time.perf_counter()
        function()
        synchronize()
        return comm.allreduce(time.perf_counter() - start, op=MPI.MAX)

    initialization = timed(sim.initialize_warpx)

    def host(value):
        return value.get() if hasattr(value, "get") else np.asarray(value)

    def values(name, direction=None):
        return np.squeeze(host(sim.fields.get(name, direction, level=0)[...])).copy()

    def number_density(species):
        if species in sim.fluid_species:
            return values("fluid_density_" + species.name)
        return np.squeeze(
            host(
                sim.particles.get(species.name).get_charge_density(lev=0, local=False)[
                    ...
                ]
            )
        ).copy() / (qe if species.name == "beam" else species.charge or -qe)

    def physical_number(density):
        nr, nz = args.cells
        dr, dz = rmax / nr, 2 * zmax / nz
        nodal = density.shape[0] == nr + 1
        radial = (np.arange(density.shape[0]) + (0 if nodal else 0.5)) * dr
        volume = 2 * np.pi * radial * dr * dz
        if nodal:
            volume[0] = np.pi * dr * dr * dz / 3
        return float(np.sum(density[:, :nz] * volume[:, None]))

    histories, snapshots = [], {}

    def sample(step):
        densities = {
            species.name: number_density(species)
            for species in [beam, electrons, *ions]
        }
        field_values = {
            kind + direction: values(kind + "field_fp", direction)
            for kind in ["E", "B"]
            for direction in ["r", "theta", "z"]
        }
        counts = {
            species.name: sim.particles.get(species.name).number_of_particles()
            for species in sim.species
        }
        counts.update({species.name: 0 for species in sim.fluid_species})
        row = dict(
            step=step,
            time=sim.extension.warpx.gett_new(lev=0),
            physical={
                name: physical_number(array) for name, array in densities.items()
            },
            macroparticles=counts,
        )
        row["electron_energy_J"] = sim.particles.get("electrons").sum_particle_energy(
            False
        )
        # Frozen kinetic ions retain their recoil/thermal momenta. Their energy
        # is the independently measured mechanical energy omitted by a matching
        # immobile destination (the ions never undergo further collisions here).
        row["ion_energy_J"] = {
            ion.name: 0.0
            if ion in sim.fluid_species
            else sim.particles.get(ion.name).sum_particle_energy(False)
            for ion in ions
        }
        row["pending"] = {
            collision.name: float(
                values(collision.name + "_product_weight_remainder").sum()
            )
            for collision in collisions
            if isinstance(collision, picmi.ProtonImpactIonizationCollisions)
        }
        row["primary_source_budgets"] = {
            collision.name: values(collision.name + "_source_budget")
            .reshape(-1, 4)
            .sum(axis=0)
            .tolist()
            for collision in collisions
            if args.beam == "fluid"
            and isinstance(collision, picmi.ProtonImpactIonizationCollisions)
        }
        histories.append(row)
        for name, array in densities.items():
            snapshots[f"{step}_{name}"] = array
        for name, array in field_values.items():
            snapshots[f"{step}_{name}"] = array
        if args.snapshot_particles:
            for species in sim.species:
                parts = [
                    np.column_stack(
                        [
                            host(tile[key])
                            for key in ["r", "theta", "z", "w", "ux", "uy", "uz"]
                        ]
                    )
                    for tile in sim.particles.get(species.name).iterator(level=0)
                ]
                local = np.concatenate(parts) if parts else np.empty((0, 7))
                particles = np.concatenate(comm.allgather(local))
                snapshots[f"{step}_particles_{species.name}"] = particles[
                    np.lexsort(particles.T[::-1])
                ]
            for collision in collisions:
                if not isinstance(collision, picmi.ProtonImpactIonizationCollisions):
                    continue
                suffixes = ["product_weight_remainder"]
                if args.beam == "fluid":
                    suffixes += ["source_budget", "sampling_counter"]
                for suffix in suffixes:
                    name = collision.name + "_" + suffix
                    snapshots[f"{step}_{name}"] = values(name)

    start_step = sim.extension.warpx.getistep(lev=0)
    sample(start_step)
    if args.check_restored:
        assert args.restart and args.snapshot_particles
        reference = np.load(args.check_restored)
        for name, actual in snapshots.items():
            if "particles_" in name or name.endswith("sampling_counter"):
                np.testing.assert_array_equal(actual, reference[name], err_msg=name)
            else:
                scale = np.max(np.abs(reference[name]), initial=0)
                np.testing.assert_allclose(
                    actual,
                    reference[name],
                    rtol=3e-13,
                    atol=3e-13 * scale,
                    err_msg=name,
                )
        if comm.rank == 0:
            print(
                "PASS: all saved fields, particles and source state restored before stepping"
            )
    timings = []
    for step in range(start_step + 1, args.steps + 1):
        timings.append(timed(lambda: sim.step(1)))
        if step in [args.steps // 2, args.steps]:
            sample(step)
    # Weighted electron energy histogram, reduced without gathering particle records.
    edges = np.concatenate(([0.0], np.geomspace(1e-5, 1e8, 161), [np.inf]))
    histogram = np.zeros(edges.size - 1)
    for tile in sim.particles.get("electrons").iterator(level=0):
        u2 = sum(host(tile[axis]) ** 2 for axis in ["ux", "uy", "uz"])
        kinetic = me * u2 / (np.sqrt(1 + u2 / c**2) + 1) / qe
        histogram += np.histogram(kinetic, bins=edges, weights=host(tile["w"]))[0]
    histogram = comm.allreduce(histogram)
    snapshots["energy_edges_eV"], snapshots["electron_histogram"] = edges, histogram
    real_bytes = 8 if sim.extension.Config.precision == "DOUBLE" else 4
    particle_real_bytes = (
        8 if sim.extension.Config.precision_particles == "DOUBLE" else 4
    )
    field_bytes = {}
    for species in sim.fluid_species:
        mf = sim.fields.get("fluid_density_" + species.name, level=0)
        local = sum(mfi.fabbox().num_pts for mfi in mf) * mf.num_comp * real_bytes
        field_bytes[species.name] = comm.allreduce(local)
    particle_bytes = {}
    for species in sim.species:
        pc = sim.particles.get(species.name)
        particle_bytes[species.name] = pc.number_of_particles() * (
            particle_real_bytes * pc.num_real_comps + 4 * pc.num_int_comps + 8
        )
    # Reserved device memory includes the AMReX arenas and CUDA/HIP context;
    # logical payload sizes above isolate species storage from those pools.
    device_bytes = None
    if backend in ["CUDA", "HIP"]:
        free, total = cupy.cuda.runtime.memGetInfo()
        device_bytes = comm.allreduce(total - free)
    if comm.rank == 0:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        checkpoint_bytes = (
            sum(
                path.stat().st_size
                for path in Path(f"diags/chk{args.steps:06d}").rglob("*")
                if path.is_file()
            )
            if args.checkpoint
            else None
        )
        result = dict(
            parameters={
                key: str(value) if isinstance(value, Path) else value
                for key, value in vars(args).items()
            },
            revision=subprocess.check_output(
                ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
            ).strip(),
            warpx_version=sim.extension.__version__,
            benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            platform=platform.platform(),
            ranks=comm.size,
            backend=backend or "CPU",
            initialization_s=initialization,
            timestep_s=timings,
            ordinary_step_median_s=float(
                np.median(timings[2:-1] if args.checkpoint else timings[2:])
            ),
            checkpoint_step_s=timings[-1] if args.checkpoint else None,
            checkpoint_bytes=checkpoint_bytes,
            fluid_density_bytes=field_bytes,
            particle_payload_bytes=particle_bytes,
            device_reserved_bytes=device_bytes,
            history=histories,
        )
        args.output.write_text(json.dumps(result, indent=2) + "\n")
        np.savez_compressed(args.output.with_suffix(".npz"), **snapshots)
    sim.finalize()


if __name__ == "__main__":
    main()
