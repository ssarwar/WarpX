#!/usr/bin/env python3
"""Compare fluid ion increments with the same frozen kinetic charge footprints."""

import argparse
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import amr, picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument(
    "--kind", choices=["proton", "ionization", "attachment"], required=True
)
parser.add_argument("--shape", type=int, default=3)
parser.add_argument(
    "--solver",
    choices=["Yee", "PSATD", "semi_implicit_em", "semi_implicit_mm"],
    default="Yee",
)
parser.add_argument("--walls", action="store_true")
parser.add_argument("--restart")
parser.add_argument("--small-tiles", action="store_true")
args = parser.parse_args()

grid = picmi.CylindricalGrid(
    number_of_cells=[16, 32],
    lower_bound=[0.0, 0.0],
    upper_bound=[1.0, 2.0],
    lower_boundary_conditions=["none", "none" if args.walls else "periodic"],
    upper_boundary_conditions=["none", "none" if args.walls else "periodic"],
    lower_boundary_conditions_particles=[
        "none",
        "absorbing" if args.walls else "periodic",
    ],
    upper_boundary_conditions_particles=[
        "absorbing",
        "absorbing" if args.walls else "periodic",
    ],
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
qe, me, mp, c = (
    picmi.constants.q_e,
    picmi.constants.m_e,
    picmi.constants.m_p,
    picmi.constants.c,
)
neutral_mass = 28.0134 * 1.66053906660e-27
mi = neutral_mass + (me if args.kind == "attachment" else -me)
ion_fluid = picmi.Species(name="ion_fluid", charge=qe, mass=mi)
electrons = picmi.Species(name="electrons", particle_type="electron", **frozen)
species = [electrons]
collisions = []
if args.kind == "proton":
    gamma = 1 + 8e8 * qe / (mp * c * c)
    beam = picmi.Species(
        name="beam",
        particle_type="proton",
        **frozen,
        initial_distribution=picmi.UniformDistribution(
            density=1e8, directed_velocity=[0, 0, c * np.sqrt(gamma * gamma - 1)]
        ),
    )
    reference_e = picmi.Species(name="reference_e", particle_type="electron", **frozen)
    reference_i = picmi.Species(name="reference_i", charge=qe, mass=mi, **frozen)
    species += [beam, reference_e, reference_i]
    for name, electron, ion in [
        ("fluid", electrons, ion_fluid),
        ("kinetic", reference_e, reference_i),
    ]:
        collisions.append(
            picmi.ProtonImpactIonizationCollisions(
                name=name,
                species=beam,
                product_species=[electron, ion],
                ionization_target="N2",
                background_density=1e25,
                background_temperature=0,
                fixed_product_weight=1e5,
                max_products_per_cell=4,
            )
        )
else:
    energy = 1000.0
    gamma = 1 + energy * qe / (me * c * c)
    electrons.initial_distribution = picmi.UniformDistribution(
        density=1e8, directed_velocity=[0, 0, c * np.sqrt(gamma * gamma - 1)]
    )
    table = (
        Path(__file__).resolve().with_name(f"background_mcc_immobile_{args.kind}.txt")
    )
    if args.kind == "ionization":
        process = dict(
            cross_section=str(table),
            energy=15.58,
            species=ion_fluid,
            energy_sharing_model="RBEQ",
            scattering_angle_model="IAA",
            rbeq_target="N2",
        )
    else:
        process = dict(
            cross_section=str(table), cross_section_units="m2", species=ion_fluid
        )
    collisions.append(
        picmi.MCCCollisions(
            name="mcc",
            species=electrons,
            background_density=1e24,
            background_temperature=0,
            background_mass=neutral_mass,
            scattering_processes={args.kind: process},
        )
    )

sim = picmi.Simulation(
    solver=solver,
    time_step_size=1e-12,
    max_steps=3,
    particle_shape=args.shape,
    warpx_collisions=collisions,
    warpx_random_seed=47,
    warpx_current_deposition_algo="direct"
    if args.solver.startswith("semi_implicit")
    else None,
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.NewtonNonlinearSolver(
            relative_tolerance=1e-10,
            use_mass_matrices_jacobian=args.solver == "semi_implicit_mm",
        )
    )
    if args.solver.startswith("semi_implicit")
    else None,
    verbose=0,
)
for sp in species:
    sim.add_species(
        sp,
        layout=picmi.GriddedLayout(
            grid=grid,
            n_macroparticle_per_cell=[1, 1, 1] if args.kind == "proton" else [4, 1, 4],
        ),
    )
if args.kind == "attachment":
    sim.add_diagnostic(picmi.Checkpoint(name="chk", period=2))
sim.initialize_inputs()
if args.small_tiles:
    # Neighboring OpenMP tiles scatter into shared nodes of the same FAB.
    warpx.get_bucket("particles").do_tiling = True
    warpx.get_bucket("particles").tile_size = [4, 4]
if args.restart:
    amr.restart = args.restart
warpx.get_bucket("fluids").species_names = ["ion_fluid"]
fluid = warpx.get_bucket("ion_fluid")
fluid.model = "immobile"
fluid.charge = -qe if args.kind == "attachment" else qe
fluid.mass = mi
sim.initialize_warpx()


def host(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def rho(name):
    mf = sim.particles.get(name).get_charge_density(lev=0, local=True)
    periodicity = sim.extension.warpx.Geom(0).periodicity()
    mf.sum_boundary(0, 1, mf.n_grow_vect, mf.n_grow_vect, periodicity)
    return host(mf[0:3j, -2j:3j]).copy()


def state():
    local = []
    for tile in sim.particles.get("electrons").iterator(level=0):
        local.append(
            np.column_stack(
                [host(tile[key]) for key in ["r", "theta", "z", "w", "ux", "uy", "uz"]]
            )
        )
    local = np.concatenate(local) if local else np.empty((0, 7))
    particles = np.concatenate(MPI.COMM_WORLD.allgather(local))
    particles = particles[np.lexsort(particles.T[::-1])]
    density = host(
        sim.fields.get("fluid_density_ion_fluid", level=0)[0:3j, -2j:3j]
    ).copy()
    return dict(particles=particles, density=density)


start = sim.extension.warpx.getistep(lev=0)
if args.restart:
    assert args.kind == "attachment" and start == 2
    directory = Path(args.restart).parents[1]
    restored = state()
    with np.load(directory / "attachment_2.npz") as saved:
        for name, array in restored.items():
            np.testing.assert_allclose(array, saved[name], rtol=2e-13, atol=0)
    with np.load(directory / "attachment_initial.npz") as saved:
        initial, squared_weights = saved["rho"], float(saved["squared_weights"])
        initial_weight = float(saved["initial_weight"])
        initial_count = int(saved["initial_count"])
        attachment_rate = float(saved["attachment_rate"])
else:
    initial = rho("electrons")
    if args.kind == "attachment":
        initial_particles = state()["particles"]
        weights = initial_particles[:, 3]
        initial_weight = weights.sum()
        squared_weights = np.sum(weights**2)
        initial_count = len(weights)
        # Frozen monoenergetic electrons in a stationary, uniform gas obey
        # dN/dt = -n_gas*sigma*v*N exactly. Compute v from stored proper momentum
        # independently of the collision code and PICMI injection convention.
        proper_speed2 = np.sum(initial_particles[:, 4:] ** 2, axis=1)
        np.testing.assert_allclose(proper_speed2, proper_speed2[0], rtol=2e-15)
        speed = np.sqrt(proper_speed2[0] / (1 + proper_speed2[0] / c**2))
        attachment_rate = 1e24 * 1e-20 * speed
        if MPI.COMM_WORLD.rank == 0:
            np.savez_compressed(
                "attachment_initial.npz",
                rho=initial,
                squared_weights=squared_weights,
                initial_weight=initial_weight,
                initial_count=initial_count,
                attachment_rate=attachment_rate,
            )
previous = np.zeros_like(initial)
for step in range(start, 3):
    sim.step(1)
    if args.kind == "proton":
        expected = rho("reference_i")
    else:
        expected = initial - rho("electrons")
    actual = fluid.charge * host(
        sim.fields.get("fluid_density_ion_fluid", level=0)[0:3j, -2j:3j]
    )
    scale = np.max(np.abs(expected))
    assert scale > 0, "The fixture emitted no products"
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-14 * scale)
    assert np.max(np.abs(actual - previous)) > 0, "No new collision increment was added"
    previous = actual.copy()
    if args.kind == "attachment":
        saved_state = state()
        probability = np.exp(-attachment_rate * sim.extension.warpx.gett_new(lev=0))
        survivors = saved_state["particles"]
        # The exact Bernoulli variance is p*(1-p)*sum(w^2), including unequal
        # cylindrical particle weights. Six-sigma bounds are fixed beforehand
        # and allow stochastic CPU/GPU and checkpoint continuations.
        variance_factor = probability * (1 - probability)
        assert abs(survivors[:, 3].sum() - probability * initial_weight) < 6 * np.sqrt(
            variance_factor * squared_weights
        )
        assert abs(len(survivors) - probability * initial_count) < 6 * np.sqrt(
            variance_factor * initial_count
        )
        if not args.restart and MPI.COMM_WORLD.rank == 0:
            np.savez_compressed(f"attachment_{step + 1}.npz", **saved_state)
if args.restart:
    with np.load(directory / "attachment_3.npz") as reference:
        # The existing MCC RNG does not promise a replay after redistribution.
        # For weighted Bernoulli survivors Var(N) <= sum(w^2)/4; six standard
        # deviations of two independent continuations gives this conservative bound.
        difference = abs(
            saved_state["particles"][:, 3].sum() - reference["particles"][:, 3].sum()
        )
        assert difference <= 6 * np.sqrt(squared_weights / 2)
print(
    "PASS: immobile ion footprints match frozen kinetic events, including grid interfaces"
)
sim.finalize()
