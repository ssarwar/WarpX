#!/usr/bin/env python3
"""Compare fluid ion increments with the same frozen kinetic charge footprints."""

import argparse
from pathlib import Path

import numpy as np

from pywarpx import picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument(
    "--kind", choices=["proton", "ionization", "attachment"], required=True
)
parser.add_argument("--shape", type=int, default=3)
parser.add_argument("--solver", choices=["Yee", "PSATD"], default="Yee")
parser.add_argument("--walls", action="store_true")
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
    grid=grid, method=args.solver, stencil_order=[8, 8]
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
mi = 28.0134 * 1.66053906660e-27 - me
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
            background_mass=mi + me,
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
sim.initialize_inputs()
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


initial = rho("electrons")
previous = np.zeros_like(initial)
for step in range(3):
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
print(
    "PASS: immobile ion footprints match frozen kinetic events, including grid interfaces"
)
sim.finalize()
