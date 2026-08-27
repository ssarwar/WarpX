#!/usr/bin/env python3
"""Sample N2 and O2 RBEQ ionization energy sharing at three energies."""

import math
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

PARTICLE_COUNT = 32768
INITIAL_ION_COUNT = 1
BACKGROUND_DENSITY = 1.0e20
CROSS_SECTION = 1.0e-20
MAJORANT_OPTICAL_DEPTH = 5.0

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e
AMU = 1.66053906660e-27

CASES = [
    ("n2_near", "N2", 15.78, 15.58, 28.0134 * AMU),
    ("n2_mid", "N2", 100.0, 15.58, 28.0134 * AMU),
    ("n2_high", "N2", 1000.0, 15.58, 28.0134 * AMU),
    ("o2_near", "O2", 12.27, 12.07, 31.9988 * AMU),
    ("o2_mid", "O2", 100.0, 12.07, 31.9988 * AMU),
    ("o2_high", "O2", 1000.0, 12.07, 31.9988 * AMU),
]


def electron_gamma(energy_ev):
    return 1.0 + energy_ev * Q_E / (M_E * C**2)


def electron_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(1.0 - 1.0 / gamma**2)


def electron_proper_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(gamma**2 - 1.0)


reference_speed = electron_speed(100.0)
nu_max = BACKGROUND_DENSITY * CROSS_SECTION * reference_speed
dt = MAJORANT_OPTICAL_DEPTH / nu_max
source_dir = Path(__file__).resolve().parent

grid = picmi.Cartesian1DGrid(
    number_of_cells=[1],
    lower_bound=[0.0],
    upper_bound=[100.0],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)


def make_electrons(name, energy_ev):
    return picmi.Species(
        particle_type="electron",
        name=f"electrons_{name}",
        initial_distribution=picmi.UniformDistribution(
            density=1.0,
            # WarpX's constant momentum injector interprets this PICMI field as
            # proper velocity even though the portable PICMI name says velocity.
            directed_velocity=[0.0, 0.0, electron_proper_speed(energy_ev)],
        ),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


def make_ions(name, neutral_mass):
    return picmi.Species(
        name=f"ions_{name}",
        charge="q_e",
        mass=neutral_mass - M_E,
        initial_distribution=picmi.UniformDistribution(
            density=1.0,
            directed_velocity=[0.0, 0.0, 0.0],
        ),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )


def make_collision(name, electrons, ions, target, energy_ev, binding_energy, mass):
    process = {
        "ionization": {
            "cross_section": str(
                source_dir / f"background_mcc_rbeq_{target.lower()}.txt"
            ),
            "energy": binding_energy,
            "energy_sharing_model": "RBEQ",
            "rbeq_target": target,
            "species": ions,
        }
    }
    return picmi.MCCCollisions(
        name=f"mcc_{name}",
        species=electrons,
        background_density=BACKGROUND_DENSITY
        * reference_speed
        / electron_speed(energy_ev),
        background_temperature=0.0,
        background_mass=mass,
        scattering_processes=process,
        nu_max=nu_max,
    )


species = {}
collisions = []
for name, target, energy_ev, binding_energy, mass in CASES:
    electrons = make_electrons(name, energy_ev)
    ions = make_ions(name, mass)
    species[name] = (electrons, ions)
    collisions.append(
        make_collision(name, electrons, ions, target, energy_ev, binding_energy, mass)
    )

sim = picmi.Simulation(
    solver=solver,
    time_step_size=dt,
    max_steps=1,
    warpx_collisions=collisions,
    verbose=1,
)
for electrons, ions in species.values():
    sim.add_species(
        electrons,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid
        ),
    )
    sim.add_species(
        ions,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[INITIAL_ION_COUNT], grid=grid
        ),
    )


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def component(container, name):
    return np.concatenate([to_numpy(pti[name]) for pti in container.iterator(level=0)])


def particle_ids(container):
    return np.concatenate(
        [
            to_numpy(libwarpx.amr.unpack_ids(pti["idcpu"]))
            for pti in container.iterator(level=0)
        ]
    )


sim.initialize_inputs()
sim.initialize_warpx()
initial_ids = {
    name: particle_ids(sim.particles.get(f"electrons_{name}")) for name, *_ in CASES
}
sim.step(1)


def kinetic_energy_ev(ux, uy, uz):
    ux = np.asarray(ux, dtype=np.float64)
    uy = np.asarray(uy, dtype=np.float64)
    uz = np.asarray(uz, dtype=np.float64)
    proper_speed_sq = ux**2 + uy**2 + uz**2
    gamma = np.sqrt(1.0 + proper_speed_sq / C**2)
    return M_E * proper_speed_sq / ((gamma + 1.0) * Q_E)


results = {}
for name, _, energy_ev, _, _ in CASES:
    container = sim.particles.get(f"electrons_{name}")
    ids = particle_ids(container)
    energies = kinetic_energy_ev(
        component(container, "ux"),
        component(container, "uy"),
        component(container, "uz"),
    )
    is_secondary = ~np.isin(ids, initial_ids[name])
    secondary_energies = energies[is_secondary]
    event_count = secondary_energies.size
    sampled_binding_mean = (PARTICLE_COUNT * energy_ev - np.sum(energies)) / event_count
    results[f"{name}_event_count"] = event_count
    results[f"{name}_secondary_energies"] = secondary_energies
    results[f"{name}_binding_mean"] = sampled_binding_mean

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez("background_mcc_rbeq_results.npz", **results)

sim.finalize()
