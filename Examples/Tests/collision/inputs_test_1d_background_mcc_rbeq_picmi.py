#!/usr/bin/env python3
"""Sample N2 and O2 RBEQ ionization from threshold through 1 GeV."""

import math
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

PARTICLE_COUNT = 32768
INITIAL_ION_COUNT = 1
# A high test density preserves the desired optical depth in a very short
# physical step. Original particles and products then remain much closer than
# their initial spacing, which permits unambiguous event-by-event matching.
BACKGROUND_DENSITY = 1.0e29
CROSS_SECTION = 1.0e-20
MAJORANT_OPTICAL_DEPTH = 5.0
MAJORANT_MARGIN = 1.001

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e
AMU = 1.66053906660e-27

CASES = [
    ("n2_near", "N2", 15.78, 15.58, 28.0134 * AMU, "IAA"),
    ("n2_partial_guard", "N2", 19.50, 15.58, 28.0134 * AMU, "IAA"),
    ("n2_mid", "N2", 100.0, 15.58, 28.0134 * AMU, "IAA"),
    ("n2_high", "N2", 1000.0, 15.58, 28.0134 * AMU, "IAA"),
    ("n2_gev", "N2", 1.0e9, 15.58, 28.0134 * AMU, "IAA"),
    ("n2_forward_gev", "N2", 1.0e9, 15.58, 28.0134 * AMU, "forward"),
    ("n2_backward_gev", "N2", 1.0e9, 15.58, 28.0134 * AMU, "backward"),
    ("n2_isotropic_gev", "N2", 1.0e9, 15.58, 28.0134 * AMU, "isotropic"),
    ("o2_near", "O2", 12.27, 12.07, 31.9988 * AMU, "IAA"),
    ("o2_sdcs_guard", "O2", 18.702, 12.07, 31.9988 * AMU, "IAA"),
    ("o2_partial_guard", "O2", 21.50, 12.07, 31.9988 * AMU, "IAA"),
    ("o2_mid", "O2", 100.0, 12.07, 31.9988 * AMU, "IAA"),
    ("o2_high", "O2", 1000.0, 12.07, 31.9988 * AMU, "IAA"),
    ("o2_gev", "O2", 1.0e9, 12.07, 31.9988 * AMU, "IAA"),
]

PARTICLE_COUNTS = {
    "n2_forward_gev": 8192,
    "n2_backward_gev": 8192,
    "n2_isotropic_gev": 8192,
    "o2_sdcs_guard": 524288,
}


def electron_gamma(energy_ev):
    return 1.0 + energy_ev * Q_E / (M_E * C**2)


def electron_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(1.0 - 1.0 / gamma**2)


def electron_proper_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(gamma**2 - 1.0)


reference_speed = electron_speed(100.0)
nu_max = MAJORANT_MARGIN * BACKGROUND_DENSITY * CROSS_SECTION * reference_speed
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


def make_collision(
    name, electrons, ions, target, energy_ev, binding_energy, mass, angle_model
):
    process = {
        "ionization": {
            "cross_section": str(
                source_dir / f"background_mcc_rbeq_{target.lower()}.txt"
            ),
            "energy": binding_energy,
            "energy_sharing_model": "RBEQ",
            "rbeq_target": target,
            "scattering_angle_model": angle_model,
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
for name, target, energy_ev, binding_energy, mass, angle_model in CASES:
    electrons = make_electrons(name, energy_ev)
    ions = make_ions(name, mass)
    species[name] = (electrons, ions)
    collisions.append(
        make_collision(
            name,
            electrons,
            ions,
            target,
            energy_ev,
            binding_energy,
            mass,
            angle_model,
        )
    )

sim = picmi.Simulation(
    solver=solver,
    time_step_size=dt,
    max_steps=1,
    warpx_collisions=collisions,
    verbose=1,
)
for name, (electrons, ions) in species.items():
    sim.add_species(
        electrons,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[PARTICLE_COUNTS.get(name, PARTICLE_COUNT)],
            grid=grid,
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


def match_positions(reference_positions, target_positions):
    """Return unique nearest-reference indices for well-separated positions."""
    order = np.argsort(reference_positions)
    sorted_reference = reference_positions[order]
    right = np.clip(
        np.searchsorted(sorted_reference, target_positions),
        0,
        sorted_reference.size - 1,
    )
    left = np.maximum(right - 1, 0)
    use_left = np.abs(sorted_reference[left] - target_positions) < np.abs(
        sorted_reference[right] - target_positions
    )
    matches = np.where(use_left, left, right)
    assert np.unique(matches).size == target_positions.size
    assert np.max(np.abs(sorted_reference[matches] - target_positions)) < 1.0e-5
    return order[matches]


sim.initialize_inputs()
sim.initialize_warpx()
initial_ids = {
    name: particle_ids(sim.particles.get(f"electrons_{name}")) for name, *_ in CASES
}
initial_ion_ids = {
    name: particle_ids(sim.particles.get(f"ions_{name}")) for name, *_ in CASES
}
sim.step(1)


def kinetic_energy_ev(ux, uy, uz, mass):
    ux = np.asarray(ux, dtype=np.float64)
    uy = np.asarray(uy, dtype=np.float64)
    uz = np.asarray(uz, dtype=np.float64)
    proper_speed_sq = ux**2 + uy**2 + uz**2
    gamma = np.sqrt(1.0 + proper_speed_sq / C**2)
    return mass * proper_speed_sq / ((gamma + 1.0) * Q_E)


def cosine_to_z(ux, uy, uz):
    proper_speed = np.sqrt(ux**2 + uy**2 + uz**2)
    return np.divide(
        uz,
        proper_speed,
        out=np.zeros_like(proper_speed, dtype=np.float64),
        where=proper_speed > 0.0,
    )


results = {}
particle_real_bytes = None
for name, _, energy_ev, binding_energy, neutral_mass, _ in CASES:
    particle_count = PARTICLE_COUNTS.get(name, PARTICLE_COUNT)
    container = sim.particles.get(f"electrons_{name}")
    ids = particle_ids(container)
    electron_ux_raw = component(container, "ux")
    if particle_real_bytes is None:
        particle_real_bytes = electron_ux_raw.dtype.itemsize
    electron_ux = electron_ux_raw.astype(np.float64)
    electron_uy = component(container, "uy").astype(np.float64)
    electron_uz = component(container, "uz").astype(np.float64)
    electron_z = component(container, "z").astype(np.float64)
    energies = kinetic_energy_ev(electron_ux, electron_uy, electron_uz, M_E)
    is_original = np.isin(ids, initial_ids[name])
    is_secondary = ~is_original
    secondary_energies = energies[is_secondary]
    event_count = secondary_energies.size

    ion_container = sim.particles.get(f"ions_{name}")
    ion_ids = particle_ids(ion_container)
    ion_ux = component(ion_container, "ux").astype(np.float64)
    ion_uy = component(ion_container, "uy").astype(np.float64)
    ion_uz = component(ion_container, "uz").astype(np.float64)
    ion_z = component(ion_container, "z").astype(np.float64)
    ion_mass = neutral_mass - M_E
    ion_energies = kinetic_energy_ev(ion_ux, ion_uy, ion_uz, ion_mass)

    # Products retain the source-particle position. Match each secondary to its
    # original primary and product ion so every event independently verifies
    # the sampled discrete binding loss and three-body energy conservation.
    secondary_z = electron_z[is_secondary]
    original_indices = np.flatnonzero(is_original)
    primary_indices = original_indices[
        match_positions(electron_z[original_indices], secondary_z)
    ]

    is_product_ion = ~np.isin(ion_ids, initial_ion_ids[name])
    product_ion_indices = np.flatnonzero(is_product_ion)
    matched_ion_indices = product_ion_indices[
        match_positions(ion_z[product_ion_indices], secondary_z)
    ]

    binding_losses = (
        energy_ev
        - energies[primary_indices]
        - secondary_energies
        - ion_energies[matched_ion_indices]
    )
    sampled_binding_mean = (
        particle_count * energy_ev - np.sum(energies) - np.sum(ion_energies)
    ) / event_count

    final_momentum = np.array(
        [np.sum(electron_ux), np.sum(electron_uy), np.sum(electron_uz)]
    )
    final_momentum += (
        ion_mass / M_E * np.array([np.sum(ion_ux), np.sum(ion_uy), np.sum(ion_uz)])
    )
    initial_momentum = np.array(
        [0.0, 0.0, particle_count * electron_proper_speed(energy_ev)]
    )
    results[f"{name}_particle_count"] = particle_count
    results[f"{name}_event_count"] = event_count
    results[f"{name}_ion_count"] = ion_ux.size - INITIAL_ION_COUNT
    results[f"{name}_secondary_energies"] = secondary_energies
    results[f"{name}_binding_losses"] = binding_losses
    results[f"{name}_binding_mean"] = sampled_binding_mean
    results[f"{name}_primary_cosines"] = cosine_to_z(
        electron_ux[primary_indices],
        electron_uy[primary_indices],
        electron_uz[primary_indices],
    )
    results[f"{name}_secondary_cosines"] = cosine_to_z(
        electron_ux[is_secondary],
        electron_uy[is_secondary],
        electron_uz[is_secondary],
    )
    results[f"{name}_momentum_residual"] = (
        final_momentum - initial_momentum
    ) / np.linalg.norm(initial_momentum)

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    results["particle_real_bytes"] = particle_real_bytes
    np.savez("background_mcc_rbeq_results.npz", **results)

sim.finalize()
