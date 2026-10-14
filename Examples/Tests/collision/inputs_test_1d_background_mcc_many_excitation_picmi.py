#!/usr/bin/env python3
"""Physics and performance test for many excitation channels sharing one DCS."""

import argparse
import math
import shutil
import time
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--particle-count", type=int, default=65536)
parser.add_argument("--process-count", type=int, default=64)
parser.add_argument("--unique-dcs-files", action="store_true")
args = parser.parse_args()

PARTICLE_COUNT = args.particle_count
PROCESS_COUNT = args.process_count
BACKGROUND_DENSITY = 1.0e20
BACKGROUND_MASS = 28.0134 * 1.66053906660e-27
INCIDENT_ENERGY = 1000.0
CROSS_SECTION_PER_PROCESS = 1.0e-22
FIRST_ENERGY_LOSS = 5.0
ENERGY_LOSS_STEP = 0.1
OPTICAL_DEPTH = 0.5
ANISOTROPY = 0.6

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e
MC2_EV = M_E * C**2 / Q_E

assert PARTICLE_COUNT > 0
assert PROCESS_COUNT > 0


def electron_gamma(energy_ev):
    return 1.0 + energy_ev / MC2_EV


def electron_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(1.0 - 1.0 / gamma**2)


def electron_proper_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(gamma**2 - 1.0)


collision_speed = electron_speed(INCIDENT_ENERGY)
nu_max = (
    BACKGROUND_DENSITY * PROCESS_COUNT * CROSS_SECTION_PER_PROCESS * collision_speed
)

dcs_path = Path("background_mcc_many_excitation_dcs.txt").resolve()
dcs_energies = np.geomspace(1.0, 1.0e9, 65)
theta = np.linspace(0.0, math.pi, 361)
dcs = 1.0 + ANISOTROPY * np.cos(theta)
np.savetxt(
    dcs_path,
    np.column_stack((dcs_energies, np.tile(dcs, (dcs_energies.size, 1)))),
    header=(
        "Synthetic shared DCS with the IAA/elmolcs row layout\n"
        "COLUMNS: theta = linspace(0, 180, 361) (deg)"
    ),
)

processes = {}
energy_losses = FIRST_ENERGY_LOSS + ENERGY_LOSS_STEP * np.arange(PROCESS_COUNT)
for process_index, energy_loss in enumerate(energy_losses):
    cross_section_path = Path(
        f"background_mcc_many_excitation_{process_index}.txt"
    ).resolve()
    np.savetxt(
        cross_section_path,
        np.array(
            [
                [0.0, 0.0],
                [energy_loss, 0.0],
                [energy_loss + 1.0e-4, CROSS_SECTION_PER_PROCESS],
                [1.0e9, CROSS_SECTION_PER_PROCESS],
            ]
        ),
    )
    process_dcs_path = dcs_path
    if args.unique_dcs_files:
        process_dcs_path = Path(
            f"background_mcc_many_excitation_dcs_{process_index}.txt"
        ).resolve()
        shutil.copyfile(dcs_path, process_dcs_path)
    processes[f"excitation{process_index}"] = {
        "cross_section": str(cross_section_path),
        "energy": energy_loss,
        "scattering_angle_model": "IAA",
        "differential_cross_section": str(process_dcs_path),
    }

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
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, electron_proper_speed(INCIDENT_ENERGY)],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
collision = picmi.MCCCollisions(
    name="mcc",
    species=electrons,
    background_density=BACKGROUND_DENSITY,
    background_temperature=0.0,
    background_mass=BACKGROUND_MASS,
    scattering_processes=processes,
    nu_max=nu_max,
)
sim = picmi.Simulation(
    solver=solver,
    time_step_size=OPTICAL_DEPTH / nu_max,
    max_steps=1,
    warpx_collisions=[collision],
    verbose=1,
)
sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid),
)


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def component(container, name):
    return np.concatenate(
        [to_numpy(pti[name]) for pti in container.iterator(level=0)]
    ).astype(np.float64)


sim.initialize_inputs()
initialization_start = time.perf_counter()
sim.initialize_warpx()
initialization_elapsed = time.perf_counter() - initialization_start
step_start = time.perf_counter()
sim.step(1)
step_elapsed = time.perf_counter() - step_start

container = sim.particles.get("electrons")
ux = component(container, "ux")
uy = component(container, "uy")
uz = component(container, "uz")
proper_speed_sq = ux**2 + uy**2 + uz**2
proper_speed = np.sqrt(proper_speed_sq)
gamma = np.sqrt(1.0 + proper_speed_sq / C**2)
outgoing_energy = M_E * proper_speed_sq / ((gamma + 1.0) * Q_E)
cosine = uz / proper_speed

incident_pc = math.sqrt(INCIDENT_ENERGY * (INCIDENT_ENERGY + 2.0 * MC2_EV))
outgoing_pc = np.sqrt(outgoing_energy * (outgoing_energy + 2.0 * MC2_EV))
recoil_pc_sq = (
    incident_pc**2 + outgoing_pc**2 - 2.0 * incident_pc * outgoing_pc * cosine
)
neutral_rest_energy = BACKGROUND_MASS * C**2 / Q_E
recoil_energy = recoil_pc_sq / (
    np.sqrt(neutral_rest_energy**2 + recoil_pc_sq) + neutral_rest_energy
)
inferred_loss = INCIDENT_ENERGY - outgoing_energy - recoil_energy
collided = inferred_loss > 0.5 * FIRST_ENERGY_LOSS

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez(
        "background_mcc_many_excitation_results.npz",
        process_count=PROCESS_COUNT,
        particle_count=PARTICLE_COUNT,
        optical_depth=OPTICAL_DEPTH,
        anisotropy=ANISOTROPY,
        first_energy_loss=FIRST_ENERGY_LOSS,
        energy_loss_step=ENERGY_LOSS_STEP,
        cosines=cosine[collided],
        inferred_losses=inferred_loss[collided],
        initialization_elapsed=initialization_elapsed,
        step_elapsed=step_elapsed,
        particle_calls_per_second=PARTICLE_COUNT / step_elapsed,
    )
    print(
        "Many-excitation MCC benchmark: "
        f"processes={PROCESS_COUNT}, particles={PARTICLE_COUNT}, "
        f"events={np.count_nonzero(collided)}, "
        f"initialization={initialization_elapsed:.6f} s, "
        f"step={step_elapsed:.6f} s, "
        f"throughput={PARTICLE_COUNT / step_elapsed:.3e} particle-calls/s"
    )

sim.finalize()
