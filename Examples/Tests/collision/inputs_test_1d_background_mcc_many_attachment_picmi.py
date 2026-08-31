#!/usr/bin/env python3
"""Performance smoke test for many attachment channels sharing one product."""

import argparse
import math
import time
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--particle-count", type=int, default=32768)
parser.add_argument("--process-count", type=int, default=64)
parser.add_argument("--steps", type=int, default=1)
parser.add_argument("--subcycles", type=int, default=1)
parser.add_argument("--cell-count", type=int, default=1)
parser.add_argument("--max-grid-size", type=int)
parser.add_argument("--target-acceptance", type=float)
parser.add_argument("--mpi-timing", action="store_true")
args = parser.parse_args()

if args.mpi_timing:
    import mpi4py

    # AMReX owns MPI finalization after pywarpx initializes the simulation.
    mpi4py.rc.finalize = False
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    wall_time = MPI.Wtime
else:
    comm = None
    wall_time = time.perf_counter

PARTICLE_COUNT = args.particle_count
PROCESS_COUNT = args.process_count
STEPS = args.steps
COLLISION_SUBCYCLES = args.subcycles
GRID_CELL_COUNT = args.cell_count
MAX_GRID_SIZE = (
    args.max_grid_size if args.max_grid_size is not None else GRID_CELL_COUNT
)
TARGET_ACCEPTANCE = args.target_acceptance
INITIAL_NEGATIVE_COUNT = GRID_CELL_COUNT
BACKGROUND_DENSITY = 1.0e20
BACKGROUND_MASS = 1000.0 * picmi.constants.m_e
CROSS_SECTION_PER_PROCESS = 1.0e-24
ELECTRON_ENERGY_EV = 1000.0
MAJORANT_OPTICAL_DEPTH = 0.2
MAJORANT_MARGIN = 1.001

assert PARTICLE_COUNT > 0
assert PROCESS_COUNT > 0
assert STEPS > 0
assert COLLISION_SUBCYCLES > 0
assert GRID_CELL_COUNT > 0
assert MAX_GRID_SIZE > 0
assert PARTICLE_COUNT % GRID_CELL_COUNT == 0
if TARGET_ACCEPTANCE is not None:
    assert 0.0 < TARGET_ACCEPTANCE < 1.0 / MAJORANT_MARGIN

PARTICLES_PER_CELL = PARTICLE_COUNT // GRID_CELL_COUNT

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e

GAMMA = 1.0 + ELECTRON_ENERGY_EV * Q_E / (M_E * C**2)
ELECTRON_SPEED = C * math.sqrt(1.0 - 1.0 / GAMMA**2)
RATE_UPPER_BOUND = BACKGROUND_DENSITY * PROCESS_COUNT * CROSS_SECTION_PER_PROCESS * C
if TARGET_ACCEPTANCE is None:
    NU_MAX = None
    DT = MAJORANT_OPTICAL_DEPTH / RATE_UPPER_BOUND
else:
    physical_rate = (
        BACKGROUND_DENSITY
        * PROCESS_COUNT
        * CROSS_SECTION_PER_PROCESS
        * ELECTRON_SPEED
    )
    NU_MAX = MAJORANT_MARGIN * physical_rate
    candidate_probability = MAJORANT_MARGIN * TARGET_ACCEPTANCE
    DT = -math.log1p(-candidate_probability) / NU_MAX

cross_section_path = Path("background_mcc_many_attachment.txt").resolve()
energies = np.geomspace(1.0e-6, 2.5e6, 1025)
sigmas = np.full_like(energies, CROSS_SECTION_PER_PROCESS)
sigmas[-1] = 0.0
if comm is None or comm.rank == 0:
    np.savetxt(cross_section_path, np.column_stack((energies, sigmas)))
if comm is not None:
    comm.Barrier()
cross_section_file = str(cross_section_path)

grid = picmi.Cartesian1DGrid(
    number_of_cells=[GRID_CELL_COUNT],
    lower_bound=[0.0],
    upper_bound=[100.0 * GRID_CELL_COUNT],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=MAX_GRID_SIZE,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)

electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, ELECTRON_SPEED],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)
negative_ions = picmi.Species(
    name="negative_ions",
    charge="-q_e",
    mass=BACKGROUND_MASS + M_E,
    initial_distribution=picmi.UniformDistribution(
        density=1.0,
        directed_velocity=[0.0, 0.0, 0.0],
    ),
    warpx_do_not_deposit=True,
    warpx_do_not_gather=True,
)

processes = {
    f"attachment{index}": {
        "cross_section": cross_section_file,
        "cross_section_units": "m2",
        "species": negative_ions,
    }
    for index in range(PROCESS_COUNT)
}
collision_options = {}
if NU_MAX is not None:
    collision_options["nu_max"] = NU_MAX
collision = picmi.MCCCollisions(
    name="mcc",
    species=electrons,
    background_density=BACKGROUND_DENSITY,
    background_temperature=0.0,
    background_mass=BACKGROUND_MASS,
    scattering_processes=processes,
    ndt_subcycle=COLLISION_SUBCYCLES,
    **collision_options,
)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=DT,
    max_steps=STEPS,
    warpx_collisions=[collision],
    verbose=1,
)
sim.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[PARTICLES_PER_CELL], grid=grid
    ),
)
sim.add_species(
    negative_ions,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[INITIAL_NEGATIVE_COUNT], grid=grid
    ),
)

sim.initialize_inputs()
if comm is not None:
    comm.Barrier()
initialization_start = wall_time()
sim.initialize_warpx()
if comm is not None:
    assert libwarpx.amr.ParallelDescriptor.NProcs() == comm.size, (
        "--mpi-timing requires a WarpX build with MPI enabled"
    )
initialization_elapsed = wall_time() - initialization_start
if comm is not None:
    initialization_elapsed = comm.allreduce(initialization_elapsed, op=MPI.MAX)
    comm.Barrier()
start = wall_time()
sim.step(STEPS)
elapsed = wall_time() - start
if comm is not None:
    elapsed = comm.allreduce(elapsed, op=MPI.MAX)

electron_count = sim.particles.get("electrons").number_of_particles(only_local=False)
negative_count = sim.particles.get("negative_ions").number_of_particles(only_local=False)
attachment_events = negative_count - INITIAL_NEGATIVE_COUNT
nominal_particles_per_second = PARTICLE_COUNT * STEPS * COLLISION_SUBCYCLES / elapsed

if (comm is not None and comm.rank == 0) or (
    comm is None and libwarpx.amr.ParallelDescriptor.MyProc() == 0
):
    np.savez(
        "background_mcc_many_attachment_results.npz",
        electron_count=electron_count,
        attachment_events=attachment_events,
        particle_count=PARTICLE_COUNT,
        process_count=PROCESS_COUNT,
        cell_count=GRID_CELL_COUNT,
        max_grid_size=MAX_GRID_SIZE,
        steps=STEPS,
        subcycles=COLLISION_SUBCYCLES,
        target_acceptance=(
            TARGET_ACCEPTANCE if TARGET_ACCEPTANCE is not None else np.nan
        ),
        initialization_elapsed=initialization_elapsed,
        elapsed=elapsed,
        nominal_particles_per_second=nominal_particles_per_second,
        product_passes_per_collision_call=1,
    )
    print(
        "MCC benchmark: "
        f"processes={PROCESS_COUNT}, particles={PARTICLE_COUNT}, steps={STEPS}, "
        f"subcycles={COLLISION_SUBCYCLES}, cells={GRID_CELL_COUNT}, "
        f"max_grid_size={MAX_GRID_SIZE}, accepted={attachment_events}, "
        f"acceptance={attachment_events / PARTICLE_COUNT:.6f}, "
        f"initialization={initialization_elapsed:.6f} s, elapsed={elapsed:.6f} s, "
        f"nominal_throughput={nominal_particles_per_second:.3e} particle-calls/s, "
        "product_passes_per_call=1"
    )

# mpi4py owns MPI initialization in distributed benchmark mode. Let Python
# shutdown finalize MPI and AMReX once, as in other MPI-enabled PICMI inputs.
if comm is None:
    sim.finalize()
