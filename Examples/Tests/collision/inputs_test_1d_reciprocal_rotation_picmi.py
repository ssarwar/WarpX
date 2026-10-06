# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Particle-level reciprocal rotation checks and complete MCC timing workloads.

All physical tables are prepared offline. Initialization, collision steps and
particle measurements are timed separately; results are MPI-reduced.
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from pywarpx import amrex, picmi

try:
    import mpi4py

    mpi4py.rc.initialize = False
    mpi4py.rc.finalize = False
    from mpi4py import MPI
except ImportError:
    MPI = None

p = argparse.ArgumentParser()
p.add_argument("--data-dir", type=Path, required=True)
p.add_argument("--source-dir", type=Path, required=True)
p.add_argument("--target", choices=["N2", "O2", "air"], default="air")
p.add_argument("--temperature", type=float, default=300)
p.add_argument("--translational-temperature", type=float)
p.add_argument("--electron-temperature", type=float, default=300)
p.add_argument("--mode", choices=["thermal", "mono", "broad"], default="thermal")
p.add_argument("--energy", type=float, default=2.47)
p.add_argument("--particles", type=int, default=32768)
p.add_argument("--cells", type=int, default=128)
p.add_argument("--steps", type=int, default=128)
p.add_argument("--subcycles", type=int, default=1)
p.add_argument("--warmup", type=int, default=0)
p.add_argument("--dt", type=float, default=1e-14)
p.add_argument("--density", type=float, default=2.5e25)
p.add_argument("--seed", type=int, default=2026)
p.add_argument("--cumulative", action="store_true")
p.add_argument("--disabled", action="store_true")
p.add_argument(
    "--pic",
    action="store_true",
    help="Include gather and deposition in timestep timings",
)
p.add_argument("--checkpoint", type=Path)
p.add_argument("--restart", type=Path)
p.add_argument("--check-thermal", action="store_true")
p.add_argument("--output", type=Path, default=Path("reciprocal-results.json"))
p.add_argument(
    "--write-input",
    type=Path,
    help="Write native inputs for the complete MCC benchmark",
)
args = p.parse_args()
if args.particles % args.cells:
    raise ValueError("Particle count must be divisible by the cell count")
grid = picmi.Cartesian1DGrid(
    number_of_cells=[args.cells],
    lower_bound=[0],
    upper_bound=[10000],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=32,
    warpx_blocking_factor=1,
)
c, me, qe, kb = (
    picmi.constants.c,
    picmi.constants.m_e,
    picmi.constants.q_e,
    picmi.constants.kb,
)
rest = me * c * c / qe
distribution = dict(density=1)
if args.mode == "thermal":
    distribution["rms_velocity"] = [np.sqrt(kb * args.electron_temperature / me)] * 3
else:
    distribution["directed_velocity"] = [
        0,
        0,
        c * np.sqrt(args.energy * (args.energy + 2 * rest)) / rest,
    ]
electron = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(**distribution),
    warpx_do_not_deposit=not args.pic,
    warpx_do_not_gather=not args.pic,
)
collisions = []
targets = {"N2": 0.78084, "O2": 0.20946} if args.target == "air" else {args.target: 1.0}
for gas, fraction in targets.items():
    source = (args.source_dir / gas / "IAA").resolve()
    bundle = (
        args.data_dir
        / gas
        / "IAA"
        / f"reciprocal_hybrid_{args.temperature:g}K"
        / "thermal_rotation.rot"
    ).resolve()
    process = dict(
        cross_section=str(source / "elastic.txt"), scattering_angle_model="IAA"
    )
    if args.disabled:
        process["differential_cross_section"] = str(source / "elastic_dcs.txt")
    else:
        process.update(
            rotation_model="reciprocal_hybrid",
            rotation_file=str(bundle),
            rotational_temperature=args.temperature,
            rotation_sampling="cumulative" if args.cumulative else "alias",
        )
    collisions.append(
        picmi.MCCCollisions(
            name="mcc_" + gas,
            species=electron,
            background_density=args.density * fraction,
            background_temperature=(
                args.temperature
                if args.translational_temperature is None
                else args.translational_temperature
            ),
            scattering_processes={"elastic": process},
            ndt_subcycle=args.subcycles,
        )
    )
simulation = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9),
    time_step_size=args.dt,
    warpx_collisions=collisions,
    warpx_random_seed=args.seed,
    verbose=0,
    warpx_amr_restart=str(args.restart) if args.restart else None,
)
simulation.add_species(
    electron,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[args.particles // args.cells], grid=grid
    ),
)
if args.checkpoint:
    simulation.add_diagnostic(
        picmi.Checkpoint(name="chk", period=args.steps, write_dir=str(args.checkpoint))
    )
simulation.initialize_inputs()
amrex.the_arena_init_size = 8 * 1024 * 1024
if args.write_input:
    simulation.write_input_file(str(args.write_input))
    sys.exit(0)
start = time.perf_counter()
simulation.initialize_warpx()
startup = time.perf_counter() - start
container = simulation.particles.get("electrons")
comm = MPI.COMM_WORLD if MPI is not None and MPI.Is_initialized() else None
rank = comm.rank if comm is not None else 0
if args.mode == "broad":
    generator = np.random.default_rng(args.seed + rank)
    for tile in container.iterator(level=0):
        kinetic = np.exp(generator.uniform(np.log(0.002), np.log(3e6), len(tile["uz"])))
        stored = (c * np.sqrt(kinetic * (kinetic + 2 * rest)) / rest).astype(
            tile["uz"].dtype
        )
        if hasattr(tile["uz"], "set"):
            tile["uz"].set(stored)
        else:
            tile["uz"][:] = stored


def statistics():
    result = np.zeros(8)
    pairs = []
    for tile in container.iterator(level=0):
        u = [
            np.asarray(tile[k].get() if hasattr(tile[k], "get") else tile[k], float)
            for k in ["ux", "uy", "uz"]
        ]
        u2 = sum(v * v for v in u)
        energy = me * u2 / (qe * (1 + np.sqrt(1 + u2 / c**2)))
        assert np.isfinite(energy).all() and (energy >= 0).all()
        if args.check_thermal:
            ids = tile["idcpu"]
            ids = np.asarray(ids.get() if hasattr(ids, "get") else ids, dtype=np.uint64)
            pairs.append((ids.copy(), energy.copy()))
        result += [
            len(energy),
            energy.sum(),
            (energy**2).sum(),
            (energy**3).sum(),
            u[2].sum(),
            (u[2] ** 2).sum(),
            np.count_nonzero(energy > 1),
            np.count_nonzero(energy > 1000),
        ]
    if comm is not None:
        total = np.zeros_like(result)
        comm.Allreduce(result, total)
        result = total
    if args.check_thermal:
        if comm is not None:
            pairs = [
                pair for rank_pairs in comm.allgather(pairs) for pair in rank_pairs
            ]
        ids = np.concatenate([pair[0] for pair in pairs])
        energy = np.concatenate([pair[1] for pair in pairs])
        order = np.argsort(ids)
        return result, (ids[order], energy[order])
    return result, None


if args.warmup:
    simulation.step(args.warmup)
initial, initial_pairs = statistics()
if comm is not None:
    comm.Barrier()
start = time.perf_counter()
simulation.step(args.steps)
container.sum_particle_weight(local=True)
elapsed = time.perf_counter() - start
if comm is not None:
    elapsed = comm.allreduce(elapsed, op=MPI.MAX)
final, final_pairs = statistics()
assert final[0] == initial[0] == args.particles
mean_initial, mean_final = initial[1] / initial[0], final[1] / final[0]
error = np.sqrt(
    max(0, initial[2] / initial[0] - mean_initial**2) / initial[0]
    + max(0, final[2] / final[0] - mean_final**2) / final[0]
)
if args.check_thermal:
    assert np.array_equal(initial_pairs[0], final_pairs[0])
    delta = final_pairs[1] - initial_pairs[1]
    error = delta.std(ddof=1) / np.sqrt(len(delta))
    if args.electron_temperature < args.temperature:
        assert mean_final - mean_initial > 3 * error, (mean_initial, mean_final, error)
    elif args.electron_temperature > args.temperature:
        assert mean_initial - mean_final > 3 * error, (mean_initial, mean_final, error)
    else:
        assert abs(mean_final - mean_initial) <= 6 * error, (
            mean_initial,
            mean_final,
            error,
        )
report = dict(
    target=args.target,
    mode=args.mode,
    temperature=args.temperature,
    electron_temperature=args.electron_temperature,
    particles=args.particles,
    steps=args.steps,
    dt=args.dt,
    startup_seconds=startup,
    step_seconds=elapsed,
    mean_initial_eV=mean_initial,
    mean_final_eV=mean_final,
    mean_change_eV=mean_final - mean_initial,
    standard_error_bound_eV=error,
    initial=initial.tolist(),
    final=final.tolist(),
    ranks=comm.size if comm else 1,
    peak_host_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    * (1 if sys.platform == "darwin" else 1024),
)
if rank == 0:
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)
simulation.finalize()
