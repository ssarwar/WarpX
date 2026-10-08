# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Warmed, fixed-population MCC/PIC timing with both rotational gases resident.

Run the same command with PYTHONPATH selecting each build being compared.
Particle positions and momenta are restored outside each timed block. The bath
is uniform; this measures computational cost, not a self-consistent beam pulse.
"""

import argparse
import json
import resource
import time
from pathlib import Path

import mpi4py
import numpy as np

mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from pywarpx import amrex, picmi
from pywarpx.LoadThirdParty import load_cupy

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--data", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--spectrum",
    choices=["thermal", "resonance", "intermediate", "relativistic", "broad"],
    default="thermal",
)
parser.add_argument("--particles", type=int, default=1048576)
parser.add_argument("--cells", type=int, default=4096)
parser.add_argument("--max-grid-size", type=int, default=1024)
parser.add_argument("--subcycles", type=int, default=64)
parser.add_argument("--dt", type=float, default=6.4e-13)
parser.add_argument("--blocks", type=int, default=20)
parser.add_argument("--steps", type=int, default=8)
parser.add_argument("--warmup", type=int, default=4)
parser.add_argument("--seed", type=int, default=2026)
parser.add_argument("--cumulative", action="store_true")
parser.add_argument("--ordinary-elastic", action="store_true")
parser.add_argument("--parser-background", action="store_true")
parser.add_argument("--inactive-ionization", action="store_true")
parser.add_argument("--pic", action="store_true")
args = parser.parse_args()
assert args.particles > 0 and args.cells > 0 and args.particles % args.cells == 0
assert args.blocks > 0 and args.steps > 0 and args.subcycles > 0 and args.dt > 0
assert args.warmup >= 0

c, me, qe, kb = (
    picmi.constants.c,
    picmi.constants.m_e,
    picmi.constants.q_e,
    picmi.constants.kb,
)
rest = me * c * c / qe
grid = picmi.Cartesian1DGrid(
    number_of_cells=[args.cells],
    lower_bound=[0],
    upper_bound=[10000],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=args.max_grid_size,
    warpx_blocking_factor=1,
)
distribution = dict(density=1)
if args.spectrum == "thermal":
    distribution["rms_velocity"] = [np.sqrt(kb * 300 / me)] * 3
else:
    energy = {"resonance": 2.47, "intermediate": 50, "relativistic": 2.5e6}.get(
        args.spectrum, 2.47
    )
    distribution["directed_velocity"] = [
        0,
        0,
        c * np.sqrt(energy * (energy + 2 * rest)) / rest,
    ]
electrons = picmi.Species(
    particle_type="electron",
    name="electrons",
    initial_distribution=picmi.UniformDistribution(**distribution),
    warpx_do_not_push=not args.pic,
    warpx_do_not_deposit=not args.pic,
    warpx_do_not_gather=not args.pic,
)
collisions, ions = [], []
args.output.parent.mkdir(parents=True, exist_ok=True)
for gas, fraction, mass in [("N2", 0.78084, 28.0134), ("O2", 0.20946, 31.9988)]:
    source = (args.data / gas / "IAA").resolve()
    elastic = dict(
        cross_section=str(source / "elastic.txt"), scattering_angle_model="IAA"
    )
    if args.ordinary_elastic:
        elastic["differential_cross_section"] = str(source / "elastic_dcs.txt")
    else:
        elastic.update(
            rotation_model="reciprocal_hybrid",
            rotation_file=str(source / "reciprocal_hybrid_300K/thermal_rotation.rot"),
            rotational_temperature=300,
            rotation_sampling="cumulative" if args.cumulative else "alias",
        )
    processes = dict(elastic=elastic)
    if args.inactive_ionization:
        # Exercise product buffers at fixed particle count. This synthetic zero
        # channel is a performance control, not a physical air mixture.
        ion = picmi.Species(
            name="inactive_ions_" + gas,
            charge=qe,
            mass=mass * 1.66053906660e-27 - me,
            warpx_do_not_deposit=True,
            warpx_do_not_gather=True,
        )
        ions.append(ion)
        table = args.output.resolve().parent / ("inactive_ionization_" + gas + ".txt")
        table.write_text("0 0\n1e9 0\n")
        processes["ionization"] = dict(
            cross_section=str(table),
            energy=15.58 if gas == "N2" else 12.07,
            species=ion,
        )
    density = 2.5e25 * fraction
    collisions.append(
        picmi.MCCCollisions(
            name="mcc_" + gas,
            species=electrons,
            background_mass=mass * 1.66053906660e-27,
            background_density=(
                f"{density:.17g} + 0*x" if args.parser_background else density
            ),
            max_background_density=density if args.parser_background else None,
            background_temperature="300 + 0*x" if args.parser_background else 300,
            scattering_processes=processes,
            ndt_subcycle=args.subcycles,
        )
    )
simulation = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9),
    time_step_size=args.dt,
    warpx_collisions=collisions,
    warpx_random_seed=args.seed,
    verbose=0,
)
simulation.add_species(
    electrons,
    layout=picmi.GriddedLayout(
        n_macroparticle_per_cell=[args.particles // args.cells], grid=grid
    ),
)
for ion in ions:
    simulation.add_species(ion, layout=None)
simulation.initialize_inputs()
amrex.the_arena_init_size = 8 * 1024 * 1024
start = time.perf_counter()
simulation.initialize_warpx()
startup = time.perf_counter() - start
comm = MPI.COMM_WORLD
startup = comm.allreduce(startup, op=MPI.MAX)
container = simulation.particles.get("electrons")
xp, _ = load_cupy()
if args.spectrum == "broad":
    generator = np.random.default_rng(args.seed + comm.rank)
    for tile in container.iterator(level=0):
        energy = np.exp(generator.uniform(np.log(0.002), np.log(3e6), len(tile["uz"])))
        stored = (c * np.sqrt(energy * (energy + 2 * rest)) / rest).astype(
            tile["uz"].dtype
        )
        tile["uz"][:] = xp.asarray(stored)


def synchronize():
    if hasattr(xp, "cuda"):
        xp.cuda.runtime.deviceSynchronize()


saved = [
    {name: tile[name].copy() for name in ("z", "ux", "uy", "uz")}
    for tile in container.iterator(level=0)
]


def restore():
    count = 0
    for index, tile in enumerate(container.iterator(level=0)):
        for name, values in saved[index].items():
            assert tile[name].shape == values.shape
            tile[name][...] = values
        count += 1
    assert count == len(saved)
    synchronize()


if args.warmup:
    simulation.step(args.warmup)
times = []
for _ in range(args.blocks):
    restore()
    comm.Barrier()
    start = time.perf_counter()
    simulation.step(args.steps)
    synchronize()
    times.append(comm.allreduce(time.perf_counter() - start, op=MPI.MAX))

count = 0
for tile in container.iterator(level=0):
    count += len(tile["w"])
    for name in ("ux", "uy", "uz"):
        assert bool(xp.isfinite(tile[name]).all())
assert comm.allreduce(count, op=MPI.SUM) == args.particles
peak_rss = comm.allreduce(
    resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.MAX
)
if comm.rank == 0:
    args.output.write_text(
        json.dumps(
            dict(
                configuration=vars(args),
                ranks=comm.size,
                startup_seconds=startup,
                block_seconds=times,
                median_ms_per_pic_step=1000 * float(np.median(times)) / args.steps,
                peak_rss_native_units=peak_rss,
            ),
            default=str,
            indent=2,
        )
        + "\n"
    )
simulation.finalize()
