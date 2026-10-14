#!/usr/bin/env python3
"""Observe checkpointed fields without splitting WarpX's evolution call."""

import argparse
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import callbacks, geometry, libwarpx

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input")
parser.add_argument("--parameter", action="append", default=[])
args = parser.parse_args()
geometry.dims = 3
seen = set()


def sample():
    simulation = libwarpx.warpx
    step = simulation.getistep(lev=0)
    if step not in (5, 6, 10) or step in seen:
        return
    fields = simulation.multifab_register()
    arrays = {"time": np.array(simulation.gett_new(lev=0))}
    for kind in ("Efield_fp", "Bfield_fp", "Efield_avg_fp", "Bfield_avg_fp"):
        for direction in "xyz":
            if fields.has(kind, direction, 0):
                values = fields.get(kind, direction, level=0)[...]
                arrays[kind + "_" + direction] = (
                    values.get()
                    if hasattr(values, "get")
                    else np.asarray(values).copy()
                )
    if MPI.COMM_WORLD.rank == 0:
        np.savez(Path(f"native_{step}.npz"), **arrays)
    seen.add(step)


callbacks.installafterdiagnostics(sample)
libwarpx.initialize(["probe_restart_fields", args.input, *args.parameter])
sample()  # A restart is observed before advancing its first resumed step.
libwarpx.warpx.evolve(-1)
libwarpx.finalize()
