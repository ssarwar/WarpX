#!/usr/bin/env python3
"""Record distinct rank/GPU placement and check a CUDA-aware MPI reduction."""

import json
import os
import socket

import cupy as cp
import cupyx  # noqa: F401 -- initialize lazy version data before NERSC PyMon exits
from mpi4py import MPI

comm = MPI.COMM_WORLD
device = cp.cuda.Device()
send = cp.array([comm.rank + 1], dtype=cp.float64)
received = cp.zeros_like(send)
cp.cuda.get_current_stream().synchronize()
comm.Allreduce(send, received, op=MPI.SUM)
assert float(received.get()[0]) == comm.size * (comm.size + 1) / 2
placement = comm.allgather(
    dict(
        rank=comm.rank,
        host=socket.gethostname(),
        visible_devices=os.environ.get("CUDA_VISIBLE_DEVICES"),
        pci_bus_id=device.pci_bus_id,
    )
)
assert len({(item["host"], item["pci_bus_id"]) for item in placement}) == comm.size
if comm.rank == 0:
    print(json.dumps(dict(placement=placement, device_buffer_allreduce="PASS"), indent=2))
