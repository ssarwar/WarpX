#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check multiple-source particle counts, density scaling and histogram weights."""

import numpy as np
import openpmd_api as io
from mpi4py import MPI

from pywarpx import picmi

comm = MPI.COMM_WORLD
grid = picmi.Cartesian3DGrid(
    number_of_cells=[8, 8, 8],
    lower_bound=[0, 0, 0],
    upper_bound=[1, 1, 1],
    lower_boundary_conditions=["periodic"] * 3,
    upper_boundary_conditions=["periodic"] * 3,
    warpx_max_grid_size=4,
    warpx_blocking_factor=4,
)
gaussian = picmi.GaussianBunchDistribution(
    n_physical_particles=80,
    rms_bunch_size=[0.01] * 3,
    centroid_position=[0.5] * 3,
    centroid_velocity=[0, 0, 0],
    warpx_do_symmetrize=False,
)
uniform = picmi.UniformDistribution(density=40)
protons = picmi.Species(
    name="protons",
    particle_type="proton",
    density_scale=2.5,
    initial_distribution=[gaussian, uniform],
    warpx_do_not_deposit=True,
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    time_step_size=1e-12,
    max_steps=1,
    particle_shape=1,
    verbose=0,
)
sim.add_species(
    protons,
    layout=[
        picmi.PseudoRandomLayout(grid=grid, n_macroparticles=64),
        picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[1, 1, 1]),
    ],
)
histogram = picmi.ReducedDiagnostic(
    diag_type="ParticleHistogram2D",
    name="counts",
    species=protons,
    period=1,
    bin_number_abs=1,
    bin_min_abs=0,
    bin_max_abs=1,
    bin_number_ord=1,
    bin_min_ord=0,
    bin_max_ord=1,
    histogram_function_abs="x",
    histogram_function_ord="z",
)
histogram.openpmd_backend = "h5"
sim.add_diagnostic(histogram)
sim.initialize_inputs()

# PICMI's ParticleListDistribution uses "multipleparticles" even for one
# particle. Add a native single-particle source first to exercise the source
# loop's continuation, followed by the PICMI Gaussian and uniform sources.
protons.species.injection_sources = ["one", "dist0", "dist1"]
for name, value in dict(
    injection_style="singleparticle",
    single_particle_pos=[0.25] * 3,
    single_particle_u=[0, 0, 0],
    single_particle_weight=3,
).items():
    protons.species.add_new_group_attr("one", name, value)

sim.initialize_warpx()
sim.step(1)
local_count, local_weight = 0, 0.0
for tile in sim.particles.get("protons").iterator(level=0):
    weights = tile["w"]
    weights = weights.get() if hasattr(weights, "get") else np.asarray(weights)
    local_count += len(weights)
    local_weight += float(weights.sum())
assert comm.allreduce(local_count) == 1 + 64 + 8**3
# N_single + scale*N_Gaussian + scale*density_uniform*volume.
expected_weight = 3 + 2.5 * 80 + 2.5 * 40 * 1
np.testing.assert_allclose(comm.allreduce(local_weight), expected_weight, rtol=1e-13)
sim.finalize()

if comm.rank == 0:
    series = io.Series("diags/reducedfiles/counts/openpmd_%T.h5", io.Access.read_only)
    for step in series.iterations:
        histogram = series.iterations[step].meshes["data"][io.Record_Component.SCALAR]
        values = histogram.load_chunk()
        series.flush()
        np.testing.assert_allclose(values.sum(), expected_weight, rtol=1e-13)
    series.close()
    print("PASS: all injection sources and the histogram preserve statistical weights")
