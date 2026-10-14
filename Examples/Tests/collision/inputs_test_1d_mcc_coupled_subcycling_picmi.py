#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check hot -> cold -> attached electrons against an exact rate equation.

Two different collision objects act on the same electrons. Equal constant
rates give hot=e^-x, cold=x*e^-x and attached=1-(1+x)*e^-x. Advancing each
operator for the entire PIC step instead gives cold=(1-e^-x)*e^-x, regardless
of how finely each operator is subcycled.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--subcycles", type=int, default=64)
parser.add_argument("--attachment-subcycles", type=int)
parser.add_argument("--reverse", action="store_true")
parser.add_argument("--solver", choices=["Yee", "semi_implicit_em"], default="Yee")
args = parser.parse_args()
n_attach = args.attachment_subcycles or args.subcycles
comm = MPI.COMM_WORLD
dt, rate, density = 1e-12, 2e12, 1e25
mass = 28.0134 * 1.66053906892e-27
me, qe, c = picmi.constants.m_e, picmi.constants.q_e, picmi.constants.c
particles, cells = 16384, 32


def speed(energy):
    tau = energy * qe / (me * c**2)
    return c * np.sqrt(tau * (tau + 2)) / (tau + 1)


if comm.rank == 0:
    sigma_hot, sigma_cold = rate / (density * speed(3)), rate / (density * speed(1))
    Path("excitation.txt").write_text(f"0 0\n2 0\n3 {sigma_hot:.17g}\n3.1 0\n100 0\n")
    # The plateau absorbs the tiny finite-mass recoil shift at one eV.
    Path("attachment.txt").write_text(
        f"0 0\n0.99 {sigma_cold:.17g}\n1.01 {sigma_cold:.17g}\n1.5 0\n100 0\n"
    )
comm.Barrier()
grid = picmi.Cartesian1DGrid(
    number_of_cells=[cells],
    lower_bound=[0],
    upper_bound=[1],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    warpx_max_grid_size=8,
    warpx_blocking_factor=8,
)
electrons = picmi.Species(
    name="electrons",
    particle_type="electron",
    initial_distribution=picmi.UniformDistribution(
        density=1, directed_velocity=[0, 0, speed(3)]
    ),
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
    warpx_do_not_deposit=True,
)
negative = picmi.Species(
    name="negative",
    charge="-q_e",
    mass=mass + me,
    warpx_do_not_push=True,
    warpx_do_not_gather=True,
    warpx_do_not_deposit=True,
)
common = dict(
    species=electrons,
    background_density=density,
    background_temperature=0,
    background_mass=mass,
    start_step=2,
)
excitation = picmi.MCCCollisions(
    name="cool",
    ndt_subcycle=args.subcycles,
    scattering_processes={
        "excitation": {"cross_section": "excitation.txt", "energy": 2}
    },
    **common,
)
attachment = picmi.MCCCollisions(
    name="attach",
    ndt_subcycle=n_attach,
    scattering_processes={
        "attachment": {
            "cross_section": "attachment.txt",
            "cross_section_units": "m2",
            "species": negative,
        }
    },
    **common,
)
collisions = [attachment, excitation] if args.reverse else [excitation, attachment]
implicit = args.solver != "Yee"
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    time_step_size=dt,
    max_steps=3,
    particle_shape=1,
    verbose=0,
    warpx_collisions=collisions,
    warpx_random_seed=173,
    warpx_current_deposition_algo="direct" if implicit else None,
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.NewtonNonlinearSolver(
            relative_tolerance=1e-12,
            linear_solver=picmi.GMRESLinearSolver(relative_tolerance=1e-12),
        )
    )
    if implicit
    else None,
)
sim.add_species(
    electrons,
    picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[particles // cells]),
)
sim.add_species(negative, layout=None)


def host(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def counts():
    result = np.zeros(3, dtype=np.int64)
    for tile in sim.particles.get("electrons").iterator(level=0):
        valid = host(libwarpx.amr.unpack_ids(host(tile["idcpu"]))) > 0
        u2 = sum(host(tile[f"u{axis}"]).astype(float) ** 2 for axis in "xyz")
        energy = me * u2 / ((np.sqrt(1 + u2 / c**2) + 1) * qe)
        result[0] += np.count_nonzero(valid & (energy > 2))
        result[1] += np.count_nonzero(valid & (energy <= 2))
    for tile in sim.particles.get("negative").iterator(level=0):
        result[2] += np.count_nonzero(
            host(libwarpx.amr.unpack_ids(host(tile["idcpu"]))) > 0
        )
    return comm.allreduce(result)


sim.step(2)
np.testing.assert_array_equal(counts(), [particles, 0, 0])
sim.step(1)
observed = counts()
sim.finalize()
assert observed.sum() == particles, "Excitation/attachment did not conserve charge"
x = rate * dt
expected = np.array([np.exp(-x), x * np.exp(-x), 1 - (1 + x) * np.exp(-x)])
noise = 6 * np.sqrt(expected * (1 - expected) / particles)
# First-order operator splitting, plus a conservative finite-neutral-mass bound.
bound = 2 / min(args.subcycles, n_attach) + 2e-4
assert np.all(abs(observed / particles - expected) < noise + bound), (
    observed,
    expected,
)
if args.subcycles == n_attach:
    # The two exact Bernoulli propagators also have a closed-form discrete solution.
    h = np.exp(-x / args.subcycles)
    cold = args.subcycles * (1 - h) * np.exp(-x) / (h if args.reverse else 1)
    discrete = np.array([np.exp(-x), cold, 1 - np.exp(-x) - cold])
    assert np.all(abs(observed / particles - discrete) < noise + 2e-4), (
        observed,
        discrete,
    )
if comm.rank == 0:
    result = dict(
        observed=(observed / particles).tolist(), exact=expected.tolist(), **vars(args)
    )
    Path("coupled-subcycling.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        "PASS: coupled collision subcycles converge to the exact rate equations", result
    )
