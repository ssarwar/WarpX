# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check beam diagnostics against a discrete distribution with exact moments."""

import itertools
import re
from pathlib import Path

import numpy as np
import pytest
from conftest import Config, rtol
from helpers import make_sim

import pywarpx
from pywarpx import picmi


@pytest.mark.parametrize("empty", [False, True])
def test_beam_relevant_covariance(empty):
    """Independent symmetric signs give exact covariance and Twiss parameters.

    For each axis, x = x0 + a*s and u = u0 + b*(h*s + k*t), with s,t
    independent and equally weighted at +/-1. Thus Var(x)=a**2,
    Var(u)=b**2*(h**2+k**2), Cov(x,u)=a*b*h, and alpha=-h/abs(k).
    Scaling u to normalized momentum u/c leaves alpha dimensionless.
    """
    c = picmi.constants.c
    a, b = 2.0**-10, 2.0**18
    centers = np.array([2.0**-8, 2.0**-9, 2.0**-8])
    u0 = np.array([2.0**21, 2.0**20, 2.0**23])
    h, k = np.array([2.0, -1.0, 3.0]), np.array([3.0, 2.0, 4.0])
    signs = np.array(list(itertools.product([-1.0, 1.0], repeat=6)))
    positions = centers + a * signs[:, ::2]
    velocities = u0 + b * (h * signs[:, ::2] + k * signs[:, 1::2])

    if pywarpx.libwarpx.geometry_dim == "rz":
        # Match make_sim's on-demand allocation for repeated GPU initialization.
        pywarpx.amrex.the_arena_init_size = 0
        pywarpx.amrex.throw_exception = 1
        pywarpx.amrex.signal_handling = 0
        grid = picmi.CylindricalGrid(
            number_of_cells=[8, 16],
            lower_bound=[0, -0.01],
            upper_bound=[0.01, 0.01],
            lower_boundary_conditions=["none", "periodic"],
            upper_boundary_conditions=["none", "periodic"],
            n_azimuthal_modes=1,
            warpx_max_grid_size=8,
        )
        sim = picmi.Simulation(
            solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
            time_step_size=1.0e-15,
            max_steps=1,
            particle_shape="quadratic",
            verbose=0,
        )
    else:
        dims = int(pywarpx.libwarpx.geometry_dim[0])
        sim = make_sim(
            lower_bound=[-0.01] * dims,
            upper_bound=[0.01] * dims,
            dt=1.0e-15,
        )
    beam = picmi.Species(
        particle_type="electron",
        name="beam",
        initial_distribution=picmi.ParticleListDistribution(
            x=positions[:, 0],
            y=positions[:, 1],
            z=positions[:, 2],
            ux=velocities[:, 0],
            uy=velocities[:, 1],
            uz=velocities[:, 2],
            weight=np.ones(len(signs)),
        )
        if not empty
        else None,
        warpx_do_not_push=True,
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )
    sim.add_species(beam, layout=None)
    sim.add_diagnostic(
        picmi.ReducedDiagnostic(
            diag_type="BeamRelevant", name="beam_moments", species=beam, period=1
        )
    )
    sim.step(1)
    if Config.have_mpi:
        from mpi4py import MPI

        if MPI.COMM_WORLD.rank != 0:
            return

    path = Path("diags/reducedfiles/beam_moments.txt")
    columns = re.findall(r"\[\d+\]([^\s(]+)\(", path.read_text().splitlines()[0])
    data = np.loadtxt(path, ndmin=2)
    if empty:
        np.testing.assert_array_equal(data[:, 2:], 0.0)
        return
    for axis, center, drift, shear, spread in zip("xyz", centers, u0, h, k):
        expected = {
            f"{axis}_mean": center,
            f"{axis}_rms": a,
            f"p{axis}_mean": picmi.constants.m_e * drift,
            f"p{axis}_rms": picmi.constants.m_e * b * np.hypot(shear, spread),
            f"emittance_{axis}": a * b * abs(spread) / c,
            f"alpha_{axis}": -shear / abs(spread),
            f"beta_{axis}": a * c / (b * abs(spread)),
        }
        for name, value in expected.items():
            if name in columns:
                np.testing.assert_allclose(
                    data[:, columns.index(name)],
                    value,
                    rtol=rtol(),
                    atol=0.0,
                    err_msg=name,
                )
    np.testing.assert_allclose(
        data[:, columns.index("charge")],
        -len(signs) * picmi.constants.q_e,
        rtol=rtol(),
        atol=0.0,
    )
