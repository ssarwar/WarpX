# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check the physical charge returned by the Python deposition wrapper in RZ."""

import numpy as np
import pytest
from conftest import rtol

import pywarpx
from pywarpx import amrex, boundary, picmi
from pywarpx.particle_containers import ParticleContainerWrapper

pytestmark = pytest.mark.skipif(
    pywarpx.libwarpx.geometry_dim != "rz", reason="uses cylindrical volume elements"
)


@pytest.mark.parametrize("filtered", [False, True])
@pytest.mark.parametrize("axis_correction", [False, True])
def test_wrapper_finalizes_accumulated_charge_once(filtered, axis_correction):
    grid = picmi.CylindricalGrid(
        number_of_cells=[16, 16],
        lower_bound=[0.0, 0.0],
        upper_bound=[1.0, 1.0],
        lower_boundary_conditions=["none", "periodic"],
        upper_boundary_conditions=["none", "periodic"],
        lower_boundary_conditions_particles=["none", "periodic"],
        upper_boundary_conditions_particles=["absorbing", "periodic"],
        n_azimuthal_modes=1,
        warpx_max_grid_size=8,
        warpx_blocking_factor=8,
    )
    sim = picmi.Simulation(
        solver=picmi.ElectromagneticSolver(grid=grid, method="Yee", divE_cleaning=True),
        max_steps=1,
        particle_shape="cubic",
        warpx_use_filter=filtered,
        verbose=0,
    )
    sim.add_species(
        picmi.Species(name="electrons", particle_type="electron"), layout=None
    )
    sim.initialize_inputs()
    pywarpx.warpx.do_electrostatic = "none"
    pywarpx.warpx.filter_npass_each_dir = [0, 1]
    boundary.verboncoeur_axis_correction = axis_correction
    amrex.throw_exception = 1
    amrex.signal_handling = 0
    amrex.the_arena_init_size = 0
    sim.initialize_warpx()

    # Include a shape crossing the axis and shapes crossing box interfaces.
    weights = np.array([1e6, 2e6, 4e6])
    sim.particles.get("electrons").add_particles(
        x=np.array([0.01, 0.49, 0.6]),
        y=np.zeros(3),
        z=np.array([0.01, 0.49, 0.7]),
        ux=np.zeros(3),
        uy=np.zeros(3),
        uz=np.zeros(3),
        w=weights,
        unique_particles=False,
    )
    expected_charge = -picmi.constants.q_e * weights.sum()
    field = sim.fields.get("rho_fp", level=0)

    def host(value):
        return value.get() if hasattr(value, "get") else np.asarray(value)

    def check_charge(copies):
        density = host(field[...])[..., 0]
        radius = host(field.mesh("r"))
        dr = dz = 1.0 / 16
        volume = 2 * np.pi * radius * dr * dz
        volume[0] = np.pi * dr**2 * dz / (3 if axis_correction else 4)
        # Exclude the duplicate periodic axial endpoint.
        actual = np.sum(density[:, :16] * volume[:, None])
        np.testing.assert_allclose(
            actual, copies * expected_charge, rtol=rtol(), atol=0
        )

    with pytest.warns(UserWarning, match="deprecated"):
        wrapper = ParticleContainerWrapper("electrons")
    wrapper.deposit_charge_density(0, clear_rho=True, sync_rho=False)
    wrapper.deposit_charge_density(0, clear_rho=False, sync_rho=True)
    check_charge(2)
    wrapper.deposit_charge_density(0, clear_rho=True, sync_rho=True)
    check_charge(1)
