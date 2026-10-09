#!/usr/bin/env python3
# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check exact magnetic-field save/restore with real/imaginary RZ modes.

This exercises production buffer allocation and the same MultiFab copy
operations as the semi-implicit retry hooks. It does not evolve the fields:
the existing implicit vector algebra only operates on the first component,
so this test makes no claim about full multimode implicit plasma dynamics.
"""

import numpy as np

from pywarpx import picmi

modes = 3
grid = picmi.CylindricalGrid(
    number_of_cells=[8, 16],
    lower_bound=[0, 0],
    upper_bound=[1, 1],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    n_azimuthal_modes=modes,
    warpx_max_grid_size=8,
    warpx_blocking_factor=8,
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    warpx_evolve_scheme=picmi.SemiImplicitEMEvolveScheme(
        nonlinear_solver=picmi.PicardNonlinearSolver()
    ),
    warpx_current_deposition_algo="direct",
    time_step_size=1e-12,
    max_steps=1,
    particle_shape=1,
    verbose=0,
)
sim.add_species(picmi.Species(name="electrons", particle_type="electron"), layout=None)
sim.initialize_inputs()
sim.initialize_warpx()

for direction in ["r", "theta", "z"]:
    current = sim.fields.get("Bfield_fp", direction, 0)
    saved = sim.fields.get("B_old", direction, 0)
    # Fail before an invalid copy if a future change undersizes the buffer.
    assert saved.n_comp == current.n_comp == 2 * modes - 1
    assert saved.n_grow_vect == current.n_grow_vect
    values = current[()]
    radial = current.imesh(0, include_ghosts=True)
    axial = current.imesh(1, include_ghosts=True)
    for component in range(2 * modes - 1):
        values[..., component] = 1e-7 * (
            1 + component + radial[:, None] / 16 + axial[None, :] / 32
        )
    current[()] = values
    reference = current[()].copy()

    # These are the save/restore copy extents in CopyVectorField. Include
    # every ghost cell; tests on valid cells alone miss an undersized buffer.
    saved.copymf(current, 0, 0, saved.n_comp, saved.n_grow_vect)
    current.set_val(0.0)
    current.copymf(saved, 0, 0, current.n_comp, current.n_grow_vect)
    np.testing.assert_array_equal(current[()], reference)

sim.finalize()
print("PASS: exact three-mode RZ magnetic rollback, including ghost cells")
