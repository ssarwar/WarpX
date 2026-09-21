#!/usr/bin/env python3
"""Check native fluid storage and restore its state before advancing a restart."""

import argparse
from pathlib import Path

import numpy as np

from pywarpx import algo, picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument("--solver", default="Yee", choices=["Yee", "PSATD"])
parser.add_argument("--restart")
args = parser.parse_args()

grid = picmi.CylindricalGrid(
    number_of_cells=[16, 32],
    lower_bound=[0.0, -0.01],
    upper_bound=[0.01, 0.01],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    warpx_max_grid_size=16,
    warpx_blocking_factor=8,
    n_azimuthal_modes=1,
)
solver = picmi.ElectromagneticSolver(
    grid=grid, method=args.solver, stencil_order=[8, 8]
)
sim = picmi.Simulation(
    solver=solver,
    time_step_size=1.0e-13,
    max_steps=4,
    particle_shape="cubic",
    warpx_amr_restart=args.restart,
    verbose=0,
)
sim.add_diagnostic(picmi.Checkpoint(name="chk", period=2, write_dir="diags"))
sim.initialize_inputs()
algo.particle_shape = 3
warpx.get_bucket("fluids").species_names = ["positive", "negative"]
for name, sign in [("positive", 1), ("negative", -1)]:
    species = warpx.get_bucket(name)
    species.model = "immobile"
    species.mass = 28 * 1.66053906660e-27
    species.charge = sign * picmi.constants.q_e
    species.profile = "parse_density_function"
    species.__setattr__(
        "density_function(x,y,z)",
        "0" if args.restart else "1e12*(1+10*x+0.01*cos(100*pi*z))",
    )
sim.initialize_warpx()


def host(values):
    return values.get() if hasattr(values, "get") else np.asarray(values)


def check():
    for name in ("positive", "negative"):
        density = sim.fields.get(f"fluid_density_{name}", level=0)
        r = host(density.mesh("r"))
        z = host(density.mesh("z"))
        expected = 1e12 * (
            1 + 10 * r[:, None] + 0.01 * np.cos(100 * np.pi * z[None, :])
        )
        np.testing.assert_allclose(np.squeeze(host(density[...])), expected, rtol=2e-14)
        assert not any(
            f"fluid_momentum_density_{name}" in key for key in sim.fields.list()
        )


check()  # The restart check deliberately precedes the first resumed timestep.
start = sim.extension.warpx.getistep(lev=0)
assert start == (2 if args.restart else 0)
sim.step(4 - start)
check()
assert np.isclose(sim.extension.warpx.gett_new(lev=0), 4e-13, rtol=2e-15)
if args.restart:
    assert Path("diags/chk000004/FluidModels").is_file()
print("PASS: native immobile density, zero momentum storage and checkpoint restoration")
sim.finalize()
