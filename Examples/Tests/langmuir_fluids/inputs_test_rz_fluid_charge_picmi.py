#!/usr/bin/env python3
"""Physical charge includes runtime ionization levels and immobile densities."""

from pathlib import Path

import numpy as np

from pywarpx import picmi

grid = picmi.CylindricalGrid(
    number_of_cells=[8, 16],
    lower_bound=[0, 0],
    upper_bound=[0.008, 0.016],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    n_azimuthal_modes=1,
    warpx_max_grid_size=8,
    warpx_blocking_factor=8,
)
helium = picmi.Species(
    name="helium",
    particle_type="He",
    charge_state=2,
    initial_distribution=picmi.UniformDistribution(density=1e10),
    warpx_do_not_deposit=True,
)
electrons = picmi.Species(name="electrons", particle_type="electron")
ions = picmi.FluidSpecies(
    name="ions", model="immobile", particle_type="proton", initial_density=1e10
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    max_steps=1,
    time_step_size=1e-15,
    particle_shape=2,
    verbose=0,
)
sim.add_species(
    helium, picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=[2, 1, 2])
)
sim.add_species(electrons, layout=None)
sim.add_fluid_species(ions)
sim.add_interaction(
    picmi.FieldIonization(
        model="ADK", ionized_species=helium, product_species=electrons
    )
)
for kind in ["ParticleNumber", "ParticleCharge"]:
    sim.add_diagnostic(picmi.ReducedDiagnostic(diag_type=kind, name=kind, period=1))
sim.step()


def reduced(kind):
    path = Path("diags/reducedfiles") / (kind + ".txt")
    names = [key.split("]", 1)[1] for key in path.read_text().splitlines()[0].split()]
    return dict(zip(names, np.loadtxt(path, ndmin=2)[-1]))


number, charge = reduced("ParticleNumber"), reduced("ParticleCharge")
qe = picmi.constants.q_e
for name, z in [("helium", 2), ("ions", 1), ("electrons", -1)]:
    np.testing.assert_allclose(
        charge[name + "(C)"], z * qe * number[name + "_weight()"], rtol=3e-13
    )
np.testing.assert_allclose(
    charge["total(C)"], charge["helium(C)"] + charge["ions(C)"], rtol=3e-13
)
sim.finalize()
