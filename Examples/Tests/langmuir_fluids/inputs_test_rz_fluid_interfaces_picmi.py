#!/usr/bin/env python3
"""Run the same fluid configuration through PICMI and a written text input."""

import argparse
import os
import subprocess
from pathlib import Path

import numpy as np

from pywarpx import fluids, new_fluid_species, picmi, warpx

parser = argparse.ArgumentParser()
parser.add_argument("--executable", required=True)
args = parser.parse_args()

beam_parameters = dict(
    name="beam",
    model="rigid_beam",
    particle_type="proton",
    kinetic_energy=800e6,
    sigma_r=0.002,
    sigma_t=25e-12,
    first_pulse_time=0.0,
    pulse_period=1 / 201.25e6,
    pulse_count=2,
    bunch_charge=0.6 * np.sqrt(2 * np.pi) * 25e-12,
)
for changes in [
    {"velocity_z": 1e8},
    {"r_rms": 0.002},
    {"peak_density": 1e12},
    {"pulse_times": [0.0]},
    {"initial_density": 1.0},
    {"sigma_r": None},
]:
    try:
        picmi.FluidSpecies(**(beam_parameters | changes))
    except ValueError:
        pass
    else:
        raise AssertionError(f"Invalid fluid input accepted: {changes}")

grid = picmi.CylindricalGrid(
    number_of_cells=[16, 32],
    lower_bound=[0.0, -0.04],
    upper_bound=[0.016, 0.04],
    lower_boundary_conditions=["none", "periodic"],
    upper_boundary_conditions=["none", "periodic"],
    n_azimuthal_modes=1,
    warpx_max_grid_size=16,
    warpx_blocking_factor=8,
)
sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    time_step_size=1e-13,
    max_steps=2,
    particle_shape=2,
    verbose=0,
)
beam = picmi.FluidSpecies(**beam_parameters)
positive = picmi.FluidSpecies(
    name="positive",
    model="immobile",
    charge=picmi.constants.q_e,
    mass=28e-27,
    initial_density="1e10*(1+x)",
)
cold = picmi.FluidSpecies(
    name="cold",
    particle_type="electron",
    initial_density=1e10,
)
for species in [beam, positive, cold]:
    sim.add_fluid_species(species)
try:
    sim.add_fluid_species(positive)
except ValueError:
    pass
else:
    raise AssertionError("Duplicate species name accepted")

fields = [
    "rho",
    "Er",
    "Ez",
    "Bt",
    "fluid_density_beam",
    "fluid_density_positive",
    "fluid_density_cold",
    "fluid_current_beamz",
    "part_per_cell",
    "part_per_cell_beam",
    "part_per_cell_positive",
    "part_per_cell_cold",
]
sim.add_diagnostic(
    picmi.FieldDiagnostic(name="fields", grid=grid, period=2, data_list=fields)
)
sim.write_input_file("inputs")
try:
    sim.add_fluid_species(picmi.FluidSpecies(name="late", particle_type="electron"))
except RuntimeError:
    pass
else:
    raise AssertionError("Fluid registration accepted after input initialization")
sim.step()
sim.finalize()
assert fluids.species_names == []
new_fluid_species("fresh", model="immobile", species_type="proton")
assert fluids.species_names == ["fresh"]
assert not any("beam." in value for value in warpx.create_argv_list())
warpx.finalize()

text_run = Path("text_run")
text_run.mkdir(exist_ok=True)
environment = os.environ.copy()
environment.pop("AMREX_INPUTS_FILE_PREFIX", None)
for key in list(environment):
    if key.startswith(("PMI_", "PMIX_", "OMPI_")):
        environment.pop(key)
with (text_run / "output.txt").open("w") as output:
    subprocess.run(
        [args.executable, str(Path("inputs").resolve())],
        cwd=text_run,
        env=environment,
        stdout=output,
        stderr=subprocess.STDOUT,
        check=True,
        timeout=30,
    )

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd().parent / ".mpl-cache"))
import yt

for step in [0, 2]:
    name = f"diags/fields{step:06d}"
    python_data = yt.load(name)
    text_data = yt.load(str(text_run / name))
    python_grid = python_data.covering_grid(
        0, python_data.domain_left_edge, python_data.domain_dimensions
    )
    text_grid = text_data.covering_grid(
        0, text_data.domain_left_edge, text_data.domain_dimensions
    )
    for field in fields:
        if field.startswith("part_per_cell"):
            np.testing.assert_array_equal(python_grid["boxlib", field].v, 0)
        np.testing.assert_allclose(
            python_grid["boxlib", field].v,
            text_grid["boxlib", field].v,
            rtol=3e-14,
            err_msg=field,
        )
print("PASS: PICMI/text fluid evolution, validation and low-level input reset")
