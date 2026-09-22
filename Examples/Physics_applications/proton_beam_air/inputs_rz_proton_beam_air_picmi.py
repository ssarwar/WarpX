#!/usr/bin/env python3
"""Illustrative 800 MeV proton bunch in prescribed dry air; see README.md.

Supply evaluated cross sections using cross_sections.example.json as the schema.
The numerical settings illustrate setup and have not been converged for a BPM.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from pywarpx import picmi

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--cross-sections", type=Path, required=True)
parser.add_argument("--pressure-pa", type=float, default=101325.0)
parser.add_argument("--temperature-k", type=float, default=293.15)
parser.add_argument("--dt", type=float, default=1.0e-13)
parser.add_argument("--steps", type=int, default=200)
parser.add_argument("--mcc-subcycles", type=int, default=1)
parser.add_argument("--product-weight", type=float, default=1.0e6)
parser.add_argument("--write-input", type=Path)
args = parser.parse_args()

for name in ("pressure_pa", "temperature_k", "dt", "product_weight"):
    value = getattr(args, name)
    if not np.isfinite(value) or value <= 0:
        parser.error(f"{name} must be finite and positive")
if args.steps < 1 or args.mcc_subcycles < 1:
    parser.error("steps and mcc-subcycles must be positive")

C = picmi.constants.c
QE = picmi.constants.q_e
ME = picmi.constants.m_e
MP = picmi.constants.m_p
MU = 1.66053906660e-27
KB = 1.380649e-23
gas_density = args.pressure_pa / (KB * args.temperature_k)
# Dry-air partial densities; omitted argon/trace gases are not renormalized away.
densities = {"N2": 0.78084 * gas_density, "O2": 0.20946 * gas_density}
masses = {"N2": 28.0134 * MU, "O2": 31.9988 * MU}

energy_ev = 800.0e6
gamma = 1.0 + energy_ev * QE / (MP * C**2)
proper_speed = C * np.sqrt((gamma - 1.0) * (gamma + 1.0))
beam_speed = proper_speed / gamma
beam_radius = 0.005
peak_current = 0.6
sigma_t = 25.0e-12  # The reported 100 ps is 4 sigma, not sigma or FWHM.
sigma_z = beam_speed * sigma_t
peak_density = peak_current / (QE * beam_speed * np.pi * beam_radius**2)

grid = picmi.CylindricalGrid(
    number_of_cells=[64, 256],
    n_azimuthal_modes=1,
    lower_bound=[0.0, 0.0],
    upper_bound=[0.02, 0.10],
    # Conducting boundaries support the initial beam self-field Poisson solve.
    lower_boundary_conditions=["none", "dirichlet"],
    upper_boundary_conditions=["dirichlet", "dirichlet"],
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    warpx_max_grid_size=64,
    warpx_blocking_factor=8,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee")
beam = picmi.Species(
    name="protons",
    particle_type="proton",
    initial_distribution=picmi.AnalyticDistribution(
        density_expression=(
            f"{peak_density:.17g}*exp(-0.5*((z-0.03)/{sigma_z:.17g})^2)"
            f"*(x*x+y*y < {beam_radius**2:.17g})"
        ),
        directed_velocity=[0.0, 0.0, proper_speed],
    ),
    # Ballistic bunch: deposit its charge/current and move at its prescribed speed.
    warpx_do_not_gather=True,
)
electrons = picmi.Species(name="electrons", particle_type="electron")
ions = {
    "N2_plus": picmi.Species(name="N2_plus", charge=QE, mass=masses["N2"] - ME),
    "O2_plus": picmi.Species(name="O2_plus", charge=QE, mass=masses["O2"] - ME),
    "O_minus": picmi.Species(name="O_minus", charge=-QE, mass=15.9994 * MU + ME),
    "O2_minus": picmi.Species(name="O2_minus", charge=-QE, mass=masses["O2"] + ME),
}

manifest_path = args.cross_sections.resolve()
manifest = json.loads(manifest_path.read_text())
collisions = []
for target in ("N2", "O2"):
    collisions.append(
        picmi.ProtonImpactIonizationCollisions(
            name=f"p_{target}",
            species=beam,
            product_species=[electrons, ions[f"{target}_plus"]],
            ionization_target=target,
            background_density=densities[target],
            background_temperature=args.temperature_k,
            fixed_product_weight=args.product_weight,
            max_products_per_cell=8,
        )
    )

for target in ("N2", "O2"):
    processes = {}
    for process_name, supplied in manifest[target].items():
        process = dict(supplied)
        if process.get("cross_section_units") == "m5":
            third_density = process.get("third_body_density")
            if (
                third_density is None
                or not np.isfinite(third_density)
                or third_density <= 0
            ):
                raise ValueError(
                    f"{target}/{process_name}: specify the third-body density in m^-3 "
                    "for the collider/mixture documented by your m^5 dataset"
                )
        for key in ("cross_section", "differential_cross_section"):
            if key in process:
                path = (manifest_path.parent / process[key]).resolve()
                if not path.is_file():
                    raise FileNotFoundError(path)
                process[key] = str(path)
        if "species" in process:
            process["species"] = ions[process["species"]]
        processes[process_name] = process
    collisions.append(
        picmi.MCCCollisions(
            name=f"e_{target}",
            species=electrons,
            background_density=densities[target],
            background_temperature=args.temperature_k,
            background_mass=masses[target],
            scattering_processes=processes,
            ndt_subcycle=args.mcc_subcycles,
        )
    )

sim = picmi.Simulation(
    solver=solver,
    time_step_size=args.dt,
    max_steps=args.steps,
    particle_shape="linear",
    warpx_collisions=collisions,
    warpx_random_seed=42,
    verbose=1,
)
sim.add_species(
    beam,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2, 4, 2], grid=grid),
    initialize_self_field=True,
)
for species in (electrons, *ions.values()):
    sim.add_species(species, layout=None)  # Initially empty kinetic products.
sim.add_diagnostic(
    picmi.FieldDiagnostic(
        name="fields", grid=grid, period=100, data_list=["E", "B", "J", "rho"]
    )
)
sim.add_diagnostic(
    picmi.ParticleDiagnostic(
        name="particles",
        period=100,
        species=[beam, electrons, *ions.values()],
        data_list=["position", "momentum", "weighting"],
    )
)
if args.write_input:
    sim.write_input_file(file_name=str(args.write_input))
else:
    sim.step()
    sim.finalize()
