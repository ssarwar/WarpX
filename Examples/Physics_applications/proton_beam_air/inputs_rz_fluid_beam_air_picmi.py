#!/usr/bin/env python3
"""800 MeV Gaussian pulse train, kinetic electrons and immobile air ions.

Supply evaluated electron tables with --cross-sections, or select --proton-only
for the source-only example matching inputs_rz_fluid_beam_air. See README.md.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from pywarpx import picmi

parser = argparse.ArgumentParser(description=__doc__)
chemistry = parser.add_mutually_exclusive_group(required=True)
chemistry.add_argument("--cross-sections", type=Path)
chemistry.add_argument("--proton-only", action="store_true")
parser.add_argument(
    "--sigma-r", type=float, required=True, help="Cartesian transverse RMS width [m]"
)
parser.add_argument(
    "--frequency", type=float, default=201.25e6, help="Pulse frequency [Hz]"
)
parser.add_argument("--pulses", type=int, default=2)
parser.add_argument("--pressure-pa", type=float, default=101325.0)
parser.add_argument("--temperature-k", type=float, default=293.15)
parser.add_argument("--dt", type=float, default=1e-13)
parser.add_argument("--steps", type=int, default=200)
parser.add_argument("--mcc-subcycles", type=int, default=1)
parser.add_argument("--source-subcycles", type=int, default=1)
parser.add_argument("--product-weight", type=float, default=1e6)
parser.add_argument("--write-input", type=Path)
parser.add_argument("--restart", type=Path)
args = parser.parse_args()
for name in [
    "sigma_r",
    "frequency",
    "pressure_pa",
    "temperature_k",
    "dt",
    "product_weight",
]:
    if not np.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
        parser.error(f"{name} must be finite and positive")
if min(args.pulses, args.steps, args.mcc_subcycles, args.source_subcycles) < 1:
    parser.error("Pulse, step and subcycle counts must be positive")

qe, me = picmi.constants.q_e, picmi.constants.m_e
amu = 1.66053906660e-27
gas_density = args.pressure_pa / (picmi.constants.kb * args.temperature_k)
densities = {"N2": 0.78084 * gas_density, "O2": 0.20946 * gas_density}
masses = {"N2": 28.0134 * amu, "O2": 31.9988 * amu}
grid = picmi.CylindricalGrid(
    number_of_cells=[64, 256],
    n_azimuthal_modes=1,
    lower_bound=[0.0, 0.0],
    upper_bound=[0.02, 0.10],
    lower_boundary_conditions=["none", "dirichlet"],
    upper_boundary_conditions=["dirichlet", "dirichlet"],
    lower_boundary_conditions_particles=["none", "absorbing"],
    upper_boundary_conditions_particles=["absorbing", "absorbing"],
    warpx_max_grid_size=64,
    warpx_blocking_factor=8,
)
beam = picmi.FluidSpecies(
    name="beam",
    model="rigid_beam",
    particle_type="proton",
    kinetic_energy=800e6,
    sigma_r=args.sigma_r,
    sigma_t=25e-12,
    peak_current=0.6,
    z_reference=0.03,
    first_pulse_time=0.0,
    pulse_period=1 / args.frequency,
    pulse_count=args.pulses,
    initialize_self_fields=True,
)
electrons = picmi.Species(name="electrons", particle_type="electron")
ions = {
    name: picmi.FluidSpecies(name=name, model="immobile", mass=mass, charge=charge)
    for name, mass, charge in [
        ("N2_plus", masses["N2"] - me, qe),
        ("O2_plus", masses["O2"] - me, qe),
        ("O_minus", 15.9994 * amu + me, -qe),
        ("O2_minus", masses["O2"] + me, -qe),
    ]
}
collisions = [
    picmi.ProtonImpactIonizationCollisions(
        name="p_" + target,
        species=beam,
        product_species=[electrons, ions[target + "_plus"]],
        ionization_target=target,
        background_density=densities[target],
        background_temperature=args.temperature_k,
        fixed_product_weight=args.product_weight,
        max_products_per_cell=8,
        ndt_subcycle=args.source_subcycles,
    )
    for target in ["N2", "O2"]
]
if args.cross_sections:
    manifest_path = args.cross_sections.resolve()
    manifest = json.loads(manifest_path.read_text())
    for target in ["N2", "O2"]:
        processes = {}
        for name, supplied in manifest[target].items():
            process = dict(supplied)
            if process.get("cross_section_units") == "m5":
                density = process.get("third_body_density")
                if density is None or not np.isfinite(density) or density <= 0:
                    raise ValueError(
                        f"{target}/{name}: supply the table's third-body density [m^-3]"
                    )
            for key in ["cross_section", "differential_cross_section"]:
                if key in process:
                    path = (manifest_path.parent / process[key]).resolve()
                    if not path.is_file():
                        raise FileNotFoundError(path)
                    process[key] = str(path)
            if "species" in process:
                process["species"] = ions[process["species"]]
            processes[name] = process
        collisions.append(
            picmi.MCCCollisions(
                name="e_" + target,
                species=electrons,
                background_density=densities[target],
                background_temperature=args.temperature_k,
                background_mass=masses[target],
                scattering_processes=processes,
                ndt_subcycle=args.mcc_subcycles,
            )
        )

sim = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    time_step_size=args.dt,
    max_steps=args.steps,
    particle_shape="linear",
    warpx_collisions=collisions,
    warpx_random_seed=42,
    warpx_amr_restart=str(args.restart) if args.restart else None,
    verbose=1,
)
sim.add_fluid_species(beam)
for ion in ions.values():
    sim.add_fluid_species(ion)
sim.add_species(electrons, layout=None)
sim.add_diagnostic(
    picmi.FieldDiagnostic(
        name="fields",
        grid=grid,
        period=100,
        data_list=[
            "E",
            "B",
            "J",
            "rho",
            "rho_beam",
            "rho_electrons",
            "fluid_density_beam",
            "fluid_current_beamz",
        ]
        + ["fluid_density_" + name for name in ions],
    )
)
sim.add_diagnostic(
    picmi.ParticleDiagnostic(
        name="particles",
        period=100,
        species=[electrons],
        data_list=["position", "momentum", "weighting"],
    )
)
sim.add_diagnostic(picmi.Checkpoint(name="checkpoint", period=100))
for kind in [
    "ParticleNumber",
    "ParticleCharge",
    "ParticleEnergy",
    "ParticleMomentum",
    "PrescribedSourceBudget",
]:
    sim.add_diagnostic(picmi.ReducedDiagnostic(name=kind, diag_type=kind, period=100))
if args.write_input:
    sim.write_input_file(str(args.write_input))
else:
    sim.step()
sim.finalize()
