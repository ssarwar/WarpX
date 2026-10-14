# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Preserve collision/fluid inputs through PICMI validation and serialization."""

import tempfile
from pathlib import Path

import numpy as np

import pywarpx
from pywarpx import picmi

beam = picmi.FluidSpecies(
    name="beam",
    particle_type="proton",
    model="rigid_beam",
    kinetic_energy=800e6,
    sigma_r=0.005,
    sigma_t=25e-12,
    pulse_times=np.array([0.0, 5e-9]),
    peak_current=0.6,
)
ions = picmi.FluidSpecies(name="ions", model="immobile", mass=4.65e-26, charge="q_e")
electrons = picmi.Species(name="electrons", particle_type="electron")
protons = picmi.ProtonImpactIonizationCollisions(
    name="protons",
    species=beam,
    product_species=[electrons, ions],
    ionization_target="N2",
    background_density="1e20*(1+z)",
    background_temperature=300.0,
    fixed_product_weight=7.0,
    source_sampling_points=12,
    gas_quadrature_points=4,
    sampling_seed=17,
    ndt_subcycle=8,
    start_step=2,
)
mcc = picmi.MCCCollisions(
    name="mcc",
    species=electrons,
    background_density=1e20,
    background_temperature="300+z",
    background_mass=4.65e-26,
    nu_max=2e10,
    ndt_subcycle=8,
    scattering_processes={
        "ionization": {
            "cross_section": "ionization.txt",
            "species": ions,
            "energy_sharing_model": "RBEQ",
            "rbeq_target": "N2",
        },
        "elastic": {
            "cross_section": "elastic.txt",
            "rotation_model": "reciprocal_hybrid",
            "rotation_file": "thermal_rotation.rot",
            "rotational_temperature": 300.0,
        },
    },
)
grid = picmi.Cartesian1DGrid(
    number_of_cells=[8],
    lower_bound=[0.0],
    upper_bound=[1.0],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
)
simulation = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(grid=grid, method="Yee"),
    warpx_collisions=[protons, mcc],
    max_steps=0,
    particle_shape="linear",
)
simulation.add_species(electrons, layout=None)
simulation.add_fluid_species(beam)
simulation.add_fluid_species(ions)
for name, kind in [("charge", "ParticleCharge"), ("budget", "PrescribedSourceBudget")]:
    simulation.add_diagnostic(
        picmi.ReducedDiagnostic(name=name, diag_type=kind, period=1)
    )

histogram = picmi.ReducedDiagnostic(
    name="histogram",
    diag_type="ParticleHistogram2D",
    species=electrons,
    period=1,
    bin_number_abs=4,
    bin_number_ord=8,
    bin_min_abs=-1,
    bin_max_abs=1,
    bin_min_ord=-1,
    bin_max_ord=1,
    histogram_function_abs="ux",
    histogram_function_ord="uz",
    value_function="w",
)
histogram.openpmd_backend = "h5"
simulation.add_diagnostic(histogram)

# A nested collision must retain its concrete type and all physical parameters.
restored = picmi.Simulation.model_validate_json(simulation.model_dump_json())
assert isinstance(restored.collisions[0], picmi.ProtonImpactIonizationCollisions)
assert isinstance(restored.collisions[1], picmi.MCCCollisions)
assert isinstance(restored.collisions[0].species, picmi.FluidSpecies)
assert isinstance(
    restored.collisions[1].scattering_processes["ionization"]["species"],
    picmi.FluidSpecies,
)
assert len(restored.fluid_species) == 2

with tempfile.TemporaryDirectory() as directory:
    inputs = Path(directory) / "inputs"
    restored.write_input_file(inputs)
    parameters = {}
    for line in inputs.read_text().splitlines():
        if " = " in line:
            key, value = line.split(" = ", 1)
            parameters[key] = value.strip('"')
    expected = {
        "fluids.species_names": "beam ions",
        "particles.species_names": "electrons",
        "beam.species_type": "proton",
        "beam.model": "rigid_beam",
        "beam.kinetic_energy": "800000000.0",
        "ions.model": "immobile",
        "ions.charge": "q_e",
        "protons.type": "proton_impact_ionization",
        "protons.species": "beam",
        "protons.product_species": "electrons ions",
        "protons.background_density(x,y,z,t)": "1e20*(1+z)",
        "protons.fixed_product_weight": "7.0",
        "protons.source_sampling_points": "12",
        "protons.gas_quadrature_points": "4",
        "protons.sampling_seed": "17",
        "protons.ndt_subcycle": "8",
        "protons.start_step": "2",
        "mcc.nu_max": "20000000000.0",
        "mcc.background_temperature(x,y,z,t)": "300+z",
        "mcc.ionization_species": "ions",
        "mcc.ionization_energy_sharing_model": "RBEQ",
        "mcc.elastic_rotation_model": "reciprocal_hybrid",
        "mcc.elastic_rotation_file": "thermal_rotation.rot",
        "charge.type": "ParticleCharge",
        "budget.type": "PrescribedSourceBudget",
        "histogram.type": "ParticleHistogram2D",
        "histogram.openpmd_backend": "h5",
        "histogram.value_function(t,x,y,z,ux,uy,uz,w)": "w",
    }
    for key, value in expected.items():
        assert parameters[key] == value, (key, parameters[key], value)

try:
    restored.add_fluid_species(picmi.FluidSpecies(name="late", particle_type="proton"))
except RuntimeError:
    pass
else:
    raise AssertionError("Fluid registration after initialization must fail")
restored.finalize()
assert not pywarpx.fluids.species_names
assert not picmi.Simulation().fluid_species
print("PICMI collision and fluid validation, serialization and inputs passed")
