#!/usr/bin/env python3
"""Check invalid Background MCC inputs in isolated subprocesses."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from pywarpx import amrex, picmi

CASES = {
    "missing_cross_section_units": (
        "Every attachment process must specify",
        "attachment",
        {},
    ),
    "missing_third_body_density": (
        "Attachment cross sections in m5 require a positive",
        "attachment",
        {"cross_section_units": "m5"},
    ),
    "invalid_cross_section_units": (
        "Attachment cross_section_units must be either m2 or m5.",
        "attachment",
        {"cross_section_units": "barn"},
    ),
    "density_with_m2": (
        "Attachment third_body_density is only valid",
        "attachment",
        {"cross_section_units": "m2", "third_body_density": 1.0e25},
    ),
    "malformed_cross_section": (
        "Each cross-section table row must contain",
        "attachment",
        {"cross_section_units": "m2"},
    ),
    "wrong_product_charge": (
        "attachment product species must have charge -q_e",
        "attachment",
        {"cross_section_units": "m2"},
    ),
    "iaa_elastic_missing_dcs": (
        "IAA elastic or excitation scattering requires a",
        "elastic",
        {"scattering_angle_model": "IAA"},
    ),
    "elastic_dcs_without_iaa": (
        "<process>_differential_cross_section requires",
        "elastic",
        {"differential_cross_section": "valid"},
    ),
    "inconsistent_elastic_dcs": (
        "Elastic differential-cross-section rows must have the same number",
        "elastic",
        {"scattering_angle_model": "IAA", "differential_cross_section": "inconsistent"},
    ),
    "zero_integral_elastic_dcs": (
        "Every elastic differential-cross-section energy row must have a",
        "elastic",
        {"scattering_angle_model": "IAA", "differential_cross_section": "zero"},
    ),
}


def run_invalid_case(case):
    amrex.throw_exception = 1
    amrex.signal_handling = 0

    _, process_type, process_options = CASES[case]
    process_options = process_options.copy()

    source_dir = Path(__file__).resolve().parent
    cross_section = source_dir / "background_mcc_attachment_m2.txt"
    if process_type == "elastic":
        cross_section = source_dir / "background_mcc_relativistic_elastic.txt"
    if case == "malformed_cross_section":
        cross_section = Path("background_mcc_malformed_cross_section.txt").resolve()
        cross_section.write_text("0.0 1.0e-22\nnot-a-number 1.0e-22\n")

    dcs_kind = process_options.pop("differential_cross_section", None)
    if dcs_kind is not None:
        differential_cross_section = Path(
            "background_mcc_invalid_elastic_dcs.txt"
        ).resolve()
        if dcs_kind == "valid":
            differential_cross_section.write_text("10 1 1 1\n100 1 1 1\n")
        elif dcs_kind == "inconsistent":
            differential_cross_section.write_text("10 1 1 1\n100 1 1 1 1\n")
        else:
            assert dcs_kind == "zero"
            differential_cross_section.write_text("10 0 0 0\n100 0 0 0\n")
        process_options["differential_cross_section"] = str(differential_cross_section)

    grid = picmi.Cartesian1DGrid(
        number_of_cells=[1],
        lower_bound=[0.0],
        upper_bound=[1.0],
        lower_boundary_conditions=["periodic"],
        upper_boundary_conditions=["periodic"],
        lower_boundary_conditions_particles=["periodic"],
        upper_boundary_conditions_particles=["periodic"],
        warpx_max_grid_size=1,
        warpx_blocking_factor=1,
    )
    solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)

    electrons = picmi.Species(
        particle_type="electron",
        name="electrons",
        initial_distribution=picmi.UniformDistribution(density=1.0),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )
    product_charge = -0.995 * picmi.constants.q_e
    if case != "wrong_product_charge":
        product_charge = -picmi.constants.q_e
    negative_ions = picmi.Species(
        name="negative_ions",
        charge=product_charge,
        mass=32.0 * picmi.constants.m_p + picmi.constants.m_e,
        initial_distribution=picmi.UniformDistribution(density=1.0),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )

    process = {"cross_section": str(cross_section), **process_options}
    if process_type == "attachment":
        process["species"] = negative_ions
    collision = picmi.MCCCollisions(
        name="mcc",
        species=electrons,
        background_density=1.0e20,
        background_temperature=0.0,
        background_mass=32.0 * picmi.constants.m_p,
        scattering_processes={process_type: process},
        nu_max=1.0e6,
    )
    sim = picmi.Simulation(
        solver=solver,
        time_step_size=1.0e-9,
        max_steps=1,
        warpx_collisions=[collision],
        verbose=0,
    )
    sim.add_species(
        electrons,
        layout=picmi.GriddedLayout(n_macroparticle_per_cell=[1], grid=grid),
    )
    sim.add_species(
        negative_ions,
        layout=picmi.GriddedLayout(n_macroparticle_per_cell=[1], grid=grid),
    )
    sim.initialize_inputs()
    sim.initialize_warpx()
    sim.step(1)


def check_invalid_cases():
    for case, (expected_message, _, _) in CASES.items():
        environment = os.environ.copy()
        environment["PYTHONFAULTHANDLER"] = "0"
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--case", case],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=60,
            env=environment,
        )
        assert result.returncode != 0, f"Invalid case {case!r} unexpectedly succeeded"
        assert expected_message in result.stdout, (
            f"Invalid case {case!r} did not report {expected_message!r}.\n"
            f"Output:\n{result.stdout[-8000:]}"
        )
        print(f"{case}: rejected as expected")


parser = argparse.ArgumentParser()
parser.add_argument("--case", choices=CASES)
args = parser.parse_args()
if args.case is None:
    check_invalid_cases()
else:
    run_invalid_case(args.case)
