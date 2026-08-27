#!/usr/bin/env python3
"""Sample a tabulated DCS for elastic and excitation MCC processes."""

import argparse
import math
import time
from pathlib import Path

import numpy as np

from pywarpx import libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--particle-count", type=int, default=65536)
parser.add_argument(
    "--scattering-angle-model", choices=("IAA", "isotropic"), default="IAA"
)
parser.add_argument("--n2-dcs", type=Path)
parser.add_argument("--o2-dcs", type=Path)
args = parser.parse_args()

PARTICLE_COUNT = args.particle_count
SCATTERING_ANGLE_MODEL = args.scattering_angle_model
CROSS_SECTION = 1.0e-20
NU_MAX = 1.0e8
# This one-step sampling test deliberately makes nearly every particle collide.
# It is not a time-converged physical simulation and therefore does not use the
# usual nu_max*dt <= 0.1 production guideline.
MAJORANT_OPTICAL_DEPTH = 12.0
RATE_FRACTION = 0.999

C = picmi.constants.c
Q_E = picmi.constants.q_e
M_E = picmi.constants.m_e
AMU = 1.66053906660e-27

assert PARTICLE_COUNT > 0
assert (args.n2_dcs is None) == (args.o2_dcs is None)


def electron_gamma(energy_ev):
    return 1.0 + energy_ev * Q_E / (M_E * C**2)


def electron_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(1.0 - 1.0 / gamma**2)


def electron_proper_speed(energy_ev):
    gamma = electron_gamma(energy_ev)
    return C * math.sqrt(gamma**2 - 1.0)


def read_dcs_row(path, energy_ev):
    rows = []
    for line in path.read_text().splitlines():
        try:
            values = np.array([float(value) for value in line.split()])
        except ValueError:
            continue
        if values.size >= 4:
            rows.append(values)
    matches = [row for row in rows if math.isclose(row[0], energy_ev, rel_tol=1.0e-12)]
    assert len(matches) == 1
    return matches[0][1:]


def dcs_statistics(dcs):
    theta = np.linspace(0.0, math.pi, dcs.size)
    density = dcs * np.sin(theta)
    density[[0, -1]] = 0.0
    normalization = np.trapezoid(density, theta)
    cosine = np.cos(theta)
    moments = np.array(
        [
            np.trapezoid(density * cosine**order, theta) / normalization
            for order in range(1, 9)
        ]
    )
    backward_probability = (
        np.trapezoid(density[theta >= 0.5 * math.pi], theta[theta >= 0.5 * math.pi])
        / normalization
    )
    return moments, backward_probability


def analytic_dcs_statistics(anisotropy):
    moments = np.array(
        [
            1.0 / (order + 1) if order % 2 == 0 else anisotropy / (order + 2)
            for order in range(1, 9)
        ]
    )
    return moments, 0.5 - 0.25 * anisotropy


def inverse_dcs_cdf(dcs, probabilities):
    theta = np.linspace(0.0, math.pi, dcs.size)
    density = dcs * np.sin(theta)
    density[[0, -1]] = 0.0
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(theta)))
    )
    cumulative /= cumulative[-1]
    return np.cos(np.interp(probabilities, cumulative, theta))


def interpolated_dcs_statistics(dcs_lo, dcs_hi, energy_fraction):
    probabilities = np.linspace(0.0, 1.0, 200001)
    cosine = (1.0 - energy_fraction) * inverse_dcs_cdf(
        dcs_lo, probabilities
    ) + energy_fraction * inverse_dcs_cdf(dcs_hi, probabilities)
    moments = np.array(
        [np.trapezoid(cosine**order, probabilities) for order in range(1, 9)]
    )
    zero_probability = np.interp(0.0, -cosine, probabilities)
    return moments, 1.0 - zero_probability


if args.n2_dcs is None:
    dcs_path = Path("background_mcc_iaa_elastic_dcs.txt").resolve()
    table_energies = np.array([10.0, 100.0, 1000.0, 10000.0, 1.0e9])
    anisotropies = np.array([-0.6, 0.0, 0.0, 0.9, 0.3])
    theta = np.linspace(0.0, math.pi, 361)
    table = np.column_stack(
        (table_energies, np.array([1.0 + a * np.cos(theta) for a in anisotropies]))
    )
    np.savetxt(
        dcs_path,
        table,
        header=(
            "Synthetic DCS with the IAA/elmolcs row layout\n"
            "COLUMNS: theta = linspace(0, 180, 361) (deg)"
        ),
    )
    dcs_rows = table[:, 1:]
    sampled_energies = (10.0, 100.0, 550.0, 5050.0, 10000.0, 1.0e9)
    sampled_statistics = (
        analytic_dcs_statistics(-0.6),
        analytic_dcs_statistics(0.0),
        analytic_dcs_statistics(0.0),
        interpolated_dcs_statistics(dcs_rows[2], dcs_rows[3], 0.45),
        analytic_dcs_statistics(0.9),
        analytic_dcs_statistics(0.3),
    )
    cases = [
        (
            f"elastic_analytic_{index}",
            energy,
            28.0134 * AMU,
            dcs_path,
            statistics,
            "elastic",
            0.0,
        )
        for index, (energy, statistics) in enumerate(
            zip(sampled_energies, sampled_statistics)
        )
    ]
    excitation_energy = 7.0
    cases += [
        (
            "excitation_threshold",
            excitation_energy,
            28.0134 * AMU,
            dcs_path,
            analytic_dcs_statistics(-0.6),
            "excitation",
            excitation_energy,
        ),
        (
            "excitation_low",
            10.0,
            28.0134 * AMU,
            dcs_path,
            analytic_dcs_statistics(-0.6),
            "excitation",
            excitation_energy,
        ),
        (
            "excitation_interpolated",
            550.0,
            28.0134 * AMU,
            dcs_path,
            analytic_dcs_statistics(0.0),
            "excitation",
            excitation_energy,
        ),
        (
            "excitation_gev",
            1.0e9,
            28.0134 * AMU,
            dcs_path,
            analytic_dcs_statistics(0.3),
            "excitation",
            excitation_energy,
        ),
    ]
else:
    cases = []
    for target, mass, excitation_energy, path in (
        ("n2", 28.0134 * AMU, 6.0, args.n2_dcs.resolve()),
        ("o2", 31.9988 * AMU, 1.0, args.o2_dcs.resolve()),
    ):
        for energy in (0.1, 1000.0, 1.0e6, 1.0e9):
            cases.append(
                (
                    f"elastic_{target}_{energy:g}",
                    energy,
                    mass,
                    path,
                    None,
                    "elastic",
                    0.0,
                )
            )
        for energy in (1000.0, 1.0e9):
            cases.append(
                (
                    f"excitation_{target}_{energy:g}",
                    energy,
                    mass,
                    path,
                    None,
                    "excitation",
                    excitation_energy,
                )
            )

case_names = []
case_energies = []
case_masses = []
case_process_types = []
case_energy_losses = []
case_rate_fractions = []
expected_moments = []
expected_backward_probabilities = []
species = []
collisions = []

grid = picmi.Cartesian1DGrid(
    number_of_cells=[1],
    lower_bound=[0.0],
    upper_bound=[100.0],
    lower_boundary_conditions=["periodic"],
    upper_boundary_conditions=["periodic"],
    lower_boundary_conditions_particles=["periodic"],
    upper_boundary_conditions_particles=["periodic"],
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)

for (
    name,
    energy_ev,
    neutral_mass,
    dcs_file,
    statistics,
    process_type,
    energy_loss,
) in cases:
    electrons = picmi.Species(
        particle_type="electron",
        name=f"electrons_{name}",
        initial_distribution=picmi.UniformDistribution(
            density=1.0,
            # WarpX's constant momentum injector interprets this PICMI field as
            # proper velocity even though the portable PICMI name says velocity.
            directed_velocity=[0.0, 0.0, electron_proper_speed(energy_ev)],
        ),
        warpx_do_not_deposit=True,
        warpx_do_not_gather=True,
    )
    cross_section_path = Path(f"background_mcc_iaa_{name}_cross_section.txt").resolve()
    if process_type == "elastic":
        cross_section_table = np.array([[0.0, CROSS_SECTION], [1.0e9, CROSS_SECTION]])
    else:
        cross_section_table = np.array(
            [
                [0.0, 0.0],
                [energy_loss, 0.0],
                [energy_loss + 1.0e-3, 0.0],
                [energy_loss + 2.0e-3, CROSS_SECTION],
                [1.0e9, CROSS_SECTION],
            ]
        )
    np.savetxt(cross_section_path, cross_section_table)
    local_cross_section = np.interp(energy_ev, *cross_section_table.T)
    case_rate_fraction = RATE_FRACTION * local_cross_section / CROSS_SECTION

    process = {
        "cross_section": str(cross_section_path),
        "scattering_angle_model": SCATTERING_ANGLE_MODEL,
    }
    if process_type == "excitation":
        process["energy"] = energy_loss
    if SCATTERING_ANGLE_MODEL == "IAA":
        process["differential_cross_section"] = str(dcs_file)
    collision = picmi.MCCCollisions(
        name=f"mcc_{name}",
        species=electrons,
        background_density=(
            RATE_FRACTION * NU_MAX / (CROSS_SECTION * electron_speed(energy_ev))
        ),
        background_temperature=0.0,
        background_mass=neutral_mass,
        scattering_processes={process_type: process},
        nu_max=NU_MAX,
    )
    species.append(electrons)
    collisions.append(collision)
    case_names.append(name)
    case_energies.append(energy_ev)
    case_masses.append(neutral_mass)
    case_process_types.append(process_type)
    case_energy_losses.append(energy_loss)
    case_rate_fractions.append(case_rate_fraction)

    if SCATTERING_ANGLE_MODEL == "isotropic":
        reference_moments, reference_backward = analytic_dcs_statistics(0.0)
    elif statistics is not None:
        reference_moments, reference_backward = statistics
    else:
        reference_moments, reference_backward = dcs_statistics(
            read_dcs_row(dcs_file, energy_ev)
        )
    expected_moments.append(reference_moments)
    expected_backward_probabilities.append(reference_backward)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=MAJORANT_OPTICAL_DEPTH / NU_MAX,
    max_steps=1,
    warpx_collisions=collisions,
    verbose=1,
)
for electrons in species:
    sim.add_species(
        electrons,
        layout=picmi.GriddedLayout(
            n_macroparticle_per_cell=[PARTICLE_COUNT], grid=grid
        ),
    )


def to_numpy(array):
    return array.get() if hasattr(array, "get") else np.asarray(array)


def component(container, name):
    return np.concatenate(
        [to_numpy(pti[name]) for pti in container.iterator(level=0)]
    ).astype(np.float64)


sim.initialize_inputs()
initialization_start = time.perf_counter()
sim.initialize_warpx()
initialization_elapsed = time.perf_counter() - initialization_start
step_start = time.perf_counter()
sim.step(1)
step_elapsed = time.perf_counter() - step_start

results = {
    "case_names": np.array(case_names),
    "case_energies": np.array(case_energies),
    "case_masses": np.array(case_masses),
    "case_process_types": np.array(case_process_types),
    "case_energy_losses": np.array(case_energy_losses),
    "case_rate_fractions": np.array(case_rate_fractions),
    "expected_moments": np.array(expected_moments),
    "expected_backward_probabilities": np.array(expected_backward_probabilities),
    "particle_count": PARTICLE_COUNT,
    "rate_fraction": RATE_FRACTION,
    "scattering_angle_model": SCATTERING_ANGLE_MODEL,
    "initialization_elapsed": initialization_elapsed,
    "step_elapsed": step_elapsed,
    "particle_calls_per_second": PARTICLE_COUNT * len(cases) / step_elapsed,
}
for name, energy_ev, _, _, _, _, _ in cases:
    container = sim.particles.get(f"electrons_{name}")
    ux = component(container, "ux")
    uy = component(container, "uy")
    uz = component(container, "uz")
    proper_speed_sq = ux**2 + uy**2 + uz**2
    proper_speed = np.sqrt(proper_speed_sq)
    gamma = np.sqrt(1.0 + proper_speed_sq / C**2)
    scattered = np.logical_or(ux != 0.0, uy != 0.0)
    results[f"{name}_cosines"] = uz[scattered] / proper_speed[scattered]
    results[f"{name}_azimuth_x"] = ux[scattered] / proper_speed[scattered]
    results[f"{name}_azimuth_y"] = uy[scattered] / proper_speed[scattered]
    results[f"{name}_energies"] = (
        M_E * proper_speed_sq[scattered] / ((gamma[scattered] + 1.0) * Q_E)
    )

if libwarpx.amr.ParallelDescriptor.MyProc() == 0:
    np.savez("background_mcc_iaa_elastic_results.npz", **results)
    print(
        "IAA differential-scattering benchmark: "
        f"model={SCATTERING_ANGLE_MODEL}, cases={len(cases)}, "
        f"particles/case={PARTICLE_COUNT}, "
        f"initialization={initialization_elapsed:.6f} s, "
        f"step={step_elapsed:.6f} s, "
        f"throughput={results['particle_calls_per_second']:.3e} particle-calls/s"
    )

sim.finalize()
