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


def dcs_screening_radius(path):
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[:3] == ["SPECIES:", "e", "/"]:
            return {"N2": 0.6052, "O2": 0.5677}.get(fields[3], 0.0)
    return 0.0


def screened_rutherford_statistics(energy_ev, screening_radius):
    fine_structure = 7.2973525693e-3
    tau = energy_ev / (M_E * C**2 / Q_E)
    k_sq = tau * (tau + 2.0) / fine_structure**2
    eta = 1.0 / (4.0 * screening_radius**2 * k_sq)
    probabilities = np.linspace(0.0, 1.0, 200001)
    deflection = 2.0 * eta * probabilities / (1.0 - probabilities + eta)
    cosine = 1.0 - deflection
    moments = np.array(
        [np.trapezoid(cosine**order, probabilities) for order in range(1, 9)]
    )
    backward_probability = eta / (1.0 + 2.0 * eta)
    return moments, backward_probability


def read_dcs_table(path):
    rows = []
    for line in path.read_text().splitlines():
        try:
            values = np.array([float(value) for value in line.split()])
        except ValueError:
            continue
        if values.size >= 4:
            rows.append(values)
    return np.array(rows)


GAUSS_NODES, GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(12)


def dcs_integrals(dcs):
    theta = np.linspace(0.0, math.pi, dcs.size)
    half_sines = np.sin(0.5 * theta)
    integrals = np.zeros(9)
    backward_integral = 0.0
    backward_boundary = math.sqrt(0.5)

    def integrate_segment(y_lo, y_hi, dcs_lo, slope):
        midpoint = 0.5 * (y_lo + y_hi)
        half_width = 0.5 * (y_hi - y_lo)
        y = midpoint + half_width * GAUSS_NODES
        density = y * (dcs_lo + slope * (y - y_lo))
        cosine = 1.0 - 2.0 * y**2
        values = np.vstack(
            (density, *(density * cosine**order for order in range(1, 9)))
        )
        return half_width * (values @ GAUSS_WEIGHTS)

    for index in range(dcs.size - 1):
        y_lo = half_sines[index]
        y_hi = half_sines[index + 1]
        slope = (dcs[index + 1] - dcs[index]) / (y_hi - y_lo)
        integrals += integrate_segment(y_lo, y_hi, dcs[index], slope)
        if y_hi > backward_boundary:
            clipped_lo = max(y_lo, backward_boundary)
            clipped_dcs = dcs[index] + slope * (clipped_lo - y_lo)
            backward_integral += integrate_segment(
                clipped_lo, y_hi, clipped_dcs, slope
            )[0]
    return integrals[0], integrals[1:], backward_integral


def dcs_statistics(dcs):
    normalization, moment_integrals, backward_integral = dcs_integrals(dcs)
    return moment_integrals / normalization, backward_integral / normalization


def analytic_dcs_statistics(anisotropy):
    moments = np.array(
        [
            1.0 / (order + 1) if order % 2 == 0 else anisotropy / (order + 2)
            for order in range(1, 9)
        ]
    )
    return moments, 0.5 - 0.25 * anisotropy


def interpolated_dcs_statistics(dcs_lo, dcs_hi, energy_fraction):
    integral_lo, moments_lo, backward_lo = dcs_integrals(dcs_lo)
    integral_hi, moments_hi, backward_hi = dcs_integrals(dcs_hi)
    normalization = (
        (1.0 - energy_fraction) * integral_lo + energy_fraction * integral_hi
    )
    moments = (
        (1.0 - energy_fraction) * moments_lo + energy_fraction * moments_hi
    ) / normalization
    backward_probability = (
        (1.0 - energy_fraction) * backward_lo + energy_fraction * backward_hi
    ) / normalization
    return moments, backward_probability


def dcs_file_statistics(path, energy_ev):
    table = read_dcs_table(path)
    energies = table[:, 0]
    if energy_ev <= energies[0]:
        return dcs_statistics(table[0, 1:])
    if energy_ev >= energies[-1]:
        return dcs_statistics(table[-1, 1:])
    index = np.searchsorted(energies, energy_ev) - 1
    energy_fraction = math.log(energy_ev / energies[index]) / math.log(
        energies[index + 1] / energies[index]
    )
    return interpolated_dcs_statistics(
        table[index, 1:], table[index + 1, 1:], energy_fraction
    )


if args.n2_dcs is None:
    dcs_path = Path("background_mcc_iaa_elastic_dcs.txt").resolve()
    o2_dcs_path = Path("background_mcc_iaa_elastic_o2_dcs.txt").resolve()
    table_energies = np.array([10.0, 100.0, 1000.0, 8000.0, 10000.0, 1.0e9])
    anisotropies = np.array([-0.6, 0.0, 0.0, 0.0, 0.9, 0.3])
    theta = np.linspace(0.0, math.pi, 361)
    table = np.column_stack(
        (table_energies, np.array([1.0 + a * np.cos(theta) for a in anisotropies]))
    )
    # Put most of the 8 keV probability inside the first 0.5-degree cell. This
    # catches both loss of the theta=0 endpoint and incorrect normalization
    # while interpolating between rows with very different angular integrals.
    table[3, 1:] = 1.0e-3
    table[3, 1] = 1.0e6
    # Deliberately make the unresolved 1 GeV row inconsistent with the analytic
    # continuation. Recognized N2/O2 metadata must select screened Rutherford.
    table[-1, 1:] = np.exp(-0.5 * (theta / 0.004) ** 2) + 1.0e-8
    np.savetxt(
        dcs_path,
        table,
        header=(
            "Synthetic DCS with the IAA/elmolcs row layout\n"
            "SPECIES: e / N2\n"
            "COLUMNS: theta = linspace(0, 180, 361) (deg)"
        ),
        comments="",
    )
    np.savetxt(
        o2_dcs_path,
        table,
        header=(
            "Synthetic DCS with the IAA/elmolcs row layout\n"
            "SPECIES: e / O2\n"
            "COLUMNS: theta = linspace(0, 180, 361) (deg)"
        ),
        comments="",
    )
    dcs_rows = table[:, 1:]
    sampled_energies = (10.0, 100.0, 550.0, 5050.0, 8000.0, 10100.0, 1.0e9)
    sampled_statistics = (
        analytic_dcs_statistics(-0.6),
        analytic_dcs_statistics(0.0),
        analytic_dcs_statistics(0.0),
        interpolated_dcs_statistics(
            dcs_rows[2],
            dcs_rows[3],
            math.log(5050.0 / 1000.0) / math.log(8000.0 / 1000.0),
        ),
        dcs_statistics(dcs_rows[3]),
        screened_rutherford_statistics(1.01e4, 0.6052),
        screened_rutherford_statistics(1.0e9, 0.6052),
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
            screened_rutherford_statistics(1.0e9, 0.6052),
            "excitation",
            excitation_energy,
        ),
        (
            "elastic_o2_10p1kev",
            1.01e4,
            31.9988 * AMU,
            o2_dcs_path,
            screened_rutherford_statistics(1.01e4, 0.5677),
            "elastic",
            0.0,
        ),
        (
            "excitation_o2_gev",
            1.0e9,
            31.9988 * AMU,
            o2_dcs_path,
            screened_rutherford_statistics(1.0e9, 0.5677),
            "excitation",
            1.0,
        ),
    ]
else:
    cases = []
    for target, mass, excitation_energy, path in (
        ("n2", 28.0134 * AMU, 6.0, args.n2_dcs.resolve()),
        ("o2", 31.9988 * AMU, 1.0, args.o2_dcs.resolve()),
    ):
        for energy in (0.1, 1000.0, 9990.0, 10010.0, 1.0e6, 1.0e9):
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
case_screening_radii = []
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
    screening_radius = dcs_screening_radius(dcs_file)
    case_screening_radii.append(screening_radius)

    if SCATTERING_ANGLE_MODEL == "isotropic":
        reference_moments, reference_backward = analytic_dcs_statistics(0.0)
    elif statistics is not None:
        reference_moments, reference_backward = statistics
    elif screening_radius > 0.0 and energy_ev >= 1.0e4:
        reference_moments, reference_backward = screened_rutherford_statistics(
            energy_ev, screening_radius
        )
    else:
        reference_moments, reference_backward = dcs_file_statistics(
            dcs_file, energy_ev
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
    "case_screening_radii": np.array(case_screening_radii),
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
    transverse_fraction = np.clip(
        (ux[scattered] ** 2 + uy[scattered] ** 2) / proper_speed_sq[scattered],
        0.0,
        1.0,
    )
    absolute_cosine = np.sqrt(1.0 - transverse_fraction)
    results[f"{name}_deflections"] = np.where(
        uz[scattered] >= 0.0,
        transverse_fraction / (1.0 + absolute_cosine),
        1.0 + absolute_cosine,
    )
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
