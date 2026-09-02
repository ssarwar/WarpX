#!/usr/bin/env python3
"""Independent checks of corrected PJG ionization and product kinematics."""

import math

import numpy as np

ELECTRON_REST_ENERGY = 510_998.95069  # eV
CLASSICAL_ELECTRON_RADIUS = 2.817_940_3205e-15  # m
PI_E4 = math.pi * (CLASSICAL_ELECTRON_RADIUS * ELECTRON_REST_ENERGY) ** 2
C = 299_792_458.0
K_B = 1.380_649e-23
Q_E = 1.602_176_634e-19

PARAMETERS = {
    "N2": {
        "electrons": 14,
        "b0": 0.029,
        "b1": 1.035,
        "e0": 8239.0,
        "t1": 53.3,
        "gamma1_fixed": 115.0,
        "k": 7.58e-20,
        "gamma_s": 11.1,
        "gamma_numerator": 1.27e4,
        "gamma_denominator": 1.81e3,
        "t_s": 4.0,
        "t_numerator": 2.03e4,
        "t_denominator": 1.97e3,
        "j": 3.39,
        "nu": -1.93e-1,
        "delta": 84.0,
        "thresholds": np.array([15.58, 16.73, 18.75, 22.0, 23.6, 40.0]),
        "fractions": np.array([0.456, 0.2, 0.104, 0.07, 0.07, 0.1]),
        "bethe": np.array([2.48, 2.66, 2.99, 3.50, 3.76, 6.37]),
    },
    "O2": {
        "electrons": 16,
        "b0": 0.030,
        "b1": 1.035,
        "e0": 8239.0,
        "t1": 68.3,
        "gamma1_fixed": 189.1,
        "k": 6.55e-20,
        "gamma_s": 13.1,
        "gamma_numerator": 5.0e5,
        "gamma_denominator": 7.60e4,
        "t_s": 6.34,
        "t_numerator": 2.52e3,
        "t_denominator": 1.28e2,
        "j": 40.3,
        "nu": 3.14e-1,
        "delta": 132.1,
        "thresholds": np.array([12.1, 16.1, 16.9, 18.2, 20.3, 23.0, 37.0]),
        "fractions": np.array([0.08, 0.19, 0.19, 0.17, 0.11, 0.16, 0.1]),
        "bethe": np.array([1.93, 2.56, 2.69, 2.90, 3.23, 3.66, 5.89]),
    },
}

# Values obtained from a separate analytic implementation of the corrected
# equations. Units are m^2 and cover the low-energy peak through 800 MeV.
REFERENCE_TOTALS = {
    "N2": np.array(
        [
            1.2048977e-19,
            5.2962138e-20,
            2.8050342e-20,
            1.6588984e-20,
            2.4482458e-21,
            3.5145421e-22,
            1.0512378e-22,
        ]
    ),
    "O2": np.array(
        [
            4.6814019e-20,
            4.8664983e-20,
            3.3626409e-20,
            2.2053017e-20,
            3.5518735e-21,
            4.7091225e-22,
            1.2730101e-22,
        ]
    ),
}
REFERENCE_ENERGIES = np.array([50e3, 200e3, 500e3, 1e6, 10e6, 100e6, 800e6])


def projectile_state(energy, rest_energy):
    gamma = 1.0 + energy / rest_energy
    beta_squared = 1.0 - gamma**-2
    mass_ratio = ELECTRON_REST_ENERGY / rest_energy
    maximum_transfer = (
        2.0
        * ELECTRON_REST_ENERGY
        * beta_squared
        * gamma**2
        / (1.0 + 2.0 * gamma * mass_ratio + mass_ratio**2)
    )
    return (
        gamma,
        beta_squared,
        0.5 * ELECTRON_REST_ENERGY * beta_squared,
        maximum_transfer,
    )


def printed_pjg_maximum_transfer(energy, rest_energy):
    return (
        energy
        * (energy + 2.0 * rest_energy)
        / (
            energy
            + ELECTRON_REST_ENERGY
            + rest_energy / ELECTRON_REST_ENERGY * (energy + rest_energy)
        )
    )


def continuum_sdcs(target, energy, secondary_energy, rest_energy):
    p = PARAMETERS[target]
    _, beta_squared, equivalent_energy, maximum_transfer = projectile_state(
        energy, rest_energy
    )
    total_energy = energy + rest_energy
    gamma_width = (
        p["gamma_numerator"] / (equivalent_energy + p["gamma_denominator"])
        + p["gamma_s"]
    )
    peak_energy = p["t_s"] - p["t_numerator"] / (equivalent_energy + p["t_denominator"])
    log_energy = math.log(equivalent_energy / p["e0"])
    reduction = p["b0"] * (log_energy**2 + p["b1"])
    power = p["nu"] + 1.0
    distortion = equivalent_energy**power / (p["j"] ** power + equivalent_energy**power)

    secondary_energy = np.asarray(secondary_energy)
    line_shape = 1.0 / (
        (secondary_energy - peak_energy) ** 2 + gamma_width**2
    ) - reduction / ((secondary_energy - p["t1"]) ** 2 + p["gamma1_fixed"] ** 2)
    contributions = []
    for threshold, fraction, bethe_constant in zip(
        p["thresholds"], p["fractions"], p["bethe"], strict=True
    ):
        bethe_factor = (
            math.log(
                4.0
                * equivalent_energy
                * bethe_constant
                / (threshold * (1.0 - beta_squared))
                + math.e
            )
            - beta_squared
        )
        soft = p["k"] * gamma_width**2 * bethe_factor * line_shape
        hard = (
            p["electrons"]
            * PI_E4
            * (
                1.0 / (2.0 * total_energy**2)
                - beta_squared
                / (
                    (maximum_transfer + threshold + p["delta"])
                    * (secondary_energy + threshold)
                )
            )
        )
        contributions.append(fraction * distortion / equivalent_energy * (soft + hard))
    return np.asarray(contributions)


def total_cross_section(target, energy, rest_energy):
    p = PARAMETERS[target]
    _, beta_squared, equivalent_energy, maximum_transfer = projectile_state(
        energy, rest_energy
    )
    total_energy = energy + rest_energy
    gamma_width = (
        p["gamma_numerator"] / (equivalent_energy + p["gamma_denominator"])
        + p["gamma_s"]
    )
    peak_energy = p["t_s"] - p["t_numerator"] / (equivalent_energy + p["t_denominator"])
    log_energy = math.log(equivalent_energy / p["e0"])
    reduction = p["b0"] * (log_energy**2 + p["b1"])
    power = p["nu"] + 1.0
    distortion = equivalent_energy**power / (p["j"] ** power + equivalent_energy**power)

    total = 0.0
    for threshold, fraction, bethe_constant in zip(
        p["thresholds"], p["fractions"], p["bethe"], strict=True
    ):
        line_integral = (
            math.atan((maximum_transfer - peak_energy) / gamma_width)
            - math.atan(-peak_energy / gamma_width)
        ) / gamma_width - reduction * (
            math.atan((maximum_transfer - p["t1"]) / p["gamma1_fixed"])
            - math.atan(-p["t1"] / p["gamma1_fixed"])
        ) / p["gamma1_fixed"]
        bethe_factor = (
            math.log(
                4.0
                * equivalent_energy
                * bethe_constant
                / (threshold * (1.0 - beta_squared))
                + math.e
            )
            - beta_squared
        )
        soft = p["k"] * gamma_width**2 * bethe_factor * line_integral
        hard = (
            p["electrons"]
            * PI_E4
            * (
                maximum_transfer / (2.0 * total_energy**2)
                - beta_squared
                / (maximum_transfer + threshold + p["delta"])
                * math.log1p(maximum_transfer / threshold)
            )
        )
        total += fraction * distortion / equivalent_energy * (soft + hard)
    return total


def kinetic_energy(ux, uy, uz, mass):
    proper_speed_squared = ux**2 + uy**2 + uz**2
    gamma = np.sqrt(1.0 + proper_speed_squared / C**2)
    return mass * proper_speed_squared / ((gamma + 1.0) * Q_E)


def reference_cdf(target, energy, rest_energy):
    maximum_transfer = projectile_state(energy, rest_energy)[3]
    log_energy = np.linspace(0.0, np.log1p(maximum_transfer), 500_001)
    secondary_energy = np.expm1(log_energy)
    sdcs = np.sum(continuum_sdcs(target, energy, secondary_energy, rest_energy), axis=0)
    assert np.min(sdcs) >= 0.0
    density_in_log_energy = sdcs * (secondary_energy + 1.0)
    cumulative = np.zeros_like(secondary_energy)
    cumulative[1:] = np.cumsum(
        0.5
        * (density_in_log_energy[:-1] + density_in_log_energy[1:])
        * np.diff(log_energy)
    )
    numerical_total = cumulative[-1]
    analytic_total = total_cross_section(target, energy, rest_energy)
    assert np.isclose(numerical_total, analytic_total, rtol=2.0e-7)
    return secondary_energy, cumulative / numerical_total


data = np.load("proton_impact_ionization_results.npz")
projectile_energy = float(data["projectile_energy"])
projectile_rest_energy = float(data["projectile_rest_energy"])

# The exact two-body result must recover the nonrelativistic limit, and the
# printed PJG expression must exhibit the diagnosed 800 MeV underestimate.
_, _, _, maximum_transfer = projectile_state(projectile_energy, projectile_rest_energy)
assert np.isclose(maximum_transfer, 2_480_739.6170, rtol=2.0e-10)
printed_transfer = printed_pjg_maximum_transfer(
    projectile_energy, projectile_rest_energy
)
assert maximum_transfer / printed_transfer > 3.69
nonrelativistic_transfer = projectile_state(1.0e3, projectile_rest_energy)[3]
nonrelativistic_mass_ratio = ELECTRON_REST_ENERGY / projectile_rest_energy
assert np.isclose(
    nonrelativistic_transfer,
    4.0 * nonrelativistic_mass_ratio / (1.0 + nonrelativistic_mass_ratio) ** 2 * 1.0e3,
    rtol=3.0e-6,
)

for target in ["N2", "O2"]:
    calculated = np.array(
        [
            total_cross_section(target, energy, projectile_rest_energy)
            for energy in REFERENCE_ENERGIES
        ]
    )
    assert np.allclose(calculated, REFERENCE_TOTALS[target], rtol=6.0e-7)

    # Independently bound the error introduced by the 256-point runtime table,
    # and verify that its entire default range defines a positive probability.
    table_energy = np.geomspace(1.0e3, 1.0e9, 256)
    table_total = np.array(
        [
            total_cross_section(target, energy, projectile_rest_energy)
            for energy in table_energy
        ]
    )
    midpoint_energy = np.sqrt(table_energy[:-1] * table_energy[1:])
    midpoint_exact = np.array(
        [
            total_cross_section(target, energy, projectile_rest_energy)
            for energy in midpoint_energy
        ]
    )
    midpoint_interpolated = 0.5 * (table_total[:-1] + table_total[1:])
    maximum_table_error = np.max(np.abs(midpoint_interpolated / midpoint_exact - 1.0))
    assert maximum_table_error < 1.0e-3

    for scan_energy in np.geomspace(1.0e3, 1.0e9, 25):
        scan_maximum = projectile_state(scan_energy, projectile_rest_energy)[3]
        scan_secondary = np.expm1(np.linspace(0.0, np.log1p(scan_maximum), 2049))
        scan_sdcs = np.sum(
            continuum_sdcs(target, scan_energy, scan_secondary, projectile_rest_energy),
            axis=0,
        )
        assert np.min(scan_sdcs) >= 0.0

    for component in ["x", "y", "z", "ux", "uy", "uz", "w", "id"]:
        assert np.array_equal(
            data[f"{target}_beam_initial_{component}"],
            data[f"{target}_beam_{component}"],
        )

    electron_weight = data[f"{target}_electrons_w"]
    ion_weight = data[f"{target}_ions_w"]
    assert electron_weight.size == ion_weight.size
    assert 10000 < electron_weight.size <= 20000
    assert np.array_equal(electron_weight, ion_weight)
    assert np.all(electron_weight == float(data["fixed_product_weight"]))
    for position in ["x", "y", "z"]:
        assert np.array_equal(
            data[f"{target}_electrons_{position}"],
            data[f"{target}_ions_{position}"],
        )
    for species in ["electrons", "ions"]:
        particle_ids = data[f"{target}_{species}_id"]
        assert np.unique(particle_ids).size == particle_ids.size

    gamma, beta_squared, _, maximum_transfer = projectile_state(
        projectile_energy, projectile_rest_energy
    )
    projectile_speed = C * math.sqrt(beta_squared)
    expected_weight = (
        float(data["projectile_density"])
        * float(data["background_density"])
        * total_cross_section(target, projectile_energy, projectile_rest_energy)
        * projectile_speed
        * float(data["time_step"])
    )
    observed_weight = np.sum(electron_weight)
    assert np.isclose(observed_weight, expected_weight, rtol=2.0e-3)

    electron_energy = kinetic_energy(
        data[f"{target}_electrons_ux"],
        data[f"{target}_electrons_uy"],
        data[f"{target}_electrons_uz"],
        9.109_383_7139e-31,
    )
    assert np.all(electron_energy >= 0.0)
    assert np.max(electron_energy) <= maximum_transfer * (1.0 + 2.0e-12)

    reference_energy, reference_probability = reference_cdf(
        target, projectile_energy, projectile_rest_energy
    )
    sorted_energy = np.sort(electron_energy)
    sampled_probability = np.interp(
        sorted_energy, reference_energy, reference_probability
    )
    empirical_probability = (
        np.arange(sorted_energy.size, dtype=float) + 0.5
    ) / sorted_energy.size
    assert np.max(np.abs(sampled_probability - empirical_probability)) < 4.0e-3

    quantiles = np.array([0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
    expected_quantiles = np.interp(quantiles, reference_probability, reference_energy)
    observed_quantiles = np.quantile(electron_energy, quantiles)
    print(
        f"{target} energy quantiles expected={expected_quantiles} "
        f"observed={observed_quantiles}"
    )
    assert np.allclose(
        observed_quantiles[:-1], expected_quantiles[:-1], rtol=2.0e-2, atol=0.03
    )
    assert np.isclose(observed_quantiles[-1], expected_quantiles[-1], rtol=8.0e-2)

    direction = data[f"{target}_direction"]
    electron_u = np.column_stack(
        [
            data[f"{target}_electrons_ux"],
            data[f"{target}_electrons_uy"],
            data[f"{target}_electrons_uz"],
        ]
    )
    electron_u_magnitude = np.linalg.norm(electron_u, axis=1)
    cosine = electron_u @ direction / electron_u_magnitude
    assert np.all(np.abs(cosine) <= 1.0 + 2.0e-14)

    contributions = np.maximum(
        continuum_sdcs(
            target, projectile_energy, electron_energy, projectile_rest_energy
        ),
        0.0,
    )
    binding_energy = (PARAMETERS[target]["thresholds"] @ contributions) / np.sum(
        contributions, axis=0
    )
    free_cosine = np.sqrt(
        electron_energy
        * (maximum_transfer + 2.0 * ELECTRON_REST_ENERGY)
        / (maximum_transfer * (electron_energy + 2.0 * ELECTRON_REST_ENERGY))
    )
    expected_cosine = (
        free_cosine
        * (electron_energy + 0.5 * binding_energy)
        / (electron_energy + binding_energy)
    )
    assert abs(np.mean(cosine) - np.mean(expected_cosine)) < 1.5e-2
    assert np.mean(cosine[electron_energy >= np.quantile(electron_energy, 0.9)]) > (
        np.mean(cosine[electron_energy <= np.quantile(electron_energy, 0.5)])
    )
    mean_direction = np.mean(electron_u / electron_u_magnitude[:, None], axis=0)
    transverse_mean = mean_direction - np.dot(mean_direction, direction) * direction
    assert np.linalg.norm(transverse_mean) < 1.5e-2

    ion_mass = float(data[f"{target}_neutral_mass"]) - 9.109_383_7139e-31
    thermal_speed = math.sqrt(K_B * float(data[f"{target}_temperature"]) / ion_mass)
    ion_u = np.column_stack(
        [
            data[f"{target}_ions_ux"],
            data[f"{target}_ions_uy"],
            data[f"{target}_ions_uz"],
        ]
    )
    assert np.all(np.abs(np.mean(ion_u, axis=0)) < 0.03 * thermal_speed)
    assert np.allclose(np.std(ion_u, axis=0), thermal_speed, rtol=0.035)

    print(
        f"{target}: pairs={electron_weight.size}, sigma={calculated[-1]:.9e} m^2, "
        f"Tmax={maximum_transfer:.6f} eV, KS="
        f"{np.max(np.abs(sampled_probability - empirical_probability)):.3e}, "
        f"table-error={maximum_table_error:.3e}"
    )

assert float(data["initialization_elapsed"]) < 30.0
assert float(data["step_elapsed"]) < 30.0
