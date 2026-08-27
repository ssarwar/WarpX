#!/usr/bin/env python3
"""Compare sampled N2/O2 ionization energies with independent RBEQ moments."""

import math

import numpy as np

PARTICLE_COUNT = 32768
MC2_EV = 510998.95069
LN2 = math.log(2.0)

CASES = [
    ("n2_near", "N2", 15.78),
    ("n2_mid", "N2", 100.0),
    ("n2_high", "N2", 1000.0),
    ("o2_near", "O2", 12.27),
    ("o2_mid", "O2", 100.0),
    ("o2_high", "O2", 1000.0),
]

SHELLS = {
    "N2": np.array(
        [
            [409.50, 603.30, 4.0, 1.000],
            [37.30, 71.13, 2.0, 0.760],
            [18.72, 63.18, 2.0, 1.000],
            [16.74, 44.30, 4.0, 0.938],
            [15.58, 54.91, 2.0, 0.792],
        ]
    ),
    "O2": np.array(
        [
            [543.80, 796.20, 4.0, 1.0000],
            [40.33, 79.73, 2.0, 0.9600],
            [27.05, 90.92, 2.0, 1.0000],
            [20.30, 71.84, 2.0, 1.0000],
            [17.08, 59.89, 4.0, 1.0000],
            [12.07, 84.88, 2.0, 0.9314],
        ]
    ),
}


def rbeq_terms(energy, shell):
    binding, kinetic, occupation, q = shell
    if energy <= binding:
        return None
    t = energy / MC2_EV
    b = binding / MC2_EV
    u = kinetic / MC2_EV
    ratio = energy / binding
    gamma_tilde = 1.0 + t + u + b
    beta_tilde_sq = 1.0 - 1.0 / gamma_tilde**2
    correction = -(1.0 + q - (5.0 - 3.0 * q) * LN2) / q
    prefactor = occupation / (2.0 * b * beta_tilde_sq)
    binding_sq = (b / gamma_tilde) ** 2
    exchange = (2.0 * gamma_tilde - 1.0) / ((1.0 + ratio) * gamma_tilde**2)
    beta_sq = t * (t + 2.0) / (1.0 + t) ** 2
    bethe_log = math.log(t * (t + 2.0)) - beta_sq - math.log(2.0 * b)
    bethe_log += correction
    total = prefactor * (
        0.5 * q * bethe_log * (1.0 - 1.0 / ratio**2)
        + (2.0 - q)
        * (
            1.0
            - 1.0 / ratio
            - math.log(ratio) * exchange
            + 0.5 * binding_sq * (ratio - 1.0)
        )
    )
    return prefactor, ratio, binding_sq, exchange, bethe_log, max(total, 0.0)


def conditional_cdf(secondary_energy, incident_energy, shell, terms):
    binding, _, _, q = shell
    prefactor, ratio, binding_sq, exchange, bethe_log, total = terms
    if total <= 0.0:
        return 0.0
    maximum = 0.5 * (incident_energy - binding)
    if secondary_energy >= maximum:
        return 1.0
    w = max(secondary_energy, 0.0) / binding
    c1 = (
        0.5
        * bethe_log
        * q
        * (1.0 - 1.0 / (w + 1.0) ** 2 + 1.0 / (ratio - w) ** 2 - 1.0 / ratio**2)
    )
    c2 = 1.0 - 1.0 / (w + 1.0) + 1.0 / (ratio - w) - 1.0 / ratio + w * binding_sq
    c3 = math.log((w + 1.0) * ratio / (ratio - w)) * exchange
    return np.clip(prefactor * (c1 + (2.0 - q) * (c2 - c3)) / total, 0.0, 1.0)


def expected_statistics(target, incident_energy, cdf_probe):
    shells = SHELLS[target]
    terms = [rbeq_terms(incident_energy, shell) for shell in shells]
    partials = np.array([term[-1] if term is not None else 0.0 for term in terms])
    probabilities = partials / np.sum(partials)

    binding_mean = np.sum(probabilities * shells[:, 0])
    binding_variance = np.sum(probabilities * (shells[:, 0] - binding_mean) ** 2)
    secondary_mean = 0.0
    secondary_second_moment = 0.0
    cdf_at_probe = 0.0
    primary_cosine_mean = 0.0
    primary_cosine_second_moment = 0.0
    secondary_cosine_mean = 0.0
    secondary_cosine_second_moment = 0.0
    for probability, shell, term in zip(probabilities, shells, terms):
        if probability == 0.0:
            continue
        maximum = 0.5 * (incident_energy - shell[0])
        grid = np.linspace(0.0, maximum, 20001)
        cdf = np.array(
            [conditional_cdf(value, incident_energy, shell, term) for value in grid]
        )
        secondary_mean += probability * np.trapezoid(1.0 - cdf, grid)
        secondary_second_moment += probability * np.trapezoid(
            2.0 * grid * (1.0 - cdf), grid
        )
        cdf_at_probe += probability * conditional_cdf(
            cdf_probe, incident_energy, shell, term
        )
        pdf = np.maximum(np.gradient(cdf, grid, edge_order=2), 0.0)
        pdf /= np.trapezoid(pdf, grid)
        available = incident_energy - shell[0]
        primary_energy = available - grid
        primary_cosine = np.sqrt(
            primary_energy
            * (available + 2.0 * MC2_EV)
            / (available * (primary_energy + 2.0 * MC2_EV))
        )
        secondary_kinematic_cosine = np.sqrt(
            grid * (available + 2.0 * MC2_EV) / (available * (grid + 2.0 * MC2_EV))
        )
        free_weight = grid / (grid + shell[0])
        secondary_conditional_mean = free_weight * secondary_kinematic_cosine
        secondary_conditional_second_moment = (
            secondary_conditional_mean**2 + (1.0 - free_weight) ** 2 / 3.0
        )
        primary_cosine_mean += probability * np.trapezoid(primary_cosine * pdf, grid)
        primary_cosine_second_moment += probability * np.trapezoid(
            primary_cosine**2 * pdf, grid
        )
        secondary_cosine_mean += probability * np.trapezoid(
            secondary_conditional_mean * pdf, grid
        )
        secondary_cosine_second_moment += probability * np.trapezoid(
            secondary_conditional_second_moment * pdf, grid
        )
    secondary_variance = secondary_second_moment - secondary_mean**2
    primary_cosine_variance = primary_cosine_second_moment - primary_cosine_mean**2
    secondary_cosine_variance = (
        secondary_cosine_second_moment - secondary_cosine_mean**2
    )
    return (
        binding_mean,
        binding_variance,
        secondary_mean,
        secondary_variance,
        cdf_at_probe,
        primary_cosine_mean,
        primary_cosine_variance,
        secondary_cosine_mean,
        secondary_cosine_variance,
    )


data = np.load("background_mcc_rbeq_results.npz")
for name, target, incident_energy in CASES:
    event_count = int(data[f"{name}_event_count"])
    secondary = data[f"{name}_secondary_energies"]
    primary_cosines = data[f"{name}_primary_cosines"]
    secondary_cosines = data[f"{name}_secondary_cosines"]
    momentum_residual = data[f"{name}_momentum_residual"]
    observed_binding_mean = float(data[f"{name}_binding_mean"])
    outer_binding = np.min(SHELLS[target][:, 0])
    cdf_probe = 0.1 * (incident_energy - outer_binding)
    (
        expected_binding_mean,
        expected_binding_variance,
        expected_secondary_mean,
        expected_secondary_variance,
        expected_cdf,
        expected_primary_cosine,
        expected_primary_cosine_variance,
        expected_secondary_cosine,
        expected_secondary_cosine_variance,
    ) = expected_statistics(target, incident_energy, cdf_probe)

    expected_event_fraction = -math.expm1(-5.0)
    event_standard_error = math.sqrt(
        expected_event_fraction * (1.0 - expected_event_fraction) / PARTICLE_COUNT
    )
    assert abs(event_count / PARTICLE_COUNT - expected_event_fraction) < (
        7.0 * event_standard_error
    )
    assert int(data[f"{name}_ion_count"]) == event_count
    assert primary_cosines.size == event_count
    assert secondary_cosines.size == event_count
    assert np.isfinite(secondary).all()
    assert np.all(secondary >= 0.0)
    assert np.max(secondary) <= 0.5 * (incident_energy - outer_binding)

    binding_standard_error = math.sqrt(expected_binding_variance / event_count)
    assert abs(observed_binding_mean - expected_binding_mean) < (
        7.0 * binding_standard_error + 5.0e-3 * max(1.0, expected_binding_mean)
    )
    secondary_standard_error = math.sqrt(expected_secondary_variance / event_count)
    assert abs(float(np.mean(secondary)) - expected_secondary_mean) < (
        7.0 * secondary_standard_error + 5.0e-3 * max(1.0, expected_secondary_mean)
    )
    observed_cdf = float(np.mean(secondary <= cdf_probe))
    cdf_standard_error = math.sqrt(expected_cdf * (1.0 - expected_cdf) / event_count)
    assert abs(observed_cdf - expected_cdf) < 7.0 * cdf_standard_error + 4.0e-3

    assert np.isfinite(primary_cosines).all()
    assert np.isfinite(secondary_cosines).all()
    assert np.all(np.abs(primary_cosines) <= 1.0)
    assert np.all(np.abs(secondary_cosines) <= 1.0)
    primary_cosine_standard_error = math.sqrt(
        expected_primary_cosine_variance / event_count
    )
    secondary_cosine_standard_error = math.sqrt(
        expected_secondary_cosine_variance / event_count
    )
    assert abs(np.mean(primary_cosines) - expected_primary_cosine) < (
        7.0 * primary_cosine_standard_error + 5.0e-3
    )
    assert abs(np.mean(secondary_cosines) - expected_secondary_cosine) < (
        7.0 * secondary_cosine_standard_error + 5.0e-3
    )
    assert np.linalg.norm(momentum_residual) < 2.0e-5

    print(
        f"{name}: events={event_count}, binding mean="
        f"{observed_binding_mean:.5f}/{expected_binding_mean:.5f} eV, "
        f"secondary mean={np.mean(secondary):.5f}/{expected_secondary_mean:.5f} eV, "
        f"P(W<={cdf_probe:.3g} eV)={observed_cdf:.5f}/{expected_cdf:.5f}, "
        f"<cos1>={np.mean(primary_cosines):.5f}/{expected_primary_cosine:.5f}, "
        f"<cos2>={np.mean(secondary_cosines):.5f}/{expected_secondary_cosine:.5f}, "
        f"|dp|/p0={np.linalg.norm(momentum_residual):.3e}"
    )
