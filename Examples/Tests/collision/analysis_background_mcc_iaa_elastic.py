#!/usr/bin/env python3
"""Verify tabulated elastic angular moments and relativistic molecular recoil."""

import math

import numpy as np
from scipy.constants import c, electron_volt, m_e

data = np.load("background_mcc_iaa_elastic_results.npz")
case_names = data["case_names"]
case_energies = data["case_energies"]
case_masses = data["case_masses"]
expected_moments = data["expected_moments"]
expected_backward_probabilities = data["expected_backward_probabilities"]
particle_count = int(data["particle_count"])
rate_fraction = float(data["rate_fraction"])

electron_rest_energy = m_e * c**2 / electron_volt


def expected_final_energy(incident_energy, neutral_mass, cosine):
    neutral_rest_energy = neutral_mass * c**2 / electron_volt
    sine_sq = np.maximum(1.0 - cosine**2, 0.0)
    recoil_root = np.sqrt(
        np.maximum(
            neutral_rest_energy**2 - electron_rest_energy**2 * sine_sq,
            0.0,
        )
    )
    recoil_factor = (
        (incident_energy + electron_rest_energy) * sine_sq
        + neutral_rest_energy * (1.0 - cosine)
        + cosine
        * electron_rest_energy**2
        * sine_sq
        / (neutral_rest_energy + recoil_root)
    )
    incident_pc_sq = incident_energy * (incident_energy + 2.0 * electron_rest_energy)
    total_energy = incident_energy + electron_rest_energy + neutral_rest_energy
    energy_loss = (
        recoil_factor * incident_pc_sq / (total_energy**2 - incident_pc_sq * cosine**2)
    )
    return incident_energy - energy_loss


assert float(data["initialization_elapsed"]) < 30.0
assert float(data["step_elapsed"]) < 30.0
assert float(data["particle_calls_per_second"]) > 0.0

for index, raw_name in enumerate(case_names):
    name = str(raw_name)
    incident_energy = float(case_energies[index])
    neutral_mass = float(case_masses[index])
    cosines = data[f"{name}_cosines"]
    azimuth_x = data[f"{name}_azimuth_x"]
    azimuth_y = data[f"{name}_azimuth_y"]
    energies = data[f"{name}_energies"]

    expected_event_fraction = rate_fraction * (-math.expm1(-12.0))
    event_standard_error = math.sqrt(
        expected_event_fraction * (1.0 - expected_event_fraction) / particle_count
    )
    assert abs(cosines.size / particle_count - expected_event_fraction) < (
        7.0 * event_standard_error
    )
    assert np.isfinite(cosines).all()
    assert np.isfinite(energies).all()
    assert np.all(np.abs(cosines) <= 1.0 + 2.0e-7)
    assert np.all(energies >= 0.0)

    reference_moments = expected_moments[index]
    observed_moments = np.array([np.mean(cosines**order) for order in range(1, 5)])
    for order, observed_moment in enumerate(observed_moments, start=1):
        expected_moment = float(reference_moments[order - 1])
        expected_square = float(reference_moments[2 * order - 1])
        moment_variance = max(expected_square - expected_moment**2, 0.0)
        moment_tolerance = 7.0 * math.sqrt(moment_variance / cosines.size) + 1.0e-4
        assert abs(float(observed_moment) - expected_moment) < moment_tolerance

    expected_backward = float(expected_backward_probabilities[index])
    observed_backward = float(np.mean(cosines <= 0.0))
    backward_standard_error = math.sqrt(
        expected_backward * (1.0 - expected_backward) / cosines.size
    )
    assert abs(observed_backward - expected_backward) < (
        7.0 * backward_standard_error + 5.0e-5
    )
    assert abs(float(np.mean(azimuth_x))) < 4.0 / math.sqrt(cosines.size)
    assert abs(float(np.mean(azimuth_y))) < 4.0 / math.sqrt(cosines.size)

    recoil_reference = expected_final_energy(incident_energy, neutral_mass, cosines)
    relative_recoil_error = np.abs(energies - recoil_reference) / max(
        incident_energy, 1.0e-12
    )
    assert np.quantile(relative_recoil_error, 0.999) < 2.0e-5
    assert np.max(relative_recoil_error) < 2.0e-4

    incident_pc = math.sqrt(
        incident_energy * (incident_energy + 2.0 * electron_rest_energy)
    )
    outgoing_pc = np.sqrt(energies * (energies + 2.0 * electron_rest_energy))
    recoil_pc_sq = (
        incident_pc**2 + outgoing_pc**2 - 2.0 * incident_pc * outgoing_pc * cosines
    )
    neutral_rest_energy = neutral_mass * c**2 / electron_volt
    neutral_recoil = recoil_pc_sq / (
        np.sqrt(neutral_rest_energy**2 + recoil_pc_sq) + neutral_rest_energy
    )
    conservation_residual = np.abs(incident_energy - energies - neutral_recoil) / max(
        incident_energy, 1.0e-12
    )
    assert np.quantile(conservation_residual, 0.999) < 2.0e-5
    assert np.max(conservation_residual) < 2.0e-4

    print(
        f"{name}: moments={np.array2string(observed_moments, precision=6)}/"
        f"{np.array2string(reference_moments[:4], precision=6)}, "
        f"P(cos(theta)<=0)={observed_backward:.6f}/{expected_backward:.6f}, "
        f"conservation q99.9={np.quantile(conservation_residual, 0.999):.3e}"
    )

print(
    "IAA elastic performance: "
    f"initialization={float(data['initialization_elapsed']):.6f} s, "
    f"step={float(data['step_elapsed']):.6f} s, "
    f"throughput={float(data['particle_calls_per_second']):.3e} particle-calls/s"
)
