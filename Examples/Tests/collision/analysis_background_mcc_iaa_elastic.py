#!/usr/bin/env python3
"""Verify tabulated elastic/excitation angles, losses and molecular recoil."""

import math

import numpy as np
from scipy.constants import c, electron_volt, m_e

data = np.load("background_mcc_iaa_elastic_results.npz")
case_names = data["case_names"]
case_energies = data["case_energies"]
case_masses = data["case_masses"]
case_process_types = data["case_process_types"]
case_energy_losses = data["case_energy_losses"]
case_rate_fractions = data["case_rate_fractions"]
case_screening_radii = data["case_screening_radii"]
expected_moments = data["expected_moments"]
expected_backward_probabilities = data["expected_backward_probabilities"]
particle_count = int(data["particle_count"])
scattering_angle_model = str(data["scattering_angle_model"])

electron_rest_energy = m_e * c**2 / electron_volt


def expected_final_energy(incident_energy, neutral_mass, energy_loss, cosine):
    neutral_rest_energy = neutral_mass * c**2 / electron_volt
    incident_total_energy = incident_energy + electron_rest_energy
    incident_pc_sq = incident_energy * (incident_energy + 2.0 * electron_rest_energy)
    total_energy = incident_total_energy + neutral_rest_energy
    phase_space_1 = np.maximum(
        2.0 * neutral_rest_energy * (incident_energy - energy_loss)
        - 2.0 * electron_rest_energy * energy_loss
        - energy_loss**2,
        0.0,
    )
    phase_space_2 = np.maximum(
        4.0 * neutral_rest_energy * electron_rest_energy
        + 2.0 * neutral_rest_energy * incident_energy
        - 2.0 * energy_loss * (neutral_rest_energy - electron_rest_energy)
        - energy_loss**2,
        0.0,
    )
    kallen_quarter = 0.25 * phase_space_1 * phase_space_2
    angular_scale = incident_pc_sq * electron_rest_energy**2
    radicand = np.maximum(kallen_quarter - angular_scale * (1.0 - cosine**2), 0.0)
    coefficient = total_energy**2 - incident_pc_sq * cosine**2
    linear_term = (
        electron_rest_energy**2
        + neutral_rest_energy * (incident_total_energy - energy_loss)
        - 0.5 * energy_loss**2
    )
    outgoing_total_energy = (
        total_energy * linear_term
        + np.sqrt(incident_pc_sq) * cosine * np.sqrt(radicand)
    ) / coefficient
    return np.maximum(outgoing_total_energy - electron_rest_energy, 0.0)


assert float(data["initialization_elapsed"]) < 30.0
assert float(data["step_elapsed"]) < 30.0
assert float(data["particle_calls_per_second"]) > 0.0

for index, raw_name in enumerate(case_names):
    name = str(raw_name)
    incident_energy = float(case_energies[index])
    neutral_mass = float(case_masses[index])
    process_type = str(case_process_types[index])
    energy_loss = float(case_energy_losses[index])
    rate_fraction = float(case_rate_fractions[index])
    screening_radius = float(case_screening_radii[index])
    cosines = data[f"{name}_cosines"]
    deflections = data[f"{name}_deflections"]
    azimuth_x = data[f"{name}_azimuth_x"]
    azimuth_y = data[f"{name}_azimuth_y"]
    energies = data[f"{name}_energies"]

    expected_event_fraction = rate_fraction * (-math.expm1(-12.0))
    event_standard_error = math.sqrt(
        expected_event_fraction * (1.0 - expected_event_fraction) / particle_count
    )
    assert abs(cosines.size / particle_count - expected_event_fraction) < (
        7.0 * event_standard_error + 1.0 / particle_count
    )
    if expected_event_fraction == 0.0:
        assert cosines.size == 0
        print(f"{name}: no events below the recoil-shifted threshold")
        continue
    assert np.isfinite(cosines).all()
    assert np.isfinite(energies).all()
    assert np.all(np.abs(cosines) <= 1.0 + 2.0e-7)
    assert np.all(energies >= 0.0)

    if (
        scattering_angle_model == "IAA"
        and screening_radius > 0.0
        and incident_energy >= 1.0e4
    ):
        fine_structure = 7.2973525693e-3
        tau = incident_energy / electron_rest_energy
        k_sq = tau * (tau + 2.0) / fine_structure**2
        eta = 1.0 / (4.0 * screening_radius**2 * k_sq)
        reconstructed_draws = deflections * (1.0 + eta) / (
            deflections + 2.0 * eta
        )
        assert np.isfinite(reconstructed_draws).all()
        assert np.all(reconstructed_draws >= -2.0e-6)
        assert np.all(reconstructed_draws <= 1.0 + 2.0e-6)
        sample_size = reconstructed_draws.size
        mean_standard_error = math.sqrt(1.0 / (12.0 * sample_size))
        second_moment_standard_error = math.sqrt((4.0 / 45.0) / sample_size)
        assert abs(float(np.mean(reconstructed_draws)) - 0.5) < (
            7.0 * mean_standard_error + 2.0e-5
        )
        assert abs(float(np.mean(reconstructed_draws**2)) - 1.0 / 3.0) < (
            7.0 * second_moment_standard_error + 2.0e-5
        )
        ordered_draws = np.sort(np.clip(reconstructed_draws, 0.0, 1.0))
        ranks = np.arange(1, sample_size + 1, dtype=np.float64) / sample_size
        ks_distance = max(
            float(np.max(ranks - ordered_draws)),
            float(np.max(ordered_draws - (ranks - 1.0 / sample_size))),
        )
        assert ks_distance < 4.0 / math.sqrt(sample_size)

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

    recoil_reference = expected_final_energy(
        incident_energy, neutral_mass, energy_loss, cosines
    )
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
    excited_neutral_rest_energy = neutral_mass * c**2 / electron_volt + energy_loss
    neutral_recoil = recoil_pc_sq / (
        np.sqrt(excited_neutral_rest_energy**2 + recoil_pc_sq)
        + excited_neutral_rest_energy
    )
    inferred_discrete_loss = incident_energy - energies - neutral_recoil
    conservation_residual = np.abs(inferred_discrete_loss - energy_loss) / max(
        incident_energy, energy_loss, 1.0e-12
    )
    assert np.quantile(conservation_residual, 0.999) < 2.0e-5
    assert np.max(conservation_residual) < 2.0e-4
    if process_type == "excitation" and incident_energy < 1.0e6:
        absolute_loss_error = np.abs(inferred_discrete_loss - energy_loss)
        assert np.quantile(absolute_loss_error, 0.999) < 2.0e-4
        assert np.max(absolute_loss_error) < 2.0e-3

    print(
        f"{name}: moments={np.array2string(observed_moments, precision=6)}/"
        f"{np.array2string(reference_moments[:4], precision=6)}, "
        f"P(cos(theta)<=0)={observed_backward:.6f}/{expected_backward:.6f}, "
        f"loss={energy_loss:g} eV, "
        f"conservation q99.9={np.quantile(conservation_residual, 0.999):.3e}"
    )

print(
    "IAA differential-scattering performance: "
    f"initialization={float(data['initialization_elapsed']):.6f} s, "
    f"step={float(data['step_elapsed']):.6f} s, "
    f"throughput={float(data['particle_calls_per_second']):.3e} particle-calls/s"
)
