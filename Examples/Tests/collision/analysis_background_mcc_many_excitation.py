#!/usr/bin/env python3
"""Verify many-channel excitation losses and shared-DCS angular moments."""

import math

import numpy as np

data = np.load("background_mcc_many_excitation_results.npz")
process_count = int(data["process_count"])
particle_count = int(data["particle_count"])
optical_depth = float(data["optical_depth"])
anisotropy = float(data["anisotropy"])
first_energy_loss = float(data["first_energy_loss"])
energy_loss_step = float(data["energy_loss_step"])
cosines = data["cosines"]
inferred_losses = data["inferred_losses"]

assert process_count == 64
expected_event_fraction = -math.expm1(-optical_depth)
event_standard_error = math.sqrt(
    expected_event_fraction * (1.0 - expected_event_fraction) / particle_count
)
assert abs(cosines.size / particle_count - expected_event_fraction) < (
    7.0 * event_standard_error
)
assert cosines.size == inferred_losses.size
assert np.isfinite(cosines).all()
assert np.isfinite(inferred_losses).all()

expected_moments = np.array([anisotropy / 3.0, 1.0 / 3.0, anisotropy / 5.0, 1.0 / 5.0])
expected_higher_moments = np.array([1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0])
observed_moments = np.array([np.mean(cosines**order) for order in range(1, 5)])
for order, observed in enumerate(observed_moments):
    variance = expected_higher_moments[order] - expected_moments[order] ** 2
    tolerance = 7.0 * math.sqrt(variance / cosines.size) + 1.0e-4
    assert abs(observed - expected_moments[order]) < tolerance

channel_coordinate = (inferred_losses - first_energy_loss) / energy_loss_step
channel_index = np.rint(channel_coordinate).astype(int)
quantization_error = np.abs(channel_coordinate - channel_index) * energy_loss_step
assert np.quantile(quantization_error, 0.999) < 2.0e-4
assert np.max(quantization_error) < 2.0e-3
assert np.all((channel_index >= 0) & (channel_index < process_count))

channel_counts = np.bincount(channel_index, minlength=process_count)
expected_channel_count = expected_event_fraction * particle_count / process_count
channel_probability = expected_event_fraction / process_count
channel_standard_error = math.sqrt(
    particle_count * channel_probability * (1.0 - channel_probability)
)
assert np.max(np.abs(channel_counts - expected_channel_count)) < (
    7.0 * channel_standard_error
)

assert float(data["initialization_elapsed"]) < 30.0
assert float(data["step_elapsed"]) < 30.0
assert float(data["particle_calls_per_second"]) > 0.0

print(
    f"{process_count} excitation channels: events={cosines.size}, "
    f"moments={np.array2string(observed_moments, precision=6)}/"
    f"{np.array2string(expected_moments, precision=6)}, "
    f"maximum loss error={np.max(quantization_error):.3e} eV, "
    f"initialization={float(data['initialization_elapsed']):.6f} s, "
    f"step={float(data['step_elapsed']):.6f} s, "
    f"throughput={float(data['particle_calls_per_second']):.3e} particle-calls/s"
)
