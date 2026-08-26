#!/usr/bin/env python3
"""Validate thermal-target rates and collision-weighted neutral velocities."""

import math

import numpy as np

PARTICLE_COUNT = 131072
MAJORANT_THERMAL_STD = 8.0
MAJORANT_OPTICAL_DEPTH = 0.5
COLLISION_SUBCYCLES = 16

results = np.load("background_mcc_thermal_attachment_results.npz")
electron_count = int(results["electron_count"])
process_events = np.array(
    [int(results["negative_a_events"]), int(results["negative_b_events"])]
)
attachment_events = int(process_events.sum())

assert electron_count + attachment_events == PARTICLE_COUNT

# For a stationary electron and a Maxwellian target, <g> = sqrt(8/pi)*v_th.
# Each MCC object accepts the following fraction of its current electrons per
# substep. The second object then acts on the survivors from the first object.
mean_rate_ratio = math.sqrt(8.0 / math.pi) / MAJORANT_THERMAL_STD
substep_fraction = (
    -math.expm1(-MAJORANT_OPTICAL_DEPTH / COLLISION_SUBCYCLES)
    * math.sqrt(8.0 / math.pi)
    / MAJORANT_THERMAL_STD
)
survival_per_object = (1.0 - substep_fraction) ** COLLISION_SUBCYCLES
expected_process_fractions = np.array(
    [1.0 - survival_per_object, survival_per_object * (1.0 - survival_per_object)]
)
expected_fraction = float(expected_process_fractions.sum())
observed_fraction = attachment_events / PARTICLE_COUNT
standard_error = math.sqrt(
    expected_fraction * (1.0 - expected_fraction) / PARTICLE_COUNT
)
assert abs(observed_fraction - expected_fraction) < 7.0 * standard_error
for count, expected_process_fraction in zip(process_events, expected_process_fractions):
    process_error = math.sqrt(
        PARTICLE_COUNT * expected_process_fraction * (1.0 - expected_process_fraction)
    )
    assert abs(count - PARTICLE_COUNT * expected_process_fraction) < 7.0 * process_error

# Sixteen subcycles are close to the continuous-time summed-rate limit, unlike
# a single one-event-per-call update at this optical depth.
continuous_fraction = -math.expm1(-2.0 * MAJORANT_OPTICAL_DEPTH * mean_rate_ratio)
assert abs(expected_fraction - continuous_fraction) < 0.005

# Collision selection weights the Maxwellian by speed. Isotropy gives
# <v_i^2>_collision = (1/3)<v^3>/<v> = (4/3)*v_th^2.
thermal_std = float(results["neutral_velocity_std"])
expected_component_std = math.sqrt(4.0 / 3.0) * thermal_std
for species, event_count in zip(["negative_a", "negative_b"], process_events):
    for direction in ["ux", "uy", "uz"]:
        values = results[f"{species}_{direction}"]
        assert values.size == event_count
        assert np.all(np.isfinite(values))
        assert abs(np.mean(values)) < 0.1 * thermal_std
        assert abs(np.std(values) / expected_component_std - 1.0) < 0.08

print(
    f"thermal attachment events={process_events.tolist()}, "
    f"fraction={observed_fraction:.6f}, expected={expected_fraction:.6f}"
)
