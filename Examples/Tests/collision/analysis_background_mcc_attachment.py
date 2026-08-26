#!/usr/bin/env python3
"""Validate multiple ionization and two- and three-body attachment channels."""

import math

import numpy as np

PARTICLE_COUNT = 65536
OPTICAL_DEPTH = 0.4
MAJORANT_SAFETY = 1.01

results = np.load("background_mcc_attachment_results.npz")
electron_count = int(results["electron_count"])
process_events = np.array(
    [
        int(results["ion_a_events"]),
        int(results["ion_b_events"]),
        int(results["negative_a_events"]),
        int(results["negative_b_events"]),
    ]
)

ionization_events = int(process_events[:2].sum())
attachment_events = int(process_events[2:].sum())
total_events = int(process_events.sum())

# One extra electron is created by ionization, while attachment replaces one
# electron by one negative ion. Each original electron can select at most one
# process during this MCC substep.
assert electron_count == PARTICLE_COUNT + ionization_events - attachment_events
assert total_events <= PARTICLE_COUNT

# The user majorant is 1% above the actual rate. Neutral thermal motion changes
# the 1 keV relative speed by far less than the statistical uncertainty here.
expected_total_fraction = -math.expm1(-OPTICAL_DEPTH) / MAJORANT_SAFETY
observed_total_fraction = total_events / PARTICLE_COUNT
standard_error = math.sqrt(
    expected_total_fraction
    * (1.0 - expected_total_fraction)
    / PARTICLE_COUNT
)
assert abs(observed_total_fraction - expected_total_fraction) < 7.0 * standard_error

expected_each = 0.25 * expected_total_fraction * PARTICLE_COUNT
process_probability = 0.25 * expected_total_fraction
process_standard_error = math.sqrt(
    PARTICLE_COUNT * process_probability * (1.0 - process_probability)
)
for count in process_events:
    assert abs(count - expected_each) < 7.0 * process_standard_error

# The m^5 table multiplied by the supplied third-body density is deliberately
# equivalent to the m^2 table. This also catches premature single-precision
# underflow of the unscaled m^5 values.
difference_error = math.sqrt(process_events[2] + process_events[3])
assert abs(process_events[2] - process_events[3]) < 7.0 * difference_error

# Product ions inherit the sampled neutral velocity. At this temperature the
# sampled neutral is non-relativistic, so proper and physical velocity agree to
# far beyond this statistical test's precision.
expected_std = float(results["neutral_velocity_std"])
for species in ["ion_a", "ion_b", "negative_a", "negative_b"]:
    for direction in ["ux", "uy", "uz"]:
        values = results[f"{species}_{direction}"]
        assert values.size > 1000
        assert abs(np.mean(values)) < 0.08 * expected_std
        assert abs(np.std(values) / expected_std - 1.0) < 0.08

print(
    "events: "
    f"ion_a={process_events[0]}, ion_b={process_events[1]}, "
    f"attachment_m2={process_events[2]}, attachment_m5={process_events[3]}, "
    f"total_fraction={observed_total_fraction:.6f}"
)
