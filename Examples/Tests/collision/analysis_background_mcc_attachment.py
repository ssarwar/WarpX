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

# Every created macroparticle must retain the incident electron's weight. Check
# the weight balance directly rather than inferring it only from particle counts.
initial_electron_w = results["initial_electron_w"]
electron_w = results["electron_w"]
assert initial_electron_w.size == PARTICLE_COUNT
assert electron_w.size == electron_count
assert np.all(np.isfinite(initial_electron_w))
assert np.all(initial_electron_w > 0.0)
source_weight = initial_electron_w[0]
assert np.all(initial_electron_w == source_weight)
assert np.all(electron_w == source_weight)

product_weight_keys = [
    "ion_a_w",
    "ion_b_w",
    "negative_a_w",
    "negative_b_w",
]
for key, event_count in zip(product_weight_keys, process_events):
    weights = results[key]
    assert weights.size == event_count
    assert np.all(weights == source_weight)

ion_weight = sum(
    np.sum(results[key], dtype=np.float64) for key in product_weight_keys[:2]
)
negative_weight = sum(
    np.sum(results[key], dtype=np.float64) for key in product_weight_keys[2:]
)
initial_electron_weight = np.sum(initial_electron_w, dtype=np.float64)
final_electron_weight = np.sum(electron_w, dtype=np.float64)
assert negative_weight == float(source_weight) * attachment_events

# Charge is expressed in units of q_e. Initial product particles are excluded
# from both sides because they are unchanged by the MCC step.
charge_residual = (
    -final_electron_weight + ion_weight - negative_weight + initial_electron_weight
)
weight_epsilon = np.finfo(initial_electron_w.dtype).eps
assert abs(charge_residual) < 64.0 * weight_epsilon * initial_electron_weight

# Newly created electrons and ions must receive valid, globally distinct IDs;
# attachment must remove exactly the selected original electron IDs.
initial_electron_ids = results["initial_electron_ids"]
electron_ids = results["electron_ids"]
assert initial_electron_ids.size == PARTICLE_COUNT
assert electron_ids.size == electron_count
assert np.unique(initial_electron_ids).size == initial_electron_ids.size
assert np.unique(electron_ids).size == electron_ids.size

new_electron_ids = np.setdiff1d(electron_ids, initial_electron_ids)
removed_electron_ids = np.setdiff1d(initial_electron_ids, electron_ids)
assert new_electron_ids.size == ionization_events
assert removed_electron_ids.size == attachment_events

created_product_ids = np.concatenate(
    [
        results["ion_a_ids"],
        results["ion_b_ids"],
        results["negative_a_ids"],
        results["negative_b_ids"],
    ]
)
all_created_ids = np.concatenate([new_electron_ids, created_product_ids])
assert np.all(all_created_ids > 0)
assert np.unique(all_created_ids).size == all_created_ids.size
assert np.intersect1d(all_created_ids, initial_electron_ids).size == 0

# The user majorant is 1% above the actual rate. Neutral thermal motion changes
# the 1 keV relative speed by far less than the statistical uncertainty here.
expected_total_fraction = -math.expm1(-OPTICAL_DEPTH) / MAJORANT_SAFETY
observed_total_fraction = total_events / PARTICLE_COUNT
standard_error = math.sqrt(
    expected_total_fraction * (1.0 - expected_total_fraction) / PARTICLE_COUNT
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
