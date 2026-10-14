#!/usr/bin/env python3
"""Validate competing elastic and ionization events from one MCC draw."""

import math

import numpy as np

PARTICLE_COUNT = 131072
OPTICAL_DEPTH = 0.5
ELASTIC_CONDITIONAL_FRACTION = 0.25
IONIZATION_CONDITIONAL_FRACTION = 0.75

results = np.load("background_mcc_selector_results.npz")
electron_count = int(results["electron_count"])
elastic_events = int(results["elastic_events"])
ionization_events = int(results["ionization_events"])
reordered_electron_count = int(results["reordered_electron_count"])
reordered_elastic_events = int(results["reordered_elastic_events"])
reordered_ionization_events = int(results["reordered_ionization_events"])

expected_total_fraction = -math.expm1(-OPTICAL_DEPTH)
standard_error = math.sqrt(
    expected_total_fraction * (1.0 - expected_total_fraction) / PARTICLE_COUNT
)
expected_elastic = (
    ELASTIC_CONDITIONAL_FRACTION * expected_total_fraction * PARTICLE_COUNT
)
expected_ionization = (
    IONIZATION_CONDITIONAL_FRACTION * expected_total_fraction * PARTICLE_COUNT
)
elastic_standard_error = math.sqrt(
    PARTICLE_COUNT
    * ELASTIC_CONDITIONAL_FRACTION
    * expected_total_fraction
    * (1.0 - ELASTIC_CONDITIONAL_FRACTION * expected_total_fraction)
)
ionization_standard_error = math.sqrt(
    PARTICLE_COUNT
    * IONIZATION_CONDITIONAL_FRACTION
    * expected_total_fraction
    * (1.0 - IONIZATION_CONDITIONAL_FRACTION * expected_total_fraction)
)


def validate_sample(electron_count, elastic_events, ionization_events):
    total_events = elastic_events + ionization_events
    observed_total_fraction = total_events / PARTICLE_COUNT

    # One ion and one electron are created for each selected ionization event.
    assert electron_count == PARTICLE_COUNT + ionization_events
    assert total_events <= PARTICLE_COUNT
    assert abs(observed_total_fraction - expected_total_fraction) < 7.0 * standard_error
    assert abs(elastic_events - expected_elastic) < 7.0 * elastic_standard_error
    assert (
        abs(ionization_events - expected_ionization) < 7.0 * ionization_standard_error
    )
    return total_events, observed_total_fraction


total_events, observed_total_fraction = validate_sample(
    electron_count, elastic_events, ionization_events
)
reordered_total_events, reordered_total_fraction = validate_sample(
    reordered_electron_count,
    reordered_elastic_events,
    reordered_ionization_events,
)

# Reversing the process dictionary must not change total or per-process
# statistics beyond the uncertainty of two independent samples.
assert abs(total_events - reordered_total_events) < 7.0 * math.sqrt(
    2.0 * PARTICLE_COUNT * expected_total_fraction * (1.0 - expected_total_fraction)
)
assert (
    abs(elastic_events - reordered_elastic_events)
    < 7.0 * math.sqrt(2.0) * elastic_standard_error
)
assert (
    abs(ionization_events - reordered_ionization_events)
    < 7.0 * math.sqrt(2.0) * ionization_standard_error
)

print(
    f"events: elastic={elastic_events}, ionization={ionization_events}, "
    f"total_fraction={observed_total_fraction:.6f}; reordered: "
    f"elastic={reordered_elastic_events}, ionization={reordered_ionization_events}, "
    f"total_fraction={reordered_total_fraction:.6f}"
)
