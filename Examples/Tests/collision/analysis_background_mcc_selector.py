#!/usr/bin/env python3
"""Validate competing elastic and ionization events from one MCC draw."""

import math

import numpy as np

PARTICLE_COUNT = 131072
OPTICAL_DEPTH = 0.5

results = np.load("background_mcc_selector_results.npz")
electron_count = int(results["electron_count"])
elastic_events = int(results["elastic_events"])
ionization_events = int(results["ionization_events"])

total_events = elastic_events + ionization_events
expected_total_fraction = -math.expm1(-OPTICAL_DEPTH)
observed_total_fraction = total_events/PARTICLE_COUNT

# One ion and one electron are created for each selected ionization event.
assert electron_count == PARTICLE_COUNT + ionization_events
assert total_events <= PARTICLE_COUNT

standard_error = math.sqrt(
    expected_total_fraction*(1.0 - expected_total_fraction)/PARTICLE_COUNT
)
assert abs(observed_total_fraction - expected_total_fraction) < 7.0*standard_error

expected_each = 0.5*expected_total_fraction*PARTICLE_COUNT
process_standard_error = math.sqrt(
    PARTICLE_COUNT*0.5*expected_total_fraction
    * (1.0 - 0.5*expected_total_fraction)
)
assert abs(elastic_events - expected_each) < 7.0*process_standard_error
assert abs(ionization_events - expected_each) < 7.0*process_standard_error

print(
    f"events: elastic={elastic_events}, ionization={ionization_events}, "
    f"total_fraction={observed_total_fraction:.6f}"
)
