#!/usr/bin/env python3
"""Validate the reduced-mass dependence of the automatic MCC majorant."""

import math

import numpy as np

data = np.load("background_mcc_nonrelativistic_majorant_results.npz")
particle_count = int(data["particle_count"])
optical_depth = float(data["optical_depth"])
event_count = int(data["event_count"])

expected_probability = -math.expm1(-optical_depth)
expected_count = particle_count * expected_probability
standard_deviation = math.sqrt(
    particle_count * expected_probability * (1.0 - expected_probability)
)
assert abs(event_count - expected_count) < 7.0 * standard_deviation

# Equal-mass backward elastic scattering transfers the projectile momentum to
# the stationary target. This also verifies the event classifier used above.
assert float(data["maximum_collided_speed"]) < (
    2.0e-5 * float(data["relative_proper_speed"])
)

print(
    f"nonrelativistic automatic majorant: events={event_count}, "
    f"expected={expected_count:.1f}, sigma={standard_deviation:.1f}"
)
