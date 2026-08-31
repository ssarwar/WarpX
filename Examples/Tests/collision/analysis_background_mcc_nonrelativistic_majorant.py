#!/usr/bin/env python3
"""Validate interior and high-energy-tail automatic MCC majorants."""

import math

import numpy as np

data = np.load("background_mcc_nonrelativistic_majorant_results.npz")
particle_count = int(data["particle_count"])
optical_depth = float(data["optical_depth"])
peak_event_count = int(data["peak_event_count"])

expected_probability = -math.expm1(-optical_depth)
expected_count = particle_count * expected_probability
standard_deviation = math.sqrt(
    particle_count * expected_probability * (1.0 - expected_probability)
)
assert abs(peak_event_count - expected_count) < 7.0 * standard_deviation

# Equal-mass backward elastic scattering transfers the projectile momentum to
# the stationary target. This also verifies the event classifier used above.
assert float(data["peak_maximum_collided_speed"]) < (
    2.0e-5 * float(data["peak_relative_proper_speed"])
)

# ScatteringProcess holds the last tabulated cross section constant. The
# automatic majorant must therefore include its speed-of-light supremum even
# for a non-electron projectile. At u=100c, the accepted event probability is
# the majorant-event probability times the ordinary-speed ratio v/c.
tail_probability = expected_probability * float(data["tail_collision_speed_over_c"])
tail_expected_count = particle_count * tail_probability
tail_standard_deviation = math.sqrt(
    particle_count * tail_probability * (1.0 - tail_probability)
)
tail_event_count = int(data["tail_event_count"])
assert abs(tail_event_count - tail_expected_count) < 7.0 * tail_standard_deviation
assert float(data["tail_maximum_collided_speed"]) < (
    2.0e-5 * float(data["tail_relative_proper_speed"])
)

print(
    f"non-electron automatic majorant: peak_events={peak_event_count}, "
    f"peak_expected={expected_count:.1f}, tail_events={tail_event_count}, "
    f"tail_expected={tail_expected_count:.1f}"
)
