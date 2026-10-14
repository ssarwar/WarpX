# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent event probabilities and signed energy-transfer moments."""

import numpy as np

data = np.load("background_mcc_rotation_results.npz")
count = int(data["particles"])
assert int(data["steps"]) == 1
for key in data:
    if not key.startswith("e_"):
        continue
    mean, second, reference, reference_second, events, probability = data[key]
    error = np.sqrt(max(second, reference_second) / count)
    assert abs(mean - reference) < 7 * error + 2e-7, (key, mean, reference)
    assert (
        abs(events - probability)
        < 7 * np.sqrt(probability * (1 - probability) / count) + 1e-5
    )
    if "_0_5" not in key:
        assert mean < 0, (key, "cold electrons must gain energy")
    else:
        assert mean > 0, (key, "hot electrons must lose energy")
print(
    "Combined MCC: event rates, zero-speed gain, signed moments and ordinary channels passed"
)
print("Startup seconds:", float(data["startup_seconds"]))
print("Complete timestep seconds:", float(data["step_seconds"]))
