#!/usr/bin/env python3
"""Guard against catastrophic many-process or broad-grid MCC regressions."""

import numpy as np

PARTICLE_COUNT = 32768

results = np.load("background_mcc_many_attachment_results.npz")
electron_count = int(results["electron_count"])
attachment_events = int(results["attachment_events"])
elapsed = float(results["elapsed"])

assert electron_count + attachment_events == PARTICLE_COUNT
assert 0.004 < attachment_events / PARTICLE_COUNT < 0.025

# This is deliberately loose. The purpose is to catch the former scan whose
# cost scaled as (E_max-E_min)/min(dE), not to benchmark shared CI hardware.
assert elapsed < 30.0

print(
    f"64-process attachment step: events={attachment_events}, "
    f"elapsed={elapsed:.6f} s"
)
