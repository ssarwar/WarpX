#!/usr/bin/env python3
"""Guard against catastrophic many-process or broad-grid MCC regressions."""

import numpy as np

results = np.load("background_mcc_many_attachment_results.npz")
electron_count = int(results["electron_count"])
attachment_events = int(results["attachment_events"])
particle_count = int(results["particle_count"])
process_count = int(results["process_count"])
steps = int(results["steps"])
subcycles = int(results["subcycles"])
product_passes_per_collision_call = int(results["product_passes_per_collision_call"])
initialization_elapsed = float(results["initialization_elapsed"])
elapsed = float(results["elapsed"])
nominal_particles_per_second = float(results["nominal_particles_per_second"])

assert process_count == 64
assert steps == 1
assert subcycles == 1
assert product_passes_per_collision_call == 1
assert electron_count + attachment_events == particle_count
assert 0.004 < attachment_events / particle_count < 0.025

# This is deliberately loose. The purpose is to catch the former scan whose
# cost scaled as (E_max-E_min)/min(dE), not to benchmark shared CI hardware.
assert elapsed < 30.0

print(
    f"{process_count}-process attachment step: events={attachment_events}, "
    f"initialization={initialization_elapsed:.6f} s, elapsed={elapsed:.6f} s, "
    f"nominal throughput={nominal_particles_per_second:.3e} particle-calls/s"
)
