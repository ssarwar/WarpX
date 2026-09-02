#!/usr/bin/env python3
"""Guard bounded creation and catastrophic source-kernel regressions."""

import numpy as np

data = np.load("proton_impact_ionization_performance_results.npz")
cell_count = int(data["cell_count"])
particle_count = int(data["particle_count"])
steps = int(data["steps"])
max_products_per_cell = int(data["max_products_per_cell"])
electron_count = int(data["electron_count"])
ion_count = int(data["ion_count"])
initialization_elapsed = float(data["initialization_elapsed"])
step_elapsed = float(data["step_elapsed"])
throughput = float(data["particle_calls_per_second"])

assert cell_count == 512
assert particle_count == 32768
assert steps == 5
assert max_products_per_cell == 8
assert electron_count == ion_count
assert electron_count == cell_count * max_products_per_cell * steps
assert initialization_elapsed < 30.0
assert 0.0 < step_elapsed < 30.0
assert throughput > 0.0

print(
    f"bounded PJG source: pairs={electron_count}, "
    f"initialization={initialization_elapsed:.6f} s, step={step_elapsed:.6f} s, "
    f"throughput={throughput:.3e} particle-calls/s"
)
