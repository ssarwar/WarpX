# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Regress the implicit reverse contribution using real molecular kernels."""

import numpy as np
from hybrid_reference import Hybrid
from reference_refinement import inserted_moments

for target, temperature, maximum, bounds in [
    ("N2", 100, 20, (1, 1.25)),
    ("O2", 300, 1000, (200, 220)),
]:
    model = Hybrid(target, temperature, maximum, 1, 48, angular_resolution=4).solve()
    indices = np.flatnonzero((model.energy > bounds[0]) & (model.energy < bounds[1]))
    self_coupled = 0
    worst = 0.0
    for i in indices:
        energy = model.energy[i]
        self_coupled += np.count_nonzero(model.future_weights(energy, i)[2]) > 0
        independently_solved, _ = inserted_moments(model, energy)
        current = model.angular_moments(energy)
        exact = independently_solved @ model.measure
        error = abs(current @ model.measure - exact) / np.maximum(abs(exact), 1e-100)
        worst = max(worst, error.max())
    assert self_coupled > 0, "The regression must exercise implicit reverse terms"
    assert worst < 1e-10, (target, worst)
    print(target, "self-coupled rows", self_coupled, "maximum error", worst, flush=True)
