# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Bounded rigid-rotor basis for offline high-energy quadratures."""

import numpy as np
from hybrid_reference import cg_array
from rotation_reference import ROTATION, converged_j, populations, weights


def make_basis(target, temperature=300):
    b = ROTATION[target]
    dissociation = {"N2": 9.759, "O2": 5.116}[target]
    ground = 0 if target == "N2" else 1
    bound = int(
        np.floor((np.sqrt(1 + 4 * (dissociation / b + ground * (ground + 1))) - 1) / 2)
    )
    if target == "O2":
        bound -= bound % 2 == 0
    maximum = converged_j(target, temperature)
    if maximum > bound:
        raise ValueError("Rotational bath extends beyond the bounded rigid-rotor model")
    population = populations(target, temperature, maximum)[0]
    initial_states = np.flatnonzero(weights(target, np.arange(maximum + 1)))
    initial, final = np.array(
        [(i, f) for i in initial_states for f in range(i % 2, bound + 1, 2)]
    ).T
    ranks = np.arange(0, bound + maximum + 1, 2)
    coefficients = cg_array(initial, final, ranks) * population[initial, None]
    loss = b * (final * (final + 1) - initial * (initial + 1))
    return {
        "initial": initial,
        "final": final,
        "ranks": ranks,
        "coefficients": coefficients,
        "loss": loss,
        "bound": bound,
    }
