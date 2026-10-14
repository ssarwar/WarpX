# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Canonical levels, degeneracies, and phase-space helpers for offline calculations."""

import numpy as np
from scipy.special import gammaln

KB = 8.617333262145e-5
REST = 510998.95069
C = 299792458.0
A0 = 5.29177210544e-11
HARTREE = 27.211386245988
MASSES = {"N2": 28.0134 * 1.66053906660e-27, "O2": 31.9988 * 1.66053906660e-27}
QE = 1.602176634e-19
ROTATION = {"N2": 0.0002477204284695341, "O2": 0.00017828927734694198}


def weights(target, j):
    j = np.asarray(j)
    spin = np.where(j % 2 == 0, 6, 3) if target == "N2" else j % 2
    return (2 * j + 1) * spin


def populations(target, temperature, maximum_j):
    if temperature < 0 or not np.isfinite(temperature):
        raise ValueError("Rotational temperature must be finite and nonnegative")
    j = np.arange(maximum_j + 1)
    ground = 0 if target == "N2" else 1
    p = np.zeros(j.size)
    if temperature == 0:
        p[ground] = 1
        return p, 0.0
    a = ROTATION[target] / (KB * temperature)
    unnormalized = weights(target, j) * np.exp(
        -a * (j * (j + 1) - ground * (ground + 1))
    )
    partition = unnormalized.sum()
    # Upper bound from the first omitted term plus the decreasing continuum
    # integral. The factor six also bounds the O2 spin weight.
    first = maximum_j + 1
    tail = (
        6
        * (2 * first + 1 + 1 / a)
        * np.exp(-a * (first * (first + 1) - ground * (ground + 1)))
    )
    if a * (2 * first + 1) ** 2 < 2:
        tail = np.inf
    return unnormalized / partition, tail / partition


def converged_j(target, temperature, tolerance=1e-10):
    for maximum in range(8, 4097, 8):
        if populations(target, temperature, maximum)[1] < tolerance:
            return maximum
    raise ValueError("Rotational population did not converge below J=4096")


def cg_squared(initial, rank, final):
    """Squared C(initial,0;rank,0|final,0), from the factorial 3j identity."""
    if abs(initial - rank) > final or initial + rank < final:
        return 0.0
    total = initial + rank + final
    if total % 2:
        return 0.0
    half = total // 2
    differences = np.array([half - initial, half - rank, half - final])
    return (2 * final + 1) * np.exp(
        2 * (gammaln(half + 1) - gammaln(differences + 1).sum())
        + gammaln(2 * differences + 1).sum()
        - gammaln(total + 2)
    )


def speed(energy):
    energy = np.asarray(energy)
    return C * np.sqrt(energy * (energy + 2 * REST)) / (energy + REST)


def threshold(target, initial, final):
    """Canonical internal-energy threshold; molecular recoil shifts are neglected."""
    return ROTATION[target] * (final * (final + 1) - initial * (initial + 1))


def phase_rate(energy, loss):
    """v(E)*p_out/p_in in the heavy-target limit, including finite gain at E=0."""
    outgoing = np.maximum(np.asarray(energy) - loss, 0)
    return C * np.sqrt(outgoing * (outgoing + 2 * REST)) / (energy + REST)


def thermal_rates(target, temperature, maximum_j, transitions, rates):
    p, tail = populations(target, temperature, maximum_j)
    if tail > 1e-10:
        raise ValueError(f"Unresolved rotational population: {tail}")
    factors = np.array([1.0, *(p[initial] for initial, _ in transitions)])
    return rates * factors[None, None, :]
