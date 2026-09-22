# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Printed proton PJG with the 1977 erratum, for comparison only.

Negative values are preserved, not silently clipped into a different model.
"""

import numpy as np
from pjg_model import kinematics
from reference import BETHE_CONSTANT, ELECTRON_REST_ENERGY, PROTON_REST_ENERGY
from scipy.integrate import simpson
from target_parameters import TARGETS


def legacy_cutoff(energy, variant="printed"):
    if variant != "printed":
        raise ValueError("Only printed PJG with the 1977 erratum is retained")
    energy = np.asarray(energy)
    m, mass = ELECTRON_REST_ENERGY, PROTON_REST_ENERGY
    return energy * (energy + 2 * mass) / (energy + m + mass / m * (energy + mass))


def legacy_sdcs(target, energy, secondary, variant="printed"):
    """Original proton Eq. (16), with the erratum's ordinary addition."""
    if variant != "printed":
        raise ValueError("Unknown legacy model")
    energy, t = np.broadcast_arrays(energy, secondary)
    p = TARGETS[target]
    mass = PROTON_REST_ENERGY
    gamma, beta2, ee, maximum = kinematics(energy)
    cutoff = legacy_cutoff(energy, variant)
    width = p.gamma_s + p.gamma_numerator / (ee + p.gamma_denominator)
    center = p.center_s - p.center_numerator / (ee + p.center_denominator)
    b0, broad_center, delta = {"N2": (0.029, 53.3, 84), "O2": (0.030, 68.3, 132.1)}[
        target
    ]
    reduction = b0 * (np.log(ee / 8239) ** 2 + 1.035)
    j, amplitude = p.j, p.k
    distortion = 1 / (1 + (j / ee) ** p.power)
    logarithm = np.sum(
        np.array(p.fractions)
        * (
            np.log(
                4
                * ee[..., None]
                * gamma[..., None] ** 2
                * np.array(p.bethe_constants)
                / p.thresholds
                + np.e
            )
            - beta2[..., None]
        ),
        axis=-1,
    )
    soft = (
        amplitude
        * width**2
        * logarithm
        * (
            1 / ((t - center) ** 2 + width**2)
            - reduction / ((t - broad_center) ** 2 + p.broad_width**2)
        )
    )
    loss = t[..., None] + p.thresholds
    remainder = 1 / (4 * (energy[..., None] + 2 * mass) ** 2)
    remainder = remainder - 1 / ((cutoff[..., None] + p.thresholds + delta) * loss)
    hard = p.electrons * BETHE_CONSTANT * np.sum(p.fractions * remainder, axis=-1)
    return np.where((t >= 0) & (t <= cutoff), distortion / ee * (soft + hard), 0)


def moments(target, energies, variant="printed", points=4097):
    energies = np.atleast_1d(energies)
    x = np.linspace(0, 1, points)[None, :] * np.log1p(
        legacy_cutoff(energies, variant)[:, None]
    )
    t = np.expm1(x)
    # Roundoff in expm1(log1p(cutoff)) must not remove the final endpoint.
    t = np.minimum(t, legacy_cutoff(energies, variant)[:, None])
    weighted = legacy_sdcs(target, energies[:, None], t, variant) * np.exp(x)
    total = simpson(weighted, x=x, axis=-1)
    return total, simpson(t * weighted, x=x, axis=-1) / total


def bhabha(target, energy, secondary):
    """Independent, unfactored free stationary-electron cross section."""
    _, beta2, ee, maximum = kinematics(energy)
    t = np.asarray(secondary)
    bracket = 1 - beta2 * t / maximum + t**2 / (2 * (energy + PROTON_REST_ENERGY) ** 2)
    return TARGETS[target].electrons * BETHE_CONSTANT / ee / t**2 * bracket
