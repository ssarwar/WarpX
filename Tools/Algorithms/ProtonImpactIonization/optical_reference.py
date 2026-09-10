# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Audit evaluated optical data without equating absorption with an SDCS.

The NIFS data include non-ionizing absorption and core excitation. Their
integral moments test the target response, not an inclusive electron-yield
sum rule. The optional PJG energy-loss response is a diagnostic of its
existing fixed channel allocation; it is not an experimental photoelectron
spectrum. No coefficients are fitted or production tables modified here.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from pjg_model import optical_coefficient
from reference import RYDBERG
from target_parameters import TARGETS


def continuum_integral(data, function, lower=0, upper=np.inf, interpolation="linear"):
    """Integrate each tabulated interval, excluding separately tabulated bands.

    Linear interpolation is the primary result. The log-log alternative
    estimates interpolation sensitivity on the coarse high-energy mesh; it
    is not an experimental uncertainty. Neither extrapolates the spectrum.
    """
    table = np.asarray(data["continuum"])
    x, y = table[:, :2].T
    left = np.maximum(x[:-1], lower)
    right = np.minimum(x[1:], upper)
    valid = right > left
    for start, end in data["continuum_gaps_eV"]:
        valid &= ~((x[:-1] < end) & (x[1:] > start))
    nodes, weights = leggauss(16)
    left, right = left[valid], right[valid]
    energy = left[:, None] + (right - left)[:, None] * (nodes + 1) / 2
    x0, x1 = x[:-1][valid, None], x[1:][valid, None]
    y0, y1 = y[:-1][valid, None], y[1:][valid, None]
    density = y0 + (y1 - y0) * (energy - x0) / (x1 - x0)
    if interpolation == "log":
        positive = (y0[:, 0] > 0) & (y1[:, 0] > 0)
        power = np.log(y1[positive] / y0[positive]) / np.log(
            x1[positive] / x0[positive]
        )
        density[positive] = y0[positive] * (energy[positive] / x0[positive]) ** power
    elif interpolation != "linear":
        raise ValueError("Expected linear or log interpolation")
    return float(
        np.sum(density * function(energy) * weights * (right - left)[:, None] / 2)
    )


def response_moments(data, interpolation="linear"):
    lines = np.asarray(data["lines"])

    def moment(function, lower=0, upper=np.inf):
        value = continuum_integral(data, function, lower, upper, interpolation)
        selected = (lines[:, 0] >= lower) & (lines[:, 0] < upper)
        value += np.sum(lines[selected, 1] * function(lines[selected, 0]))
        # Only integrated band strength is supplied. Use its midpoint for
        # moments other than S(0), and retain endpoint bounds for comparison.
        for start, end, strength in data["bands"]:
            if start >= lower and end <= upper:
                value += strength * function(np.array((start + end) / 2))
            elif start < upper and end > lower:
                raise ValueError("Moment boundary cuts through an unresolved band")
        return float(value)

    strength = moment(np.ones_like)
    logarithmic = moment(np.log)
    inverse_second = moment(lambda energy: (2 * RYDBERG / energy) ** 2)
    bands_alpha_bounds = [0.0, 0.0]
    for start, end, strength_band in data["bands"]:
        midpoint = strength_band * (2 * RYDBERG / ((start + end) / 2)) ** 2
        bands_alpha_bounds[0] += strength_band * (2 * RYDBERG / end) ** 2 - midpoint
        bands_alpha_bounds[1] += strength_band * (2 * RYDBERG / start) ** 2 - midpoint
    return {
        "strength_through_100keV": strength,
        "strength_over_electron_count": strength / data["electrons"],
        "log_mean_excitation_eV": float(np.exp(logarithmic / strength)),
        "polarizability_a0_cubed": inverse_second,
        "unresolved_band_polarizability_bounds": [
            inverse_second + shift for shift in bands_alpha_bounds
        ],
        "sub_ionization_strength": moment(
            np.ones_like, upper=data["first_threshold_eV"]
        ),
        "continuum_25_100eV_strength": continuum_integral(
            data, np.ones_like, 25, 100, interpolation
        ),
        "above_100eV_strength": moment(np.ones_like, lower=100),
        "notes": "Finite tabulated range only; no above-100-keV tail added. "
        "Sum rules influenced source selection. This is not independent validation.",
    }


def pjg_loss_density(target, loss, parameters, center_s=None):
    """Implied df/dW under the existing PJG soft-channel weights.

    If A_j(T)=f_j*A(T), photoionization gives df_j/dW=W*A_j(W-I_j).
    This reconstructs a *model* energy-loss response without shifting an
    experimental total spectrum by a single assumed binding energy.
    It omits non-ionizing transitions and any explicit inner-shell response.
    """
    loss = np.asarray(loss, dtype=float)
    result = np.zeros_like(loss)
    target_data = TARGETS[target]
    for fraction, threshold in zip(
        target_data.fractions, target_data.thresholds, strict=True
    ):
        opened = loss >= threshold
        result[opened] += (
            fraction
            * loss[opened]
            * optical_coefficient(
                target, loss[opened] - threshold, parameters, center_s=center_s
            )
        )
    return result
