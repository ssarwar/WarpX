# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Constrain the existing PJG shape with an N2 optical reference.

This is a non-production, nonrelativistic diagnostic. First fit K, Gamma_s,
and Lambda to the photoelectron Bethe coefficient inferred from Rudd (1992),
Table II. Freeze those values, then fit J, p, and a common multiplier of the
printed C_j to the recommended proton totals and spectra. By default remove
gamma_1, hence also its unused denominator gamma_2. Compare the retained
printed width correction with --retain-width-correction.

The photon parameterization is not a new independent measurement, nor a
complete inclusive electron-production response including all decay channels.
No optical high-energy continuation or new target-response function is fitted.
"""

import argparse
from dataclasses import replace

import numpy as np
from fit_pjg_repair import moments, report, spectral_grid
from pjg_repair import (
    initial_parameters,
    nonrelativistic_sdcs,
    optical_coefficient,
)
from reference import RYDBERG, rudd_total
from scipy.optimize import least_squares


def n2_photoelectron_coefficient(secondary):
    """Return sum_j [(df_j/dW)/W] at W=T+I_j, in eV^-2.

    Rudd et al., Rev. Mod. Phys. 64, 441 (1992), Eqs. (32)--(33), Table II,
    https://doi.org/10.1103/RevModPhys.64.441. The tabulated dimensionless
    quantity is W*df/dW, so divide by W twice, not once. Each row already
    includes its orbital population. There is no extra occupancy multiplier.

    Restrict this diagnostic to 0 <= T <= 100 eV. The table's footnote warns
    that its high-energy asymptote requires further fitting. This restriction
    is a chosen audit window, not a claim of experimental coverage/accuracy.
    """
    t = np.asarray(secondary, dtype=float)
    if np.any(~np.isfinite(t)) or np.any((t < 0) | (t > 100)):
        raise ValueError("The optical diagnostic requires 0 <= T <= 100 eV")
    thresholds = (15.59, 28.8, 410)
    gaussians = (
        (
            (5.5, 0.465, 0.155, 0),
            (0.8, 0.59, 0.035, 0),
            (3.4, 0.71, 0.14, 0),
            (1.3, 0.81, 0.09, 0),
            (1.6, 0.26, 0.08, 0),
            (10, 0.11, 0.1, 1.1),
        ),
        (
            (1.6, 0.26, 0.06, 0),
            (0.7, 0.36, 0.055, 0),
            (0.5, 0.525, 0.02, 2),
            (0.5, 0.47, 0.08, 1),
            (0.6, 0.12, 0.08, 0.1),
            (0.4, 0.18, 0.05, 1),
        ),
    )
    result = np.zeros_like(t)
    for threshold, coefficients in zip(thresholds[:2], gaussians, strict=True):
        loss = t + threshold
        r = RYDBERG / loss
        strength = sum(
            a * np.exp(-(((r - b) / c) ** 2)) * r**d for a, b, c, d in coefficients
        )
        result += strength / loss**2
    loss = t + thresholds[2]
    r = RYDBERG / loss
    result += (2.99 * r + 5813 * r**2 + 30506 * r**3 - 1987000 * r**4) / loss**2
    return result


def fit_optical_shape():
    """Fit only three existing shape parameters on a fixed 2--100 eV grid."""
    t = np.geomspace(2, 100, 40)
    expected = n2_photoelectron_coefficient(t)
    defaults = initial_parameters("N2")

    def unpack(values):
        return replace(
            defaults, amplitude=values[0], gamma_s=values[1], broad_excess=values[2]
        )

    fit = least_squares(
        lambda values: np.log(
            optical_coefficient("N2", t, unpack(values), relativistic=False) / expected
        ),
        [1, defaults.gamma_s, defaults.broad_excess],
        bounds=([0.01, 0.1, 1], [20, 100, 1000]),
        max_nfev=600,
        ftol=1e-9,
        xtol=1e-9,
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return unpack(fit.x)


def fit_proton_dependence(parameters, retain_width_correction=False):
    """Freeze the optical shape and refit only existing proton-energy terms."""
    energies = np.geomspace(5e3, 4e6, 31)
    expected = rudd_total("N2", energies)
    se, st, sy = spectral_grid("N2")

    def unpack(values):
        return replace(
            parameters,
            j=values[0],
            power=values[1],
            bethe_scale=values[2],
            gamma_numerator_scale=float(retain_width_correction),
        )

    def residual(values):
        p = unpack(values)
        total, _ = moments("N2", energies, p)
        return np.concatenate(
            (
                np.log(total / expected) / np.sqrt(len(energies)),
                np.log(nonrelativistic_sdcs("N2", se, st, p) / sy) / np.sqrt(len(sy)),
            )
        )

    fit = least_squares(
        residual,
        [30, 0.8, 1],
        bounds=([0.01, 0.1, 1e-4], [1e5, 4, 100]),
        max_nfev=600,
        ftol=1e-9,
        xtol=1e-9,
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return unpack(fit.x)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retain-width-correction", action="store_true")
    args = parser.parse_args()
    shape = fit_optical_shape()
    t = np.geomspace(2, 100, 40)
    ratio = optical_coefficient("N2", t, shape, relativistic=False)
    ratio /= n2_photoelectron_coefficient(t)
    print(
        "Optical ratio percentiles 0/10/50/90/100:",
        np.percentile(ratio, [0, 10, 50, 90, 100]),
    )
    print(
        "Unfitted T=0 optical ratio:",
        optical_coefficient("N2", 0, shape, relativistic=False)
        / n2_photoelectron_coefficient(0),
    )
    report("N2", fit_proton_dependence(shape, args.retain_width_correction))


if __name__ == "__main__":
    main()
