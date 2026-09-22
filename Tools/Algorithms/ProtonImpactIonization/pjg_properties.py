# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Distribution and moment diagnostics for the final energy model."""

import numpy as np
from pjg_model import kinematics
from pjg_moments import LossGrid, free_second_moment, mass_stopping
from scipy.integrate import cumulative_trapezoid, simpson
from target_parameters import TARGETS

ENERGIES = np.array(
    [
        5e3,
        7e3,
        7.5e3,
        10e3,
        15e3,
        20e3,
        30e3,
        50e3,
        70e3,
        100e3,
        150e3,
        200e3,
        300e3,
        500e3,
        700e3,
        1e6,
        1.5e6,
        2e6,
        3e6,
        4e6,
        10e6,
        100e6,
        800e6,
        1e9,
        1e10,
    ]
)


def power_law_moments(x, y, orders=(0, 1, 2)):
    """Integrate a positive piecewise log-log interpolant, without extrapolation.

    For y=C*T**a each cell contributes C*(b**r-a**r)/r,
    r=a+order+1. The expm1 form retains the logarithmic r -> 0 limit.
    Sparse measured points do not resolve resonances or unmeasured tails.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if (
        x.ndim != 1
        or y.shape != x.shape
        or len(x) < 2
        or np.any(~np.isfinite(x))
        or np.any(~np.isfinite(y))
        or np.any(x <= 0)
        or np.any(y <= 0)
        or np.any(np.diff(x) <= 0)
    ):
        raise ValueError("Require ordered positive finite energies and SDCS values")
    log_step = np.log(x[1:] / x[:-1])
    slope = np.log(y[1:] / y[:-1]) / log_step
    result = []
    for order in orders:
        exponent = slope + order + 1
        factor = np.divide(
            np.expm1(exponent * log_step),
            exponent,
            out=log_step.copy(),
            where=np.abs(exponent) > 1e-14,
        )
        result.append(np.sum(y[:-1] * x[:-1] ** (order + 1) * factor))
    return np.asarray(result)


def mapped_cdf(secondary, weighted_density, dx):
    """Return a positive, monotone CDF across contiguous mapped segments.

    Trapezoidal integration preserves positivity. Moments elsewhere use
    Simpson integration, so the CDF normalization difference is recorded.
    This research calculation is not the GPU sampling implementation.
    """
    secondary = np.asarray(secondary)
    weighted_density = np.asarray(weighted_density)
    if (
        secondary.ndim != 2
        or secondary.shape[0] < 1
        or secondary.shape[1] < 2
        or secondary.shape != weighted_density.shape
        or np.any(~np.isfinite(weighted_density))
        or np.any(weighted_density < 0)
        or dx <= 0
        or not np.isfinite(dx)
        or np.any(~np.isfinite(secondary))
        or np.any(np.diff(secondary, axis=1) < 0)
        or not np.allclose(secondary[1:, 0], secondary[:-1, -1], rtol=1e-14, atol=0)
    ):
        raise ValueError("Require contiguous ordered segments and nonnegative density")
    cumulative = cumulative_trapezoid(weighted_density, dx=dx, axis=-1, initial=0)
    offset = np.r_[0, np.cumsum(cumulative[:-1, -1])]
    cumulative += offset[:, None]
    normalization = cumulative[-1, -1]
    if normalization <= 0 or not np.isfinite(normalization):
        raise ValueError("The distribution must have a positive finite integral")
    x = np.r_[secondary[0], secondary[1:, 1:].ravel()]
    cdf = np.r_[cumulative[0], cumulative[1:, 1:].ravel()] / normalization
    return x, cdf, float(normalization)


def distribution_properties(target, parameters, energies=ENERGIES, points=4097):
    grid = LossGrid(target, energies, points)
    moment = grid(parameters)
    soft, hard = grid.grid.parts(parameters)
    density = soft + hard
    hard_total = simpson(hard * grid.jacobian, dx=grid.dx, axis=-1).sum(axis=1)
    hard_loss = simpson(hard * grid.jacobian * grid.secondary, dx=grid.dx, axis=-1).sum(
        axis=1
    )
    mean = moment.kinetic / moment.total
    second = moment.kinetic_second / moment.total
    variance = second - mean**2
    _, beta2, _, maximum = kinematics(grid.energies)
    records = []
    for index, energy in enumerate(grid.energies):
        x, cdf, area = mapped_cdf(
            grid.secondary[index], density[index] * grid.jacobian[index], grid.dx
        )
        # Keep the left side of a CDF plateau. The requested quantiles are
        # interior, so an underflowed far-tail plateau cannot move an endpoint.
        keep = np.r_[True, np.diff(cdf) > 0]
        quantiles = np.interp([0.5, 0.9, 0.99], cdf[keep], x[keep])
        thresholds = [2, min(TARGETS[target].thresholds), 100, 1000, 10000]
        probability_below = np.interp(thresholds, x, cdf)
        sigma = moment.total[index]
        records.append(
            {
                "energy_eV": float(energy),
                "sigma_e_cm2": float(sigma),
                "mean_T_eV": float(mean[index]),
                "mean_T2_eV2": float(second[index]),
                "standard_deviation_T_eV": float(np.sqrt(variance[index])),
                "quantile_50_90_99_eV": quantiles.tolist(),
                "thresholds_eV": thresholds,
                "probability_below_threshold": probability_below.tolist(),
                "free_Tmax_eV": float(maximum[index]),
                "molecular_Tmax_eV": float(grid.secondary[index, -1, -1]),
                "above_free_count_energy_second_fraction": [
                    float(moment.tail_total[index] / sigma),
                    float(moment.tail_kinetic[index] / moment.kinetic[index]),
                    float(
                        moment.tail_kinetic_second[index] / moment.kinetic_second[index]
                    ),
                ],
                "effective_binding_eV": float(moment.binding[index] / sigma),
                "effective_pair_cost_eV": float(moment.ionization[index] / sigma),
                "kinetic_mass_stopping": float(
                    mass_stopping(target, moment.kinetic[index])
                ),
                "binding_mass_stopping": float(
                    mass_stopping(target, moment.binding[index])
                ),
                "pair_mass_stopping": float(
                    mass_stopping(target, moment.ionization[index])
                ),
                "hard_count_fraction": float(hard_total[index] / sigma),
                "hard_kinetic_fraction": float(
                    hard_loss[index] / moment.kinetic[index]
                ),
                "second_moment_to_free": float(
                    moment.kinetic_second[index] / free_second_moment(target, energy)
                ),
                # Convert sigma from cm^2 to m^2; c is exact in SI units.
                "rate_coefficient_m3_per_s": float(
                    sigma * 1e-4 * 299792458 * np.sqrt(beta2[index])
                ),
                "iid_energy_relative_variance_times_N": float(
                    variance[index] / mean[index] ** 2
                ),
                "compound_poisson_relative_variance_times_expected_N": float(
                    second[index] / mean[index] ** 2
                ),
                "cdf_to_simpson_relative_difference": float(area / sigma - 1),
                "nonnegative": bool(np.all(density[index] >= 0)),
                "cdf_monotone": bool(np.all(np.diff(cdf) >= 0)),
                "bounded_second_moment": bool(second[index] <= x[-1] * mean[index]),
            }
        )
    return records


def property_convergence(target, parameters):
    """Check quantiles independently of the high-order moment integration."""
    energy = [5e3, 50e3, 1e6, 800e6, 1e10]
    coarse = distribution_properties(target, parameters, energy, points=2049)
    fine = distribution_properties(target, parameters, energy, points=8193)
    return {
        "energies_eV": energy,
        "relative_change_2049_to_8193": {
            key: float(
                np.max(
                    np.abs(
                        np.array([r[key] for r in coarse])
                        / np.array([r[key] for r in fine])
                        - 1
                    )
                )
            )
            for key in (
                "sigma_e_cm2",
                "mean_T_eV",
                "mean_T2_eV2",
                "quantile_50_90_99_eV",
                "pair_mass_stopping",
            )
        },
    }
