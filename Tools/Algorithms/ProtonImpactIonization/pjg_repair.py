# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""PJG-preserving research candidates, not production cross sections.

Keep the PJG widths, peak-energy dependence, continuum fractions, logarithm,
and low-velocity distortion. Replace its unconstrained subtraction by a
positive difference of common-center Lorentzians and a normalized hard term.
The two functions below deliberately separate a relativistic free-domain
limit test from a nonrelativistic data-fit experiment with a molecular tail.
Neither is an accepted full relativistic molecular SDCS. See PJG_REPAIR.md.
"""

from dataclasses import dataclass

import numpy as np
from reference import (
    BETHE_CONSTANT,
    ELECTRON_REST_ENERGY,
    PROTON_REST_ENERGY,
    RYDBERG,
    free_maximum_transfer,
)
from scipy.special import expit


@dataclass(frozen=True)
class TargetParameters:
    electrons: int
    k: float
    gamma_s: float
    gamma_numerator: float
    gamma_denominator: float
    center_s: float
    center_numerator: float
    center_denominator: float
    j: float
    power: float
    broad_width: float
    thresholds: tuple
    fractions: tuple
    bethe_constants: tuple


# PJG (1976), Table III, proton parameters. The 1977 erratum makes the
# apparent superscript '+' in Eq. (16) an ordinary addition, not clipping.
# It does not change these proton parameters. The printed cutoff shift and
# empirical B(E) parameters are intentionally absent from the repair.
TARGETS = {
    "N2": TargetParameters(
        14,
        7.58e-16,
        11.1,
        1.27e4,
        1.81e3,
        4,
        2.03e4,
        1.97e3,
        3.39,
        0.807,
        115,
        (15.58, 16.73, 18.75, 22, 23.6, 40),
        (0.456, 0.2, 0.104, 0.07, 0.07, 0.1),
        (2.48, 2.66, 2.99, 3.50, 3.76, 6.37),
    ),
    "O2": TargetParameters(
        16,
        6.55e-16,
        13.1,
        5e5,
        7.60e4,
        6.34,
        2.52e3,
        1.28e2,
        40.3,
        1.314,
        189.1,
        (12.1, 16.1, 16.9, 18.2, 20.3, 23, 37),
        (0.08, 0.19, 0.19, 0.17, 0.11, 0.16, 0.1),
        (1.93, 2.56, 2.69, 2.90, 3.23, 3.66, 5.89),
    ),
}


@dataclass(frozen=True)
class RepairParameters:
    j: float
    power: float
    amplitude: float
    gamma_s: float
    broad_excess: float

    def __post_init__(self):
        values = (self.j, self.power, self.amplitude, self.gamma_s, self.broad_excess)
        if not all(np.isfinite(value) and value > 0 for value in values):
            raise ValueError("All repair parameters must be finite and positive")


def initial_parameters(target):
    p = TARGETS[target]
    return RepairParameters(p.j, p.power, 1, p.gamma_s, p.broad_width)


# Five-parameter fits to recommended model curves, not new experimental
# parameter determinations. No production coefficient is taken from here.
DIAGNOSTIC_FITS = {
    "N2": RepairParameters(
        30.80263212989071,
        0.7378262198393548,
        1.5821956151244818,
        4.488907513359014,
        80.48665474362566,
    ),
    "O2": RepairParameters(
        6.89268267689334,
        1.351401623673788,
        0.8847004941692884,
        8.471382305977615,
        319.9835070937648,
    ),
}


def lorentzian_pair(secondary, width, center, broad_excess):
    """Return the broad Lorentzian and positive narrow-minus-broad difference.

    Gamma_b^2 = Gamma^2 + Lambda^2. The factored difference retains its
    T^-4 tail when directly subtracting the two T^-2 functions would cancel.
    """
    t, width, center, extra = np.broadcast_arrays(
        secondary, width, center, broad_excess
    )
    if any(np.any(~np.isfinite(value)) for value in (t, width, center, extra)):
        raise ValueError("The Lorentzian arguments must be finite")
    if np.any(width <= 0) or np.any(extra <= 0):
        raise ValueError("Both widths must be positive")
    denominator = (t - center) ** 2 + width**2
    broad = 1 / (denominator + extra**2)
    difference = (extra**2 / denominator) * broad
    return broad, difference


def _shape(target, equivalent_energy, secondary, parameters):
    p = TARGETS[target]
    width = parameters.gamma_s + p.gamma_numerator / (
        equivalent_energy + p.gamma_denominator
    )
    center = p.center_s - p.center_numerator / (
        equivalent_energy + p.center_denominator
    )
    broad, difference = lorentzian_pair(
        secondary, width, center, parameters.broad_excess
    )
    return width, broad, difference


def distortion(kinetic_energy, parameters):
    """PJG-shaped suppression with an unbounded energy scale, so D tends to 1."""
    equivalent = np.asarray(kinetic_energy) * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY
    if np.any(~np.isfinite(equivalent)) or np.any(equivalent <= 0):
        raise ValueError("The incident kinetic energy must be positive")
    return expit(parameters.power * np.log(equivalent / parameters.j))


def relativistic_core_parts(target, kinetic_energy, secondary, parameters):
    """Return soft and hard parts inside the free domain for asymptotic tests.

    This is not a bound-tail prescription or a molecular threshold model.
    D multiplies both parts, as in PJG. Thus the complete molecular candidate
    approaches Bhabha in the joint fast-projectile/hard-secondary limit; it
    is not exactly a free-electron cross section at finite D < 1.
    """
    energy, t = np.broadcast_arrays(kinetic_energy, secondary)
    maximum = free_maximum_transfer(energy)
    if (
        np.any(~np.isfinite(energy))
        or np.any(~np.isfinite(t))
        or np.any(energy <= 0)
        or np.any(t < 0)
        or np.any(t > maximum)
    ):
        raise ValueError("Require positive E and 0 <= T <= free Tmax")
    p = TARGETS[target]
    mass = PROTON_REST_ENERGY
    gamma = 1 + energy / mass
    beta_squared = (energy / (energy + mass)) * ((energy + 2 * mass) / (energy + mass))
    equivalent = ELECTRON_REST_ENERGY * beta_squared / 2
    width, broad, difference = _shape(target, equivalent, t, parameters)
    logarithm = sum(
        f * (np.log(4 * equivalent * c * gamma**2 / i + np.e) - beta_squared)
        for f, c, i in zip(p.fractions, p.bethe_constants, p.thresholds, strict=True)
    )
    x = t / maximum
    factor = (1 - x) + x / gamma**2 + 0.5 * (t / (energy + mass)) ** 2
    scale = distortion(energy, parameters) / equivalent
    soft = scale * parameters.amplitude * p.k * width**2 * logarithm * difference
    hard = scale * p.electrons * BETHE_CONSTANT * broad * factor
    return soft, hard


def nonrelativistic_sdcs(target, kinetic_energy, secondary, parameters):
    """PJG-shaped fit diagnostic in cm^2/eV, restricted to 5--4000 keV.

    Multiply the two-Lorentzian repair by Rudd's Eqs. (41)--(42) cutoff,
    separately for each retained PJG continuum. Alpha is fixed to Rudd's
    published molecular value (N2: 0.70, O2: 0.59), not fitted. The cutoff
    is empirical; its center is not a redefinition of exact free Tmax.
    This nonrelativistic experiment omits the Bhabha relativistic terms.
    It must not be promoted by extending that polynomial above free Tmax.
    """
    energy, t = np.broadcast_arrays(kinetic_energy, secondary)
    if np.any(~np.isfinite(energy)) or np.any((energy < 5000) | (energy > 4e6)):
        raise ValueError("The nonrelativistic fit diagnostic requires 5--4000 keV")
    if np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError("The secondary energy must be finite and nonnegative")
    p = TARGETS[target]
    equivalent = energy * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY
    width, broad, difference = _shape(target, equivalent, t, parameters)
    result = np.zeros_like(t, dtype=float)
    alpha = {"N2": 0.70, "O2": 0.59}[target]
    for fraction, threshold, constant in zip(
        p.fractions, p.thresholds, p.bethe_constants, strict=True
    ):
        v = np.sqrt(equivalent / threshold)
        center = 4 * equivalent - 2 * v * threshold - RYDBERG / 4
        cutoff = expit(alpha * (center - t) / (v * threshold))
        logarithm = np.log(4 * equivalent * constant / threshold + np.e)
        soft = parameters.amplitude * p.k * width**2 * logarithm * difference
        hard = p.electrons * BETHE_CONSTANT * broad
        # Necessary energy bound only; this is not an exact three-body closure.
        result += np.where(
            t + threshold <= energy, fraction * cutoff * (soft + hard), 0
        )
    return distortion(energy, parameters) / equivalent * result


def optical_coefficient(target, secondary, parameters, relativistic=True):
    """Implied Bethe-log coefficient a(T), not a measured optical spectrum.

    Relativistic PJG width functions saturate at E_e = m_e c^2/2. The formal
    nonrelativistic high-velocity limit instead has width Gamma_s and center
    T_s. Distinguish these two tests rather than interchanging their limits.
    """
    equivalent = ELECTRON_REST_ENERGY / 2 if relativistic else np.inf
    width, _, difference = _shape(target, equivalent, secondary, parameters)
    return (
        parameters.amplitude
        * TARGETS[target].k
        * width**2
        / BETHE_CONSTANT
        * difference
    )
