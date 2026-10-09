# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Original PJG Table III constants and stable common-center line shapes."""

from dataclasses import dataclass

import numpy as np


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
# empirical B(E) parameters are specified in original_pjg.py.
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
