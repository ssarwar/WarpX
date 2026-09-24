# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent, dimensional RBEQ reference from thesis Eqs. 11.119--11.122.

The parameters are transcribed from Table 11.12 and the archived fitted blocks.
This reference does not read the implementation's parameter header or CDF.
"""

from functools import lru_cache

import numpy as np
from scipy.integrate import quad

REST = 510998.95069
RYDBERG = 13.605693122994
NORM = 4 * np.pi * (5.29177210544e-11) ** 2 * (7.2973525643e-3) ** 4
THESIS = {
    "N2": [
        [409.5, 603.3, 4, 1],
        [37.3, 71.13, 2, 0.760],
        [18.72, 63.18, 2, 1],
        [16.74, 44.30, 4, 0.938],
        [15.58, 54.91, 2, 0.792],
    ],
    "O2": [
        [543.8, 796.2, 4, 1],
        [40.33, 79.73, 2, 0.9600],
        [27.05, 90.92, 2, 1],
        [20.30, 71.84, 2, 1],
        [17.08, 59.89, 4, 1],
        [12.07, 84.88, 2, 0.9314],
    ],
}
SNAPSHOT = {
    "N2": [
        [409.5, 603.3, 4, 1],
        [37.8, 71.13, 2, 0.760],
        [23.6, 63.18, 1, 1],
        [18.746, 63.18, 1, 1],
        [16.716, 44.30, 4, 0.938],
        [15.58, 54.91, 2, 0.792],
    ],
    "O2": [[531, 796.2, 4, 1], *THESIS["O2"][1:]],
}


@lru_cache(maxsize=None)
def dipole(binding, q, snapshot):
    if not snapshot:
        return -(1 + q - (5 - 3 * q) * np.log(2)) / q
    # Direct quadrature of the finite-energy logarithmic moment, independent
    # of the algebraic antiderivative used by the production exporter.
    upper = (1e6 - binding) / 2

    def weight(w):
        return 1 / (w + binding) ** 2 + (w + binding) / (1e6 - w) ** 3

    denominator = quad(weight, 0, upper, epsabs=1e-18, epsrel=1e-11)[0]
    moment = (
        quad(
            lambda w: weight(w) * np.log((w + binding) / RYDBERG),
            0,
            upper,
            epsabs=1e-18,
            epsrel=1e-11,
        )[0]
        / denominator
    )
    ratio = 1e6 / binding
    normalization = 1 + 0.5 / ratio * (1 - 1 / ratio) - 2 / (ratio + 1)
    # The archived lnBm antiderivative has ln(t)/(4*(t+1)), whereas direct
    # integration gives ln(t)/(2*(t+1)). Reproduce the snapshot literally.
    snapshot_correction = 0.5 * np.log(ratio) / ((ratio + 1) * normalization)
    return snapshot_correction + (
        2 * (np.log(binding / RYDBERG) - moment) + 1 - 1 / q + (5 / q - 3) * np.log(2)
    )


def differential(energy, secondary, shell, snapshot=False):
    b, u, n, q = shell
    if energy <= b:
        return 0.0
    gamma = 1 + (energy + u + b) / REST
    beta2 = (gamma - 1) * (gamma + 1) / gamma**2
    t, w = energy / b, secondary / b
    incident = energy / REST
    log_term = (
        np.log(incident * (incident + 2))
        - incident * (incident + 2) / (1 + incident) ** 2
        - np.log(2 * b / REST)
        + dipole(b, q, snapshot)
    )
    exchange = (2 * gamma - 1) / ((1 + t) * gamma**2)
    return (
        NORM
        * n
        / (2 * (b / REST) * beta2 * b)
        * (
            q * log_term * ((w + 1) ** -3 + (t - w) ** -3)
            + (2 - q)
            * (
                (w + 1) ** -2
                + (t - w) ** -2
                - exchange * ((w + 1) ** -1 + (t - w) ** -1)
                + (b / (REST * gamma)) ** 2
            )
        )
    )


def partial(energy, shell, snapshot=False):
    if energy <= shell[0]:
        return 0.0
    # Integrate in log(W+B) to resolve the strongly asymmetric high-energy tail.
    binding = shell[0]
    upper = np.log((energy + binding) / (2 * binding))
    return (
        quad(
            lambda x: (
                differential(energy, binding * np.expm1(x), shell, snapshot)
                * binding
                * np.exp(x)
                / 1e-20
            ),
            0,
            upper,
            epsabs=1e-12,
            epsrel=1e-10,
        )[0]
        * 1e-20
    )
