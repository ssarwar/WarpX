# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent ionization references; not a production molecular SDCS.

Energies are in eV, lengths in cm. Rudd's recommended functions are empirical
model curves, not independent experimental measurements. The PWBA quadrature
implements the longitudinal and transverse kernel of Salvat and Heredia,
Nucl. Instrum. Methods B 546, 165157 (2024), Eqs. (68)--(71). Its common
longitudinal/transverse GOS approximation is appropriate to the dipole limit;
it does not replace the full free-electron Bhabha kernel for hard collisions.
"""

from functools import lru_cache

import numpy as np
from scipy.special import expit

ELECTRON_REST_ENERGY = 510998.95069
PROTON_REST_ENERGY = 938272088.16
RYDBERG = 13.605693122994
BOHR_RADIUS = 5.29177210544e-9
BETHE_CONSTANT = 4 * np.pi * BOHR_RADIUS**2 * RYDBERG**2

# Rudd et al., Rev. Mod. Phys. 64, 441 (1992), Table I(a): orbital ionization
# energies and occupancies. These are not the PJG continuum-fraction weights.
ORBITALS = {
    "N2": ((15.59, 16.96, 18.78, 37.3, 409.9, 409.9), (2, 4, 2, 2, 2, 2)),
    "O2": (
        (12.07, 16.42, 18.88, 25.69, 40.3, 543.5, 543.5),
        (2, 2, 4, 2, 2, 2, 2),
    ),
}

# Rudd (1992), Table V: A1, B1, C1, D1, E1, A2, B2, C2, D2, alpha.
RUDD_PARAMETERS = {
    "N2": (1.05, 12, 0.74, -0.39, 0.80, 0.95, 1.20, 1, 1.30, 0.70),
    "O2": (1.02, 50, 0.40, 0.12, 0.30, 1, 5, 0.55, 0, 0.59),
    "inner": (1.25, 0.5, 1, 1, 3, 1.10, 1.30, 1, 0, 0.66),
}


def free_maximum_transfer(kinetic_energy, rest_energy=PROTON_REST_ENERGY):
    """Independent invariant-energy expression, not a molecular cutoff."""
    energy = np.asarray(kinetic_energy, dtype=float)
    mass = rest_energy
    electron = ELECTRON_REST_ENERGY
    return (
        2
        * electron
        * energy
        * (energy + 2 * mass)
        / ((mass + electron) ** 2 + 2 * electron * energy)
    )


def rudd_total(target, kinetic_energy):
    """Recommended proton electron-production total in cm^2, Rudd (1985).

    The fit range is 5--4000 keV. Extrapolation is not a measurement and does
    not include a relativistic correction.
    """
    a, b, c, d = {"N2": (3.82, 2.78, 1.8, 0.70), "O2": (4.77, 0, 1.76, 0.93)}[target]
    u = np.asarray(kinetic_energy) * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY / RYDBERG
    if np.any(u <= 0):
        raise ValueError("The projectile kinetic energy must be positive")
    return 4 * np.pi * BOHR_RADIUS**2 / (u / (a * np.log1p(u) + b) + 1 / (c * u**d))


def rudd_sdcs(target, kinetic_energy, secondary_energy):
    """Recommended nonrelativistic SDCS in cm^2/eV, Rudd (1992), Eqs. 41--48.

    N2 data span 5--1700 keV; O2 differential data span 7.5--300 keV. Values
    outside those ranges are extrapolations of this recommended model. Its
    empirical exponential tail is not truncated at the free-electron Tmax.
    """
    energy, secondary = np.broadcast_arrays(kinetic_energy, secondary_energy)
    if np.any(energy <= 0) or np.any(secondary < 0):
        raise ValueError("Require positive projectile and nonnegative secondary energy")
    thresholds, occupation = ORBITALS[target]
    result = np.zeros_like(secondary, dtype=float)
    for threshold, number in zip(thresholds, occupation, strict=True):
        parameters = RUDD_PARAMETERS[
            "inner" if threshold > 2 * thresholds[0] else target
        ]
        a1, b1, c1, d1, e1, a2, b2, c2, d2, alpha = parameters
        v = np.sqrt(energy * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY / threshold)
        w = secondary / threshold
        h1 = a1 * np.log1p(v**2) / (v**2 + b1 / v**2)
        l1 = c1 * v**d1 / (1 + e1 * v ** (d1 + 4))
        h2 = a2 / v**2 + b2 / v**4
        l2 = c2 * v**d2
        cutoff = 4 * v**2 - 2 * v - RYDBERG / (4 * threshold)
        result += (
            BETHE_CONSTANT
            * number
            / threshold**3
            * (h1 + l1 + h2 * l2 / (h2 + l2) * w)
            / (1 + w) ** 3
            * expit(-alpha * (w - cutoff) / v)
        )
    return result


@lru_cache(maxsize=8)
def _quadrature(points):
    if points < 2:
        raise ValueError("At least two quadrature points are required")
    return np.polynomial.legendre.leggauss(points)


def pwba_momentum_grid(
    kinetic_energy, energy_loss, rest_energy=PROTON_REST_ENERGY, points=192
):
    """Return Q and separate, positive longitudinal/transverse integration weights.

    With G(Q,W)=df(Q,W)/dW, sum(weights*G)/W times Z^2*C_B/E_e
    approximates d(sigma)/dW. The weights integrate W times the kernel in
    Eq. (68). Energy loss W includes binding; it is not electron kinetic energy.
    The function contains no target-response or empirical fitting parameters.

    Sampling the logarithm of y=(cq)^2-W^2 resolves both the distant Coulomb
    logarithm and the transverse near-light-cone peak. Rationalized kinematics
    avoid subtracting nearly equal quantities for slow or relativistic beams.
    """
    energy, loss = np.broadcast_arrays(
        np.asarray(kinetic_energy, dtype=float), np.asarray(energy_loss, dtype=float)
    )
    if rest_energy <= 0 or np.any(loss <= 0) or np.any(loss >= energy):
        raise ValueError(
            "Require positive rest energy and 0 < energy loss < kinetic energy"
        )
    mass = rest_energy
    electron = ELECTRON_REST_ENERGY
    momentum_squared = energy * (energy + 2 * mass)
    momentum = np.sqrt(momentum_squared)
    final_momentum = np.sqrt((energy - loss) * (energy - loss + 2 * mass))
    total = energy + mass
    q_minus = loss * (2 * total - loss) / (momentum + final_momentum)
    y_min = (
        (mass / (total + momentum))
        * (mass / (total - loss + final_momentum))
        * (q_minus + loss) ** 2
    )
    y_max = (momentum + final_momentum - loss) * (momentum + final_momentum + loss)
    nodes, weights = _quadrature(points)
    span = np.log(y_max / y_min)
    y = y_min[..., None] * np.exp(span[..., None] * (nodes + 1) / 2)
    q_squared = y + loss[..., None] ** 2
    recoil = q_squared / (np.sqrt(q_squared + electron**2) + electron)
    # Factored 1-cos(theta_r)^2 is nonnegative throughout the allowed interval.
    sin_squared = ((y - y_min[..., None]) / q_squared) * (
        (y_max[..., None] - y) / (4 * momentum_squared[..., None])
    )
    beta_squared = (energy / total) * ((energy + 2 * mass) / total)
    jacobian = (
        electron / np.sqrt(q_squared + electron**2) * y * span[..., None] * weights / 2
    )
    longitudinal = jacobian / q_squared
    transverse = (
        jacobian * beta_squared[..., None] * sin_squared * (loss[..., None] / y) ** 2
    )
    return recoil, longitudinal, transverse
