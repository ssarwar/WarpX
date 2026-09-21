#!/usr/bin/env python3
"""Independent continuum field quadrature for a translating 3D Gaussian.

Integrate the Coulomb Green function in the beam rest frame, then Lorentz
transform the fields. The integration variable is a Gaussian variance (m^2).
This unbounded, untruncated solution requires separate domain and mesh studies
when compared with WarpX's finite-domain Poisson solve.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import c, epsilon_0


def gaussian_fields(r, z, sigma_r, sigma_z, charge, velocity, points=192):
    gamma = 1 / np.sqrt(1 - (velocity / c) ** 2)
    r, z = np.broadcast_arrays(r, z)
    rest_z = gamma * z
    nodes, weights = leggauss(points)
    u = (nodes + 1) / 2
    # A squared rational map removes the endpoint square-root behavior at
    # infinity, keeping the quadrature smooth at both limits.
    tau = sigma_r**2 * (u / (1 - u)) ** 2
    jacobian = 2 * sigma_r**2 * u / (1 - u) ** 3
    er = np.zeros_like(r, dtype=float)
    ez = np.zeros_like(r, dtype=float)
    for variance, weight in zip(tau, weights * jacobian / 2):
        radial = sigma_r**2 + variance
        longitudinal = (gamma * sigma_z) ** 2 + variance
        kernel = np.exp(-(r**2) / (2 * radial) - rest_z**2 / (2 * longitudinal))
        kernel *= weight / (radial * np.sqrt(longitudinal))
        er += kernel * r / radial
        ez += kernel * rest_z / longitudinal
    factor = charge / (4 * np.pi * epsilon_0 * np.sqrt(2 * np.pi))
    er *= gamma * factor
    ez *= factor
    return er, ez, velocity * er / c**2


def check():
    """Check the quadrature against the closed spherical Gaussian solution."""
    from scipy.special import erf

    sigma, charge = 0.003, 1e-10
    r = sigma * np.geomspace(0.01, 20, 200)
    exact = (
        charge
        / (4 * np.pi * epsilon_0 * r**2)
        * (
            erf(r / (np.sqrt(2) * sigma))
            - np.sqrt(2 / np.pi) * r / sigma * np.exp(-(r**2) / (2 * sigma**2))
        )
    )
    er, ez, bt = gaussian_fields(r, 0, sigma, sigma, charge, 0)
    np.testing.assert_allclose(er, exact, rtol=3e-7)
    np.testing.assert_array_equal(ez, 0)
    np.testing.assert_array_equal(bt, 0)
    # Integrating the Green function at two orders bounds quadrature error.
    for speed in [0.8 * c, -0.8 * c]:
        low = gaussian_fields(r, r / 2, sigma, 4 * sigma, charge, speed, 192)
        high = gaussian_fields(r, r / 2, sigma, 4 * sigma, charge, speed, 384)
        for a, b in zip(low, high):
            np.testing.assert_allclose(a, b, rtol=3e-7)
    print("PASS: continuum Gaussian field quadrature and Lorentz transform")


if __name__ == "__main__":
    check()
