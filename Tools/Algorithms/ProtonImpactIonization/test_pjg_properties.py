# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent checks of property calculations, not validation of fitted physics."""

import unittest

import numpy as np
from pjg_model import PARAMETERS, SpectrumGrid
from pjg_properties import (
    distribution_properties,
    mapped_cdf,
    power_law_moments,
)
from scipy.integrate import quad


class PropertyChecks(unittest.TestCase):
    def test_exact_power_law_window_integrals(self):
        x = np.array([2, 3, 20, 50, 100], dtype=float)
        a, b, c = x[0], x[-1], 7.3
        expected = c * np.array([1 / a - 1 / b, np.log(b / a), b - a])
        np.testing.assert_allclose(
            power_law_moments(x, c / x**2), expected, rtol=2e-15, atol=0
        )
        expected = c * np.array([np.log(b / a), b - a, (b**2 - a**2) / 2])
        np.testing.assert_allclose(
            power_law_moments(x, c / x), expected, rtol=2e-15, atol=0
        )

    def test_window_mean_is_independent_of_absolute_normalization(self):
        x, y = np.array([3, 7, 50, 100]), np.array([100, 40, 8, 2])
        original = power_law_moments(x, y)
        rescaled = power_law_moments(x, y * 1e-20)
        np.testing.assert_allclose(
            original / original[0], rescaled / rescaled[0], rtol=2e-15, atol=0
        )

    def test_mapped_uniform_cdf_across_segment_boundary(self):
        u = np.linspace(0, 1, 33)
        secondary = np.stack([u, 1 + 2 * u])
        jacobian = np.stack([np.ones_like(u), 2 * np.ones_like(u)])
        x, cdf, total = mapped_cdf(secondary, jacobian, 1 / 32)
        self.assertEqual(total, 3)
        np.testing.assert_allclose(cdf, x / 3, rtol=1e-15, atol=0)
        self.assertEqual(np.interp(0.5, cdf, x), 1.5)

    def test_invalid_measurements_or_density_are_rejected(self):
        for x, y in (([1, 1], [2, 3]), ([0, 2], [2, 3]), ([1, 2], [2, -3])):
            with self.assertRaises(ValueError):
                power_law_moments(x, y)
        with self.assertRaises(ValueError):
            mapped_cdf([[0, 1], [2, 3]], [[1, 1], [1, 1]], 1)
        with self.assertRaises(ValueError):
            mapped_cdf([[0, 1]], [[1, -1]], 1)

    def test_quantiles_against_independent_adaptive_quadrature(self):
        for target, p in PARAMETERS.items():
            result = distribution_properties(target, p, energies=[50e3], points=4097)[0]
            for q, t in zip(
                (0.5, 0.9, 0.99), result["quantile_50_90_99_eV"], strict=True
            ):
                integral, _ = quad(
                    lambda x: (
                        float(SpectrumGrid(target, 50e3, np.expm1(x))(p))
                        * np.exp(x)
                        * 1e16
                    ),
                    0,
                    np.log1p(t),
                    epsabs=1e-10,
                    epsrel=1e-10,
                )
                self.assertAlmostEqual(
                    integral / (result["sigma_e_cm2"] * 1e16), q, delta=2e-6
                )

    def test_moment_noise_identities_and_positivity(self):
        for target, p in PARAMETERS.items():
            results = distribution_properties(
                target, p, energies=[5e3, 1e6, 1e10], points=2049
            )
            for row in results:
                self.assertTrue(row["nonnegative"])
                self.assertTrue(row["cdf_monotone"])
                self.assertTrue(row["bounded_second_moment"])
                self.assertAlmostEqual(
                    row["compound_poisson_relative_variance_times_expected_N"]
                    - row["iid_energy_relative_variance_times_N"],
                    1,
                    delta=1e-10,
                )
                self.assertGreater(row["effective_pair_cost_eV"], row["mean_T_eV"])
                self.assertTrue(np.all(np.diff(row["quantile_50_90_99_eV"]) > 0))


if __name__ == "__main__":
    unittest.main()
