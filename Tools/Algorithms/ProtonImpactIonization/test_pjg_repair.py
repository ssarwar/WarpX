# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent algebra, quadrature, and limit tests; not fit acceptance tests."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext

import numpy as np
from fit_pjg_optical import fit_optical_shape, n2_photoelectron_coefficient
from pjg_repair import (
    DIAGNOSTIC_FITS,
    TARGETS,
    RepairParameters,
    distortion,
    lorentzian_pair,
    nonrelativistic_sdcs,
    optical_coefficient,
    relativistic_core_parts,
)
from reference import (
    BETHE_CONSTANT,
    ELECTRON_REST_ENERGY,
    PROTON_REST_ENERGY,
    free_maximum_transfer,
)
from scipy.integrate import quad, simpson


class PJGRepairChecks(unittest.TestCase):
    def test_continuum_fractions_and_bethe_scales(self):
        for p in TARGETS.values():
            self.assertAlmostEqual(sum(p.fractions), 1)
            self.assertEqual(len(p.thresholds), len(p.fractions))
            self.assertEqual(len(p.thresholds), len(p.bethe_constants))
            self.assertTrue(all(value > 0 for value in p.fractions))
            # The printed C_j/I_j ratios are almost constant, but preserve
            # their original rounding rather than silently replacing them.
            ratios = np.array(p.bethe_constants) / p.thresholds
            self.assertLess(np.max(ratios) / np.min(ratios), 1.004)

    def test_positive_difference_against_direct_subtraction(self):
        for center in (-13, 0, 6):
            for width, extra in ((1, 1), (18, 80), (100, 3)):
                t = np.linspace(0, 1000, 1001)
                broad, difference = lorentzian_pair(t, width, center, extra)
                broad_direct = 1 / ((t - center) ** 2 + width**2 + extra**2)
                np.testing.assert_allclose(broad, broad_direct, rtol=1.0e-14)
                # Direct double-precision subtraction loses relative accuracy
                # in precisely the tail this repair addresses. Evaluate the
                # independent, unfactored reference at higher precision.
                with localcontext() as context:
                    context.prec = 60
                    expected = []
                    for value in t:
                        d = Decimal(int(value) - center) ** 2 + Decimal(width) ** 2
                        expected.append(float(1 / d - 1 / (d + Decimal(extra) ** 2)))
                np.testing.assert_allclose(difference, expected, rtol=1.0e-10)
                self.assertTrue(np.all(difference > 0))

    def test_separate_hard_and_soft_power_laws(self):
        t, width, center, extra = 1.0e12, 18, -6, 100
        broad, difference = lorentzian_pair(t, width, center, extra)
        self.assertAlmostEqual(t**2 * broad, 1, delta=2.0e-11)
        self.assertAlmostEqual(t**4 * difference / extra**2, 1, delta=3.0e-11)

    def test_finite_oscillator_strength_integral(self):
        # Integrate (T+I)*(L_n-L_b). The analytic primitive is independent
        # of the implementation's factored rational expression.
        for width, extra, center, threshold in ((18, 100, -6, 16), (4, 80, 4, 40)):
            broad_width = np.hypot(width, extra)
            integral0 = (np.pi / 2 + np.arctan(center / width)) / width - (
                np.pi / 2 + np.arctan(center / broad_width)
            ) / broad_width
            expected = 0.5 * np.log1p(extra**2 / (center**2 + width**2))
            expected += (center + threshold) * integral0
            numeric, _ = quad(
                lambda t: (t + threshold) * lorentzian_pair(t, width, center, extra)[1],
                0,
                np.inf,
                epsabs=1.0e-11,
                epsrel=1.0e-11,
            )
            self.assertAlmostEqual(numeric / expected, 1, delta=1.0e-11)

    def test_joint_fast_projectile_hard_secondary_limit(self):
        energy = 1.0e15
        mass = PROTON_REST_ENERGY
        beta2 = energy * (energy + 2 * mass) / (energy + mass) ** 2
        equivalent = ELECTRON_REST_ENERGY * beta2 / 2
        maximum = free_maximum_transfer(energy)
        for target, parameters in DIAGNOSTIC_FITS.items():
            for fraction in (0.01, 0.5, 1):
                t = fraction * maximum
                soft, hard = relativistic_core_parts(target, energy, t, parameters)
                factor = 1 - beta2 * t / maximum + t**2 / (2 * (energy + mass) ** 2)
                expected = (
                    TARGETS[target].electrons
                    * BETHE_CONSTANT
                    / equivalent
                    / t**2
                    * factor
                )
                self.assertAlmostEqual((soft + hard) / expected, 1, delta=1.0e-6)

    def test_distortion_does_not_freeze_below_unity(self):
        for p in DIAGNOSTIC_FITS.values():
            values = distortion(np.array([5.0e3, 1.0e8, 1.0e15]), p)
            self.assertTrue(np.all(np.diff(values) > 0))
            self.assertGreater(values[-1], 1 - 1.0e-6)

    def test_relativistic_bethe_coefficient(self):
        gamma = np.array([1.0e5, 1.0e7])
        energy = (gamma - 1) * PROTON_REST_ENERGY
        equivalent = ELECTRON_REST_ENERGY / 2 * (1 - 1 / gamma**2)
        for target, parameters in DIAGNOSTIC_FITS.items():
            for t in (0, 2, 20, 100):
                soft, _ = relativistic_core_parts(target, energy, t, parameters)
                slope = np.diff(equivalent * soft / BETHE_CONSTANT) / np.diff(
                    np.log(gamma**2)
                )
                expected = optical_coefficient(target, t, parameters)
                np.testing.assert_allclose(slope, expected, rtol=1.0e-5)

    def test_nonrelativistic_quadrature_and_positivity(self):
        for target, parameters in DIAGNOSTIC_FITS.items():
            for energy in (5.0e3, 1.0e5, 4.0e6):
                x = np.linspace(0, np.log1p(energy), 4097)
                t = np.expm1(x)
                y = nonrelativistic_sdcs(target, energy, t, parameters) * 1.0e16
                self.assertTrue(np.all(np.isfinite(y)))
                self.assertTrue(np.all(y >= 0))
                numeric = simpson(y * np.exp(x), x=x)
                adaptive, _ = quad(
                    lambda u: (
                        1.0e16
                        * np.exp(u)
                        * nonrelativistic_sdcs(target, energy, np.expm1(u), parameters)
                    ),
                    0,
                    np.log1p(energy),
                    epsabs=1.0e-9,
                    epsrel=1.0e-9,
                )
                self.assertAlmostEqual(numeric / adaptive, 1, delta=1.0e-8)
                self.assertEqual(
                    nonrelativistic_sdcs(target, energy, energy, parameters), 0
                )

    def test_no_accidental_relativistic_tail_extrapolation(self):
        p = DIAGNOSTIC_FITS["N2"]
        for energy in (100, 1.0e8, np.nan):
            with self.assertRaises(ValueError):
                nonrelativistic_sdcs("N2", energy, 10, p)
        with self.assertRaises(ValueError):
            relativistic_core_parts("N2", 5000, 2 * free_maximum_transfer(5000), p)
        with self.assertRaises(ValueError):
            RepairParameters(0, 1, 1, 1, 1)

    def test_optical_diagnostic_window(self):
        t = np.linspace(0, 100, 101)
        values = n2_photoelectron_coefficient(t)
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertTrue(np.all(values > 0))
        for t in (-1, 101, np.nan):
            with self.assertRaises(ValueError):
                n2_photoelectron_coefficient(t)

    def test_optical_shape_fit_is_reproducible(self):
        # Check the stated finite-window approximation, not physical acceptance
        # of the optical input or a continuation to unmeasured energies.
        p = fit_optical_shape()
        t = np.geomspace(2, 100, 40)
        ratio = optical_coefficient("N2", t, p, relativistic=False)
        ratio /= n2_photoelectron_coefficient(t)
        self.assertLess(np.max(np.abs(ratio - 1)), 0.15)
        q = replace(p, gamma_numerator_scale=0)
        ratio = optical_coefficient("N2", t, q) / n2_photoelectron_coefficient(t)
        self.assertLess(np.max(np.abs(ratio - 1)), 0.15)

    def test_bethe_scale_is_not_an_optical_amplitude(self):
        p = DIAGNOSTIC_FITS["N2"]
        q = replace(p, bethe_scale=3)
        t = np.array([0, 2, 20, 100])
        np.testing.assert_array_equal(
            optical_coefficient("N2", t, p), optical_coefficient("N2", t, q)
        )
        soft_p, hard_p = relativistic_core_parts("N2", 1e6, t, p)
        soft_q, hard_q = relativistic_core_parts("N2", 1e6, t, q)
        np.testing.assert_array_equal(hard_p, hard_q)
        self.assertTrue(np.all(soft_q > soft_p))
        for scale in (-1, np.nan):
            with self.assertRaises(ValueError):
                replace(p, gamma_numerator_scale=scale)
        q = replace(p, gamma_numerator_scale=0)
        self.assertTrue(np.all(nonrelativistic_sdcs("N2", 5000, t, q) > 0))


if __name__ == "__main__":
    unittest.main()
