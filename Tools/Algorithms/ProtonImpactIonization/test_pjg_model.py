# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent physical and numerical checks for the matched PJG reference.

These test algebra and numerical accuracy, not experimental acceptance or
accelerator execution. No external optical table or fit is needed.
"""

import unittest

import numpy as np
from fit_support import MomentGrid
from original_pjg import bhabha
from pjg_model import (
    MOLECULAR_REST_ENERGIES,
    PARAMETERS,
    SpectrumGrid,
    endpoint_momentum_broadening,
    hard_factor,
    molecular_endpoint,
    optical_coefficient,
    sdcs,
)
from reference import (
    BETHE_CONSTANT,
    ELECTRON_REST_ENERGY,
    PROTON_REST_ENERGY,
    free_maximum_transfer,
)
from scipy.integrate import quad
from target_parameters import TARGETS


class MatchedPJGChecks(unittest.TestCase):
    def test_complete_bhabha_bracket(self):
        for energy in (1e3, 1e6, 8e8, 1e12, 1e15):
            maximum = free_maximum_transfer(energy)
            t = maximum * np.linspace(0, 1, 201)
            mass = PROTON_REST_ENERGY
            beta2 = energy * (energy + 2 * mass) / (energy + mass) ** 2
            expected = 1 - beta2 * t / maximum + t**2 / (2 * (energy + mass) ** 2)
            np.testing.assert_allclose(hard_factor(energy, t), expected, rtol=2e-12)

    def test_tail_value_derivative_and_positivity(self):
        for energy in (1e6, 1e8, 8e8, 1e10):
            maximum = free_maximum_transfer(energy)
            delta = maximum * 1e-5
            left2, left, middle, right, right2 = hard_factor(
                energy, maximum + delta * np.array([-2, -1, 0, 1, 2])
            )
            # One-sided second-order derivatives avoid confusing the finite
            # curvature of a C1 join with a derivative discontinuity.
            left_slope = 3 * middle - 4 * left + left2
            right_slope = -3 * middle + 4 * right - right2
            self.assertAlmostEqual(left_slope / right_slope, 1, delta=2e-4)
            t = maximum * np.geomspace(1, 100, 201)
            values = hard_factor(energy, t)
            self.assertTrue(np.all(np.isfinite(values)))
            self.assertTrue(np.all(values >= 0))
            self.assertTrue(np.all(np.diff(values) <= 0))

    def test_independent_moving_electron_broadening(self):
        # Differentiate a moving-target collinear two-body solution, not the
        # closed expression implemented by endpoint_momentum_broadening.
        m, mass = np.longdouble(ELECTRON_REST_ENERGY), np.longdouble(PROTON_REST_ENERGY)
        ratio = m / mass
        delta = m * np.longdouble("1e-4")
        for value in (5e3, 1e6, 8e8, 1e10, 1e15):
            energy = np.longdouble(value)
            total = energy + mass
            momentum = np.sqrt(energy * (energy + 2 * mass))

            def moving_endpoint(k):
                electron_energy = np.hypot(m, k)
                invariant = (
                    mass**2 + m**2 + 2 * (total * electron_energy - momentum * k)
                )
                return k**2 / (electron_energy + m) + (
                    2
                    * (momentum + k)
                    * (momentum * electron_energy - total * k)
                    / invariant
                )

            derivative = (moving_endpoint(delta) - moving_endpoint(-delta)) / (
                2 * delta
            )
            nr_derivative = -2 * momentum / total * (1 - ratio) / (1 + ratio) ** 2
            self.assertAlmostEqual(
                float(derivative / nr_derivative) / endpoint_momentum_broadening(value),
                1,
                delta=2e-8,
            )

    def test_bhabha_recovery_in_a_relevant_hard_window(self):
        # Binding effects need not disappear near the broadened binary edge.
        # This checks a genuine hard window well inside free support.
        for target, parameters in PARAMETERS.items():
            t = np.geomspace(1e4, 2e6, 101)
            ratio = sdcs(target, 8e8, t, parameters) / bhabha(target, 8e8, t)
            self.assertLess(np.max(np.abs(ratio - 1)), 0.002)
            ratio = sdcs(target, 8e8, 1e5, parameters) / bhabha(target, 8e8, 1e5)
            self.assertLess(abs(ratio - 1), 0.0002)

    def test_relativistic_logarithmic_coefficient(self):
        gamma = np.array([1e5, 1e7])
        energy = (gamma - 1) * PROTON_REST_ENERGY
        equivalent = ELECTRON_REST_ENERGY / 2 * (1 - 1 / gamma**2)
        for target, parameters in PARAMETERS.items():
            for t in (0, 2, 20, 100):
                soft, _ = SpectrumGrid(target, energy, t).parts(parameters)
                slope = np.diff(equivalent * soft / BETHE_CONSTANT) / np.diff(
                    np.log(gamma**2)
                )
                np.testing.assert_allclose(
                    slope, optical_coefficient(target, t, parameters), rtol=1e-6
                )

    def test_positivity_threshold_support_and_zero_secondary_limit(self):
        for target, parameters in PARAMETERS.items():
            binding = min(TARGETS[target].thresholds)
            threshold = binding * (
                1 + PROTON_REST_ENERGY / MOLECULAR_REST_ENERGIES[target]
            )
            self.assertEqual(sdcs(target, 0, 0, parameters), 0)
            self.assertEqual(sdcs(target, 0.99 * threshold, 0, parameters), 0)
            for energy in np.geomspace(20, 1e15, 65):
                maximum = molecular_endpoint(target, energy, binding)
                t = np.expm1(np.linspace(0, np.log1p(maximum), 513))
                values = sdcs(target, energy, t, parameters)
                self.assertTrue(np.all(np.isfinite(values)))
                self.assertTrue(np.all(values >= 0))
                self.assertEqual(sdcs(target, energy, 1.01 * maximum, parameters), 0)
                np.testing.assert_allclose(
                    sdcs(target, energy, 1e-9, parameters),
                    values[0],
                    rtol=1e-8,
                    atol=1e-35,
                )

    def test_quadrature_against_adaptive_integration(self):
        for target, parameters in PARAMETERS.items():
            for energy in (5e3, 1e5, 8e8):
                maximum = molecular_endpoint(
                    target, energy, min(TARGETS[target].thresholds)
                )
                upper = np.log1p(maximum)
                split = np.log1p(free_maximum_transfer(energy))
                expected = []
                for power in (0, 1):

                    def integrand(x):
                        t = np.expm1(x)
                        return (
                            float(sdcs(target, energy, t, parameters))
                            * 1e16
                            * np.exp(x)
                            * t**power
                        )

                    integral, _ = quad(
                        integrand,
                        0,
                        upper,
                        points=[split],
                        epsabs=1e-9,
                        epsrel=1e-9,
                        limit=200,
                    )
                    expected.append(integral)
                # Resolve the narrow relativistic binary edge as well as
                # the low-energy peak before testing the first moment.
                total, mean = MomentGrid(target, [energy], 8193)(parameters)
                self.assertAlmostEqual(total[0] * 1e16 / expected[0], 1, delta=2e-7)
                self.assertAlmostEqual(
                    mean[0] / (expected[1] / expected[0]), 1, delta=2e-7
                )

    def test_molecular_endpoint_against_direct_cm_solution(self):
        m = np.longdouble(ELECTRON_REST_ENERGY)
        for target, neutral_float in MOLECULAR_REST_ENERGIES.items():
            neutral = np.longdouble(neutral_float)
            for mass_float in (PROTON_REST_ENERGY, 4 * PROTON_REST_ENERGY):
                mass = np.longdouble(mass_float)
                for energy_float in (5e3, 1e5, 8e8):
                    energy = np.longdouble(energy_float)
                    binding = np.longdouble(20)
                    invariant = (mass + neutral) ** 2 + 2 * neutral * energy
                    root = np.sqrt(invariant)
                    residual = mass + neutral - m + binding
                    star_energy = (invariant + m**2 - residual**2) / (2 * root)
                    star_momentum = np.sqrt((star_energy - m) * (star_energy + m))
                    momentum = np.sqrt(energy * (energy + 2 * mass))
                    expected = (
                        (energy + mass + neutral) * star_energy
                        + momentum * star_momentum
                    ) / root - m
                    actual = molecular_endpoint(target, energy_float, 20, mass_float)
                    self.assertAlmostEqual(actual / float(expected), 1, delta=2e-9)

    def test_effective_binding_and_invalid_inputs(self):
        for target, parameters in PARAMETERS.items():
            t = np.geomspace(1e-6, 1e5, 101)
            values = SpectrumGrid(target, 1e5, t).effective_binding(parameters)
            self.assertTrue(np.all(values >= min(TARGETS[target].thresholds)))
            self.assertTrue(np.all(values <= max(TARGETS[target].thresholds)))
            for energy, secondary in ((-1, 0), (np.nan, 0), (1e5, -1), (1e5, np.inf)):
                with self.assertRaises(ValueError):
                    sdcs(target, energy, secondary, parameters)


if __name__ == "__main__":
    unittest.main()
