# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent quadrature and asymptotic checks for the research references.

Run with python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization.
No fit, simulation output, private table, or digitized plot is needed.
"""

import unittest

import numpy as np
from reference import (
    ELECTRON_REST_ENERGY,
    ORBITALS,
    PROTON_REST_ENERGY,
    free_maximum_transfer,
    pwba_momentum_grid,
    rudd_sdcs,
    rudd_total,
)
from scipy.integrate import quad


def integrated_rudd(target, energy, lower=0):
    # Cross-section scaling makes the absolute quadrature tolerance meaningful.
    result, _ = quad(
        lambda x: 1.0e16 * rudd_sdcs(target, energy, np.expm1(x)) * np.exp(x),
        np.log1p(lower),
        np.log1p(energy),
        epsabs=1.0e-10,
        epsrel=1.0e-10,
    )
    return result * 1.0e-16


class ReferenceChecks(unittest.TestCase):
    def test_orbital_electron_counts(self):
        self.assertEqual(sum(ORBITALS["N2"][1]), 14)
        self.assertEqual(sum(ORBITALS["O2"][1]), 16)

    def test_rudd_bound_tail(self):
        # These are regression values for the recommended model, not measured
        # tail fractions. They establish why its SDCS cannot be capped at Tmax.
        for target, fraction in (("N2", 0.14073635), ("O2", 0.22586266)):
            energy = 5000.0
            maximum = free_maximum_transfer(energy)
            total = integrated_rudd(target, energy)
            tail = integrated_rudd(target, energy, maximum)
            self.assertAlmostEqual(tail / total, fraction, delta=1.0e-7)
            self.assertGreater(rudd_sdcs(target, energy, 2 * maximum), 0)

    def test_reference_totals_are_not_identical_fits(self):
        # The two recommendations have independent parameterizations. A joint
        # fit must allow their stated uncertainties, not force both to agree.
        for target, ratio in (("N2", 0.9704850425), ("O2", 1.1197288922)):
            integral_ratio = integrated_rudd(target, 5000) / rudd_total(target, 5000)
            self.assertAlmostEqual(integral_ratio, ratio, delta=1.0e-8)

    def test_longitudinal_against_analytic_integral(self):
        electron = ELECTRON_REST_ENERGY
        mass = PROTON_REST_ENERGY
        for energy in (1.0e3, 5.0e3, 1.0e6, 8.0e8, 1.0e12):
            for loss in (0.1, 20.0, 500.0):
                recoil, longitudinal, transverse = pwba_momentum_grid(energy, loss)
                p = np.sqrt(energy * (energy + 2 * mass))
                pf = np.sqrt((energy - loss) * (energy - loss + 2 * mass))
                z_min = (loss * (2 * (energy + mass) - loss) / (p + pf)) ** 2
                z_max = (p + pf) ** 2
                q_min = z_min / (np.sqrt(z_min + electron**2) + electron)
                q_max = z_max / (np.sqrt(z_max + electron**2) + electron)
                analytic = np.log1p(2 * electron / q_min) - np.log1p(
                    2 * electron / q_max
                )
                self.assertAlmostEqual(longitudinal.sum() / analytic, 1, delta=2.0e-12)
                self.assertTrue(np.all(recoil > 0))
                self.assertTrue(np.all(longitudinal > 0))
                self.assertTrue(np.all(transverse > 0))

    def test_transverse_against_independent_q_quadrature(self):
        electron = ELECTRON_REST_ENERGY
        mass = PROTON_REST_ENERGY
        for energy in (5.0e3, 1.0e6, 8.0e8):
            for loss in (20.0, 200.0):
                recoil, _, transverse = pwba_momentum_grid(energy, loss)
                result = np.sum(transverse / (1 + recoil / 50) ** 2)
                p = np.sqrt(energy * (energy + 2 * mass))
                pf = np.sqrt((energy - loss) * (energy - loss + 2 * mass))
                q_momentum_min = loss * (2 * (energy + mass) - loss) / (p + pf)
                z_min, z_max = q_momentum_min**2, (p + pf) ** 2
                q_min = z_min / (np.sqrt(z_min + electron**2) + electron)
                q_max = z_max / (np.sqrt(z_max + electron**2) + electron)
                beta_squared = (p / (energy + mass)) ** 2

                def integrand(log_q):
                    q = np.exp(log_q)
                    z = q * (q + 2 * electron)
                    cosine = (z + loss * (2 * (energy + mass) - loss)) / (
                        2 * p * np.sqrt(z)
                    )
                    return (
                        beta_squared
                        * (1 - cosine**2)
                        * loss**2
                        * 2
                        * electron
                        * q
                        / (z - loss**2) ** 2
                        / (1 + q / 50) ** 2
                    )

                expected, _ = quad(
                    integrand,
                    np.log(q_min),
                    np.log(q_max),
                    epsabs=1.0e-16,
                    epsrel=1.0e-10,
                )
                self.assertAlmostEqual(result / expected, 1, delta=1.0e-8)

    def test_relativistic_dipole_logarithm(self):
        # At fixed small energy loss the transverse term supplies ln(gamma^2).
        # This tests its coefficient, not a target's optical oscillator strength.
        gamma_values = (1.0e3, 1.0e4, 1.0e5)
        integrals = []
        for gamma in gamma_values:
            energy = (gamma - 1) * PROTON_REST_ENERGY
            _, _, transverse = pwba_momentum_grid(energy, 20.0, points=384)
            integrals.append(transverse.sum())
        slopes = np.diff(integrals) / np.diff(np.log(np.asarray(gamma_values) ** 2))
        np.testing.assert_allclose(slopes, 1, rtol=2.0e-5, atol=0)

    def test_transverse_bethe_subtraction(self):
        # In the dipole limit the complete transverse integral tends to
        # ln(gamma^2)-beta^2, not just the logarithm. A small test energy loss
        # isolates this kernel limit; it is not a molecular ionization energy.
        for gamma in (1.1, 2.0, 10.0, 100.0):
            energy = (gamma - 1) * PROTON_REST_ENERGY
            _, _, transverse = pwba_momentum_grid(energy, 0.01, points=384)
            expected = 2 * np.log(gamma) - (1 - 1 / gamma**2)
            self.assertAlmostEqual(transverse.sum(), expected, delta=1.0e-10)

    def test_bound_response_can_extend_above_free_endpoint(self):
        energy, loss = 5000.0, 40.0
        self.assertGreater(loss - 15.59, free_maximum_transfer(energy))
        recoil, longitudinal, transverse = pwba_momentum_grid(energy, loss)
        # Example smooth positive finite-q response, not an N2 fit.
        response = (1 + recoil / 50) ** -2
        self.assertGreater(np.sum((longitudinal + transverse) * response), 0)

    def test_invalid_kinematics_rejected(self):
        for loss in (0, -1, 5000, 6000):
            with self.assertRaises(ValueError):
                pwba_momentum_grid(5000, loss)


if __name__ == "__main__":
    unittest.main()
