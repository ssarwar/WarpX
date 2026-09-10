# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent integration and regression checks of the frozen energy model.

An agreement tolerance in these tests is a numerical or regression bound,
not an experimental uncertainty or proof of a complete molecular response.
"""

import unittest

import numpy as np
from calibrated_pjg import PARAMETERS, sdcs, total_cross_section
from fit_pjg import centered_shape_residual
from original_pjg import bhabha, legacy_cutoff, legacy_sdcs, moments
from pjg_model import kinematics, molecular_endpoint
from pjg_moments import LossGrid, load_pstar, mass_stopping
from reference import BETHE_CONSTANT, PROTON_REST_ENERGY, free_maximum_transfer
from scipy.integrate import quad
from source_datasets import pjg_1976_figure5, rudd_1979_figure7
from target_parameters import TARGETS


class CalibrationChecks(unittest.TestCase):
    def test_original_total_against_analytic_integral(self):
        # Integrate the printed, signed Eq. (16) independently. Its hard
        # remainder is logarithmic, not the positive final Bhabha kernel.
        for target, p in TARGETS.items():
            for energy in (5e3, 5e4, 8e8):
                gamma, beta2, ee, _ = kinematics(energy)
                upper = legacy_cutoff(energy)
                width = p.gamma_s + p.gamma_numerator / (ee + p.gamma_denominator)
                center = p.center_s - p.center_numerator / (ee + p.center_denominator)
                b0, broad_center, delta = {
                    "N2": (0.029, 53.3, 84),
                    "O2": (0.030, 68.3, 132.1),
                }[target]
                reduction = b0 * (np.log(ee / 8239) ** 2 + 1.035)

                def primitive(c, w):
                    return (np.arctan((upper - c) / w) + np.arctan(c / w)) / w

                thresholds = np.asarray(p.thresholds)
                logarithm = np.dot(
                    p.fractions,
                    np.log(
                        4 * ee * gamma**2 * np.asarray(p.bethe_constants) / thresholds
                        + np.e
                    )
                    - beta2,
                )
                soft = (
                    p.k
                    * width**2
                    * logarithm
                    * (
                        primitive(center, width)
                        - reduction * primitive(broad_center, p.broad_width)
                    )
                )
                hard = (
                    p.electrons
                    * BETHE_CONSTANT
                    * np.dot(
                        p.fractions,
                        upper / (4 * (energy + 2 * PROTON_REST_ENERGY) ** 2)
                        - np.log1p(upper / thresholds) / (upper + thresholds + delta),
                    )
                )
                analytic = (soft + hard) / ee / (1 + (p.j / ee) ** p.power)
                numeric = moments(target, [energy])[0][0]
                self.assertAlmostEqual(numeric / analytic, 1, delta=2e-9)
                # The support remains the printed edge, not free Tmax.
                self.assertEqual(legacy_sdcs(target, energy, upper * 1.01), 0)

    def test_parameter_reductions(self):
        for target, p in PARAMETERS.items():
            self.assertEqual(p.power, TARGETS[target].power)
            self.assertEqual(p.low_width_ratio, 1)
        self.assertEqual(PARAMETERS["O2"].edge_scale, 1)
        with self.assertRaises(TypeError):
            PARAMETERS["N2"] = PARAMETERS["O2"]

    def test_incident_range_and_secondary_validation(self):
        for target in PARAMETERS:
            for energy in (0, 4999, 1e10 + 1, np.nan, np.inf):
                with self.assertRaises(ValueError):
                    sdcs(target, energy, 0)
                with self.assertRaises(ValueError):
                    total_cross_section(target, energy)
            for secondary in (-1, np.nan, np.inf):
                with self.assertRaises(ValueError):
                    sdcs(target, 5e3, secondary)

    def test_total_against_independent_adaptive_quadrature(self):
        for target in PARAMETERS:
            for energy in (5e3, 50e3, 1e6):
                endpoint = molecular_endpoint(
                    target, energy, min(TARGETS[target].thresholds)
                )
                # Integrate directly in log(1+T), with no use of LossGrid's
                # segments or Jacobians. Scaling makes quad's absolute
                # tolerance meaningful for molecular cross sections.
                integral, _ = quad(
                    lambda x: (
                        float(sdcs(target, energy, np.expm1(x))) * np.exp(x) * 1e16
                    ),
                    0,
                    np.log1p(endpoint),
                    epsabs=1e-9,
                    epsrel=1e-10,
                    limit=200,
                )
                self.assertAlmostEqual(
                    integral / (total_cross_section(target, energy) * 1e16),
                    1,
                    delta=2e-8,
                )

    def test_shape_preserving_total_array_interface(self):
        energy = np.array([[5e3, 1e6], [1e6, 5e3]])
        for target in PARAMETERS:
            total = total_cross_section(target, energy)
            self.assertEqual(total.shape, energy.shape)
            self.assertEqual(total[0, 0], total[1, 1])
            self.assertEqual(total[0, 1], total[1, 0])

    def test_soft_limit_bound_tail_and_molecular_support(self):
        for target in PARAMETERS:
            for energy in (5e3, 50e3, 1e10):
                zero = sdcs(target, energy, 0)
                self.assertGreater(zero, 0)
                self.assertTrue(np.isfinite(zero))
                self.assertGreater(
                    sdcs(target, energy, free_maximum_transfer(energy) * 1.001), 0
                )
                endpoint = molecular_endpoint(
                    target, energy, min(TARGETS[target].thresholds)
                )
                np.testing.assert_array_equal(
                    sdcs(target, energy, endpoint * np.array([1, 1.01, 2])), 0
                )

    def test_hard_secondary_limit(self):
        for target in PARAMETERS:
            t = np.array([1e4, 1e5])
            ratio = sdcs(target, 800e6, t) / bhabha(target, 800e6, t)
            self.assertLess(abs(ratio[0] - 1), 0.002)
            self.assertLess(abs(ratio[1] - 1), 0.0002)

    def test_moment_inequalities_and_independent_stopping_budget(self):
        for target, p in PARAMETERS.items():
            table = load_pstar(target)
            keep = table.energies >= 5e3
            energy = table.energies[keep]
            m = LossGrid(target, energy, points=1025)(p)
            endpoint = molecular_endpoint(
                target, energy, min(TARGETS[target].thresholds)
            )
            self.assertTrue(np.all(m.kinetic**2 <= m.total * m.kinetic_second))
            self.assertTrue(np.all(m.kinetic_second <= endpoint * m.kinetic))
            self.assertTrue(np.all((m.tail_total > 0) & (m.tail_total < m.total)))
            self.assertTrue(
                np.all(mass_stopping(target, m.ionization) < table.electronic[keep])
            )

    def test_profiled_shape_is_normalization_invariant(self):
        ratio, weights = np.array([0.5, 1, 1.7]), np.array([0.5, 1, 1])
        np.testing.assert_allclose(
            centered_shape_residual(ratio, weights),
            centered_shape_residual(7.1 * ratio, weights),
            rtol=1e-14,
            atol=1e-15,
        )

    def test_measured_nitrogen_digitization_provenance_and_units(self):
        pjg = pjg_1976_figure5()
        rudd = rudd_1979_figure7()
        self.assertEqual(pjg.shape, (49, 4))
        self.assertEqual(rudd.shape, (18, 3))
        np.testing.assert_array_equal(np.unique(pjg[:, 0]), [5e4, 1e5, 3e5, 1e6])
        np.testing.assert_array_equal(np.unique(rudd[:, 0]), [5e3, 2e4, 7e4])
        for data in (pjg, rudd):
            self.assertTrue(np.all(data[:, :3] > 0))
            self.assertTrue(np.all((data[:, 2] > 1e-23) & (data[:, 2] < 1e-16)))
        # A decade in PJG's ordinate is 184 pixels. In Rudd's m^2/eV
        # ordinate it is 838/4 pixels; each reader returns cm^2/eV.
        np.testing.assert_allclose(
            pjg_1976_figure5((0, 184))[:, 2], pjg[:, 2] / 10, rtol=2e-14
        )
        np.testing.assert_allclose(
            rudd_1979_figure7((0, 838 / 4))[:, 2], rudd[:, 2] / 10, rtol=2e-14
        )


if __name__ == "__main__":
    unittest.main()
