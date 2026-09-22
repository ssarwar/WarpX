# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent unit, quadrature, support, and source-semantics checks.

These tests use small analytic fixtures; they require no network, downloaded
PDFs, fitted coefficients, or production GPU kernels.
"""

import unittest
from io import StringIO

import numpy as np
from extract_nifs_oscillators import rows
from optical_reference import continuum_integral, response_moments
from reference import BOHR_RADIUS, RYDBERG
from source_datasets import (
    HC_EV_NM,
    PHOTO_FACTOR,
    RUDD_1983_ENERGIES,
    RUDD_1983_FITTED_TOTALS,
    OpticalTable,
    cheng_1989_figure1,
    load_leiden,
    rudd_1983_uncertainty,
)


class SourceTests(unittest.TestCase):
    def test_photoabsorption_units(self):
        factor = 4 * np.pi**2 * (1 / 137.035999177) * BOHR_RADIUS**2 * RYDBERG * 1e18
        self.assertAlmostEqual(PHOTO_FACTOR / factor, 1, delta=2e-6)

    def test_cross_section_axis_change_has_no_jacobian(self):
        table = load_leiden(StringIO("10 3e-18 1e-18 2e-18\n20 6e-18 2e-18 4e-18\n"))
        np.testing.assert_allclose(
            table["ionization"].energy, HC_EV_NM / np.array([20, 10])
        )
        np.testing.assert_allclose(table["ionization"].cross_section_Mb, [4, 2])

    def test_subtraction_roundoff_is_not_physical_negativity(self):
        table = load_leiden(StringIO("10 3e-18 -1e-33 3e-18\n20 6e-18 0 6e-18\n"))
        np.testing.assert_allclose(table["ionization"].cross_section_Mb, [6, 3])
        with self.assertRaises(ValueError):
            load_leiden(StringIO("10 3e-18 -1e-19 3e-18\n20 6e-18 0 6e-18\n"))

    def test_absorption_and_channel_checks(self):
        data = "10 3e-18 0 2e-18\n20 6e-18 0 4e-18\n"
        with self.assertRaises(ValueError):
            load_leiden(StringIO(data))
        table = load_leiden(StringIO(data), check_closure=False)
        self.assertEqual(table["ionization"].cross_section_Mb[0], 4)
        with self.assertRaises(ValueError):
            load_leiden(
                StringIO("10 1e-18 0 2e-18\n20 2e-18 0 4e-18\n"), check_closure=False
            )

    def test_optical_interpolation_does_not_extrapolate(self):
        table = OpticalTable([10, 20], [2, 4])
        self.assertEqual(table(15), 3)
        with self.assertRaises(ValueError):
            table(21)
        with self.assertRaises(ValueError):
            table.integrate(10, 21)
        np.testing.assert_array_equal(table([9, 21], outside="zero"), [0, 0])
        with self.assertRaises(ValueError):
            OpticalTable([10, 10], [2, 4])
        with self.assertRaises(ValueError):
            OpticalTable([10, 20], [2, -1])

    def test_optical_integral_against_analytic_moments(self):
        table = OpticalTable([1, 2, 5, 10], [2, 3, 6, 11])
        a, b = 1.5, 8.3
        self.assertAlmostEqual(
            table.integrate(a, b), b - a + (b**2 - a**2) / 2, places=12
        )
        self.assertAlmostEqual(
            table.integrate(a, b, -1), np.log(b / a) + b - a, places=10
        )

    def test_nifs_scientific_tokens_and_column_guard(self):
        self.assertEqual(
            rows("1.000E+012.000E-01\n2.000E+013.000E-01", 2), [[10, 0.2], [20, 0.3]]
        )
        with self.assertRaises(ValueError):
            rows("1.000E+012.000E-01", 3)

    def test_band_gap_is_not_bridged_by_interpolation(self):
        data = {
            "continuum": [[1, 2], [2, 2], [4, 2], [5, 2]],
            "continuum_gaps_eV": [[2, 4]],
        }
        self.assertEqual(continuum_integral(data, np.ones_like), 4)
        self.assertEqual(continuum_integral(data, np.ones_like, 2, 4), 0)

    def test_log_interpolation_power_law(self):
        data = {"continuum": [[1, 1], [10, 0.01]], "continuum_gaps_eV": []}
        self.assertAlmostEqual(
            continuum_integral(data, np.ones_like, interpolation="log"), 0.9, places=7
        )

    def test_oscillator_lines_and_bands_are_not_density_samples(self):
        data = {
            "continuum": [[1, 0], [2, 0], [4, 0], [5, 0]],
            "continuum_gaps_eV": [[2, 4]],
            "lines": [[1.5, 2]],
            "bands": [[2, 4, 3]],
            "electrons": 5,
            "first_threshold_eV": 4,
        }
        result = response_moments(data)
        self.assertEqual(result["strength_through_100keV"], 5)
        expected = np.exp((2 * np.log(1.5) + 3 * np.log(3)) / 5)
        self.assertAlmostEqual(result["log_mean_excitation_eV"], expected)
        self.assertEqual(result["sub_ionization_strength"], 5)

    def test_rudd_table_is_in_range_and_not_capture(self):
        self.assertEqual(len(RUDD_1983_ENERGIES), 18)
        self.assertLessEqual(RUDD_1983_ENERGIES[-1], 4e6)
        self.assertEqual(RUDD_1983_FITTED_TOTALS["N2"][0], 2.11e-16)
        np.testing.assert_allclose(
            rudd_1983_uncertainty([5e3, 10e3, 25e3, 100e3, 1e6]),
            [0.25, 0.20, 0.15, 0.10, 0.08],
        )
        with self.assertRaises(ValueError):
            rudd_1983_uncertainty(5e6)

    def test_cheng_digitized_points_cover_three_spectra(self):
        data = cheng_1989_figure1()
        self.assertEqual(data.shape, (30, 3))
        self.assertTrue(np.all(data > 0))
        np.testing.assert_array_equal(np.unique(data[:, 0]), [7500, 50000, 150000])
        for energy in np.unique(data[:, 0]):
            group = data[data[:, 0] == energy]
            self.assertTrue(np.all(np.diff(group[:, 1]) > 0))


if __name__ == "__main__":
    unittest.main()
