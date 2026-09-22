# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check experimental-table units, alignment, and independent provenance."""

import unittest

import numpy as np
from experimental import (
    CROOKS_RUDD_1971_ENERGIES,
    CROOKS_RUDD_1971_TOTALS,
    RUDD_1979_N2_ENERGIES,
    RUDD_1979_N2_MEAN_ENERGIES,
    RUDD_1979_N2_TOTALS,
)
from reference import rudd_total


class ExperimentalTableChecks(unittest.TestCase):
    def test_table_alignment_and_units(self):
        for values in CROOKS_RUDD_1971_TOTALS.values():
            self.assertEqual(values.shape, CROOKS_RUDD_1971_ENERGIES.shape)
            self.assertTrue(np.all(values > 0))
        self.assertEqual(RUDD_1979_N2_TOTALS.shape, RUDD_1979_N2_ENERGIES.shape)
        self.assertEqual(RUDD_1979_N2_MEAN_ENERGIES.shape, RUDD_1979_N2_ENERGIES.shape)
        self.assertEqual(CROOKS_RUDD_1971_ENERGIES[1], 100000)
        # Convert the original SI values independently, not via the table scale.
        self.assertAlmostEqual(CROOKS_RUDD_1971_TOTALS["O2"][1] / 5.44e-20, 10000)
        self.assertAlmostEqual(RUDD_1979_N2_TOTALS[0] / 2.46e-20, 10000)
        self.assertEqual(RUDD_1979_N2_MEAN_ENERGIES[0], 5.42)

    def test_original_values_are_not_renormalized_recommendations(self):
        # The low-energy disagreement is information to retain, not a reason
        # to change the transcription or silently overwrite a measurement.
        ratio = RUDD_1979_N2_TOTALS / rudd_total("N2", RUDD_1979_N2_ENERGIES)
        np.testing.assert_allclose(ratio[0], 1.2378418011637426, rtol=1.0e-12)
        self.assertNotEqual(RUDD_1979_N2_TOTALS[0], rudd_total("N2", 5000))


if __name__ == "__main__":
    unittest.main()
