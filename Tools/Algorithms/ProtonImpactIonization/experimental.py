# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Published electron-production measurements, not recommended model curves.

These small tables were transcribed from the typeset Table I in each paper.
Energies are in eV and cross sections in cm^2, as in reference.py. Preserve
the original values: do not silently rescale them to a later recommendation.
"""

import numpy as np

# Crooks and Rudd, Phys. Rev. A 3, 1628 (1971), Table I, p. 1632.
# https://doi.org/10.1103/PhysRevA.3.1628
# Printed cross-section unit: 10^-20 m^2 = 10^-16 cm^2.
# Their total uncertainty is about 17%, not an independent error on each
# point. Detector efficiency and gas pressure give shared systematic errors.
CROOKS_RUDD_1971_ENERGIES = 1000 * np.array([50, 100, 150, 200, 250, 300])
CROOKS_RUDD_1971_TOTALS = {
    "N2": 1.0e-16 * np.array([5.53, 5.37, 4.69, 4.13, 3.63, 3.24]),
    "O2": 1.0e-16 * np.array([5.18, 5.44, 4.91, 4.34, 3.87, 3.47]),
}

# Rudd, Phys. Rev. A 20, 787 (1979), Table I, p. 791.
# https://doi.org/10.1103/PhysRevA.20.787
# Nitrogen's first printed value is 2.46e-20 m^2; the hydrogen column has
# a different exponent and is deliberately not transcribed into this table.
# The N2 data were normalized at each angle to Crooks and Rudd's 50-keV
# measurements. They are therefore not an independent absolute calibration.
# Some lowest-energy measurements have uncertainties as large as a factor
# of two. The paper does not supply pointwise covariance or standard errors.
RUDD_1979_N2_ENERGIES = 1000 * np.array([5, 7, 10, 15, 20, 30, 50, 70])
RUDD_1979_N2_TOTALS = 1.0e-16 * np.array(
    [2.46, 2.79, 3.22, 3.78, 4.22, 4.82, 5.33, 5.33]
)
RUDD_1979_N2_MEAN_ENERGIES = np.array([5.42, 6.67, 7.74, 9.99, 11.9, 15.2, 20.4, 24.4])
