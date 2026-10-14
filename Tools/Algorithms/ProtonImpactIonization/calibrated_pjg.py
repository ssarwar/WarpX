# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Frozen N2/O2 proton energy-model calibration of September 2026.

See the proton_impact_ionization.rst theory documentation for the formulation,
sources, uncertainty treatment and limitations. This is an offline reference,
not the production GPU sampler.
The observable is inclusive electron yield represented by effective pairs.
Incident energies above 4 MeV are extrapolations; the interface is limited
to the 5 keV--10 GeV numerical audit range. No secondary-energy cut is used.
"""

import numpy as np
from pjg_model import PARAMETERS
from pjg_model import sdcs as reference_sdcs
from pjg_moments import LossGrid


def checked_energies(energy):
    """Reject incident energies outside the numerical audit range."""
    energy = np.asarray(energy, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any((energy < 5e3) | (energy > 1e10)):
        raise ValueError("The calibrated reference requires 5 keV <= E <= 10 GeV")
    return energy


def sdcs(target, energy, secondary):
    """Inclusive electron SDCS, cm^2/eV per molecule; both energies in eV."""
    return reference_sdcs(
        target, checked_energies(energy), secondary, PARAMETERS[target]
    )


def total_cross_section(target, energy, points=2049):
    """Integrate this same SDCS, in cm^2; no independent normalization fit.

    The segmented quadrature is a host-side reference for table generation.
    Calling it separately for every PIC collision would be inappropriate.
    """
    energy = checked_energies(energy)
    result = LossGrid(target, energy.ravel(), points)(PARAMETERS[target]).total
    return result.reshape(energy.shape)
