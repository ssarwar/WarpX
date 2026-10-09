# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Reproduce the inferred N2 rank fractions from the readable Jung digitization."""

import json
import os
from pathlib import Path

import numpy as np
from rotation_reference import ROTATION, cg_squared, populations
from scipy.optimize import nnls


def infer(data):
    population, _ = populations("N2", data["temperature_K"], 96)
    branches = dict(zip(data["integral_branches_delta_J"], data["integral_branches"]))
    result = {}
    for energy, values in data["dcs"].items():
        energy = float(energy)
        values = np.asarray(values)
        response = np.zeros((2, 2))
        for row, delta_j in enumerate([2, 4]):
            for column, rank in enumerate([2, 4]):
                for initial, weight in enumerate(population):
                    final = initial + delta_j
                    loss = ROTATION["N2"] * (
                        final * (final + 1) - initial * (initial + 1)
                    )
                    response[row, column] += (
                        weight
                        * cg_squared(initial, rank, final)
                        * np.sqrt(max(1 - loss / energy, 0))
                    )
        strengths = np.column_stack([nnls(response, row)[0] for row in values[1:].T])
        assert (
            np.max(abs(response @ strengths - values[1:]))
            < data["minimum_manual_reading_uncertainty_in_dcs_units"]
        )
        inclusive = values[0] + sum(
            (1 + branches[-j] / branches[j]) * values[k]
            for k, j in enumerate([2, 4], 1)
        )
        unchanged = inclusive - strengths.sum(axis=0)
        assert np.all(unchanged > 0)
        result[str(energy)] = (np.vstack((unchanged, strengths)) / inclusive).tolist()
    return result


if __name__ == "__main__":
    directory = (
        Path(os.environ["WARPX_CROSS_SECTION_DATA"]) / "N2/IAA/reciprocal_sources"
    )
    inferred = infer(json.loads((directory / "n2_jung_digitization.json").read_text()))
    reference = json.loads((directory / "n2_jung.json").read_text())
    for energy, values in inferred.items():
        np.testing.assert_allclose(values, reference[energy], rtol=1e-12, atol=1e-14)
    print("PASS: Jung rank fractions reproduced from the digitized branches")
