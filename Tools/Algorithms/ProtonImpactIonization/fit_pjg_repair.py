# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Reproduce PJG-preserving fit diagnostics, without modifying production data.

Default: report the recorded five-parameter diagnostic. With --refit, fit
J alone, then J/p/K, then J/p/K/Gamma_s/Lambda. A dense recommended curve is
not a set of independent measurements. The original experimental tables are
reported separately and are not included in this diagnostic objective.
"""

import argparse
from dataclasses import astuple

import numpy as np
from experimental import (
    CROOKS_RUDD_1971_ENERGIES,
    CROOKS_RUDD_1971_TOTALS,
    RUDD_1979_N2_ENERGIES,
    RUDD_1979_N2_MEAN_ENERGIES,
    RUDD_1979_N2_TOTALS,
)
from pjg_repair import (
    DIAGNOSTIC_FITS,
    RepairParameters,
    initial_parameters,
    nonrelativistic_sdcs,
)
from reference import ELECTRON_REST_ENERGY, PROTON_REST_ENERGY, rudd_sdcs, rudd_total
from scipy.integrate import simpson
from scipy.optimize import least_squares


def integration_grid(energies, points=2049):
    x = np.linspace(0, 1, points)[None, :] * np.log1p(np.asarray(energies)[:, None])
    return x, np.expm1(x)


def moments(target, energies, parameters):
    x, t = integration_grid(energies)
    weighted = nonrelativistic_sdcs(
        target, np.asarray(energies)[:, None], t, parameters
    ) * np.exp(x)
    total = simpson(weighted, x=x, axis=-1)
    mean = simpson(t * weighted, x=x, axis=-1) / total
    return total, mean


def spectral_grid(target):
    energies = (
        [5e3, 10e3, 30e3, 50e3, 100e3, 300e3, 1e6]
        if target == "N2"
        else [7.5e3, 10e3, 30e3, 50e3, 100e3, 300e3]
    )
    se, st, sy = [], [], []
    for energy in energies:
        t = np.geomspace(
            2, min(8 * energy * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY, 3000), 48
        )
        y = rudd_sdcs(target, energy, t)
        keep = y > y[0] * 1.0e-4
        se.extend([energy] * np.count_nonzero(keep))
        st.extend(t[keep])
        sy.extend(y[keep])
    return np.array(se), np.array(st), np.array(sy)


def report(target, parameters):
    energies = np.geomspace(5e3, 4e6, 31)
    total, _ = moments(target, energies, parameters)
    se, st, sy = spectral_grid(target)
    ratio = nonrelativistic_sdcs(target, se, st, parameters) / sy
    print(target, parameters)
    print(
        "  recommended total ratio min/max:",
        np.min(total / rudd_total(target, energies)),
        np.max(total / rudd_total(target, energies)),
    )
    print(
        "  recommended SDCS ratio percentiles 0/10/50/90/100:",
        np.percentile(ratio, [0, 10, 50, 90, 100]),
    )
    experimental, _ = moments(target, CROOKS_RUDD_1971_ENERGIES, parameters)
    print(
        "  original 1971 total ratios:", experimental / CROOKS_RUDD_1971_TOTALS[target]
    )
    if target == "N2":
        experimental, mean = moments(target, RUDD_1979_N2_ENERGIES, parameters)
        print("  original 1979 total ratios:", experimental / RUDD_1979_N2_TOTALS)
        print("  original 1979 mean-energy ratios:", mean / RUDD_1979_N2_MEAN_ENERGIES)
    # Resolve the free endpoint as an integration boundary, rather than masking
    # cells that happen to straddle it in a full-spectrum quadrature.
    x, t = integration_grid([5000])
    maximum = 4 * 5000 * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY
    xtail = np.linspace(np.log1p(maximum), x[0, -1], 2049)
    tail = simpson(
        nonrelativistic_sdcs(target, 5000, np.expm1(xtail), parameters) * np.exp(xtail),
        x=xtail,
    )
    print(
        "  nonrelativistic free-endpoint tail at 5 keV:",
        tail / moments(target, [5000], parameters)[0][0],
    )


def refit(target, count):
    defaults = np.array(astuple(initial_parameters(target)))
    energies = np.geomspace(5e3, 4e6, 31)
    expected = rudd_total(target, energies)
    se, st, sy = spectral_grid(target)

    def unpack(values):
        parameters = defaults.copy()
        parameters[:count] = values * defaults[:count]
        return RepairParameters(*parameters)

    def residual(values):
        parameters = unpack(values)
        total, _ = moments(target, energies, parameters)
        totals = np.log(total / expected) / np.sqrt(len(energies))
        spectra = np.log(
            nonrelativistic_sdcs(target, se, st, parameters) / sy
        ) / np.sqrt(len(sy))
        return np.concatenate((totals, spectra))

    lower = np.array([0.01, 0.1, 0.01, 0.1, 1])[:count] / defaults[:count]
    upper = np.array([1e5, 4, 20, 100, 1000])[:count] / defaults[:count]
    fit = least_squares(
        residual,
        np.ones(count),
        bounds=(lower, upper),
        max_nfev=600,
        ftol=1e-9,
        xtol=1e-9,
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return unpack(fit.x)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refit", action="store_true")
    parser.add_argument("--target", choices=("N2", "O2"))
    args = parser.parse_args()
    for target in (args.target,) if args.target else ("N2", "O2"):
        if args.refit:
            for count in (1, 3, 5):
                print("Number of fitted parameters:", count, flush=True)
                report(target, refit(target, count))
        else:
            report(target, DIAGNOSTIC_FITS[target])


if __name__ == "__main__":
    main()
