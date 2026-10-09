# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Shared grids and log-parameter minimization for the final joint refit."""

from dataclasses import asdict, astuple

import numpy as np
from experimental import (
    CROOKS_RUDD_1971_ENERGIES,
    CROOKS_RUDD_1971_TOTALS,
    RUDD_1979_N2_ENERGIES,
    RUDD_1979_N2_MEAN_ENERGIES,
    RUDD_1979_N2_TOTALS,
)
from optical_reference import pjg_loss_density
from pjg_model import (
    MatchedParameters,
    SpectrumGrid,
    hard_factor,
    initial_parameters,
    kinematics,
    molecular_endpoint,
)
from pjg_moments import free_second_moment, mass_stopping
from reference import (
    BETHE_CONSTANT,
    ELECTRON_REST_ENERGY,
    PROTON_REST_ENERGY,
    rudd_sdcs,
    rudd_total,
)
from scipy.integrate import simpson
from scipy.optimize import least_squares
from source_datasets import (
    PHOTO_FACTOR,
    RUDD_1983_ENERGIES,
    RUDD_1983_FITTED_TOTALS,
    cheng_1989_figure1,
    rudd_1983_uncertainty,
)
from target_parameters import TARGETS


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


class MomentGrid:
    def __init__(self, target, energies, points=1025, continuation="exponential"):
        self.energy = np.asarray(energies, dtype=float)
        maximum = molecular_endpoint(
            target, self.energy, min(TARGETS[target].thresholds)
        )
        self.x = np.linspace(0, 1, points)[None, :] * np.log1p(maximum[:, None])
        self.t = np.expm1(self.x)
        self.grid = SpectrumGrid(target, self.energy[:, None], self.t, continuation)

    def __call__(self, parameters, hard_distortion="fading", center_s=None):
        weighted = self.grid(parameters, hard_distortion, center_s) * np.exp(self.x)
        total = simpson(weighted, x=self.x, axis=-1)
        mean = simpson(self.t * weighted, x=self.x, axis=-1) / total
        return total, mean


class FitProblem:
    def fit(self, active, initial=None):
        defaults = np.array(astuple(initial or initial_parameters(self.target)))
        active = np.asarray(active)
        lower = np.array([0.01, 0.1, 0.01, 0.1, 1, 1e-4, 1e-4, 0.02, 0.25])[active]
        upper = np.array([1e5, 4, 10, 100, 1000, 100, 5, 5, 4])[active]

        def unpack(values):
            parameters = defaults.copy()
            parameters[active] = np.exp(values)
            return MatchedParameters(*parameters)

        fit = least_squares(
            lambda values: self.residual(unpack(values)),
            np.log(defaults[active]),
            bounds=(np.log(lower), np.log(upper)),
            max_nfev=600,
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-9,
        )
        if not fit.success:
            raise RuntimeError(fit.message)
        return unpack(fit.x), fit


def ratio_summary(ratio):
    ratio = np.asarray(ratio)
    return {
        "minimum": float(ratio.min()),
        "maximum": float(ratio.max()),
        "log_rms": float(np.sqrt(np.mean(np.log(ratio) ** 2))),
    }


class SourceFitProblem(FitProblem):
    """Group-weighted residuals, not independent measurements or chi-squared."""

    def __init__(self, target, optical, optical_weight=1, mean_weight=0.5):
        self.target = target
        self.hard_distortion = "fading"
        self.optical_weight = optical_weight
        # A residual-block weight, not an additional SDCS parameter or an
        # experimental inverse uncertainty.
        self.mean_weight = mean_weight
        self.energies = RUDD_1983_ENERGIES
        self.total_reference = RUDD_1983_FITTED_TOTALS[target]
        self.total_uncertainty = rudd_1983_uncertainty(self.energies)
        self.moments = MomentGrid(target, self.energies)
        se, st, self.spectral_reference = spectral_grid(target)
        self.spectra = SpectrumGrid(target, se, st)
        self.loss = np.geomspace(25, 100, 60)
        self.optical_reference = optical(self.loss) / PHOTO_FACTOR
        self.crooks = MomentGrid(target, CROOKS_RUDD_1971_ENERGIES)
        self.original_moments = (
            MomentGrid(target, RUDD_1979_N2_ENERGIES) if target == "N2" else None
        )
        self.cheng = []
        if target == "O2":
            points = cheng_1989_figure1()
            for energy in np.unique(points[:, 0]):
                group = points[points[:, 0] == energy]
                self.cheng.append(
                    (SpectrumGrid(target, energy, group[:, 1]), group[:, 2])
                )


def summarize(problem, parameters, loss_grid, pstar):
    total, _ = problem.moments(parameters)
    losses = loss_grid(parameters)
    stopping = mass_stopping(problem.target, losses.ionization) / pstar
    maximum = int(np.argmax(stopping))
    crooks, _ = problem.crooks(parameters)
    energy = 800e6
    t = np.array([1e3, 1e4, 1e5])
    _, _, equivalent, _ = kinematics(energy)
    bhabha = TARGETS[problem.target].electrons * BETHE_CONSTANT / equivalent
    bhabha = bhabha * hard_factor(energy, t) / t**2
    result = {
        "parameters": asdict(parameters),
        "objective_norm": float(np.linalg.norm(problem.residual(parameters))),
        "rudd_1983_fitted_total_ratio": ratio_summary(total / problem.total_reference),
        "rudd_1985_recommended_total_ratio": ratio_summary(
            total / rudd_total(problem.target, problem.energies)
        ),
        "crooks_1971_measured_total_ratio": ratio_summary(
            crooks / CROOKS_RUDD_1971_TOTALS[problem.target]
        ),
        "rudd_1992_formula_sdcs_ratio": ratio_summary(
            problem.spectra(parameters) / problem.spectral_reference
        ),
        "optical_fixed_channel_loss_density_ratio_25_100eV": ratio_summary(
            pjg_loss_density(problem.target, problem.loss, parameters)
            / problem.optical_reference
        ),
        "maximum_pair_to_pstar": float(stopping[maximum]),
        "maximum_pair_to_pstar_energy_eV": float(loss_grid.energies[maximum]),
        "pair_to_pstar_rows": np.c_[loss_grid.energies, stopping].tolist(),
        "hard_800MeV_T_1_10_100keV_ratio": (
            SpectrumGrid(problem.target, energy, t)(parameters) / bhabha
        ).tolist(),
        "second_moment_to_free_at_10GeV": float(
            losses.kinetic_second[-1]
            / free_second_moment(problem.target, loss_grid.energies[-1])
        ),
    }
    if problem.original_moments:
        total, mean = problem.original_moments(parameters)
        result["rudd_1979_measured_total_ratio"] = ratio_summary(
            total / RUDD_1979_N2_TOTALS
        )
        result["rudd_1979_measured_mean_energy_ratio"] = ratio_summary(
            mean / RUDD_1979_N2_MEAN_ENERGIES
        )
    result["cheng_1989_rescaled_data_ratios"] = [
        {
            "energy_eV": float(grid.energy[0]),
            "ratio": ratio_summary(grid(parameters) / reference),
            "shape_log_rms": float(np.std(np.log(grid(parameters) / reference))),
        }
        for grid, reference in problem.cheng
    ]
    return result
