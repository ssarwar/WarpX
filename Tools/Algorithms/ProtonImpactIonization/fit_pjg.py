# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Finalize the empirical proton PJG energy model, independently of angles.

Measured spectral shapes take priority over Rudd's recommended formula.
Shared normalizations and uncertain low-energy electron measurements are
not counted as independent precise constraints. The optical forward map
is approximate; the objective is a documented compromise, not chi-squared.
No production table, commit, or remote state is changed by this script.
"""

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from calibrated_pjg import PARAMETERS
from experimental import (
    CROOKS_RUDD_1971_ENERGIES,
    CROOKS_RUDD_1971_TOTALS,
    RUDD_1979_N2_ENERGIES,
    RUDD_1979_N2_MEAN_ENERGIES,
    RUDD_1979_N2_TOTALS,
)
from fit_support import SourceFitProblem, ratio_summary, summarize
from optical_reference import pjg_loss_density
from original_pjg import bhabha, legacy_cutoff, legacy_sdcs, moments
from pjg_model import MatchedParameters, SpectrumGrid
from pjg_moments import LossGrid, load_pstar, mass_stopping
from pjg_properties import distribution_properties, property_convergence
from source_datasets import (
    PHOTO_FACTOR,
    RUDD_1983_ENERGIES,
    RUDD_1983_FITTED_TOTALS,
    OpticalTable,
    cheng_1989_figure1,
    load_leiden,
    pjg_1976_figure5,
    rudd_1979_figure7,
    rudd_1983_uncertainty,
    source_digest,
)


def centered_shape_residual(ratio, weights):
    """Profile a shared log normalization, without an SDCS fit coefficient."""
    residual = np.log(ratio)
    return weights * (residual - np.average(residual, weights=weights**2))


class FinalFitProblem(SourceFitProblem):
    def __init__(self, target, optical):
        super().__init__(target, optical, optical_weight={"N2": 1, "O2": 2}[target])
        self.measured_spectra = []
        if target == "N2":
            data = pjg_1976_figure5()
            keys = np.unique(data[:, [0, 3]], axis=0)
            groups = [
                data[(data[:, 0] == e) & (data[:, 3] == source)] for e, source in keys
            ]
        else:
            data = cheng_1989_figure1()
            groups = [data[data[:, 0] == e] for e in np.unique(data[:, 0])]
        for group in groups:
            self.measured_spectra.append(
                (
                    SpectrumGrid(target, group[:, 0], group[:, 1]),
                    group[:, 2],
                    np.where(group[:, 1] < 10, 0.5, 1.0),
                )
            )

    def residual(self, parameters, center_s=None):
        total, _ = self.moments(parameters, center_s=center_s)
        blocks = [
            np.log(total / self.total_reference)
            * 0.15
            / self.total_uncertainty
            / np.sqrt(len(total)),
            # A smooth literature recommendation is a weak guide, not a
            # dense independent experiment that can overwhelm measured data.
            0.25
            * np.log(
                self.spectra(parameters, center_s=center_s) / self.spectral_reference
            )
            / np.sqrt(len(self.spectral_reference)),
            self.optical_weight
            * np.log(
                pjg_loss_density(self.target, self.loss, parameters, center_s)
                / self.optical_reference
            )
            / np.sqrt(len(self.loss)),
        ]
        crooks, _ = self.crooks(parameters, center_s=center_s)
        blocks.append(
            0.5 * np.log(crooks / CROOKS_RUDD_1971_TOTALS[self.target]) / np.sqrt(6)
        )
        if self.original_moments:
            original, mean = self.original_moments(parameters, center_s=center_s)
            # The source explicitly warns of larger uncertainties at low
            # proton energy. These are conservative relative block weights,
            # not fabricated pointwise experimental standard deviations.
            confidence = np.interp(
                np.log(RUDD_1979_N2_ENERGIES),
                np.log([5e3, 10e3, 30e3]),
                [0.5, 2 / 3, 1],
            )
            blocks += [
                0.4 * confidence * np.log(original / RUDD_1979_N2_TOTALS) / np.sqrt(8),
                confidence * np.log(mean / RUDD_1979_N2_MEAN_ENERGIES) / np.sqrt(8),
            ]
        count = len(self.measured_spectra)
        for grid, measured, confidence in self.measured_spectra:
            ratio = grid(parameters, center_s=center_s) / measured
            blocks.append(
                centered_shape_residual(ratio, confidence)
                / np.sqrt(count * len(measured))
            )
        return np.concatenate(blocks)


def measured_checks(target, parameters):
    sets = (
        {
            "pjg_figure5_measurements": pjg_1976_figure5(),
            "rudd1979_held_out_sdcs": rudd_1979_figure7(),
        }
        if target == "N2"
        else {"cheng1989_measurements": cheng_1989_figure1()}
    )
    output = {}
    for name, data in sets.items():
        records = []
        for energy in np.unique(data[:, 0]):
            group = data[data[:, 0] == energy]
            ratio = (
                SpectrumGrid(target, group[:, 0], group[:, 1])(parameters) / group[:, 2]
            )
            records.append(
                {
                    "energy_eV": float(energy),
                    "points": len(group),
                    "ratio": ratio_summary(ratio),
                    "shape_log_rms": float(np.std(np.log(ratio))),
                    "points_eV_cm2_per_eV_and_model_ratio": np.c_[
                        group[:, 1:3], ratio
                    ].tolist(),
                }
            )
        output[name] = records
    return output


def fit_models(optical, frozen=False):
    result = {
        "status": "Calibrated production energy model; angular closure outside the fit.",
        "targets": {},
    }
    for target in ("N2", "O2"):
        problem = FinalFitProblem(target, optical[target])
        pstar = load_pstar(target)
        keep = pstar.energies >= 5e3
        losses = LossGrid(target, pstar.energies[keep])
        initial = PARAMETERS[target]
        records = {}
        for name, active, start in (
            ("six", [0, 2, 3, 4, 6, 8], initial),
            (
                "five_zero_center_numerator",
                [0, 2, 3, 4, 8],
                replace(initial, center_scale=0),
            ),
            ("five_fixed_edge", [0, 2, 3, 4, 6], replace(initial, edge_scale=1)),
        ):
            if frozen and name != {"N2": "six", "O2": "five_fixed_edge"}[target]:
                continue
            p, fit = (initial, None) if frozen else problem.fit(active, start)
            record = summarize(problem, p, losses, pstar.electronic[keep])
            record.update(
                optimizer_success=bool(fit.success) if fit is not None else None,
                evaluations=int(fit.nfev) if fit is not None else 0,
                frozen=frozen,
                active_parameter_count=len(active),
                measured_checks=measured_checks(target, p),
            )
            records[name] = record
            print(
                target,
                name,
                "objective",
                record["objective_norm"],
                "max pair/NIST",
                record["maximum_pair_to_pstar"],
                "parameters",
                record["parameters"],
                flush=True,
            )
        result["targets"][target] = records
    return result


def add_validation(result, optical):
    for target, records in result["targets"].items():
        # O2's unrestricted edge multiplier is 1.00073: fixing it to the
        # literature value removes an unneeded coefficient. N2's 1.19114
        # multiplier has a material effect on its measured spectral shape.
        selected = {"N2": "six", "O2": "five_fixed_edge"}[target]
        records["selected"] = selected
        record = records[selected]
        p = MatchedParameters(**record["parameters"])
        problem = FinalFitProblem(target, optical[target])
        alternate, fit = problem.fit(
            [0, 2, 3, 4, 6, 8] if selected == "six" else [0, 2, 3, 4, 6],
            replace(
                p, j=p.j * 1.5, width=p.width * 0.8, broad_excess=p.broad_excess * 1.2
            ),
        )
        record["independent_start"] = {
            "success": bool(fit.success),
            "max_spectral_relative_change": float(
                np.max(np.abs(problem.spectra(alternate) / problem.spectra(p) - 1))
            ),
        }
        record["properties"] = distribution_properties(target, p)
        record["numerical_convergence"] = property_convergence(target, p)
        energy = 800e6
        t = np.array([1e3, 1e4, 1e5])
        record["hard_tail_ratio"] = (
            SpectrumGrid(target, energy, t)(p) / bhabha(target, energy, t)
        ).tolist()
        table = load_pstar(target)
        energy = table.energies[table.energies >= 5000]
        m = LossGrid(target, energy, points=4097)(p)
        record["stopping_rows"] = np.c_[
            energy,
            mass_stopping(target, m.kinetic),
            mass_stopping(target, m.binding),
            mass_stopping(target, m.ionization),
            table.electronic[table.energies >= 5000],
        ].tolist()
        if target == "N2":
            total, mean = problem.original_moments(p)
            record["measured_mean_rows"] = np.c_[
                RUDD_1979_N2_ENERGIES, RUDD_1979_N2_MEAN_ENERGIES, mean, total
            ].tolist()


STYLES = (
    ("printed", "Original PJG + 1977 erratum", "#777777", ":"),
    ("final", "Final matched PJG", "#0072B2", "-"),
)


def selected_parameters(result, target):
    records = result["targets"][target]
    return MatchedParameters(**records[records["selected"]]["parameters"])


def plot_all(result, optical, directory):
    import matplotlib.pyplot as plt

    directory.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.framealpha": 0.96,
        }
    )
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12, 8.2),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.5, 1]},
    )
    e = np.geomspace(5e3, 4e6, 180)
    for col, target in enumerate(("N2", "O2")):
        p = selected_parameters(result, target)
        for variant, label, color, style in STYLES:
            total = (
                LossGrid(target, e, 1025)(p).total
                if variant == "final"
                else moments(target, e, variant, points=2049)[0]
            )
            axes[0, col].loglog(
                e / 1000,
                total * 1e16,
                color=color,
                linestyle=style,
                linewidth=2,
                label=label,
            )
        ref = RUDD_1983_FITTED_TOTALS[target]
        uncertainty = rudd_1983_uncertainty(RUDD_1983_ENERGIES)
        axes[0, col].fill_between(
            RUDD_1983_ENERGIES / 1000,
            ref * (1 - uncertainty) * 1e16,
            ref * (1 + uncertainty) * 1e16,
            color="#009E73",
            alpha=0.12,
            label="Rudd 1983 fit: quoted uncertainty",
        )
        axes[0, col].plot(
            RUDD_1983_ENERGIES / 1000,
            ref * 1e16,
            color="#009E73",
            linestyle="-.",
            linewidth=1,
            label="Rudd 1983 authors' fit (not raw points)",
        )
        axes[0, col].scatter(
            CROOKS_RUDD_1971_ENERGIES / 1000,
            CROOKS_RUDD_1971_TOTALS[target] * 1e16,
            facecolor="white",
            edgecolor="black",
            marker="o",
            s=34,
            label="Crooks & Rudd 1971 measurements",
            zorder=5,
        )
        if target == "N2":
            axes[0, col].scatter(
                RUDD_1979_N2_ENERGIES / 1000,
                RUDD_1979_N2_TOTALS * 1e16,
                facecolor="white",
                edgecolor="#D55E00",
                marker="s",
                s=30,
                label="Rudd 1979 measurements",
                zorder=5,
            )
        axes[0, col].set(
            title=target,
            ylabel=r"Electron-production cross section ($10^{-16}$ cm$^2$)",
            xlim=(5, 4000),
        )
        for variant, label, color, style in STYLES:
            total = (
                LossGrid(target, RUDD_1983_ENERGIES, 1025)(p).total
                if variant == "final"
                else moments(target, RUDD_1983_ENERGIES, variant)[0]
            )
            axes[1, col].semilogx(
                RUDD_1983_ENERGIES / 1000, total / ref, color=color, linestyle=style
            )
        axes[1, col].fill_between(
            RUDD_1983_ENERGIES / 1000,
            1 - uncertainty,
            1 + uncertainty,
            color="#009E73",
            alpha=0.12,
        )
        axes[1, col].axhline(1, color="black", linewidth=0.7)
        axes[1, col].set(
            xlabel="Proton kinetic energy (keV)",
            ylabel="Model / Rudd 1983 fit",
            ylim=(0, 2.1),
            xlim=(5, 4000),
        )
        for row in range(2):
            axes[row, col].grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7.7)
    fig.suptitle(
        "N₂ and O₂: final energy model versus the original PJG and measured yields"
    )
    fig.savefig(directory / "total_cross_sections.png", dpi=170)
    plt.close(fig)

    for target, energies in (
        ("N2", [5e3, 20e3, 50e3, 70e3, 300e3, 1e6]),
        ("O2", [7.5e3, 50e3, 150e3]),
    ):
        fig, axes = plt.subplots(
            2 if target == "N2" else 1,
            3,
            figsize=(13, 8 if target == "N2" else 4.5),
            constrained_layout=True,
        )
        axes = np.asarray(axes).ravel()
        p = selected_parameters(result, target)
        for axis, energy in zip(axes, energies, strict=True):
            observed = (
                np.concatenate([pjg_1976_figure5()[:, :3], rudd_1979_figure7()])
                if target == "N2"
                else cheng_1989_figure1()
            )
            observed = observed[observed[:, 0] == energy]
            upper = (
                130
                if energy <= 7500
                else 250
                if energy <= 70e3
                else 1000
                if energy <= 300e3
                else 1500
            )
            upper = max(upper, 1.3 * observed[:, 1].max())
            t = np.geomspace(0.7, upper, 1300)
            model_maximum = 0
            for variant, label, color, style in STYLES:
                y = (
                    SpectrumGrid(target, energy, t)(p)
                    if variant == "final"
                    else legacy_sdcs(target, energy, t, variant)
                )
                axis.loglog(
                    t,
                    np.where(y > 0, y, np.nan),
                    label=label,
                    color=color,
                    linestyle=style,
                    linewidth=1.8,
                )
                model_maximum = max(model_maximum, float(np.max(y)))
                if variant != "final":
                    cutoff = float(legacy_cutoff(energy, variant))
                    if cutoff < upper:
                        axis.axvline(
                            cutoff,
                            color=color,
                            linestyle=style,
                            linewidth=0.6,
                            alpha=0.45,
                        )
            if target == "N2":
                data = pjg_1976_figure5()
                for source, marker, label in (
                    (0, "o", "Crooks/Rudd: markers in PJG Fig. 5"),
                    (1, "s", "Toburen: markers in PJG Fig. 5"),
                ):
                    group = data[(data[:, 0] == energy) & (data[:, 3] == source)]
                    if len(group):
                        axis.scatter(
                            group[:, 1],
                            group[:, 2],
                            marker=marker,
                            facecolor="white",
                            edgecolor="black",
                            s=31,
                            linewidth=1,
                            label=label,
                            zorder=5,
                        )
                data = rudd_1979_figure7()
                group = data[data[:, 0] == energy]
                if len(group):
                    axis.scatter(
                        group[:, 1],
                        group[:, 2],
                        marker="D",
                        facecolor="white",
                        edgecolor="#D55E00",
                        s=28,
                        label="Rudd 1979 measured markers (held out)",
                        zorder=5,
                    )
            else:
                data = cheng_1989_figure1()
                group = data[data[:, 0] == energy]
                axis.scatter(
                    group[:, 1],
                    group[:, 2],
                    facecolor="white",
                    edgecolor="black",
                    s=31,
                    label="Cheng 1989 measured markers (rescaled)",
                    zorder=5,
                )
            measured = observed[:, 2]
            # Retain every measured marker and the large original-PJG
            # excess. A data-only vertical range would hide its 5-keV curve.
            axis.set(
                title=f"{target}, {energy / 1000:g} keV protons",
                xlabel="Secondary electron energy T (eV)",
                ylabel=r"SDCS (cm$^2$/eV per molecule)",
                xlim=(0.7, upper),
                ylim=(
                    max(1e-24, measured.min() / 6),
                    max(measured.max() * 5, model_maximum * 1.4),
                ),
            )
            axis.grid(alpha=0.2, which="both")
            axis.legend(fontsize=6.8, loc="lower left")
        fig.suptitle(
            "Measured SDCS markers; no recommended-formula curves used as experimental points\nOld curves end at their cutoffs; zero/negative values cannot appear on logarithmic axes",
            fontsize=11,
        )
        fig.savefig(directory / f"{target.lower()}_sdcs.png", dpi=170)
        plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), constrained_layout=True)
    for target, color in (("N2", "#0072B2"), ("O2", "#D55E00")):
        records = result["targets"][target]
        record = records[records["selected"]]
        p = selected_parameters(result, target)
        w = np.geomspace(25, 100, 601)
        axes[0].semilogx(
            w,
            PHOTO_FACTOR * pjg_loss_density(target, w, p) / optical[target](w),
            color=color,
            label=target,
        )
        rows = np.asarray(record["stopping_rows"])
        axes[1].semilogx(
            rows[:, 0] / 1e6, rows[:, 3] / rows[:, 4], color=color, label=target
        )
        if target == "N2":
            rows = np.asarray(record["measured_mean_rows"])
            axes[2].semilogx(
                rows[:, 0] / 1000, rows[:, 2], color=color, label="Final N₂ mean"
            )
            axes[2].scatter(
                rows[:, 0] / 1000,
                rows[:, 1],
                facecolor="white",
                edgecolor="black",
                label="Rudd 1979 measured mean",
            )
    axes[0].set(
        xlabel="Photon / energy-loss W (eV)",
        ylabel="Model optical response / evaluated data",
    )
    axes[1].set(
        xlabel="Proton energy (MeV)",
        ylabel="Effective-pair loss / NIST electronic",
        ylim=(0, 1.08),
    )
    axes[2].set(xlabel="Proton energy (keV)", ylabel="Mean secondary energy (eV)")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    for axis in axes[:2]:
        axis.axhline(1, color="black", linewidth=0.7, linestyle=":")
    fig.suptitle(
        "Final calibration: optical shape, independent stopping budget, and measured N₂ mean"
    )
    fig.savefig(directory / "physics_validation.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nifs", type=Path, required=True)
    parser.add_argument("--o2-leiden", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    parser.add_argument("--read-results", action="store_true")
    parser.add_argument(
        "--frozen",
        action="store_true",
        help="Evaluate and plot the committed calibration without refitting it",
    )
    args = parser.parse_args()
    table = np.asarray(json.loads(args.nifs.read_text())["N2"]["continuum"])
    optical = {
        "N2": OpticalTable(table[:, 0], table[:, 2]),
        "O2": load_leiden(args.o2_leiden)["ionization"],
    }
    if args.read_results:
        result = json.loads(args.output.read_text())
    else:
        result = fit_models(optical, frozen=args.frozen)
        add_validation(result, optical)
        result["source_sha256"] = {
            name: source_digest(getattr(args, name))
            for name in ("nifs", "o2_leiden")
            if getattr(args, name) is not None
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    plot_all(result, optical, args.figures)


if __name__ == "__main__":
    main()
