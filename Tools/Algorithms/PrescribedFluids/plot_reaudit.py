#!/usr/bin/env python3
"""Plot selected archived studies without mixing builds or mesh resolutions."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter


def summaries(report, jobs, phase):
    rows = []
    for job in jobs:
        for study in report[str(job)]["studies"]:
            if not study["directory"].startswith(phase + "-"):
                continue
            for name, data in study["comparisons"].items():
                if name.endswith("summary.json"):
                    rows.append(data)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument(
        "--noise-job",
        "--noise-jobs",
        dest="noise_jobs",
        type=int,
        nargs="+",
        required=True,
    )
    parser.add_argument("--scaling-jobs", type=int, nargs=2, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    libraries = {
        study["provenance"]["python-library-sha256.txt"]
        for job in [*args.noise_jobs, *args.scaling_jobs]
        for study in report[str(job)]["studies"]
    }
    assert len(libraries) == 1, "Selected studies must use the same compiled libraries"
    noise_summaries = summaries(report, args.noise_jobs, "noise")
    noise = {summary["settings"]["mode"]: summary for summary in noise_summaries}
    assert len(noise) == len(noise_summaries), "Select only one study per noise family"
    scaling = summaries(report, args.scaling_jobs, "scaling")
    assert all(len(summary["ranks"]) == 1 for summary in scaling)
    assert noise["fields"]["settings"]["cells"] == noise["coupled"]["settings"]["cells"]
    assert noise["fields"]["settings"]["cells"] == [64, 256]
    assert noise["coupled"]["settings"]["mcc_all"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.3), constrained_layout=True)
    for beam, label, color in [
        ("fluid", "Fluid beam", "C0"),
        ("quiet", "Quiet particle beam", "C1"),
        ("random", "Random particle beam", "C2"),
    ]:
        for axis, mode, metric in [
            (axes[0, 0], "fields", "beam_relative_rms_error"),
            (axes[0, 1], "coupled", "electrons_relative_seed_noise"),
        ]:
            groups = sorted(
                (
                    row
                    for row in noise[mode]["groups"]
                    if row["beam"] == beam and row["ions"] == "fluid"
                ),
                key=lambda row: row["ppc"],
            )
            assert groups and all(row["timestep_s"]["samples"] >= 6 for row in groups)
            x = [1000 * row["timestep_s"]["mean"] for row in groups]
            y = [
                100 * (row[metric]["mean"] if mode == "fields" else row[metric])
                for row in groups
            ]
            axis.errorbar(
                x,
                y,
                xerr=[1000 * row["timestep_s"]["ci99_half_width"] for row in groups],
                fmt="o-",
                color=color,
                label=label,
            )
            if beam != "fluid":
                for xx, yy, row in zip(x, y, groups):
                    offset = (3, 5)
                    if row["ppc"] == 256:
                        offset = (-22, -14) if beam == "quiet" else (-22, 9)
                    axis.annotate(
                        str(row["ppc"]),
                        (xx, yy),
                        xytext=offset,
                        textcoords="offset points",
                        fontsize=8,
                    )
            axis.set(xlabel="Timestep wall time (ms)", xscale="log")
            axis.margins(x=0.12, y=0.17)
    axes[0, 0].set(
        title="Beam sampling error: fields only",
        ylabel="Density RMS error relative to fluid (%)",
    )
    axes[0, 0].set_yscale("symlog", linthresh=0.001)
    axes[0, 0].set_ylim(bottom=-0.0002)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set(
        title="PJG + all electron MCC channels", ylabel="Plasma-density seed noise (%)"
    )
    for axis, mode in [(axes[1, 0], "push"), (axes[1, 1], "coupled")]:
        selected = [row for row in scaling if row["settings"]["mode"] == mode]
        assert len(selected) == 4 and {row["ranks"][0] for row in selected} == {
            1,
            2,
            4,
            8,
        }
        assert len({tuple(row["settings"]["cells"]) for row in selected}) == 1
        assert selected[0]["settings"]["cells"] == [256, 1024]
        for beam, ions, label, color, style in [
            ("fluid", "fluid", "Fluid beam / fluid ions", "C0", "o-"),
            ("quiet", "fluid", "Quiet beam / fluid ions", "C1", "o-"),
            ("fluid", "frozen", "Fluid beam / frozen ions", "C0", "s--"),
            ("quiet", "frozen", "Quiet beam / frozen ions", "C1", "s--"),
        ]:
            if mode == "push" and ions == "frozen":
                continue
            data = sorted(
                (row["ranks"][0], group["timestep_s"])
                for row in selected
                for group in row["groups"]
                if group["beam"] == beam and group["ions"] == ions
            )
            assert len(data) == 4
            assert all(value["samples"] >= 3 for _, value in data)
            axis.errorbar(
                [rank for rank, _ in data],
                [1000 * value["mean"] for _, value in data],
                yerr=[1000 * value["ci99_half_width"] for _, value in data],
                fmt=style,
                color=color,
                label=label,
            )
        axis.set(
            xlabel="A100 GPUs (8 GPUs = 2 nodes)",
            ylabel="Timestep wall time (ms)",
            xscale="log",
            yscale="log",
        )
        axis.set_xticks([1, 2, 4, 8], [1, 2, 4, 8])
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.yaxis.set_major_formatter(ScalarFormatter())
        axis.yaxis.set_minor_formatter(ScalarFormatter())
        axis.legend(fontsize=8)
    axes[1, 0].set_title("Strong scaling: 1,048,576 identical electrons")
    axes[1, 1].set_title("Strong scaling: coupled production and MCC")
    for axis in axes.flat:
        axis.grid(True, alpha=0.25)
    fig.suptitle(
        "Perlmutter A100, RZ: measured noise and cost\n"
        "Top: 64 × 256, 6 seeds; bottom: 256 × 1024, 3 seeds, quiet beam 64 particles/cell\n"
        "Bars: 99% timing intervals; labels: beam particles/cell; synthetic MCC regression rates",
        fontsize=11,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        fig.savefig(args.output.with_suffix("." + extension), dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
