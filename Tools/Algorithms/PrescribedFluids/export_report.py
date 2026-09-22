#!/usr/bin/env python3
"""Export a compact record and plots from the completed Perlmutter studies.

Keep fields/, push/, source/, coupled/, storage/storage-{20,100}/ and
profile/profile/ below the input directory, preserving each result.json,
result.npz, run.log and summary.json. Full particle checkpoints are not needed.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
from analyze import interval


def read(path):
    return json.loads(path.read_text())


def profile_rows(directory):
    rows = []
    pattern = re.compile(
        r"^ProtonImpactIonizationCollision::\S+\s+(\d+)"
        r"\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+[\d.]+%$",
        re.MULTILINE,
    )
    for path in sorted(directory.glob("*/result.json")):
        case = read(path)
        matches = pattern.findall(path.with_name("run.log").read_text())
        if not matches:
            raise ValueError(f"No source profiler record: {path}")
        # TinyProfiler prints exclusive and inclusive tables, sometimes again
        # by region. Inclusive time is the maximum for this one-rank study.
        rows.append(
            dict(
                case=path.parent.name,
                collision_calls=int(matches[0][0]),
                inclusive_source_time_s=max(float(row[3]) for row in matches),
                first_step_s=case["timestep_s"][0],
                ordinary_step_median_s=case["ordinary_step_median_s"],
                initialization_s=case["initialization_s"],
                beam_macroparticles=case["history"][-1]["macroparticles"]["beam"],
                electron_macroparticles=case["history"][-1]["macroparticles"][
                    "electrons"
                ],
                revision=case["revision"],
                warpx_version=case["warpx_version"],
            )
        )
    if not rows:
        raise ValueError(f"Missing source profiling study: {directory}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--output", type=Path, required=True, help="Output filename stem"
    )
    parser.add_argument("--date", required=True)
    parser.add_argument("--hardware", required=True)
    parser.add_argument("--compiler", required=True)
    args = parser.parse_args()
    report = dict(
        date=args.date,
        hardware=args.hardware,
        compiler=args.compiler,
        confidence="Two-sided 99% Student intervals across runs; timestep statistic "
        "is the mean of per-run medians after two warm-up steps.",
        ensembles={},
    )
    for phase in ["fields", "push", "source", "coupled"]:
        summary = read(args.directory / phase / "summary.json")
        expected_runs = sum(row["timestep_s"]["samples"] for row in summary["groups"])
        if len(list((args.directory / phase).glob("*/result.json"))) != expected_runs:
            raise ValueError(f"Incomplete {phase} ensemble")
        reference = next(
            row
            for row in summary["groups"]
            if row["beam"] == "fluid" and row["ions"] == "fluid"
        )
        for row in summary["groups"]:
            if "electron_cdf" in row:
                row["electron_cdf_distance_to_fluid"] = float(
                    np.max(
                        np.abs(
                            np.array(row["electron_cdf"]) - reference["electron_cdf"]
                        )
                    )
                )
        for row in summary["groups"]:
            row.pop("electron_cdf", None)
        report["ensembles"][phase] = summary

    cache = args.directory / "acceptance" / "CMakeCache.txt"
    report["build_configuration"] = {
        key: value
        for key, value in re.findall(
            r"^(WarpX_\w+):\w+=(.*)$", cache.read_text(), re.MULTILINE
        )
        if key
        in [
            "WarpX_COMPUTE",
            "WarpX_DIMS",
            "WarpX_PRECISION",
            "WarpX_PARTICLE_PRECISION",
            "WarpX_MPI",
            "WarpX_FFT",
            "WarpX_OPENPMD",
            "WarpX_PYTHON",
        ]
    }

    source = args.directory / "source"
    density_error, histogram_error, electron_budget_error, recoil_error = [], [], [], []
    for path in sorted(source.glob("fluid-fluid-*/result.json")):
        partner = source / path.parent.name.replace("-fluid-ppc", "-frozen-ppc")
        fluid, frozen = read(path), read(partner / "result.json")
        with (
            np.load(path.with_suffix(".npz")) as a,
            np.load(partner / "result.npz") as b,
        ):
            step = fluid["parameters"]["steps"]
            for species in ["electrons", "N2plus", "O2plus", "Ominus"]:
                key = f"{step}_{species}"
                density_error.append(
                    float(np.max(np.abs(a[key] - b[key])))
                    / max(1.0, float(np.max(np.abs(a[key]))))
                )
            histogram_error.append(
                float(np.max(np.abs(a["electron_histogram"] - b["electron_histogram"])))
            )
        state, frozen_state = fluid["history"][-1], frozen["history"][-1]
        budgets = state["primary_source_budgets"]
        expected = sum(row[1] for row in budgets.values())
        electron_budget_error.append(abs(state["electron_energy_J"] / expected - 1))
        for gas in ["N2", "O2"]:
            discarded = budgets["p_" + gas][3]
            recoil = frozen_state["ion_energy_J"][gas + "plus"]
            recoil_error.append(abs(discarded / recoil - 1))
    report["identical_primary_events"] = dict(
        maximum_density_difference_over_peak=max(density_error),
        maximum_histogram_weight_difference=max(histogram_error),
        maximum_electron_energy_budget_relative_error=max(electron_budget_error),
        maximum_discarded_recoil_relative_error=max(recoil_error),
    )

    report["storage"] = []
    for steps in [20, 100]:
        for ions in ["fluid", "frozen"]:
            cases = [
                read(path)
                for path in sorted(
                    (args.directory / "storage" / f"storage-{steps}").glob(
                        f"fluid-{ions}-*/result.json"
                    )
                )
            ]
            if not cases:
                raise ValueError(f"Missing storage study: {steps} steps, {ions} ions")
            report["storage"].append(
                dict(
                    steps=steps,
                    ions=ions,
                    electron_macroparticles=interval(
                        [
                            case["history"][-1]["macroparticles"]["electrons"]
                            for case in cases
                        ]
                    ),
                    ion_density_bytes=interval(
                        [
                            sum(
                                value
                                for key, value in case["fluid_density_bytes"].items()
                                if key != "beam"
                            )
                            for case in cases
                        ]
                    ),
                    ion_particle_payload_bytes=interval(
                        [
                            sum(
                                value
                                for key, value in case["particle_payload_bytes"].items()
                                if key not in ["beam", "electrons"]
                            )
                            for case in cases
                        ]
                    ),
                    checkpoint_bytes=interval(
                        [case["checkpoint_bytes"] for case in cases]
                    ),
                    checkpoint_step_s=interval(
                        [case["checkpoint_step_s"] for case in cases]
                    ),
                )
            )
    report["source_profile_including_startup"] = profile_rows(
        args.directory / "profile" / "profile"
    )
    eager = args.directory / "profile-eager" / "profile"
    if eager.exists():
        report["source_profile_eager_module_loading"] = profile_rows(eager)
        differences = {"fluid": [], "quiet": []}
        for path in sorted(
            (args.directory / "profile" / "profile").glob("*/result.json")
        ):
            partner = eager / path.parent.name / "result.json"
            original, repeated = read(path), read(partner)
            if original["parameters"] != repeated["parameters"]:
                raise ValueError(f"Different profiling configurations: {path}")
            if (
                original["history"][-1]["macroparticles"]
                != repeated["history"][-1]["macroparticles"]
            ):
                raise ValueError(f"Different profiling populations: {path}")
            with (
                np.load(path.with_suffix(".npz")) as a,
                np.load(partner.with_suffix(".npz")) as b,
            ):
                np.testing.assert_array_equal(
                    a["electron_histogram"], b["electron_histogram"]
                )
                step = original["parameters"]["steps"]
                for species in ["electrons", "N2plus", "O2plus"]:
                    key = f"{step}_{species}"
                    differences[original["parameters"]["beam"]].append(
                        float(np.max(np.abs(a[key] - b[key])))
                        / max(1.0, float(np.max(np.abs(a[key]))))
                    )
        report["module_loading_comparison"] = dict(
            identical_product_counts_and_energy_histograms=True,
            maximum_density_difference_over_peak={
                beam: max(values) for beam, values in differences.items()
            },
            particle_positions="Kinetic projectile selection depends on unordered GPU "
            "dense-bin entries. Paired spatial snapshots require statistical comparisons.",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import ticker

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for beam, label, color in [
        ("fluid", "Fluid beam", "C0"),
        ("quiet", "Quiet particle beam", "C1"),
        ("random", "Random particle beam", "C2"),
    ]:
        for axis, phase, metric in [
            (axes[0], "fields", "beam_relative_rms_error"),
            (axes[1], "coupled", "electrons_relative_seed_noise"),
        ]:
            rows = sorted(
                [
                    row
                    for row in report["ensembles"][phase]["groups"]
                    if row["beam"] == beam and row["ions"] == "fluid"
                ],
                key=lambda row: row["ppc"],
            )
            x = [1e3 * row["timestep_s"]["mean"] for row in rows]
            y = [
                100 * (row[metric]["mean"] if phase == "fields" else row[metric])
                for row in rows
            ]
            axis.errorbar(
                x,
                y,
                xerr=[1e3 * row["timestep_s"]["ci99_half_width"] for row in rows],
                fmt="o-",
                color=color,
                label=label,
            )
            if beam != "fluid":
                for xx, yy, row in zip(x, y, rows):
                    axis.annotate(
                        str(row["ppc"]),
                        (xx, yy),
                        xytext=(3, 5),
                        textcoords="offset points",
                        fontsize=8,
                    )
            axis.set(xlabel="Timestep wall time (ms)", xscale="log")
            axis.grid(True, alpha=0.25)
    axes[0].set(
        title="Beam projection: fields only",
        ylabel="Beam-density RMS error relative to fluid (%)",
    )
    axes[0].set_yscale("symlog", linthresh=0.001, linscale=0.3)
    axes[0].set_ylim(bottom=-0.0002)
    axes[0].legend(fontsize=8, loc="lower right")
    axes[1].set(
        title="Coupled PJG and electron collisions",
        ylabel="Plasma-density seed noise (%)",
    )
    for axis, phase in zip(axes, ["fields", "coupled"]):
        axis.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1, 2, 5)))
        axis.xaxis.set_major_formatter(ticker.ScalarFormatter())
        axis.xaxis.set_minor_formatter(ticker.NullFormatter())
        rows = [
            row
            for row in report["ensembles"][phase]["groups"]
            if row["ions"] == "fluid"
        ]
        axis.set_xlim(
            0.8e3 * min(row["timestep_s"]["mean"] for row in rows),
            1.2e3 * max(row["timestep_s"]["mean"] for row in rows),
        )
    nr, nz = report["ensembles"]["fields"]["settings"]["cells"]
    samples = report["ensembles"]["fields"]["groups"][0]["timestep_s"]["samples"]
    fig.suptitle(
        f"{args.hardware}, {nr} x {nz} RZ, {samples} seeds, immobile ion fluids\n"
        "Labels: beam particles/cell; horizontal bars: 99% timing intervals",
        fontsize=11,
    )
    fig.savefig(args.output.with_name(args.output.name + "-noise-cost.png"), dpi=180)
    fig.savefig(args.output.with_name(args.output.name + "-noise-cost.pdf"))
    plt.close(fig)
    print(f"Exported {args.output}")


if __name__ == "__main__":
    main()
