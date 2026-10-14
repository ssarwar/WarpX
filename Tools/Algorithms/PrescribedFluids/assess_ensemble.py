#!/usr/bin/env python3
"""Apply the documented coupled-population bounds to independent seed groups.

Spectrum differences are reported with simultaneous seed-based confidence
bands. Particle counts are not treated as independent trials: quiet sampling,
variable weights and repeated collisions correlate records within a run.
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from analyze import interval
from scipy.stats import t as student_t


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    summary = json.loads((args.directory / "summary.json").read_text())
    assert summary["settings"]["mode"] == "coupled"
    groups = summary["groups"]
    baseline = next(g for g in groups if g["beam"] == g["ions"] == "fluid")
    observables = [
        "electrons_number",
        "N2plus_number",
        "O2plus_number",
        "Ominus_number",
        "O2minus_number",
        "electron_energy_J",
    ]
    observables = [key for key in observables if key in baseline]
    spectra = defaultdict(list)
    for path in sorted(args.directory.glob("*/result.json")):
        case = json.loads(path.read_text())
        p = case["parameters"]
        key = (p["beam"], p["ions"], p["ppc"])
        for state in case["history"]:
            n = state["physical"]
            charge = (
                n["N2plus"]
                + n["O2plus"]
                - n["Ominus"]
                - n.get("O2minus", 0)
                - n["electrons"]
            )
            assert abs(charge) <= 3e-12 * max(1, n["electrons"]), (path, state)
        with np.load(path.with_suffix(".npz")) as arrays:
            counts = arrays["electron_histogram"]
            assert counts.sum() > 0, path
            spectra[key].append(np.cumsum(counts) / counts.sum())
    reference = np.asarray(spectra["fluid", "fluid", 1])
    assert len(reference) >= 6
    # Nominal 99% simultaneous Student bands: Bonferroni covers both ensembles,
    # both tails, all reported bins and all representation comparisons.
    comparisons = max(1, len(spectra) - 1) * reference.shape[1]
    quantile = 1 - 0.01 / (4 * comparisons)

    def bands(samples):
        assert len(samples) >= 6
        return (
            student_t.ppf(quantile, len(samples) - 1)
            * np.std(samples, axis=0, ddof=1)
            / np.sqrt(len(samples))
        )

    results = []
    for group in groups:
        key = (group["beam"], group["ions"], group["ppc"])
        if key == ("fluid", "fluid", 1):
            continue
        checks = {}
        for name in observables:
            ref, actual = baseline[name], group[name]
            assert min(ref["samples"], actual["samples"]) >= 6
            delta = abs(actual["mean"] - ref["mean"])
            bound = ref["ci99_half_width"] + actual["ci99_half_width"]
            bound += 3e-12 * max(abs(ref["mean"]), abs(actual["mean"]))
            checks[name] = dict(difference=delta, bound=bound, passed=delta <= bound)
        samples = np.asarray(spectra[key])
        difference = np.abs(samples.mean(axis=0) - reference.mean(axis=0))
        band = bands(samples) + bands(reference) + 3e-12
        results.append(
            dict(
                beam=key[0],
                ions=key[1],
                ppc=key[2],
                population_checks=checks,
                maximum_cdf_difference=float(difference.max()),
                maximum_cdf_band_excess=float(np.maximum(difference - band, 0).max()),
                seed_maximum_cdf_difference=interval(
                    [np.max(np.abs(row - reference.mean(axis=0))) for row in samples]
                ),
            )
        )
    output = args.directory / "acceptance.json"
    output.write_text(json.dumps(results, indent=2) + "\n")
    assert all(
        check["passed"]
        for row in results
        for check in row["population_checks"].values()
    ), output
    assert all(row["maximum_cdf_band_excess"] == 0 for row in results), output
    print("PASS: coupled populations, energy, spectra and charge balance")


if __name__ == "__main__":
    main()
