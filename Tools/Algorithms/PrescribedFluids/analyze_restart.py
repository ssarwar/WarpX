#!/usr/bin/env python3
"""Check deterministic source continuation and paired stochastic MCC ensembles."""

import argparse
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=Path)
args = parser.parse_args()
pairs = sorted(args.directory.glob("*-seed*"))
assert len(pairs) >= 6, "At least six paired independent MCC continuations are required"
observables = {}
for pair in pairs:
    records = [
        json.loads((pair / run / "result.json").read_text())
        for run in ["full", "restart"]
    ]
    final = [record["history"][-1] for record in records]
    restored = records[1]["history"][0]
    saved = records[0]["history"][1]
    assert saved["step"] == restored["step"] == 10
    assert saved["time"] == restored["time"]
    assert saved["macroparticles"] == restored["macroparticles"]
    # Complete particle/native-field comparison precedes stepping in benchmark.py.
    assert "restored before stepping" in (pair / "restart/run.log").read_text()
    for collision in ["p_N2", "p_O2"]:
        np.testing.assert_allclose(
            final[1]["primary_source_budgets"][collision],
            final[0]["primary_source_budgets"][collision],
            rtol=3e-13,
        )
        np.testing.assert_allclose(
            final[1]["pending"][collision], final[0]["pending"][collision], rtol=3e-13
        )
    for state in final:
        number = state["physical"]
        plasma_charge = (
            number["N2plus"]
            + number["O2plus"]
            - number["Ominus"]
            - number["O2minus"]
            - number["electrons"]
        )
        assert abs(plasma_charge) < 3e-12 * number["electrons"]
    for name in ["electrons", "N2plus", "O2plus", "Ominus", "O2minus"]:
        observables.setdefault(name, []).append(
            [state["physical"][name] for state in final]
        )
    observables.setdefault("electron_energy_J", []).append(
        [state["electron_energy_J"] for state in final]
    )
summary = {}
for name, values in observables.items():
    values = np.asarray(values)
    difference = values[:, 1] - values[:, 0]
    mean = float(difference.mean())
    standard_error = float(difference.std(ddof=1) / np.sqrt(len(difference)))
    # Five paired standard errors is set before measuring this ensemble.
    # The numerical floor is only roundoff; MCC is not bitwise reproducible
    # after changing the number of ranks or the device RNG decomposition.
    bound = 5 * standard_error + 3e-12 * float(np.abs(values).max())
    summary[name] = dict(
        mean_difference=mean,
        standard_error=standard_error,
        bound=bound,
        passed=abs(mean) <= bound,
    )
(args.directory / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
assert all(item["passed"] for item in summary.values()), summary
print("PASS: deterministic source budgets and stochastic MCC restart continuation")
