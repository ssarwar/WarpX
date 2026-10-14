# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Compare compact q-bin outcomes with the coupled source at 1--10 keV."""

import json

import numpy as np
from hybrid_reference import REST, Elastic
from reference_paths import OUTPUT

ROOT = OUTPUT
out = []
for target in ["N2", "O2"]:
    data = np.load(ROOT / f"binned-q-moments-{target}.npz")
    edges = data["edges"]
    mom = data["moments"]
    elastic = Elastic(target)
    R = {"N2": 2.068, "O2": 2.281}[target]
    ref = json.loads((ROOT / f"q-transition-{target}.json").read_text())
    for row in ref:
        E = row["E_eV"]
        zmax = R * np.sqrt(E * (E + 2 * REST)) / (7.2973525693e-3 * REST)
        mass = np.diff(elastic.cdf(E, np.minimum(edges / zmax, 1)))
        actual = mass @ mom
        expected = np.array(row["moments"]) / row["total_rate"]
        error = abs(actual / expected - 1)
        out.append(
            {
                "target": target,
                "E_eV": E,
                "relative_errors": error.tolist(),
                "max_relative_error": float(error.max()),
            }
        )
        print(target, E, error, flush=True)
(ROOT / "q-transition-check.json").write_text(json.dumps(out, indent=2) + "\n")
