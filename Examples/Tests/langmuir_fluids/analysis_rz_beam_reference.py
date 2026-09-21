#!/usr/bin/env python3
"""Check convergence of native ballistic particle fields to the rigid beam."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(
    0, str(Path(__file__).resolve().parents[3] / "Tools/Algorithms/PrescribedFluids")
)
from analyze import weighted_norm

parser = argparse.ArgumentParser()
parser.add_argument("fluid", type=Path)
parser.add_argument("coarse", type=Path)
parser.add_argument("fine", type=Path)
args = parser.parse_args()
paths = [args.fluid, args.coarse, args.fine]
results = [json.loads(path.read_text()) for path in paths]
arrays = [dict(np.load(path.with_suffix(".npz"))) for path in paths]
cells = results[0]["parameters"]["cells"]
step = results[0]["parameters"]["steps"]
for sample in [0, step]:
    for field in ["beam", "Er", "Ez", "Btheta"]:
        name = f"{sample}_{field}"
        scale = weighted_norm(arrays[0][name], cells)
        coarse, fine = [
            weighted_norm(array[name] - arrays[0][name], cells) / scale
            for array in arrays[1:]
        ]
        # Composite midpoint particle quadrature is second order in sampling
        # spacing. Increasing 4 -> 256 ppc refines each coordinate eightfold.
        # Allow solver/time-discretization errors to limit the final decrease.
        assert fine < 1e-3, (name, coarse, fine)
        assert fine < 0.15 * coarse, (name, coarse, fine)
        print(name, "coarse error", coarse, "fine error", fine)
for result in results:
    counts = [row["physical"]["beam"] for row in result["history"]]
    np.testing.assert_allclose(counts, counts[0], rtol=3e-13)
print("PASS: ballistic particle beam normalization, propagation and field convergence")
