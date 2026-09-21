#!/usr/bin/env python3
"""Compare fluid plotfile/openPMD records with saved native arrays and output times."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd().parent / ".mpl-cache"))

import numpy as np
import yt

parser = argparse.ArgumentParser()
parser.add_argument("--openpmd", action="store_true")
args = parser.parse_args()
qe = 1.602176634e-19


def cell_center(values):
    values = np.squeeze(values)
    for axis, cells in enumerate([16, 32]):
        if values.shape[axis] == cells + 1:
            values = (
                np.take(values, range(cells), axis=axis)
                + np.take(values, range(1, cells + 1), axis=axis)
            ) / 2
    return values


def expected(saved):
    budget = saved["N2_immobile_source_budget"]
    return {
        "fluid_density_beam": cell_center(saved["beam"]),
        "fluid_current_beamz": cell_center(saved["beam_current"]),
        "fluid_density_i_N2_immobile": cell_center(saved["i_N2_immobile"]),
        "N2_immobile_product_weight_remainder": cell_center(
            saved["N2_immobile_product_weight_remainder"]
        ),
        "N2_immobile_emitted_number": cell_center(budget[..., 0]),
        "N2_immobile_electron_energy": cell_center(budget[..., 1]),
        "N2_immobile_binding_energy": cell_center(budget[..., 2]),
    }


outputs = sorted(
    path for path in Path("diags").glob("fields[0-9]?????") if path.is_dir()
)
steps = [int(path.name[-6:]) for path in outputs]
restart = "restart" in Path.cwd().name
assert steps == ([4, 6] if restart else [0, 2, 4, 6]), steps
for path in outputs:
    step = int(path.name[-6:])
    ds = yt.load(str(path))
    np.testing.assert_allclose(float(ds.current_time), step * 1e-12, rtol=2e-15)
    assert ds.domain_dimensions.tolist() == [16, 32, 1]
    grid = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    if step == 0:
        continue
    with np.load(f"state_{step}.npz") as saved:
        for name, reference in expected(saved).items():
            np.testing.assert_allclose(
                np.squeeze(grid["boxlib", name].v), reference, rtol=3e-14, err_msg=name
            )
    ion_charge = grid["boxlib", "rho_i_N2_immobile"].v
    electron_charge = grid["boxlib", "rho_e_N2_immobile"].v
    np.testing.assert_allclose(
        ion_charge + electron_charge, 0, atol=3e-14 * np.max(np.abs(ion_charge))
    )
    # On Yee, native charge is unfiltered and can also be compared directly.
    if "PSATD" not in Path.cwd().name:
        np.testing.assert_allclose(
            ion_charge, qe * grid["boxlib", "fluid_density_i_N2_immobile"].v, rtol=3e-14
        )

for diagnostic in [
    "ParticleNumber",
    "ParticleCharge",
    "ParticleEnergy",
    "ParticleMomentum",
    "PrescribedSourceBudget",
]:
    rows = np.loadtxt(Path("diags/reducedfiles") / (diagnostic + ".txt"), ndmin=2)
    # Re-running this test may append to an existing restart diagnostic. Each
    # new continuation must still have the correct cadence and physical times.
    last = rows[-4:] if restart else rows
    np.testing.assert_array_equal(
        last[:, 0], np.arange(3, 7) if restart else np.arange(7)
    )
    np.testing.assert_allclose(last[:, 1], last[:, 0] * 1e-12, rtol=2e-15)

if args.openpmd:
    import openpmd_api as io

    paths = sorted(Path("diags/pmd").rglob("*.h5"))
    assert paths, "Missing openPMD output"
    found_steps = set()
    for path in paths:
        series = io.Series(str(path), io.Access.read_only)
        for step, iteration in series.iterations.items():
            found_steps.add(step)
            np.testing.assert_allclose(
                iteration.time * iteration.time_unit_SI, step * 1e-12, rtol=2e-15
            )
            if step == 0:
                continue
            with np.load(f"state_{step}.npz") as saved:
                references = expected(saved)
            for name, reference in references.items():
                record = iteration.meshes[name]
                component = record[io.Mesh_Record_Component.SCALAR]
                values = component.load_chunk()
                series.flush()
                np.testing.assert_allclose(
                    np.squeeze(values).T * component.unit_SI,
                    reference,
                    rtol=3e-14,
                    err_msg=name,
                )
                np.testing.assert_allclose(component.position, [0.5, 0.5])
                assert record.axis_labels == ["z", "r"]
                assert record.time_offset == 0
                dimensions = np.zeros(7)
                if name.startswith("fluid_density_"):
                    dimensions[0] = -3
                elif name.startswith("fluid_current_"):
                    dimensions[0], dimensions[3] = -2, 1
                elif name.endswith("_energy"):
                    dimensions[:3] = [2, 1, -2]
                np.testing.assert_array_equal(record.unit_dimension, dimensions)
        series.close()
    assert found_steps == set(steps), found_steps
print("PASS: fluid/source output values, mesh centering, units, timestamps and cadence")
