#!/usr/bin/env python3
"""Quantify every restart discrepancy without changing the regression tolerance."""

import argparse
import json
from pathlib import Path

import numpy as np
import yt


def compare(reference, actual, tolerance=1e-12):
    datasets = [yt.load(str(path)) for path in (reference, actual)]
    grids = []
    for ds in datasets:
        ds.force_periodicity()
        grids.append(ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions))
    assert set(datasets[0].field_list) == set(datasets[1].field_list)
    orders = {}
    identities = {}
    for species, attribute in datasets[0].field_list:
        if attribute != "particle_id":
            continue
        pairs = []
        orders[species] = []
        for grid in grids:
            ids = grid[species, "particle_id"].v.reshape(-1)
            cpus = grid[species, "particle_cpu"].v.reshape(-1)
            order = np.lexsort((ids, cpus))
            ordered = np.column_stack((cpus[order], ids[order]))
            assert not np.any(np.all(ordered[1:] == ordered[:-1], axis=1)), species
            pairs.append(ordered)
            orders[species].append(order)
        identities[species] = bool(np.array_equal(*pairs))
        assert identities[species], f"Particle identities differ: {species}"
    rows = {}
    for field in datasets[0].field_list:
        expected = grids[0][field]
        value = grids[1][field].to(expected.units).v
        units = str(expected.units)
        expected = expected.v
        if field[0] in orders:
            expected = expected.reshape(-1)[orders[field[0]][0]]
            value = value.reshape(-1)[orders[field[0]][1]]
        assert value.shape == expected.shape, field
        assert np.all(np.isfinite(expected)) and np.all(np.isfinite(value)), field
        delta = np.abs(value - expected)
        peak = float(np.max(np.abs(expected), initial=0))
        difference = float(np.max(delta, initial=0))
        error = difference / peak if peak else difference
        rows["/".join(field)] = dict(
            units=units,
            shape=list(expected.shape),
            changed_elements=int(np.count_nonzero(delta)),
            maximum_absolute_difference=difference,
            reference_peak=peak,
            component_normalized_error=error,
            passes_original_criterion=bool(error < tolerance),
        )
    # Component normalization can amplify error in a small transverse component.
    # Report vector scales as additional evidence, never as replacement criteria.
    families = [
        ("boxlib", [prefix + axis for axis in "xyz"]) for prefix in ("B", "E", "j")
    ]
    families += [
        (species, ["particle_momentum_" + axis for axis in "xyz"]) for species in orders
    ]
    for species, attributes in families:
        keys = [species + "/" + attribute for attribute in attributes]
        if not all(key in rows for key in keys):
            continue
        scale = max(rows[key]["reference_peak"] for key in keys)
        for key in keys:
            rows[key]["vector_peak"] = scale
            rows[key]["difference_over_vector_peak"] = (
                rows[key]["maximum_absolute_difference"] / scale if scale else None
            )
    metadata = [
        dict(
            time_s=float(ds.current_time.to("s")),
            dimensions=ds.domain_dimensions.tolist(),
            lower_m=ds.domain_left_edge.to("m").v.tolist(),
            upper_m=ds.domain_right_edge.to("m").v.tolist(),
        )
        for ds in datasets
    ]
    assert metadata[0] == metadata[1], "Plotfile times, meshes or domains differ"
    report = dict(
        reference=str(reference),
        actual=str(actual),
        original_tolerance=tolerance,
        identity_matches=identities,
        metadata=metadata,
        fields=rows,
    )
    report["failed_fields"] = [
        key for key, row in rows.items() if not row["passes_original_criterion"]
    ]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("actual", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = compare(args.reference.resolve(), args.actual.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {key: report["fields"][key] for key in report["failed_fields"]}, indent=2
        )
    )


if __name__ == "__main__":
    main()
