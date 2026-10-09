# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Readable probability bundles; aliases are an initialization-time representation.

Low-rank factors are a numerical encoding, not molecular transfer ranks. Every
cell retains its exact discrete support. Reconstruction is checked against the
uncompressed probabilities, including rare energy-transfer tails.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from bundle import TYPES, read_bundle

OMITTED = {"aliases", "angular_lookup", "conditional_lookup"}
BLOCK_CELLS = 256
L1_TOLERANCE = 2e-7
MOMENT_TOLERANCE = 2e-6


def spectra(arrays, first, last):
    offsets, entries = arrays["cell_offsets"], arrays["aliases"]
    ids = np.unique(entries["outcome"][offsets[first] : offsets[last]])
    result = np.zeros((len(ids), last - first))
    for j, cell in enumerate(range(first, last)):
        values = entries[offsets[cell] : offsets[cell + 1]]
        p = values["cut"].astype(float) / len(values)
        np.add.at(p, values["alias"], (1 - values["cut"].astype(float)) / len(values))
        result[np.searchsorted(ids, values["outcome"]), j] = p
    return ids, result


def factors(probability, labels):
    """Choose rank against decoded probabilities, not singular values alone."""
    support = probability > 0
    features = np.vstack(
        [
            np.ones(len(labels)),
            labels == 0,
            labels > 0,
            labels < 0,
            np.maximum(labels, 0),
            np.maximum(-labels, 0),
            labels**2,
            labels**4,
        ]
    )
    expected = features @ probability
    # Weight rare, energetic outcomes as well as the bulk distribution.
    positive = expected.max(axis=1) > 0
    scale = np.max(
        features[positive]
        / np.maximum(expected[positive].max(axis=1, keepdims=True), 1e-100),
        axis=0,
    )
    scale = np.maximum(scale, 1)
    u, singular, v = np.linalg.svd(scale[:, None] * probability, full_matrices=False)
    for rank in range(1, len(singular) + 1):
        left = u[:, :rank] * singular[:rank] / scale[:, None]
        right = v[:rank]
        # Decimal serialization must be included in the error check.
        left = np.array([float(f"{x:.10g}") for x in left.ravel()]).reshape(left.shape)
        right = np.array([float(f"{x:.10g}") for x in right.ravel()]).reshape(
            right.shape
        )
        fitted = np.where(support, np.maximum(left @ right, 0), 0)
        total = fitted.sum(axis=0)
        if np.any(total <= 0):
            continue
        fitted /= total
        l1 = float(abs(fitted - probability).sum(axis=0).max())
        moments = float(
            np.max(abs(features @ fitted - expected) / np.maximum(expected, 1e-100))
        )
        if l1 <= L1_TOLERANCE and moments <= MOMENT_TOLERANCE:
            return left, right, support, l1, moments
    raise ValueError("Probability factorization did not meet its accuracy budget")


def write_array(path, values, kind, offsets=None):
    """Exact decimal values, with explicit reuse of identical physical rows."""
    values = np.asarray(values).ravel()
    if offsets is None:
        offsets = [0, len(values)]
    rows = {}
    with path.open("w") as stream:
        stream.write(f"# {path.stem}: {kind}, {len(values)} values. See README.md.\n")
        for first, last in zip(offsets[:-1], offsets[1:]):
            row = values[first:last]
            if not len(row):
                continue
            key = row.tobytes()
            if key in rows:
                stream.write(f"copy {len(row)} {rows[key]}\n")
                continue
            rows[key] = first
            if np.all(row == row[0]) and len(row) > 3:
                stream.write(f"repeat {len(row)} {row[0]:.17g}\n")
            else:
                stream.write(f"values {len(row)}\n")
                precision = ".9g" if kind == "f32" else ".17g"
                for k in range(0, len(row), 8):
                    stream.write(
                        " ".join(format(x, precision) for x in row[k : k + 8]) + "\n"
                    )


def numeric_tokens(path):
    for line in Path(path).read_text().splitlines():
        yield from line.split("#", 1)[0].split()


def read_array(path, kind, count):
    tokens = iter(numeric_tokens(path))
    values = np.empty(count, TYPES[kind])
    position = 0
    for command in tokens:
        n = int(next(tokens))
        if n <= 0 or position + n > count:
            raise ValueError("Invalid text array length")
        if command == "values":
            values[position : position + n] = [next(tokens) for _ in range(n)]
        elif command == "repeat":
            values[position : position + n] = next(tokens)
        elif command == "copy":
            start = int(next(tokens))
            if start < 0 or start + n > position:
                raise ValueError("Invalid text array copy")
            values[position : position + n] = values[start : start + n]
        else:
            raise ValueError("Unknown text array command")
        position += n
    if position != count:
        raise ValueError("Truncated text array")
    return values


def write_text_bundle(source, directory):
    metadata, arrays = read_bundle(source)
    if metadata[0].split()[-1] != "alias":
        raise ValueError("Convert a validated alias bundle")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    descriptions = []
    for name, data in arrays.items():
        if name in OMITTED:
            continue
        kind = next(k for k, dtype in TYPES.items() if dtype == data.dtype)
        grid = (
            arrays.get("angular_offsets")
            if name in {"angular_u", "deflection", "changing"}
            else None
        )
        if name == "conditional_u":
            grid = arrays["conditional_offsets"]
        write_array(directory / f"{name}.txt", data, kind, grid)
        descriptions.append(f"{name} {kind} {len(data)} {name}.txt")
    count = len(arrays["cell_offsets"]) - 1
    maximum_l1 = maximum_moments = 0
    coefficients = 0
    with (directory / "probabilities.txt").open("w") as stream:
        stream.write(
            "# Probability encoding; columns are angular/momentum-transfer cells.\n"
        )
        stream.write(
            "# Signed factors are not cross sections. See README.md for reconstruction.\n"
        )
        stream.write(
            f"WARPX_ROTATIONAL_PROBABILITIES {count} {len(arrays['aliases'])}\n"
        )
        for first in range(0, count, BLOCK_CELLS):
            last = min(count, first + BLOCK_CELLS)
            ids, probability = spectra(arrays, first, last)
            left, right, support, l1, moments = factors(
                probability, arrays["outcomes"][2 * ids]
            )
            maximum_l1 = max(maximum_l1, l1)
            maximum_moments = max(maximum_moments, moments)
            coefficients += left.size + right.size
            stream.write(f"block {first} {last - first} {len(ids)} {left.shape[1]}\n")
            stream.write(
                "# outcome_id, support interval count, [first,end) cell intervals, left factors\n"
            )
            for outcome, mask, row in zip(ids, support, left):
                edges = np.flatnonzero(np.diff(np.r_[False, mask, False])).reshape(
                    -1, 2
                )
                stream.write(
                    f"{outcome} {len(edges)} " + " ".join(map(str, edges.ravel())) + " "
                )
                stream.write(" ".join(f"{x:.10g}" for x in row) + "\n")
            stream.write("# right factors, one row per cell\n")
            np.savetxt(stream, right.T, fmt="%.10g")
            if first % (BLOCK_CELLS * 32) == 0:
                print(first, count, maximum_l1, maximum_moments, flush=True)
    header = [
        "WARPX_RECIPROCAL_ROTATION",
        metadata[0].rsplit(" ", 1)[0] + " probabilities",
        metadata[1],
        str(len(descriptions)),
        *descriptions,
        f"probabilities alias {len(arrays['aliases'])} probabilities.txt",
    ]
    (directory / "thermal_rotation.rot").write_text("\n".join(header) + "\n")
    return {
        "bytes": sum(p.stat().st_size for p in directory.glob("*.txt"))
        + (directory / "thermal_rotation.rot").stat().st_size,
        "coefficients": coefficients,
        "maximum_l1_error": maximum_l1,
        "maximum_moment_error": maximum_moments,
    }


def read_probabilities(path, offsets):
    tokens = iter(numeric_tokens(path))
    if next(tokens) not in (
        "WARPX_ROTATIONAL_PROBABILITIES",
        "WARPX_PROBABILITY_FACTORS_V1",
    ):
        raise ValueError("Invalid probability factors")
    cells, entries = int(next(tokens)), int(next(tokens))
    if cells + 1 != len(offsets) or entries != offsets[-1]:
        raise ValueError("Probability dimensions mismatch")
    cdf = np.empty(entries, dtype=np.float64)
    outcome_ids = np.empty(entries, dtype=np.uint16)
    position = 0
    while position < cells:
        if next(tokens) != "block":
            raise ValueError("Missing probability block")
        first, count, rows, rank = [int(next(tokens)) for _ in range(4)]
        if first != position:
            raise ValueError("Probability blocks are not consecutive")
        support = np.zeros((rows, count), bool)
        ids = np.empty(rows, dtype=np.uint16)
        left = np.empty((rows, rank))
        for i in range(rows):
            ids[i] = int(next(tokens))
            intervals = int(next(tokens))
            for _ in range(intervals):
                lo, hi = int(next(tokens)), int(next(tokens))
                support[i, lo:hi] = True
            left[i] = [float(next(tokens)) for _ in range(rank)]
        right = np.array([float(next(tokens)) for _ in range(count * rank)]).reshape(
            count, rank
        )
        fitted = np.where(support, np.maximum(left @ right.T, 0), 0)
        fitted /= fitted.sum(axis=0)
        for j in range(count):
            keep = support[:, j]
            probability = fitted[keep, j]
            values = cdf[offsets[first + j] : offsets[first + j + 1]]
            if len(values) != keep.sum():
                raise ValueError("Discrete support mismatch")
            # Independent cumulative reference; the C++ reader prepares aliases.
            prefix = np.cumsum(probability, dtype=np.longdouble)
            values[:] = prefix / prefix[-1]
            outcome_ids[offsets[first + j] : offsets[first + j + 1]] = ids[keep]
        position += count
    if next(tokens, None) is not None:
        raise ValueError("Trailing probability data")
    return cdf, outcome_ids


def read_text_bundle(index):
    from lookup_index import lookup_arrays

    index = Path(index)
    lines = index.read_text().splitlines()
    arrays = {}
    for line in lines[4:-1]:
        name, kind, count, file = line.split()
        arrays[name] = read_array(index.parent / file, kind, int(count))
    arrays["cdf"], arrays["outcome_ids"] = read_probabilities(
        index.parent / "probabilities.txt", arrays["cell_offsets"]
    )
    arrays.update(
        {
            name: np.asarray(values, dtype=TYPES[kind])
            for name, kind, values in lookup_arrays(arrays)
        }
    )
    metadata = [lines[1].replace(" probabilities", " cumulative"), lines[2]]
    return metadata, arrays


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result = write_text_bundle(args.bundle, args.output)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(result, flush=True)
