# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""V6 scientific arrays, split into ordinary Git-sized binary files."""

from pathlib import Path

import numpy as np

ALIAS = np.dtype([("cut", "<f4"), ("alias", "<u2"), ("outcome", "<u2")])
TYPES = {
    "f64": np.dtype("<f8"),
    "f32": np.dtype("<f4"),
    "u32": np.dtype("<u4"),
    "u16": np.dtype("<u2"),
    "alias": ALIAS,
}
PART_BYTES = 32 * 1024 * 1024


def write_bundle(directory, metadata, arrays):
    """Write only metadata and scientific arrays, without validation records."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    descriptions = []
    parts = []
    offset = 0
    stream = None
    for name, kind, values in arrays:
        data = np.ascontiguousarray(values, dtype=TYPES[kind]).ravel()
        padding = (-offset) % 8
        payloads = [memoryview(bytes(padding)), memoryview(data).cast("B")]
        descriptions.append(f"{name} {kind} {data.size} {offset + padding}")
        for payload in payloads:
            while len(payload):
                if offset % PART_BYTES == 0:
                    if stream is not None:
                        stream.close()
                    filename = f"tables-{len(parts):03d}.bin"
                    stream = (directory / filename).open("wb")
                    parts.append(filename)
                take = min(len(payload), PART_BYTES - offset % PART_BYTES)
                stream.write(payload[:take])
                payload = payload[take:]
                offset += take
    if stream is not None:
        stream.close()
    for old in directory.glob("tables-*.bin"):
        if old.name not in parts:
            old.unlink()
    index = [
        "WARPX_THERMAL_ROTATION_V6",
        *metadata,
        f"{len(descriptions)} {len(parts)} {offset}",
        *descriptions,
    ]
    index.extend(f"{p} {(directory / p).stat().st_size}" for p in parts)
    (directory / "thermal_rotation.rot").write_text("\n".join(index) + "\n")
    return offset


def read_bundle(index):
    """Read packed data for offline verification; no source-model imports."""
    index = Path(index)
    lines = index.read_text().splitlines()
    if lines[0] != "WARPX_THERMAL_ROTATION_V6":
        raise ValueError("Not a reciprocal V6 bundle")
    count, parts, size = map(int, lines[3].split())
    if len(lines) != 4 + count + parts:
        raise ValueError("Invalid V6 index length")
    payload = bytearray()
    for line in lines[4 + count :]:
        filename, expected = line.split()
        path = index.parent / filename
        if path.parent != index.parent or path.stat().st_size != int(expected):
            raise ValueError("Invalid V6 part")
        payload.extend(path.read_bytes())
    if len(payload) != size:
        raise ValueError("Invalid V6 payload size")
    arrays = {}
    end = 0
    for line in lines[4 : 4 + count]:
        name, kind, n, offset = line.split()
        n, offset = int(n), int(offset)
        dtype = TYPES[kind]
        if name in arrays or offset < end or offset % 8:
            raise ValueError("Invalid V6 array layout")
        end = offset + n * dtype.itemsize
        if end > size:
            raise ValueError("Truncated V6 array")
        arrays[name] = np.frombuffer(payload, dtype=dtype, count=n, offset=offset)
    if end != size:
        raise ValueError("Trailing V6 bytes")
    return lines[1:3], arrays
