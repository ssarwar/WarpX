# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Verify the frozen PJG research payloads without extracting or executing them."""

import argparse
import hashlib
import json
import math
import tarfile
from pathlib import Path, PurePosixPath

ARCHIVE = Path(__file__).resolve().parent / "Research"
GUIDES = {
    "README.md",
    "DATA_SOURCES.md",
    "REPRODUCING.md",
    "VALIDATION.md",
    "FIGURES.md",
    "manifest.json",
}


def reject_constant(value):
    """Reject nonstandard JSON NaN/Infinity instead of accepting them silently."""
    raise ValueError(f"Non-finite JSON constant: {value}")


def load_json(path):
    return json.loads(
        path.read_text(), parse_constant=reject_constant, parse_float=finite_float
    )


def finite_float(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite JSON number: {value}")
    return result


def safe_path(name):
    """Require a canonical relative POSIX path, also on Windows hosts."""
    path = PurePosixPath(name)
    if (
        not name
        or not path.parts
        or path.is_absolute()
        or path.as_posix() != name
        or any(part in ("..", ".") for part in path.parts)
        or any(character in name for character in ("\\", ":", "\0", "\n", "\r"))
    ):
        raise ValueError(f"Unsafe archive path: {name!r}")
    return path


def stream_digest(stream):
    checksum = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        checksum.update(block)
    return checksum.hexdigest()


def check_record(stream, size, record):
    if size != record["bytes"] or stream_digest(stream) != record["sha256"]:
        raise ValueError(f"Size or SHA-256 mismatch: {record['path']}")


def verify_members(path, records):
    """Check exact regular-file membership; never trust tar extraction paths."""
    expected = {}
    for record in records:
        name = safe_path(record["path"]).as_posix()
        if name in expected:
            raise ValueError(f"Duplicate manifest member: {name}")
        expected[name] = record
    seen = set()
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            name = safe_path(member.name).as_posix()
            if not member.isfile() or name in seen or name not in expected:
                raise ValueError(f"Unexpected or non-regular tar member: {name}")
            seen.add(name)
            with archive.extractfile(member) as stream:
                check_record(stream, member.size, expected[name])
    if seen != expected.keys():
        raise ValueError(f"Missing tar members: {sorted(expected.keys() - seen)}")
    return len(seen)


def verify_archive(root=ARCHIVE):
    """Check payload bytes, compressed members and accidental unlisted files.

    Guides and the manifest are ordinary Git-versioned text, not self-hashed.
    Production-source changes do not silently rewrite this historical record.
    """
    root = Path(root).resolve()
    manifest = load_json(root / "manifest.json")
    if manifest["schema_version"] != 1:
        raise ValueError("Unsupported research manifest version")
    seen = set()
    total_bytes = 0
    member_count = 0
    for record in manifest["artifacts"]:
        name = safe_path(record["path"]).as_posix()
        if name in seen or name in GUIDES:
            raise ValueError(f"Duplicate or reserved payload: {name}")
        seen.add(name)
        path = root / name
        if (
            not path.is_file()
            or path.is_symlink()
            or not path.resolve().is_relative_to(root)
        ):
            raise ValueError(f"Missing or unsafe payload: {name}")
        with path.open("rb") as stream:
            check_record(stream, path.stat().st_size, record)
        total_bytes += record["bytes"]
        if "members" in record:
            member_count += verify_members(path, record["members"])
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink in research archive: {path}")
        if path.is_file():
            name = path.relative_to(root).as_posix()
            if name not in seen | GUIDES:
                raise ValueError(f"Unlisted research payload: {name}")
    return {"artifacts": len(seen), "bytes": total_bytes, "tar_members": member_count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    args = parser.parse_args()
    print(json.dumps(verify_archive(args.archive), indent=2))


if __name__ == "__main__":
    main()
