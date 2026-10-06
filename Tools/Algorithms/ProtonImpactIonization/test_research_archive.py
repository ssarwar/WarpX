# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check preservation and interpretation of the September 2026 calibration."""

import hashlib
import io
import json
import tarfile
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

import numpy as np
from calibrated_pjg import PARAMETERS
from verify_research_archive import (
    ARCHIVE,
    check_record,
    finite_float,
    load_json,
    reject_constant,
    safe_path,
    verify_archive,
    verify_members,
)


class ResearchArchiveTests(unittest.TestCase):
    def test_manifest_and_members(self):
        result = verify_archive()
        self.assertGreater(result["artifacts"], 100)
        self.assertGreater(result["tar_members"], 100)

    def test_selected_parameters_are_frozen(self):
        for filename in ("fit-results.json", "frozen-validation.json"):
            result = load_json(ARCHIVE / "final" / filename)
            for target, selected in (("N2", "six"), ("O2", "five_fixed_edge")):
                records = result["targets"][target]
                self.assertEqual(records["selected"], selected)
                self.assertEqual(
                    records[selected]["parameters"], asdict(PARAMETERS[target])
                )
        # A repeated optimizer run is evidence, not a replacement coefficient set.
        repeated = load_json(ARCHIVE / "final/production-refit.json")
        for target in ("N2", "O2"):
            records = repeated["targets"][target]
            self.assertEqual(
                records["selected"], {"N2": "six", "O2": "five_fixed_edge"}[target]
            )
            actual = records[records["selected"]]["parameters"]
            for name, expected in asdict(PARAMETERS[target]).items():
                self.assertLess(abs(actual[name] / expected - 1), 1e-4)

    def test_no_generated_metadata(self):
        manifest = load_json(ARCHIVE / "manifest.json")
        for record in manifest["artifacts"]:
            for item in [record, *record.get("members", [])]:
                path = Path(item["path"])
                self.assertNotIn("__pycache__", path.parts)
                self.assertNotIn(".DS_Store", path.parts)
                self.assertFalse(any(part.startswith("._") for part in path.parts))
                self.assertNotEqual(path.suffix, ".pyc")

    def test_source_hash_chain(self):
        validation = load_json(ARCHIVE / "validation/archive-validation.json")
        self.assertEqual(
            validation["checks"]["frozen_recomputation"]["result_sha256"],
            hashlib.sha256(
                (ARCHIVE / "final/frozen-validation.json").read_bytes()
            ).hexdigest(),
        )
        sources = {
            "nifs": "inputs/nifs-optical.json",
            "o2_leiden": "inputs/leiden-o2.txt",
            "n2_leiden": "inputs/leiden-n2-0.1nm.txt",
            "mahla": "inputs/mahla-o2-o3-supplement.zip",
            "initial_fits": "history/results/source-refresh/source-refits.json",
            "fits": "history/results/source-refresh/source-refits.json",
        }
        for filename in (
            "final/fit-results.json",
            "final/frozen-validation.json",
            "final/production-refit.json",
            "history/results/source-refresh/source-refits.json",
            "history/results/source-refresh/mean-sensitivity.json",
        ):
            record = load_json(ARCHIVE / filename)
            for name, expected in record["source_sha256"].items():
                actual = hashlib.sha256(
                    (ARCHIVE / sources[name]).read_bytes()
                ).hexdigest()
                self.assertEqual(actual, expected, (filename, name))
        properties = load_json(
            ARCHIVE / "history/results/source-refresh/properties.json"
        )
        self.assertEqual(
            properties["source_fit_sha256"],
            hashlib.sha256((ARCHIVE / sources["fits"]).read_bytes()).hexdigest(),
        )

    def test_json_and_reference_dimensions(self):
        for path in ARCHIVE.rglob("*.json"):
            load_json(path)
        nifs = load_json(ARCHIVE / "inputs/nifs-optical.json")
        pstar = load_json(ARCHIVE / "inputs/pstar_reference.json")
        for target, rows in (("N2", 239), ("O2", 729)):
            self.assertEqual(np.asarray(nifs[target]["continuum"]).shape, (rows, 5))
            self.assertEqual(np.asarray(pstar["rows"][target]).shape, (132, 7))
        with self.assertRaises(ValueError):
            json.loads('{"bad": NaN}', parse_constant=reject_constant)
        with self.assertRaises(ValueError):
            json.loads('{"bad": 1e999}', parse_float=finite_float)

    def test_legacy_table_is_numerical_data(self):
        with tarfile.open(
            ARCHIVE / "validation/legacy-numerical-data.tar.gz"
        ) as archive:
            with archive.extractfile("tmp/pjg_table_audit.bin") as stream:
                table = np.frombuffer(stream.read(), dtype="<f8").reshape(148, 2049, 10)
        self.assertTrue(np.all(np.isfinite(table)))
        np.testing.assert_array_equal(np.unique(table[:, 0, 0]), [0, 1])
        np.testing.assert_array_equal(table[:, 0, 2], 0)
        np.testing.assert_array_equal(table[:, -1, 2], 1)
        self.assertTrue(np.all(np.diff(table[:, :, 2], axis=1) > 0))

    def test_reject_unsafe_paths_and_changed_bytes(self):
        for name in (
            "",
            ".",
            "/absolute",
            "../parent",
            "x/../y",
            "./x",
            "x//y",
            "x\\y",
            "C:x",
        ):
            with self.assertRaises(ValueError, msg=name):
                safe_path(name)
        record = {
            "path": "data",
            "bytes": 4,
            "sha256": hashlib.sha256(b"data").hexdigest(),
        }
        with self.assertRaises(ValueError):
            check_record(io.BytesIO(b"edit"), 4, record)
        with self.assertRaises(ValueError):
            check_record(io.BytesIO(b"data"), 5, record)

    def test_reject_non_regular_and_duplicate_members(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.tar.gz"
            for kind in ("link", "duplicate", "missing"):
                with tarfile.open(path, "w:gz") as archive:
                    member = tarfile.TarInfo("entry")
                    if kind == "link":
                        member.type = tarfile.SYMTYPE
                        member.linkname = "outside"
                    if kind != "missing":
                        archive.addfile(member, io.BytesIO())
                    if kind == "duplicate":
                        archive.addfile(member, io.BytesIO())
                record = {
                    "path": "entry",
                    "bytes": 0,
                    "sha256": hashlib.sha256(b"").hexdigest(),
                }
                with self.assertRaises(ValueError, msg=kind):
                    verify_members(path, [record])


if __name__ == "__main__":
    unittest.main()
