#!/usr/bin/env python3
"""Check that interrupted job records remain explicit in the evidence archive."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import export_reaudit


class TestInterruptedEvidence(unittest.TestCase):
    def test_archive_keeps_valid_results_and_records_incomplete_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            study = root / "validation-42"
            nested = study / "rank-0"
            nested.mkdir(parents=True)
            (study / "topology.json").write_text("")
            (study / "comparison.json").write_text('{"passed":')
            (nested / "incomplete.xml").write_text("<testsuite>")
            (nested / "valid.xml").write_text(
                '<testsuite><testcase name="good" time="1"/>'
                '<testcase name="bad"><failure/><system-out>failed assertion</system-out>'
                "</testcase></testsuite>"
            )
            (nested / "summary.json").write_text('{"passed": false}')
            output = root / "archive.json"
            with (
                patch(
                    "sys.argv",
                    ["export", str(root), "--jobs", "42", "--output", str(output)],
                ),
                patch("subprocess.check_output", return_value="42|TIMEOUT|0:0\n"),
            ):
                export_reaudit.main()
            record = json.loads(output.read_text())["42"]["studies"][0]
            self.assertNotIn("topology", record["provenance"])
            self.assertEqual(
                {
                    case["name"]: case["passed"]
                    for case in record["tests"]["rank-0/valid.xml"]
                },
                {"good": True, "bad": False},
            )
            self.assertEqual(
                record["comparisons"], {"rank-0/summary.json": {"passed": False}}
            )
            self.assertEqual(
                record["tests"]["rank-0/valid.xml"][1]["failure_output"],
                "failed assertion",
            )
            self.assertEqual(
                set(record["unreadable_records"]),
                {"topology.json", "comparison.json", "rank-0/incomplete.xml"},
            )
            self.assertEqual(
                record["unreadable_records"]["topology.json"]["sha256"],
                hashlib.sha256(b"").hexdigest(),
            )


if __name__ == "__main__":
    unittest.main()
