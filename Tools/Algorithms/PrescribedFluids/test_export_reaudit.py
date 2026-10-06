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
    def test_archive_retains_restart_pre_step_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            study = root / "restart-expanded-42"
            for run in ["full", "restart"]:
                directory = study / "seed41" / run
                directory.mkdir(parents=True)
                (directory / "result.json").write_text('{"seed": 41}')
                (directory / "run.log").write_text(
                    "PASS: all saved fields, particles and source state restored before stepping"
                    if run == "restart"
                    else "Uninterrupted run"
                )
            output = root / "archive.json"
            with (
                patch(
                    "sys.argv",
                    ["export", str(root), "--jobs", "42", "--output", str(output)],
                ),
                patch("subprocess.check_output", return_value="42|COMPLETED|0:0\n"),
            ):
                export_reaudit.main()
            runs = json.loads(output.read_text())["42"]["studies"][0]["restart_runs"]
            self.assertEqual(len(runs), 2)
            self.assertEqual(runs["seed41/restart"]["result"], {"seed": 41})
            self.assertTrue(
                runs["seed41/restart"]["native_state_checked_before_stepping"]
            )
            self.assertFalse(
                runs["seed41/full"]["native_state_checked_before_stepping"]
            )

    def test_primary_failure_survives_shutdown_profiling(self):
        failure = "Traceback (most recent call last):\nTypeError: device array needs a host copy\n"
        output = failure + "Profiler memory usage\n" * 1000
        self.assertNotIn("TypeError", output[-6000:])
        self.assertIn(
            "TypeError: device array needs a host copy",
            export_reaudit.failure_context(output),
        )

    def test_archive_keeps_valid_results_and_records_incomplete_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            study = root / "gpu_failures-42"
            nested = study / "rank-0"
            nested.mkdir(parents=True)
            (study / "topology.json").write_text("")
            (study / "comparison.json").write_text('{"passed":')
            (nested / "incomplete.xml").write_text("<testsuite>")
            (nested / "valid.xml").write_text(
                '<testsuite><testcase name="good" time="1"/>'
                '<testcase name="bad"><failure/><system-out>failed assertion</system-out>'
                '</testcase><testcase name="geometry">'
                '<skipped message="Cartesian fixture">RZ is not configured</skipped>'
                "</testcase></testsuite>"
            )
            (nested / "run.log").write_text(
                "Traceback (most recent call last):\nValueError: original failure\n"
                + "Shutdown profiling\n" * 1000
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
                {"good": True, "bad": False, "geometry": False},
            )
            self.assertEqual(
                record["comparisons"], {"rank-0/summary.json": {"passed": False}}
            )
            self.assertEqual(
                record["tests"]["rank-0/valid.xml"][1]["failure_output"],
                "failed assertion",
            )
            self.assertEqual(
                record["tests"]["rank-0/valid.xml"][2]["skip_reason"],
                "Cartesian fixture\nRZ is not configured",
            )
            self.assertIn(
                "ValueError: original failure",
                record["control_logs"]["rank-0/run.log"]["failure_context"],
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
