#!/usr/bin/env python3
"""Run every configured regression, retaining checksum and physics results.

Checksums remain in the CTest dependency graph so restart prerequisites run
in order, including when CTest runs in parallel. Their failures are reported
separately from the physical assertions, as required by the repository policy.
An optional launcher runs the Python unit suites in an allocated MPI/GPU step.
"""

import argparse
import json
import os
import re
import shlex
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--pytest-launcher")
    args = parser.parse_args()
    args.build = args.build.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    inventory = json.loads(
        subprocess.check_output(
            ["ctest", "--test-dir", str(args.build), "--show-only=json-v1"], text=True
        )
    )
    (args.output / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    tests = {test["name"]: test for test in inventory["tests"]}
    cache = (args.build / "CMakeCache.txt").read_text()
    excluded, launched = {}, {}
    for name, test in tests.items():
        if re.search(r"^WarpX_QED:BOOL=OFF$", cache, re.MULTILINE) and name.startswith(
            ("test_3d_collider_diagnostics.", "test_3d_beam_beam_collision.")
        ):
            excluded[name] = "Requires QED photon emission; WarpX_QED=OFF"
        elif args.pytest_launcher and name.startswith("pytest.WarpX."):
            launched[name] = test
    chosen = sorted(set(tests) - set(excluded) - set(launched))
    assert chosen
    manifest = args.output / "tests.txt"
    manifest.write_text("\n".join(chosen) + "\n")
    report = dict(excluded=excluded, launched_units={}, ctest_tests=len(chosen))
    summary = args.output / "summary.json"
    summary.write_text(json.dumps(report, indent=2) + "\n")
    command = [
        "ctest",
        "--test-dir",
        str(args.build),
        "--tests-from-file",
        str(manifest),
        "-j",
        str(args.jobs),
        "--output-on-failure",
        "--timeout",
        str(args.timeout),
        "--output-junit",
        str(args.output / "ctest.xml"),
    ]
    with (args.output / "ctest.log").open("w") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    report["ctest_returncode"] = result.returncode
    cases = ET.parse(args.output / "ctest.xml").getroot().findall(".//testcase")
    assert len(cases) == len(chosen), (len(cases), len(chosen))
    failures = [
        case.attrib["name"] for case in cases if case.find("failure") is not None
    ]
    report["checksum_failures"] = [
        name for name in failures if name.endswith(".checksum")
    ]
    report["physics_failures"] = [
        name for name in failures if not name.endswith(".checksum")
    ]
    report["skipped"] = [
        case.attrib["name"] for case in cases if case.find("skipped") is not None
    ]
    summary.write_text(json.dumps(report, indent=2) + "\n")
    for name, test in launched.items():
        properties = {item["name"]: item["value"] for item in test["properties"]}
        environment = os.environ.copy()
        environment.update(
            item.split("=", 1) for item in properties.get("ENVIRONMENT", [])
        )
        command = [
            *shlex.split(args.pytest_launcher),
            *test["command"],
            "--junitxml=" + str(args.output / (name + ".xml")),
        ]
        with (args.output / (name + ".log")).open("w") as stream:
            result = subprocess.run(
                command,
                env=environment,
                cwd=properties["WORKING_DIRECTORY"],
                stdout=stream,
                stderr=subprocess.STDOUT,
            )
        report["launched_units"][name] = dict(
            command=command, returncode=result.returncode
        )
        summary.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    raise SystemExit(
        int(
            bool(report["physics_failures"] or report["skipped"])
            or any(item["returncode"] for item in report["launched_units"].values())
        )
    )


if __name__ == "__main__":
    main()
