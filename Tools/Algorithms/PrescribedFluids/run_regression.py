#!/usr/bin/env python3
"""Select the audit regressions and recursively include their prerequisites."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("build", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--list-only", action="store_true")
parser.add_argument(
    "--include", help="Restrict the audit selection before adding prerequisites"
)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
inventory = json.loads(
    subprocess.check_output(
        ["ctest", "--test-dir", str(args.build), "--show-only=json-v1"], text=True
    )
)
tests = {test["name"]: test for test in inventory["tests"]}
selection = re.compile(
    "background_mcc|proton_impact|collision|langmuir_fluid|implicit|restart|diagnostic|"
    "reduced_diags|field_probe|rigid_|immobile|fluid_|prescribed_|attachment_checkpoint|"
    "beam_reference|openpmd|plotfile"
)
chosen = {name for name in tests if selection.search(name) and "checksum" not in name}
if args.include:
    chosen = {name for name in chosen if re.search(args.include, name)}
excluded = {}
cache = (args.build / "CMakeCache.txt").read_text()
if re.search(r"^WarpX_QED:BOOL=OFF$", cache, re.MULTILINE):
    for name in sorted(chosen):
        if name.startswith(
            ("test_3d_collider_diagnostics.", "test_3d_beam_beam_collision.")
        ):
            excluded[name] = "Requires QED photon emission; WarpX_QED=OFF"
            chosen.remove(name)
pending = list(chosen)
aliases = {}
while pending:
    name = pending.pop()
    properties = {item["name"]: item["value"] for item in tests[name]["properties"]}
    for dependency in properties.get("DEPENDS", []):
        if dependency.endswith(".checksum"):
            # Checksums are platform-dependent read-only analyses. Retain the
            # physics analysis or run that actually produces the restart data.
            base = dependency.removesuffix(".checksum")
            replacement = (
                base + ".analysis" if base + ".analysis" in tests else base + ".run"
            )
            aliases[dependency] = replacement
            dependency = replacement
        if dependency not in tests and dependency.endswith(".analysis"):
            # Some stock restart tests name an analysis prerequisite even when
            # that test has analysis=OFF. Its run still creates the checkpoint.
            run = dependency.removesuffix(".analysis") + ".run"
            if run in tests:
                aliases[dependency] = run
                dependency = run
        if dependency not in chosen:
            assert dependency in tests, (name, dependency)
            chosen.add(dependency)
            pending.append(dependency)
assert chosen and not any("checksum" in name for name in chosen)
manifest = args.output.resolve() / "tests.txt"
manifest.write_text("\n".join(sorted(chosen)) + "\n")
(args.output / "test-inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
(args.output / "dependency-aliases.json").write_text(
    json.dumps(aliases, indent=2) + "\n"
)
(args.output / "excluded-tests.json").write_text(json.dumps(excluded, indent=2) + "\n")
print(f"Selected {len(chosen)} tests, including prerequisites", flush=True)
if not args.list_only:
    subprocess.run(
        [
            "ctest",
            "--test-dir",
            str(args.build),
            "--tests-from-file",
            str(manifest),
            "-j",
            "1",
            "--output-on-failure",
            "--timeout",
            "180",
            "--output-junit",
            str(args.output.resolve() / "ctest.xml"),
        ],
        check=True,
    )
