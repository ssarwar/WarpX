#!/usr/bin/env python3
"""Archive compact Perlmutter evidence without copying checkpoints or raw arrays.

Keep REAUDIT.md's interpretation alongside this machine-readable record.
Slurm completion is not physics acceptance: XML failures and numerical
comparisons are retained even when another phase in the same job passed.
"""

import argparse
import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path


def read_record(path, directory, errors, *, xml=False):
    """Retain interrupted output as evidence, without treating it as a pass."""
    contents = path.read_bytes()
    try:
        return ET.fromstring(contents) if xml else json.loads(contents)
    except (ET.ParseError, json.JSONDecodeError, UnicodeDecodeError) as error:
        errors[str(path.relative_to(directory))] = dict(
            error=str(error),
            sha256=hashlib.sha256(contents).hexdigest(),
            bytes=len(contents),
            tail=contents[-3000:].decode(errors="replace"),
        )
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, nargs="+", required=True)
    args = parser.parse_args()
    records = {}
    for job in args.jobs:
        record = {"job": job, "studies": []}
        accounting = subprocess.check_output(
            [
                "sacct",
                "-X",
                "-j",
                str(job),
                "--noheader",
                "--parsable2",
                "--format=JobID,State,ExitCode,Elapsed,AllocNodes,AllocCPUS,NodeList",
            ],
            text=True,
        ).strip()
        record["slurm"] = accounting
        for directory in sorted(args.audit.glob(f"*-{job}")):
            if not directory.is_dir():
                continue
            study = {
                "directory": directory.name,
                "provenance": {},
                "tests": {},
                "comparisons": {},
                "log_summaries": {},
                "unreadable_records": {},
            }
            for name in [
                "revision.txt",
                "stock-revision.txt",
                "validation-build.txt",
                "python-library-sha256.txt",
                "physics-build.txt",
                "physics-executable-sha256.txt",
                "library-sha256.txt",
                "compiler.txt",
                "device.txt",
                "runtime.txt",
                "modules.txt",
                "sanitizer.txt",
                "results.txt",
                "regression-driver-sha256.txt",
                "checkpoint-filesystem-layout.txt",
            ]:
                path = directory / name
                if path.is_file():
                    study["provenance"][name] = path.read_text()
            for name in ["source.patch", "source-sha256.json", "python-packages.txt"]:
                path = directory / name
                if path.is_file():
                    study["provenance"][name + ".sha256"] = hashlib.sha256(
                        path.read_bytes()
                    ).hexdigest()
            topology = directory / "topology.json"
            if topology.is_file():
                parsed = read_record(topology, directory, study["unreadable_records"])
                if parsed is not None:
                    study["provenance"]["topology"] = parsed
            for cache in sorted(directory.glob("*CMakeCache.txt")):
                selected = {}
                for line in cache.read_text().splitlines():
                    if re.match(
                        r"(?:WarpX_|AMReX_|MKL_|CMAKE_(?:BUILD_TYPE|CXX_COMPILER|CUDA_ARCHITECTURES)|MPIEXEC_)",
                        line,
                    ):
                        key, _, value = line.partition("=")
                        selected[key] = value
                study["provenance"][cache.name] = selected
            for xml in sorted(directory.rglob("*.xml")):
                parsed = read_record(
                    xml, directory, study["unreadable_records"], xml=True
                )
                if parsed is None:
                    continue
                tests = []
                for test in parsed.iter("testcase"):
                    tests.append(
                        dict(
                            name=test.attrib["name"],
                            seconds=float(test.get("time", 0)),
                            passed=(
                                test.find("failure") is None
                                and test.find("error") is None
                                and test.find("skipped") is None
                            ),
                            skipped=test.find("skipped") is not None,
                            failure_output=(test.findtext("system-out") or "")[-6000:]
                            if test.find("failure") is not None
                            or test.find("error") is not None
                            else None,
                        )
                    )
                study["tests"][str(xml.relative_to(directory))] = tests
            for name in ["python-physics.log", "gaussian.log", "analysis.log"]:
                path = directory / name
                if path.is_file():
                    study["log_summaries"][name] = path.read_text()[-3000:]
            for name in [
                "summary.json",
                "comparison.json",
                "acceptance.json",
                "excluded-tests.json",
                "*state-difference.json",
                "*load-balance-difference.json",
            ]:
                for path in sorted(directory.rglob(name)):
                    parsed = read_record(path, directory, study["unreadable_records"])
                    if parsed is not None:
                        study["comparisons"][str(path.relative_to(directory))] = parsed
            # Profiling phases intentionally keep synchronized timings separate
            # from the ordinary, uninstrumented timing ensembles. Preserve the
            # printed table headings as well as source/MCC/deposition rows so
            # minimum/average/maximum rank times cannot be mistaken for samples.
            if directory.name.startswith("profile-"):
                study["profiles"] = {}
                for path in sorted(directory.rglob("result.json")):
                    parsed = read_record(path, directory, study["unreadable_records"])
                    if parsed is None:
                        continue
                    log = path.with_name("run.log")
                    study["profiles"][str(path.parent.relative_to(directory))] = dict(
                        result=parsed,
                        log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
                        timer_rows=[
                            line
                            for line in log.read_text().splitlines()
                            if re.search(
                                r"^\s*Name\s+|(?:Call|Excl|Incl|Min|Avg|Max).*Time|"
                                r"ProtonImpactIonizationCollision::|BackgroundMCCCollision::|"
                                r"WarpXFluidContainer::CommitDensityIncrement|"
                                r"MultiFluidContainer::PrepareImmobileCharge|"
                                r"DepositCharge|DepositCurrent",
                                line,
                            )
                        ],
                    )
            if directory.name.startswith("memcheck-"):
                study["sanitizer_inputs"] = {
                    name: dict(
                        sha256=hashlib.sha256(
                            (directory / name).read_bytes()
                        ).hexdigest(),
                        text=(directory / name).read_text(),
                    )
                    for name in ["mpi-control.py", "positive.cu", "suppression.xml"]
                    if (directory / name).is_file()
                }
                study["sanitizer_logs"] = {
                    str(path.relative_to(directory)): dict(
                        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        head=path.read_text()[:12000],
                        tail=path.read_text()[-12000:],
                    )
                    for path in sorted(directory.rglob("run.log"))
                }
            record["studies"].append(study)
        # Build-only jobs can have no study directory.
        suite_status = args.audit / f"suite-{job}.status"
        if suite_status.is_file():
            record["suite_status"] = suite_status.read_text()
        records[str(job)] = record
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(records, indent=2) + "\n")
    print(f"Archived {len(records)} jobs in {args.output}")


if __name__ == "__main__":
    main()
