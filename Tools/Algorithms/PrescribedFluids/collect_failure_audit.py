#!/usr/bin/env python3
"""Archive restart-failure evidence and its Perlmutter runtime provenance."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def provenance(directory):
    files = [
        "revision.txt",
        "stock-revision.txt",
        "parameters.txt",
        "CMakeCache.txt",
        "python-library-sha256.txt",
        "source-sha256.json",
        "source.patch",
        "physics-executable-sha256.txt",
    ]
    result = {"directory": str(directory.resolve()), "files": {}}
    for name in files:
        path = directory / name
        if path.exists():
            result["files"][name] = hashlib.sha256(path.read_bytes()).hexdigest()
            if name.endswith(".txt") and name != "CMakeCache.txt":
                result[name.removesuffix(".txt")] = path.read_text().strip()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acceleration", type=Path, nargs="+", required=True)
    parser.add_argument("--psatd", type=Path, nargs="+", required=True)
    parser.add_argument("--existing", type=Path, required=True)
    parser.add_argument("--physics", type=Path, nargs="+", default=[])
    parser.add_argument("--build-job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--verify-plotfiles",
        action="store_true",
        help="Recompute the archived plotfile reports, including exact metadata checks.",
    )
    args = parser.parse_args()
    report = {"acceleration": [], "psatd": [], "existing": {}}
    jobs = [str(args.build_job)]
    for directory in args.acceleration:
        record = provenance(directory)
        record["comparisons"] = json.loads(
            (directory / "cases/comparison.json").read_text()
        )
        if args.verify_plotfiles:
            from report_plotfile_difference import compare

            for case in record["comparisons"]:
                assert "execution_error" not in case, case
                for result in case["comparisons"].values():
                    if "reference" in result:
                        assert compare(result["reference"], result["actual"]) == result
                assert all(
                    row["changed_elements"] == 0
                    for row in case["comparisons"]["restart_native_5"].values()
                ), case["test"]
            record["plotfile_reports_reproduced"] = True
        report["acceleration"].append(record)
        jobs.append(directory.name.rsplit("-", 1)[1])
    for directory in args.psatd:
        record = provenance(directory)
        record["run_status"] = dict(
            (name, int(code))
            for name, code in (
                line.split()
                for line in (directory / "status.txt").read_text().splitlines()
            )
        )
        record["topology"] = json.loads((directory / "topology.json").read_text())
        record["differences"] = {
            path.stem: json.loads(path.read_text())
            for path in sorted(directory.glob("*-state-difference.json"))
        }
        report["psatd"].append(record)
        jobs.append(directory.name.rsplit("-", 1)[1])
    for path in sorted(args.existing.glob("*.json")):
        report["existing"][path.stem] = json.loads(path.read_text())
    assert report["existing"], args.existing
    report["physics"] = []
    for directory in args.physics:
        record = provenance(directory)
        record["logs"] = {
            name: (directory / name).read_text()
            for name in (
                "fluid-tests.log",
                "pjg.log",
                "mcc-physics.log",
                "python-physics.log",
                "gaussian.log",
            )
        }
        report["physics"].append(record)
        jobs.append(directory.name.rsplit("-", 1)[1])
    report["scheduler"] = subprocess.check_output(
        [
            "sacct",
            "-j",
            ",".join(jobs),
            "--format=JobID,State,Elapsed,ExitCode",
            "-P",
        ],
        text=True,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
