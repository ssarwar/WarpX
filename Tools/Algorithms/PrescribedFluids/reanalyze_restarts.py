#!/usr/bin/env python3
"""Apply one canonical restart analysis to existing stock or branch outputs.

The CTest inventory supplies the original paths and tolerance arguments. No
simulation is rerun, no output is changed, and unsuccessful analyses are kept.
"""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("analysis", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for test in json.loads(args.inventory.read_text())["tests"]:
        command = test.get("command", [])
        if len(command) < 2:
            continue
        name = Path(command[1]).name
        if name == "analysis_default_restart.py":
            analysis = args.analysis.resolve()
        elif name == "analysis_restart_mw_start_stop.py":
            # Retain the upstream wrapper's domain-bound checks. Its local
            # import also selects the canonical particle-ID comparison.
            analysis = args.analysis.resolve().parent / "Tests" / "restart" / name
        else:
            continue
        properties = {item["name"]: item["value"] for item in test["properties"]}
        directory = Path(properties["WORKING_DIRECTORY"])
        if not directory.is_dir():
            continue
        command[1] = str(analysis)
        log = args.output / (test["name"] + ".log")
        with log.open("w") as stream:
            result = subprocess.run(
                command,
                cwd=directory,
                stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=180,
            )
        records.append(
            dict(
                test=test["name"],
                command=command,
                directory=str(directory),
                returncode=result.returncode,
                log=log.name,
                log_tail=log.read_text()[-2000:],
            )
        )
    assert records, "No existing restart outputs were selected"
    (args.output / "comparison.json").write_text(json.dumps(records, indent=2) + "\n")
    for record in records:
        print(record["test"], record["returncode"])
    raise SystemExit(int(any(record["returncode"] != 0 for record in records)))


if __name__ == "__main__":
    main()
