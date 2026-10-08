#!/usr/bin/env python3
"""Reproduce GPU suite failures and vary only the identified execution conditions.

Original assertions are retained. Four-rank controls keep the mesh fixed while
reducing boxes per GPU. ID unpacking explicitly copies to host. Reflection
probabilities zero and one give exact independent particle-count expectations.
"""

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

from discharge_convergence import run_simulation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--build-name", default="build_pm_gpu_sync")
    parser.add_argument(
        "--versions",
        nargs="+",
        choices=["branch", "stock"],
        default=["branch", "stock"],
    )
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=["original", "placement", "ids", "reflection"],
        default=["original", "placement", "ids", "reflection"],
    )
    args = parser.parse_args()
    args.root, args.output = args.root.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    records = []
    summary = args.output / "summary.json"

    def run_case(tests, version, variant, base, ranks=None, edit=None):
        directory = args.output / version / variant / base
        directory.mkdir(parents=True, exist_ok=False)
        record = dict(version=version, variant=variant, test=base, stages=[])
        for stage in ["run", "analysis"]:
            test = tests.get(base + "." + stage)
            if test is None:
                continue
            if stage == "analysis" and record["stages"][0]["returncode"]:
                record["analysis_not_run"] = (
                    "Simulation failed; do not analyze stale output"
                )
                break
            properties = {item["name"]: item["value"] for item in test["properties"]}
            environment = os.environ.copy()
            overrides = dict(
                item.split("=", 1) for item in properties.get("ENVIRONMENT", [])
            )
            environment.update(overrides)
            environment["MPLBACKEND"] = "Agg"
            command = list(test["command"])
            if stage == "run":
                if ranks is not None:
                    command[command.index("-n") + 1] = str(ranks)
                if edit is not None:
                    index = next(
                        i for i, item in enumerate(command) if item.endswith(".py")
                    )
                    source = Path(command[index]).read_text()
                    adapted = edit(source)
                    script = directory / "inputs.py"
                    script.write_text(adapted)
                    record["original_input_sha256"] = hashlib.sha256(
                        source.encode()
                    ).hexdigest()
                    record["adapted_input_sha256"] = hashlib.sha256(
                        adapted.encode()
                    ).hexdigest()
                    record["input_changed"] = adapted != source
                    command[index] = str(script)
            started = time.monotonic()
            returncode = run_simulation(
                command, directory, environment, 600, log_name=stage + ".log"
            )
            record["stages"].append(
                dict(
                    stage=stage,
                    command=command,
                    environment=overrides,
                    returncode=returncode,
                    seconds=time.monotonic() - started,
                )
            )
        records.append(record)
        summary.write_text(json.dumps(records, indent=2) + "\n")
        print(
            version,
            variant,
            base,
            [stage["returncode"] for stage in record["stages"]],
            flush=True,
        )

    def host_ids(source):
        old = '    idcpu = pti["idcpu"]'
        new = "    idcpu = pti.soa().get_idcpu_data().to_numpy(copy=True)"
        if source.count(old) == 1:
            return source.replace(old, new)
        assert source.count(new) == 1
        return source

    def exact_reflection(source, probability):
        for old in [
            'warpx_reflection_model_zhi="0.5"',
            "sim.step(max_steps)",
            "assert n == 63",
            "assert n == 67",
        ]:
            assert source.count(old) == 1, old
        source = source.replace(
            'warpx_reflection_model_zhi="0.5"',
            f'warpx_reflection_model_zhi="{probability}"',
        )
        source = source.replace(
            "sim.step(max_steps)",
            "sim.initialize_inputs()\nsim.initialize_warpx()\n"
            'initial_count = sim.particles.get("electrons").number_of_particles()\n'
            "assert initial_count > 0\nsim.step(max_steps)",
        )
        source = source.replace(
            "assert n == 63",
            "assert n == " + ("initial_count" if probability == 0 else "0"),
        )
        return source.replace(
            "assert n == 67",
            "assert n == " + ("0" if probability == 0 else "initial_count"),
        )

    for version in args.versions:
        root = args.root if version == "branch" else args.root / "build/stock-warpx"
        build = root / args.build_name
        inventory = json.loads(
            subprocess.check_output(
                ["ctest", "--test-dir", str(build), "--show-only=json-v1"], text=True
            )
        )
        (args.output / (version + "-inventory.json")).write_text(
            json.dumps(inventory, indent=2) + "\n"
        )
        tests = {test["name"]: test for test in inventory["tests"]}
        if "original" in args.variants and version == "stock":
            for base in [
                "test_2d_collisions_split_momentum_push_electromagnetic",
                "test_2d_theta_implicit_jfnk_vandb",
                "test_2d_theta_implicit_jfnk_vandb_filtered",
                "test_2d_langmuir_multi_mr_psatd",
                "test_2d_particle_reflection_picmi",
                "test_2d_id_cpu_read_picmi",
                "test_2d_pml_x_psatd",
                "test_2d_pml_x_psatd_restart",
            ]:
                run_case(tests, version, "original", base)
        if "placement" in args.variants:
            for base in [
                "test_2d_theta_implicit_jfnk_vandb",
                "test_2d_theta_implicit_jfnk_vandb_filtered",
                "test_2d_langmuir_multi_mr_psatd",
            ]:
                run_case(tests, version, "four-ranks", base, ranks=4)
        if "ids" in args.variants:
            run_case(
                tests, version, "host-ids", "test_2d_id_cpu_read_picmi", edit=host_ids
            )
        if "reflection" in args.variants:
            for probability in [0, 1]:
                run_case(
                    tests,
                    version,
                    f"reflection-{probability}",
                    "test_2d_particle_reflection_picmi",
                    edit=lambda source, probability=probability: exact_reflection(
                        source, probability
                    ),
                )
    raise SystemExit(
        int(
            any(
                stage["returncode"]
                for row in records
                if row["variant"] != "original"
                for stage in row["stages"]
            )
        )
    )


if __name__ == "__main__":
    main()
