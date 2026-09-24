# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Audit an IAA rotational bundle before permitting production export.

Run with the archived elmolcs source on PYTHONPATH. Failed decompositions and
unresolved coverage are written to JSON; no production bundle is emitted.
"""

import argparse
import hashlib
import tarfile
from pathlib import Path

import numpy as np
from elmolcs import reader
from rotation_reference import (
    A0,
    HARTREE,
    ROTATION,
    construct,
    converged_j,
    energy_grid,
    write_audit,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    digest = hashlib.sha256(args.archive.read_bytes()).hexdigest()
    if digest != "b864381086120fff37eb8c658dfcb0f150b80969cf9add626a38d89bab0a3c38":
        raise ValueError("The archive is not the inspected elmolcs snapshot")
    # Verify the actual import and data roots too; merely possessing the archive
    # must not allow an unrelated installed package to generate its namesake data.
    package = Path(reader.__file__).resolve().parent
    with tarfile.open(args.archive) as archive:
        required = {
            "src/elmolcs/reader.py": package / "reader.py",
            "src/elmolcs/tablehandler.py": package / "tablehandler.py",
            **{
                f"Data/{kind}/{gas}/{name}.e-{gas}": Path(reader.ROOT)
                / kind
                / gas
                / f"{name}.e-{gas}"
                for gas in ["N2", "O2"]
                for kind, name in [("cs", "Cross section")]
            },
        }
        for suffix, local in required.items():
            members = [
                m
                for m in archive.getmembers()
                if m.isfile() and m.name.endswith("/" + suffix)
            ]
            if (
                len(members) != 1
                or archive.extractfile(members[0]).read() != local.read_bytes()
            ):
                raise ValueError(
                    f"Source file does not match archived snapshot: {local}"
                )
    args.output.mkdir(parents=True, exist_ok=True)
    # Only integral rates enter the rotational bundle. The ordinary elastic
    # DCS remains the sole angular source in the collision operator.
    for target in ["N2", "O2"]:
        data = reader.readCS(target, db="iaa*", skip=False)
        residual = next(
            r["data"]
            for r in data
            if r["kind"] == "ELASTIC" and r["subkind"] == "RESIDUAL"
        )
        edges = np.array([-1.0, 1.0])

        def inclusive(energies):
            return np.interp(
                energies,
                np.asarray(residual.index, float),
                np.asarray(residual["CS"], float),
            )[:, None]

        elementary = {}
        observations = []
        if target == "N2":
            for rank in [2, 4, 6]:
                row = next(r for r in data if r.get("final") == f"J=0-{rank}")
                table = row["data"]
                e = np.asarray(table.index, float)
                sigma = np.asarray(table["CS"], float)
                threshold = ROTATION[target] * rank * (rank + 1)
                valid = (e > threshold) & (sigma > 0)
                e, sigma = e[valid], sigma[valid]
                amplitude = sigma / np.sqrt(1 - threshold / e)
                observations.append(
                    {
                        "rank": rank,
                        "threshold_eV": threshold,
                        "peak_energy_eV": float(e[np.argmax(sigma)]),
                        "first_positive_energy_eV": float(e[0]),
                        "maximum_eV": float(e[-1]),
                    }
                )

                def reduced(energies, e=e, amplitude=amplitude, rank=rank):
                    if np.max(energies) > e[-1]:
                        raise ValueError(
                            "Unspecified elementary high-energy continuation"
                        )
                    # Explicit threshold continuation in the reduced amplitude.
                    # Its low-energy validity remains a production gate below.
                    amplitude_e = np.exp(
                        np.interp(
                            np.log(np.maximum(energies, e[0])),
                            np.log(e),
                            np.log(amplitude),
                        )
                    )
                    return amplitude_e[:, None]

                elementary[rank] = reduced
            model = "elastic_dcs"
        else:
            # Angular integral of Eq. 11.21b at J=0, with the outgoing/incoming
            # momentum ratio removed. This supplies a Born integral-rate model;
            # it is not used to sample angles. Its validity remains sub-eV.
            def born_integral(energies):
                a = -0.29 / 3
                b = 4.93 * np.pi * np.asarray(energies) / (8 * HARTREE)
                return ((16 * np.pi / 5) * A0**2 * (a * a + 2 * a * b + 4 * b * b / 3))[
                    :, None
                ]

            elementary[2] = born_integral
            model = "elastic_dcs"
        maximum_j = converged_j(target, 1000)
        energies = energy_grid(target, maximum_j, elementary, 20, count=192)
        bundle = construct(
            target, energies, edges, elementary, inclusive, maximum_j, model=model
        )
        errors = []
        try:
            bundle.validate()
        except ValueError as error:
            errors.append(
                "Negative unchanged integral rate in the inclusive decomposition"
                if np.any(bundle.rates[:, :, 0] < 0)
                else str(error)
            )
        # A numerical positivity pass cannot establish the missing physics.
        errors.extend(
            [
                "Low-energy/resonant integral-rate continuation requires validation",
                "No bound on omitted higher-rank rotational energy-transfer moments",
                "No validated rotational continuation from 20 eV through beam electron energies",
            ]
        )
        index = np.unravel_index(
            np.argmin(bundle.rates[:, :, 0]), bundle.rates[:, :, 0].shape
        )
        report = {
            "target": target,
            "model": model,
            "archive_sha256": digest,
            "reference_temperature_K": 0,
            "reference_temperature_is_model_assumption": True,
            "maximum_temperature_K": 1000,
            "maximum_j": maximum_j,
            "elementary_transitions": observations,
            "production_ready": False,
            "gates": errors,
            "minimum_unchanged_rate_m3_s": float(bundle.rates[index[0], index[1], 0]),
            "minimum_at_energy_eV": float(energies[index[0]]),
            "angular_sampling": "existing elastic DCS",
            "conditional_model": "integral rates conditioned only on exact recoil accessibility",
        }
        write_audit(args.output / f"{target}_rotation_audit.json", report)
        print(target, "production export blocked:", *errors, sep="\n  ")


if __name__ == "__main__":
    main()
