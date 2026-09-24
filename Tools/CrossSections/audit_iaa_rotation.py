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
    spectator_bins,
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
                for kind, name in [("cs", "Cross section"), ("dcs", "DCS")]
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
    # Resolve individual angular-bin integrals using source interpolation in
    # sin(theta/2) and log(E), independently of the WarpX inverse-CDF tables.
    for target in ["N2", "O2"]:
        data = reader.readCS(target, db="iaa*", skip=False)
        residual = next(
            r["data"]
            for r in data
            if r["kind"] == "ELASTIC" and r["subkind"] == "RESIDUAL"
        )
        dcs = reader.readDCS(target, format="grid")[0]["data"]
        if dcs.index.name != "Angle" or dcs.columns.name != "Energy":
            raise ValueError("Unexpected DCS grid axes")
        source_energy = np.asarray(dcs.columns, dtype=float)
        source_angles = np.asarray(dcs.index, dtype=float)
        if source_angles[0] != 0 or source_angles[-1] != 180:
            raise ValueError("Expected source DCS angles in degrees, spanning [0,180]")
        source_y = np.sin(np.deg2rad(source_angles) / 2)
        values = np.asarray(dcs, dtype=float).T
        edges = np.linspace(-1, 1, 65)
        nodes, weights = np.polynomial.legendre.leggauss(8)
        mu = (edges[1:, None] + edges[:-1, None]) / 2 + np.diff(edges)[
            :, None
        ] * nodes / 2
        y = np.sqrt((1 - mu) / 2)
        angular_rows = np.array([np.interp(y, source_y, row) for row in values])
        bins = (angular_rows * weights * np.diff(edges)[:, None] / 2).sum(axis=-1)

        def inclusive(energies):
            shapes = np.array(
                [
                    np.interp(
                        np.log(np.maximum(energies, source_energy[0])),
                        np.log(source_energy),
                        bins[:, a],
                    )
                    for a in range(len(edges) - 1)
                ]
            ).T
            shapes /= shapes.sum(axis=1)[:, None]
            sigma = np.interp(
                energies,
                np.asarray(residual.index, float),
                np.asarray(residual["CS"], float),
            )
            return shapes * sigma[:, None]

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
                    return amplitude_e[:, None] * spectator_bins(
                        energies, rank, 2.0743, edges
                    )

                elementary[rank] = reduced
            model = "iaa_sudden_spectator"
        else:
            # Eq. 11.21b at J=0, with the outgoing/incoming momentum ratio
            # removed. The q^2 term and all constants are in atomic units.
            def born(energies):
                q2 = 4 * np.asarray(energies)[:, None, None] / HARTREE * (1 - mu)
                dcs = (4 / 5) * (-0.29 / 3 + 4.93 * np.pi * q2 / 32) ** 2 * A0**2
                return (dcs * (2 * np.pi) * weights * np.diff(edges)[:, None] / 2).sum(
                    axis=-1
                )

            elementary[2] = born
            model = "iaa_born"
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
                "Negative unchanged angular-bin rate in the inclusive decomposition"
                if np.any(bundle.rates[:, :, 0] < 0)
                else str(error)
            )
        # A numerical positivity pass cannot establish the missing physics.
        errors.extend(
            [
                "Low-energy/resonant angular approximation has no validated unified continuation",
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
            "minimum_at_cosine_interval": edges[index[1] : index[1] + 2].tolist(),
        }
        write_audit(args.output / f"{target}_rotation_audit.json", report)
        print(target, "production export blocked:", *errors, sep="\n  ")


if __name__ == "__main__":
    main()
