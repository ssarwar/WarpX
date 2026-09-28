# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Offline conversion of elmolcs elastic and rotational cross sections.

Run with elmolcs on PYTHONPATH. No elmolcs dependency is used by WarpX.
Only cross sections and loadable collision data are written to --output;
diagnostics are printed, never installed in the data repository.
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from rotation_reference import (
    REST,
    Bundle,
    C,
    cg_squared,
    converged_j,
    energy_grid,
    populations,
    refine_grid,
    speed,
    threshold,
    weights,
)


class CrossSection:
    """Positive source knots with log-log interpolation and a threshold law.

    A finite elastic intercept is joined linearly to the first positive energy.
    For excitation, interpolate sigma*p_in/p_out, preserving every positive
    source value and the square-root opening at the canonical threshold.
    The reduced amplitude is held at its first value below the first datum.
    Values beyond the source endpoint are errors, not endpoint extrapolations.
    """

    def __init__(self, table, loss=0.0):
        self.loss = loss
        energy = np.asarray(table.index, float)
        sigma = np.asarray(table["CS"], float)
        self.zero = sigma[0] if energy[0] == 0 and loss == 0 else 0.0
        valid = (energy > loss) & (sigma > 0)
        self.energy = energy[valid]
        self.sigma = sigma[valid]
        self.amplitude = self.sigma / self.phase(self.energy)
        self.knots = np.r_[loss, self.energy]
        self.maximum = self.energy[-1]

    def phase(self, energy):
        energy = np.asarray(energy)
        out = np.maximum(energy - self.loss, 0)
        return np.sqrt(out * (out + 2 * REST) / (energy * (energy + 2 * REST)))

    def reduced(self, energy):
        energy = np.asarray(energy)
        if np.any(energy > self.maximum * (1 + 4 * np.finfo(float).eps)):
            raise ValueError("Requested energy exceeds the rotational source")
        return np.exp(
            np.interp(
                np.log(np.maximum(energy, self.energy[0])),
                np.log(self.energy),
                np.log(self.amplitude),
            )
        )

    def __call__(self, energy):
        energy = np.asarray(energy)
        if self.loss:
            return self.reduced(energy) * self.phase(np.maximum(energy, self.loss))
        values = self.reduced(energy)
        return np.where(
            energy < self.energy[0],
            self.zero + energy / self.energy[0] * (self.sigma[0] - self.zero),
            values,
        )


def tabulate(function, knots, tolerance=5e-4):
    """Adapt to the ordinary WarpX linear interpolator, retaining source knots."""
    knots = np.unique(knots)
    scale = np.max(function(knots))
    fractions = np.array([0.25, 0.5, 0.75])
    for _ in range(32):
        values = function(knots)
        probes = knots[:-1, None] + np.diff(knots)[:, None] * fractions
        exact = function(probes)
        interpolated = (
            values[:-1, None] * (1 - fractions) + values[1:, None] * fractions
        )
        split = (
            np.abs(interpolated - exact) > tolerance * np.maximum(exact, 1e-6 * scale)
        ).any(axis=1)
        split &= np.diff(knots) > 4 * np.finfo(np.float32).eps * np.maximum(
            knots[:-1], 1e-30
        )
        if not split.any():
            return knots, values
        knots = np.sort(np.r_[knots, probes[split, 1]])
    raise ValueError("Cross-section interpolation did not converge")


def write_cross_section(path, function, knots, description):
    energy, sigma = tabulate(function, knots)
    np.savetxt(
        path,
        np.column_stack((energy, sigma)),
        fmt=["%.17e", "%.12e"],
        header=description + "\nEnergy (eV), cross section (m^2); linear interpolation",
    )
    print(
        f"{path.name}: {len(energy)} points, {energy[0]:g}--{energy[-1]:g} eV",
        flush=True,
    )
    return energy


class Source:
    """Source-resolved rates; all expensive state sums stay in the offline tool."""

    def __init__(self, target):
        from elmolcs import cscoll, reader

        self.target = target
        rows = reader.readCS(target, db="iaa*", skip=False)
        self.residual = CrossSection(
            next(r["data"] for r in rows if r.get("subkind") == "RESIDUAL")
        )
        self.born = cscoll.Elastic.gen_Born(target)
        self.born_scale = float(
            self.residual(self.residual.maximum) / self.born(self.residual.maximum)
        )
        self.elementary = {}
        if target == "N2":
            for rank in [2, 4, 6]:
                row = next(r for r in rows if r.get("final") == f"J=0-{rank}")
                self.elementary[rank] = CrossSection(
                    row["data"], threshold(target, 0, rank)
                )
            self.unchanged = CrossSection(
                next(r["data"] for r in rows if r.get("final") == "J=0-0")
            )
            self.ground = 0
        else:
            # Use the supplied integral table, not the previous replacement
            # polynomial for the Born DCS. The latter did not reproduce elmolcs.
            numerical = reader.readCS(target, db="iaa", skip=False)
            row = next(r for r in numerical if r.get("final") == "J=1-3")
            self.elementary[2] = CrossSection(row["data"], threshold(target, 1, 3))
            self.unchanged = None
            self.ground = 1
        self.maximum_j = converged_j(target, 1000)
        states = [j for j in range(self.maximum_j + 1) if weights(target, j)]
        self.transitions = [
            (i, f)
            for i in states
            for f in states
            if 0 < f - i <= max(self.elementary) and (f - i) % 2 == 0
        ]
        self.loss = np.array([threshold(target, i, f) for i, f in self.transitions])
        self.coefficients = np.array(
            [
                [
                    cg_squared(i, rank, f)
                    / cg_squared(self.ground, rank, self.ground + rank)
                    for rank in self.elementary
                ]
                for i, f in self.transitions
            ]
        )
        self.ratio = np.array(
            [weights(target, i) / weights(target, f) for i, f in self.transitions]
        )
        # Reverse rates need forward data at E+Delta. Reserve that margin;
        # never extrapolate or silently clamp rotational tails.
        self.maximum = min(s.maximum for s in self.elementary.values()) - max(self.loss)
        self.knots = np.unique(
            np.concatenate(
                [self.residual.knots, *(s.knots for s in self.elementary.values())]
            )
        )
        self.pairs = [pair for i, f in self.transitions for pair in [(i, f), (f, i)]]
        self.diagonal = np.array(
            [
                [
                    cg_squared(j, rank, j)
                    / cg_squared(self.ground, rank, self.ground + rank)
                    for rank in self.elementary
                ]
                for j in states
            ]
        )
        self.pairs.extend((j, j) for j in states)
        self.losses = np.r_[
            0, np.column_stack((self.loss, -self.loss)).ravel(), np.zeros(len(states))
        ]

    def elastic(self, energy):
        """Residual table followed continuously by the package's fitted Born model."""
        energy = np.asarray(energy)
        endpoint = self.residual.maximum
        return np.where(
            energy <= endpoint,
            self.residual(np.minimum(energy, endpoint)),
            self.born_scale * self.born(np.maximum(energy, endpoint)),
        )

    def ground_inelastic(self, energy):
        return sum(s(energy) for s in self.elementary.values())

    def inclusive(self, energy):
        if self.unchanged is None:
            return self.residual(energy)
        return self.unchanged(energy) + self.ground_inelastic(energy)

    def rates(self, energies):
        energies = np.asarray(energies)
        reduced = np.stack(
            [s.reduced(energies) for s in self.elementary.values()], axis=-1
        )
        outgoing = np.maximum(energies[:, None] - self.loss, 0)
        upward = (
            C
            * np.sqrt(outgoing * (outgoing + 2 * REST))
            / (energies[:, None] + REST)
            * (reduced @ self.coefficients.T)
        )
        shifted = energies[:, None] + self.loss
        reverse = sum(
            s.reduced(shifted) * self.coefficients[None, :, rank]
            for rank, s in enumerate(self.elementary.values())
        )
        downward = (
            self.ratio
            * C
            * np.sqrt(shifted * (shifted + 2 * REST))
            / (energies[:, None] + REST)
            * reverse
        )
        isotropic_rank = (
            self.unchanged(energies)
            if self.unchanged is not None
            else self.residual(energies)
            - self.ground_inelastic(energies)
            - reduced @ self.diagonal[0]
        )
        if np.any(isotropic_rank < 0):
            raise ValueError(
                "Rotational source exceeds the selected inclusive elastic rate"
            )
        unchanged = isotropic_rank[:, None] + reduced @ self.diagonal.T
        return np.column_stack(
            (
                np.zeros(len(energies)),
                np.stack((upward, downward), axis=-1).reshape(len(energies), -1),
                speed(energies)[:, None] * unchanged,
            )
        )

    def moments(self, energies):
        rates = self.rates(energies)
        result = []
        for temperature in [0, 100, 300, 1000]:
            p, _ = populations(self.target, temperature, self.maximum_j)
            weighted = rates * np.r_[1, [p[i] for i, _ in self.pairs]]
            for mask in [self.losses == 0, self.losses > 0, self.losses < 0]:
                result.extend(
                    weighted[:, mask] @ np.abs(self.losses[mask]) ** power
                    for power in [0, 1, 2]
                )
        return np.stack(result, axis=-1)

    def bundle(self):
        seed = energy_grid(
            self.target, self.maximum_j, self.elementary, self.maximum, count=192
        )
        seed = np.unique(np.r_[seed, self.knots[self.knots <= self.maximum]])
        grid = refine_grid(seed, self.moments, tolerance=5e-4, relative_floor=1e-5)
        return Bundle(
            self.target,
            "elastic_dcs",
            self.maximum_j,
            0,
            grid,
            np.array([-1.0, 1.0]),
            self.pairs,
            self.rates(grid)[:, None, :],
        )


def export(target, output, include_bundle):
    from elmolcs import reader

    output.mkdir(parents=True, exist_ok=True)
    source = Source(target)
    write_cross_section(
        output / "elastic.txt",
        source.elastic,
        np.r_[source.residual.knots, np.geomspace(6000, 1e9, 80)[1:]],
        f"e + {target}: vibrationally elastic residual; elmolcs iaa*, Schmalzried (2023), Eq. 11.10.\nAbove 6000 eV: elmolcs Elastic.gen_Born, matched continuously to the last residual datum.",
    )
    # Retain the original angular/energy resolution: the source has only
    # 75--106 rows. Above 10 keV WarpX uses its analytic angular continuation.
    dcs = reader.readDCS(target, format="grid")[0]["data"]
    np.savetxt(
        output / "elastic_dcs.txt",
        np.column_stack((np.asarray(dcs.columns, float), np.asarray(dcs).T)),
        fmt="%.10e",
        header=f"SPECIES: e / {target}\nelmolcs IAA vibrationally elastic DCS; Schmalzried (2023), Sec. 12.1.\nRows: energy (eV), followed by DCS (m^2/sr) at 0, 0.5, ..., 180 degrees.",
    )
    print(target, "DCS:", dcs.shape[1], "energies,", dcs.shape[0], "angles", flush=True)
    for rank, section in source.elementary.items():
        i, f = source.ground, source.ground + rank
        write_cross_section(
            output / f"rotation_{i}_{f}.txt",
            section,
            section.knots,
            f"e + {target}, J={i}->{f}; elmolcs {'iaa* (Itikawa/Kutz-Meyer)' if target == 'N2' else 'iaa (Takayanagi-Itikawa)'}; canonical rigid-rotor threshold.",
        )
    if source.unchanged is not None:
        write_cross_section(
            output / "rotation_0_0.txt",
            source.unchanged,
            source.unchanged.knots,
            "e + N2, J=0->0; elmolcs iaa*, Kutz and Meyer (1995), Fig. 7a.",
        )
    if include_bundle:
        bundle = source.bundle()
        bundle.write(output / "thermal_rotation.rot")
        # A common grid keeps the inclusive cross-section check and the rate
        # table consistent, without exporting hundreds of MCC processes.
        write_cross_section(
            output / "thermal_rotation_elastic.txt",
            source.inclusive,
            bundle.energies,
            f"e + {target}: ground-rotor reference inclusive rate for thermal_rotation.rot; source-matched family, not an additional elastic channel.",
        )
        print(
            target,
            "rotation:",
            len(bundle.energies),
            "energies,",
            len(bundle.transitions),
            "transitions",
            flush=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, required=True, help="warpx-data/MCC_cross_sections"
    )
    parser.add_argument(
        "--test-bundles",
        action="store_true",
        help="also prepare source-driven sampler test families in a build directory; these do not resolve the source normalization conflict",
    )
    args = parser.parse_args()
    logging.getLogger().setLevel(logging.ERROR)
    for target in ["N2", "O2"]:
        export(target, args.output / target / "IAA", args.test_bundles)


if __name__ == "__main__":
    main()
