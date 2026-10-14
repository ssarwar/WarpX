# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Source-reader classes used by the WarpX offline exporter."""

import numpy as np
from rotation_reference import (
    REST,
    C,
    cg_squared,
    converged_j,
    populations,
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
