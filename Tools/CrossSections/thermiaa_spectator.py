# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Offline, discrete Thetaermiaa rotational probabilities from real elmolcs data.

IAA Eqs. (11.24) and (11.35) supply differential transition weights. Each
initial state's open outcomes are normalized at the sampled scattering angle,
as probabilities in Eq. (2.48). Their integrals are consequently determined by
the inclusive elastic DCS; they are not constrained to the input rotational
integrals. No source model is evaluated by WarpX.
"""

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from export_elmolcs import Source
from rotation_reference import (
    REST,
    cg_squared,
    converged_j,
    populations,
    threshold,
)
from scipy.special import spherical_jn


class Spectator:
    """Finite-threshold sudden probabilities, normalized separately for each J.

    The supplied elementary set has ranks 0, 2, 4 and 6. This is a truncated
    spectator model; it does not establish convergence of higher-rank scattering.
    De-excitation follows the same incident-energy sudden kernel as excitation.
    The subsequent probability normalization does not impose detailed balance.
    """

    def __init__(self, maximum_temperature=1000.0):
        from elmolcs import reader

        self.source = Source("N2")
        self.maximum_j = converged_j("N2", maximum_temperature)
        self.ranks = np.array([0, 2, 4, 6])
        self.transitions = np.array(
            [(i, f) for i in range(self.maximum_j + 1) for f in range(i - 6, i + 7, 2)],
            dtype=np.int32,
        )
        self.initial = self.transitions[:, 0]
        self.starts = np.r_[0, np.flatnonzero(np.diff(self.initial)) + 1]
        self.diagonal = self.transitions[:, 0] == self.transitions[:, 1]
        self.losses = threshold("N2", *self.transitions.T)
        self.coefficients = np.array(
            [
                [cg_squared(i, rank, f) for rank in self.ranks]
                for i, f in self.transitions
            ]
        )
        self.population = np.array(
            [
                populations("N2", temperature, self.maximum_j)[0][self.initial]
                for temperature in [0, 100, 300, maximum_temperature]
            ]
        )
        self.maximum = min(
            self.source.unchanged.maximum,
            *(section.maximum for section in self.source.elementary.values()),
        )
        # Kutz--Meyer rigid-rotor internuclear distance, in Bohr radii.
        self.distance = 2.068
        self.mu, self.quadrature_weights = np.polynomial.legendre.leggauss(128)
        dcs = reader.readDCS("N2", format="grid")[0]["data"]
        self.dcs_energies = np.asarray(dcs.columns, float)
        self.dcs = np.asarray(dcs).T
        self.dcs_angles = np.sin(np.asarray(dcs.index, float) * np.pi / 360)
        functions = np.array(
            [
                (self.losses == 0).astype(float),
                (self.losses > 0).astype(float),
                (self.losses < 0).astype(float),
                np.maximum(self.losses, 0),
                np.maximum(-self.losses, 0),
                self.losses**2,
            ]
        )
        self.moment_weights = (self.population[:, None, :] * functions).reshape(
            -1, len(self.losses)
        )

    def differential_weights(self, energy):
        """Separable transition weights; the common incident momentum cancels."""
        if not 0 <= energy <= self.maximum:
            raise ValueError("Energy is outside the elementary source range")
        outgoing = np.maximum(energy - self.losses, 0)
        phase = np.sqrt(outgoing * (outgoing + 2 * REST))
        amplitude = np.array(
            [
                self.source.unchanged(energy),
                *(
                    section.reduced(energy)
                    for section in self.source.elementary.values()
                ),
            ]
        )
        values = self.coefficients * amplitude * phase[:, None]
        scale = values.max()
        if scale > 0:
            values /= scale
        return values.reshape(self.maximum_j + 1, 7, 4)

    def angular_basis(self, energy, half_sine):
        half_sine = np.atleast_1d(half_sine)
        k_r = (
            self.distance
            * np.sqrt(energy * (energy + 2 * REST))
            / (7.2973525693e-3 * REST)
        )
        if k_r == 0:
            # Normalized small-argument limit of j_L(k R sin(theta/2))^2.
            angular = np.array(
                [(rank + 1) / 2 * half_sine ** (2 * rank) for rank in self.ranks]
            ).T
        else:
            norm = np.array(
                [
                    self.quadrature_weights
                    @ spherical_jn(rank, k_r * np.sqrt((1 - self.mu) / 2)) ** 2
                    for rank in self.ranks
                ]
            )
            angular = (
                np.array(
                    [spherical_jn(rank, k_r * half_sine) ** 2 for rank in self.ranks]
                ).T
                / norm
            )
        return angular

    def probabilities(self, energy, half_sine):
        """Return P(J->J'; E, theta), before Boltzmann population weighting."""
        return conditional_probabilities(
            self.differential_weights(energy), self.angular_basis(energy, half_sine)
        )

    def moments(self, probabilities):
        """Heating/cooling probabilities and absolute transfer moments."""
        return probabilities @ self.moment_weights.T

    def angular_density(self, energy, half_sine):
        """IAA DCS shape, with the runtime's interpolation in log(E) and sin(theta/2)."""
        energy = max(energy, self.dcs_energies[0])
        i = np.clip(
            np.searchsorted(self.dcs_energies, energy) - 1, 0, len(self.dcs) - 2
        )
        t = np.clip(
            np.log(energy / self.dcs_energies[i])
            / np.log(self.dcs_energies[i + 1] / self.dcs_energies[i]),
            0,
            1,
        )
        return np.interp(
            half_sine, self.dcs_angles, (1 - t) * self.dcs[i] + t * self.dcs[i + 1]
        )


def conditional_probabilities(weights, angular_basis):
    score = (angular_basis @ weights.reshape(-1, 4).T).reshape(
        len(angular_basis), weights.shape[0], 7
    )
    normalization = score.sum(axis=-1)
    result = np.divide(
        score,
        normalization[:, :, None],
        out=np.zeros_like(score),
        where=normalization[:, :, None] > 0,
    )
    # The singular zero-score endpoint has zero elastic rate or angular measure.
    result[:, :, 3] += normalization == 0
    return result.reshape(len(angular_basis), -1)


@dataclass
class ConditionalBundle:
    maximum_j: int
    energies: np.ndarray
    offsets: np.ndarray
    angles: np.ndarray
    basis: np.ndarray
    weights: np.ndarray

    @property
    def losses(self):
        initial = np.arange(self.maximum_j + 1)[:, None]
        return threshold("N2", initial, initial + 2 * np.arange(-3, 4)).ravel()

    def validate(self):
        if self.energies[0] != 0 or np.any(
            np.diff(self.energies.astype(np.float32)) <= 0
        ):
            raise ValueError(
                "Conditional energies must start at zero and remain distinct"
            )
        if (
            self.offsets.shape != (len(self.energies) + 1,)
            or self.offsets[0] != 0
            or self.offsets[-1] != len(self.angles)
            or np.any(np.diff(self.offsets) < 2)
        ):
            raise ValueError("Invalid conditional angular offsets")
        if self.weights.shape != (
            len(self.energies),
            self.maximum_j + 1,
            7,
            4,
        ) or self.basis.shape != (len(self.angles), 4):
            raise ValueError("Invalid conditional weight dimensions")
        for values in [self.weights, self.basis]:
            if not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError("Invalid conditional weight")
        for row, (e, start, end) in enumerate(
            zip(self.energies, self.offsets[:-1], self.offsets[1:], strict=True)
        ):
            angular = self.angles[start:end]
            if angular[0] != 0 or angular[-1] != 1 or np.any(np.diff(angular) <= 0):
                raise ValueError("Invalid sin(theta/2) grid")
            if np.any(self.weights[row].reshape(-1, 4)[self.losses > e] != 0):
                raise ValueError("Subthreshold excitation weight")
        for j in range(min(6, self.maximum_j + 1)):
            if np.any(self.weights[:, j, np.arange(7) < 3 - j // 2] != 0):
                raise ValueError("Negative final rotational state")

    def write(self, path):
        self.validate()
        with Path(path).open("wb") as output:
            output.write(
                (
                    "WARPX_THERMAL_ROTATION_V5\n"
                    f"N2 iaa_spectator {self.maximum_j} 0\n"
                    f"{len(self.energies)} {len(self.angles)}\n"
                ).encode("ascii")
            )
            for values, dtype in [
                (self.energies, "<f8"),
                (self.offsets, "<i4"),
                (self.angles, "<f4"),
                (self.basis, "<f4"),
                (self.weights, "<f4"),
            ]:
                output.write(np.asarray(values, dtype=dtype).tobytes())

    @classmethod
    def read(cls, path):
        with Path(path).open("rb") as source:
            if source.readline() != b"WARPX_THERMAL_ROTATION_V5\n":
                raise ValueError("Expected a conditional V5 bundle")
            target, model, maximum_j, reserved = source.readline().decode().split()
            if (target, model, reserved) != ("N2", "iaa_spectator", "0"):
                raise ValueError("Unknown conditional model")
            ne, na = map(int, source.readline().split())
            energy = np.fromfile(source, "<f8", ne)
            offsets = np.fromfile(source, "<i4", ne + 1)
            angles = np.fromfile(source, "<f4", na)
            basis = np.fromfile(source, "<f4", 4 * na).reshape(na, 4)
            weights = np.fromfile(
                source, "<f4", ne * (int(maximum_j) + 1) * 28
            ).reshape(ne, int(maximum_j) + 1, 7, 4)
            if source.read(1):
                raise ValueError("Trailing conditional data")
        result = cls(int(maximum_j), energy, offsets, angles, basis, weights)
        result.validate()
        return result

    def at(self, energy, half_sine):
        index = np.clip(
            np.searchsorted(self.energies, energy, side="right") - 1,
            0,
            len(self.energies) - 2,
        )
        fraction = (energy - self.energies[index]) / (
            self.energies[index + 1] - self.energies[index]
        )
        if index == 0:
            fraction = np.sqrt(fraction)
        rows = []
        for row in [index, index + 1]:
            start, end = self.offsets[row : row + 2]
            angle = self.angles[start:end]
            j = np.clip(
                np.searchsorted(angle, half_sine, side="right") - 1, 0, len(angle) - 2
            )
            f = (half_sine - angle[j]) / (angle[j + 1] - angle[j])
            basis = (1 - f[:, None]) * self.basis[start + j] + f[:, None] * self.basis[
                start + j + 1
            ]
            rows.append(
                conditional_probabilities(self.weights[row].astype(float), basis)
            )
        return (1 - fraction) * rows[0] + fraction * rows[1]


def adaptive_angles(model, energy, tolerance):
    grid = np.linspace(0, 1, 9)
    nodes, weights = np.polynomial.legendre.leggauss(3)
    for _ in range(24):
        basis = model.angular_basis(energy, grid)
        fraction = (nodes + 1) / 2
        width = np.diff(grid)
        probes = grid[:-1, None] + width[:, None] * fraction
        reference = model.moments(model.probabilities(energy, probes.ravel())).reshape(
            len(width), len(nodes), -1
        )
        interpolated_basis = (1 - fraction[None, :, None]) * basis[
            :-1, None
        ] + fraction[None, :, None] * basis[1:, None]
        interpolated = model.moments(
            conditional_probabilities(
                model.differential_weights(energy), interpolated_basis.reshape(-1, 4)
            )
        ).reshape(len(width), len(nodes), -1)
        # dOmega is proportional to y dy. Resolve the error in the angularly
        # averaged collision moments, not pointwise probabilities in negligible cones.
        measure = (
            probes * model.angular_density(energy, probes) * width[:, None] * weights
        )
        measure /= measure.sum()
        error = np.sum(measure[:, :, None] * abs(reference - interpolated), axis=1)
        scale = np.sum(measure[:, :, None] * reference, axis=(0, 1))
        goal = tolerance * np.maximum(scale, 1e-14)
        bad = error.sum(axis=0) > goal
        split = (error[:, bad] > goal[bad] / (2 * len(width))).any(axis=1)
        if not split.any():
            return grid.astype(np.float32), basis.astype(np.float32)
        grid = np.sort(np.r_[grid, ((grid[:-1] + grid[1:]) / 2)[split]])
    raise ValueError("Conditional angular grid did not converge")


def construct(model, tolerance=5e-4):
    knots = np.unique(
        np.r_[
            0,
            np.geomspace(1e-8, model.maximum, 120),
            model.losses[model.losses > 0],
            model.source.knots[model.source.knots <= model.maximum],
            model.source.unchanged.knots,
        ]
    )
    mu, quadrature = np.polynomial.legendre.leggauss(192)
    probes = np.sqrt((1 - mu) / 2)
    cache = {}

    def feature(e):
        if e not in cache:
            cache[e] = model.moments(model.probabilities(e, probes))
        return cache[e]

    def measure(e):
        density = quadrature * model.angular_density(e, probes)
        return density / density.sum()

    def rate(e):
        return (
            model.source.residual(e)
            * 299792458
            * np.sqrt(e * (e + 2 * REST))
            / (e + REST)
        )

    peak = np.max(
        np.array([rate(e) * (measure(e) @ feature(e)) for e in knots]), axis=0
    )

    for _ in range(32):
        values = np.array([feature(e) for e in knots])
        midpoint = (knots[:-1] + knots[1:]) / 2
        midpoint[0] = knots[1] / 4
        reference = np.array([feature(e) for e in midpoint])
        interpolated = (values[:-1] + values[1:]) / 2
        mass = np.array([measure(e) for e in midpoint])
        velocity_rate = np.array([rate(e) for e in midpoint])
        error = velocity_rate[:, None] * np.sum(
            mass[:, :, None] * abs(reference - interpolated), axis=1
        )
        exact_moments = velocity_rate[:, None] * np.sum(
            mass[:, :, None] * reference, axis=1
        )
        scale = np.maximum(exact_moments, peak * 1e-5)
        split = (error > tolerance * np.maximum(scale, 1e-300)).any(axis=1)
        split &= np.diff(knots) > 4 * np.finfo(np.float32).eps * np.maximum(
            knots[:-1], 1e-30
        )
        if not split.any():
            break
        knots = np.sort(np.r_[knots, midpoint[split]])
    else:
        raise ValueError("Conditional energy grid did not converge")
    print(f"Conditional energy grid: {len(knots)} rows", flush=True)
    angles, basis, weights = [], [], []
    offsets = [0]
    for e in knots:
        x, h = adaptive_angles(model, e, tolerance)
        angles.append(x)
        basis.append(h)
        weights.append(model.differential_weights(e).astype(np.float32))
        offsets.append(offsets[-1] + len(x))
    return ConditionalBundle(
        model.maximum_j,
        knots,
        np.array(offsets, np.int32),
        np.concatenate(angles),
        np.concatenate(basis),
        np.array(weights),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum-temperature", type=float, default=1000)
    parser.add_argument("--tolerance", type=float, default=5e-4)
    args = parser.parse_args()
    logging.getLogger().setLevel(logging.ERROR)
    model = Spectator(args.maximum_temperature)
    bundle = construct(model, args.tolerance)
    bundle.write(args.output)
    print(
        f"{args.output}: {len(bundle.angles)} angular nodes, {args.output.stat().st_size / 1e6:.1f} MB",
        flush=True,
    )


if __name__ == "__main__":
    main()
