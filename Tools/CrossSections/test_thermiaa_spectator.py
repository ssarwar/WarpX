# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent source quadrature for discrete IAA spectator scattering.

Uses only the exported, real cross sections and SciPy/NumPy. Clebsch--Gordan
weights are evaluated through triple-Legendre quadrature instead of the
exporter's factorial formula. Validation output belongs in a build directory.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn
from thermiaa_spectator import ConditionalBundle

KB = 8.617333262145e-5
B = 0.0002477204284695341
REST = 510998.95069


class Reference:
    def __init__(self, directory, maximum_j):
        self.tables = [
            np.loadtxt(directory / f"rotation_0_{rank}.txt") for rank in [0, 2, 4, 6]
        ]
        self.elastic = np.loadtxt(directory / "elastic.txt")
        self.dcs = np.loadtxt(directory / "elastic_dcs.txt")
        self.initial = np.arange(maximum_j + 1)[:, None]
        self.final = self.initial + 2 * np.arange(-3, 4)
        self.loss = B * (
            self.final * (self.final + 1) - self.initial * (self.initial + 1)
        )
        x, w = np.polynomial.legendre.leggauss(128)
        legendre = np.polynomial.legendre.legvander(x, maximum_j + 6).T
        self.cg = np.zeros((maximum_j + 1, 7, 4))
        for n, rank in enumerate([0, 2, 4, 6]):
            product = (legendre * (w * legendre[rank])) @ legendre.T
            for i in range(maximum_j + 1):
                for choice, f in enumerate(self.final[i]):
                    if f >= 0 and abs(i - rank) <= f <= i + rank:
                        self.cg[i, choice, n] = (2 * f + 1) * product[i, f] / 2
        assert self.cg.min() >= 0
        self.norm_mu, self.norm_w = np.polynomial.legendre.leggauss(384)
        # Eight points on each tabulated DCS interval also integrate its
        # piecewise-linear-in-sin(theta/2) angular marginal accurately.
        angular = np.sin(np.linspace(0, np.pi / 2, self.dcs.shape[1] - 1))
        node, weight = np.polynomial.legendre.leggauss(8)
        self.y = (
            angular[:-1, None] + np.diff(angular)[:, None] * (node + 1) / 2
        ).ravel()
        self.measure = (np.diff(angular)[:, None] * weight / 2).ravel() * self.y
        self.angular = angular

    def population(self, temperature):
        j = np.arange(1024)
        p = np.zeros_like(j, dtype=float)
        if temperature == 0:
            p[0] = 1
        else:
            p = (
                (2 * j + 1)
                * np.where(j % 2, 3, 6)
                * np.exp(-B * j * (j + 1) / (KB * temperature))
            )
            p /= p.sum()
        assert p[len(self.initial) :].sum() < 1e-10
        return p[: len(self.initial)]

    def probability(self, energy, y):
        amplitude = []
        for rank, table in zip([0, 2, 4, 6], self.tables, strict=True):
            delta = B * rank * (rank + 1)
            if rank == 0:
                amplitude.append(np.interp(energy, table[:, 0], table[:, 1]))
            else:
                positive = table[:, 1] > 0
                e, sigma = table[positive][0]
                if energy <= e:
                    amplitude.append(
                        sigma
                        * np.sqrt(
                            e * (e + 2 * REST) / ((e - delta) * (e - delta + 2 * REST))
                        )
                    )
                else:
                    sigma = np.interp(energy, table[:, 0], table[:, 1])
                    amplitude.append(
                        sigma
                        * np.sqrt(
                            energy
                            * (energy + 2 * REST)
                            / ((energy - delta) * (energy - delta + 2 * REST))
                        )
                    )
        k_r = 2.068 * np.sqrt(energy * (energy + 2 * REST)) / (7.2973525693e-3 * REST)
        h = []
        for rank in [0, 2, 4, 6]:
            if k_r == 0:
                h.append((rank + 1) / 2 * y ** (2 * rank))
            else:
                norm = (
                    self.norm_w
                    @ spherical_jn(rank, k_r * np.sqrt((1 - self.norm_mu) / 2)) ** 2
                )
                h.append(spherical_jn(rank, k_r * y) ** 2 / norm)
        outgoing = np.maximum(energy - self.loss, 0)
        weights = (
            self.cg * amplitude * np.sqrt(outgoing * (outgoing + 2 * REST))[:, :, None]
        )
        score = (np.asarray(h).T @ weights.reshape(-1, 4).T).reshape(len(y), -1, 7)
        norm = score.sum(axis=2)
        p = np.divide(
            score,
            norm[:, :, None],
            out=np.zeros_like(score),
            where=norm[:, :, None] > 0,
        )
        p[:, :, 3] += norm == 0
        return p

    def angular_weights(self, energy):
        energy = max(energy, self.dcs[0, 0])
        i = np.clip(np.searchsorted(self.dcs[:, 0], energy) - 1, 0, len(self.dcs) - 2)
        t = np.clip(
            np.log(energy / self.dcs[i, 0])
            / np.log(self.dcs[i + 1, 0] / self.dcs[i, 0]),
            0,
            1,
        )
        density = np.interp(
            self.y, self.angular, (1 - t) * self.dcs[i, 1:] + t * self.dcs[i + 1, 1:]
        )
        mass = self.measure * density
        return mass / mass.sum()

    def moments(self, probability, temperature, mass):
        thermal = probability * self.population(temperature)[None, :, None]
        mu = 1 - 2 * self.y**2
        loss = self.loss.ravel()
        thermal = thermal.reshape(len(self.y), -1)
        value = np.column_stack(
            [
                thermal @ loss,
                thermal @ loss**2,
                mu,
                mu * (thermal @ loss),
                thermal[:, loss == 0].sum(axis=1),
                thermal[:, loss > 0].sum(axis=1),
                thermal[:, loss < 0].sum(axis=1),
                mu**2,
            ]
        )
        square = np.column_stack(
            [
                thermal @ loss**2,
                thermal @ loss**4,
                mu**2,
                mu**2 * (thermal @ loss**2),
                value[:, 4],
                value[:, 5],
                value[:, 6],
                mu**4,
            ]
        )
        return mass @ value, mass @ square


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    bundle = ConditionalBundle.read(args.bundle)
    reference = Reference(args.data_dir / "N2" / "IAA", bundle.maximum_j)
    rng = np.random.default_rng(932)
    probes = np.unique(
        np.r_[
            0.0008,
            0.0021,
            0.01,
            0.025,
            0.1,
            0.5,
            1.0,
            2.3,
            2.6,
            10.0,
            100.0,
            900.0,
            np.exp(rng.uniform(np.log(0.0005), np.log(950), 24)),
        ]
    )
    comparisons = []
    with (args.output / "reference.txt").open("w") as output:
        for energy in probes:
            p = reference.probability(energy, reference.y)
            tabulated = bundle.at(energy, reference.y).reshape(p.shape)
            mass = reference.angular_weights(energy)
            assert np.max(abs(tabulated.sum(axis=2) - 1)) < 5e-7
            assert np.max(tabulated[:, reference.loss > energy], initial=0) == 0
            for temperature in [0, 100, 300, 1000]:
                exact, squares = reference.moments(p, temperature, mass)
                actual, _ = reference.moments(tabulated, temperature, mass)
                # Bound interpolation relative to heating+cooling, so a small
                # net loss cannot conceal cancellation of two large powers.
                thermal = p * reference.population(temperature)[None, :, None]
                absolute_first = mass @ (
                    thermal.reshape(len(mass), -1) @ abs(reference.loss.ravel())
                )
                scales = np.array(
                    [absolute_first, exact[1], 1, absolute_first, 1, 1, 1, 1]
                )
                errors = abs(actual - exact) / np.maximum(scales, 1e-14)
                assert errors.max() < 0.002, (energy, temperature, errors)
                comparisons.append(float(errors.max()))
                output.write(
                    '"'
                    + str(args.bundle.resolve())
                    + '" '
                    + " ".join(
                        f"{x:.17g}" for x in [temperature, energy, 0, *exact, *squares]
                    )
                    + "\n"
                )
    result = {
        "cases": 4 * len(probes),
        "maximum_moment_error": max(comparisons),
        "bundle_bytes": args.bundle.stat().st_size,
        "definition": "Conditional kinetic Eq. 2.48; no assertion of equilibrium detailed balance",
    }
    (args.output / "quadrature.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
