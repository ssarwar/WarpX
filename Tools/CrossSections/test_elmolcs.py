# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Validate real elmolcs exports and prepare independent device-test references.

Requires elmolcs only for the offline source comparison. Output belongs in a
build directory, never in warpx-data. Runtime tests read the exported files.
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from export_elmolcs import Source
from rotation_reference import KB, REST, ROTATION, Bundle, cg_squared
from scipy.special import eval_legendre


def linear_rates(bundle, energy):
    energy = np.asarray(energy)
    index = np.clip(
        np.searchsorted(bundle.energies, energy, side="right") - 1,
        0,
        len(bundle.energies) - 2,
    )
    fraction = (energy - bundle.energies[index]) / (
        bundle.energies[index + 1] - bundle.energies[index]
    )
    fraction = np.where(index == 0, np.sqrt(fraction), fraction)
    rates = bundle.rates[:, 0]
    return (
        rates[index] * (1 - fraction[..., None])
        + rates[index + 1] * fraction[..., None]
    )


def population(target, temperature, maximum):
    # Independently sum a longer partition function: do not normalize the
    # truncated bundle state list or use the exporter's convergence routine.
    j = np.arange(512)
    p = np.zeros(512)
    if temperature == 0:
        p[0 if target == "N2" else 1] = 1
    else:
        spin = np.where(j % 2 == 0, 6, 3) if target == "N2" else j % 2
        p = (
            (2 * j + 1)
            * spin
            * np.exp(-ROTATION[target] * j * (j + 1) / (KB * temperature))
        )
        p /= sum(p)
    assert sum(p[maximum + 1 :]) < 1e-10
    return p[: maximum + 1]


def angular_moments(path, energy):
    if energy == 0:
        return np.array([0, 1 / 3, 0, 1 / 5])
    table = np.loadtxt(path)
    x = np.log(table[:, 0])
    e = np.log(energy)
    i = np.clip(np.searchsorted(x, e) - 1, 0, len(x) - 2)
    t = np.clip((e - x[i]) / (x[i + 1] - x[i]), 0, 1)
    dcs = (1 - t) * table[i, 1:] + t * table[i + 1, 1:]
    y = np.sin(np.linspace(0, np.pi / 2, len(dcs)))
    nodes, w = np.polynomial.legendre.leggauss(8)
    width = np.diff(y)[:, None]
    u = y[:-1, None] + width * (nodes + 1) / 2
    density = u * (dcs[:-1, None] + (dcs[1:] - dcs[:-1])[:, None] * (nodes + 1) / 2)
    measure = density * width * w / 2
    mu = 1 - 2 * u * u
    return np.array([np.sum(measure * mu**n) / np.sum(measure) for n in range(1, 5)])


def validate(target, directory, output):
    source = Source(target)
    bundle = Bundle.read(directory / "thermal_rotation.rot")
    rng = np.random.default_rng(5803)
    # Adaptive points are checked at new, randomly located probes, not just
    # the midpoints used to build the table.
    errors = []
    sections = {"elastic.txt": (source.elastic, np.r_[source.residual.knots, 1e9])}
    for rank, cs in source.elementary.items():
        sections[f"rotation_{source.ground}_{source.ground + rank}.txt"] = (
            cs,
            cs.knots,
        )
    if source.unchanged:
        sections["rotation_0_0.txt"] = (source.unchanged, source.unchanged.knots)
    for name, (function, knots) in sections.items():
        table = np.loadtxt(directory / name)
        np.testing.assert_allclose(
            np.interp(knots, table[:, 0], table[:, 1]),
            function(knots),
            rtol=1e-9,
            atol=1e-33,
        )
        lo, hi = table[:-1, 0], table[1:, 0]
        probe = lo[:, None] + (hi - lo)[:, None] * rng.uniform(size=(len(lo), 5))
        exact = function(probe)
        actual = np.interp(probe, table[:, 0], table[:, 1])
        error = abs(actual - exact) / np.maximum(exact, function(knots).max() * 1e-6)
        resolved = (hi - lo) > 4 * np.finfo(np.float32).eps * np.maximum(lo, 1e-30)
        assert error[resolved].max() < 0.002, (name, error.max())
        errors.append(error[resolved].max())

    # Independent Clebsch--Gordan reference from a triple-Legendre integral.
    mu, w = np.polynomial.legendre.leggauss(128)
    for i in [0, 1, 3, 15, 47, source.maximum_j - 6]:
        for rank in source.elementary:
            for f in range(abs(i - rank), i + rank + 1, 2):
                integral = (
                    (2 * f + 1)
                    / 2
                    * np.dot(
                        w,
                        eval_legendre(i, mu)
                        * eval_legendre(rank, mu)
                        * eval_legendre(f, mu),
                    )
                )
                assert abs(integral - cg_squared(i, rank, f)) < 3e-11

    lo, hi = bundle.energies[:-1], bundle.energies[1:]
    probe = lo + (hi - lo) * rng.uniform(0.05, 0.95, len(lo))
    exact_rates = source.rates(probe)
    approximate_rates = linear_rates(bundle, probe)
    maximum_moment_error = 0
    resolved = np.diff(bundle.energies) > 4 * np.finfo(np.float32).eps * np.maximum(
        lo, 1e-30
    )
    for temperature in [0, 50, 100, 300, 600, 1000]:
        p = population(target, temperature, source.maximum_j)
        factors = np.r_[1, [p[i] for i, _ in bundle.transitions]]
        exact = exact_rates * factors
        actual = approximate_rates * factors
        for mask in [bundle.losses == 0, bundle.losses > 0, bundle.losses < 0]:
            for power in [0, 1, 2]:
                a = exact[:, mask] @ abs(bundle.losses[mask]) ** power
                b = actual[:, mask] @ abs(bundle.losses[mask]) ** power
                scale = max(a.max(), 1e-300)
                error = abs(a - b) / np.maximum(abs(a), 1e-5 * scale)
                # No float32 grid can resolve arbitrarily close approaches to
                # a square-root threshold. Bound these intervals separately
                # by their absolute contribution, instead of relaxing the
                # relative tolerance on resolved rates or moments.
                assert (
                    np.max(abs(a[~resolved] - b[~resolved]), initial=0) < 1e-6 * scale
                )
                maximum_moment_error = max(maximum_moment_error, error[resolved].max())
                if error[resolved].max() >= 0.002:
                    index = error.argmax()
                    print(
                        "moment interpolation failure",
                        target,
                        temperature,
                        power,
                        probe[index],
                        error[index],
                        lo[index],
                        hi[index],
                        flush=True,
                    )
    assert maximum_moment_error < 0.002, maximum_moment_error

    # Integrate both powers independently at the SAME incident energy grid.
    # This tests the table's reverse interpolation and Boltzmann factors, not
    # an algebraic cancellation using a shifted copy of the forward power.
    grid = np.unique(
        np.r_[0, np.geomspace(1e-10, 4, 32000), bundle.energies[bundle.energies < 4]]
    )
    rates = linear_rates(bundle, grid)
    loss = bundle.losses
    for temperature in [100, 300, 1000]:
        p = population(target, temperature, source.maximum_j)
        factors = np.r_[1, [p[i] for i, _ in bundle.transitions]]
        phase = np.sqrt(grid * (grid + 2 * REST)) * (grid + REST)
        distribution = phase * np.exp(-grid / (KB * temperature))
        distribution /= np.trapezoid(distribution, grid)
        powers = (
            rates
            * factors
            @ np.column_stack((np.maximum(loss, 0), np.maximum(-loss, 0)))
        )
        up, down = np.trapezoid(powers * distribution[:, None], grid, axis=0)
        imbalance = abs(up - down) / (up + down)
        assert imbalance < 0.001, (target, temperature, imbalance)
        print(
            target,
            temperature,
            "K equilibrium relative power imbalance:",
            imbalance,
            flush=True,
        )

    # Source knots include the N2 resonance and the O2 Born-model minimum.
    energies = [0, 0.0008, 0.0021, 0.01, 0.1, 0.5, 2.3, 10, source.maximum * 0.9]
    for temperature in [0, 100, 300, 1000]:
        p = population(target, temperature, source.maximum_j)
        for energy in energies:
            r = (
                source.rates(np.array([energy]))[0]
                * np.r_[1, [p[i] for i, _ in bundle.transitions]]
            )
            total = sum(r)
            if total == 0:
                continue
            probability = r / total
            moments = [np.dot(probability, loss**n) for n in range(1, 5)]
            angular = angular_moments(directory / "elastic_dcs.txt", energy)
            p0, pup, pdown = [
                sum(probability[m]) for m in [loss == 0, loss > 0, loss < 0]
            ]
            expected = [
                moments[0],
                moments[1],
                angular[0],
                angular[0] * moments[0],
                p0,
                pup,
                pdown,
                angular[1],
            ]
            square = [
                moments[1],
                moments[3],
                angular[1],
                angular[1] * moments[1],
                p0,
                pup,
                pdown,
                angular[3],
            ]
            fields = [temperature, energy, total, *expected, *square]
            output.write(
                '"'
                + str((directory / "thermal_rotation.rot").resolve())
                + '" '
                + " ".join(f"{v:.17g}" for v in fields)
                + "\n"
            )
    print(target, "maximum cross-section interpolation error:", max(errors), flush=True)
    print(
        target,
        "maximum rate/transfer-moment interpolation error:",
        maximum_moment_error,
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path, required=True, help="warpx-data/MCC_cross_sections"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.getLogger().setLevel(logging.ERROR)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "reference.txt").open("w") as output:
        for target in ["N2", "O2"]:
            validate(target, args.data_dir / target / "IAA", output)


if __name__ == "__main__":
    main()
