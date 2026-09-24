# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent rigid-rotor identities and finite-mass detailed-balance checks."""

import sys
from pathlib import Path

import numpy as np
from scipy.special import eval_legendre

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "Tools/CrossSections"))
from rotation_reference import (  # noqa: E402
    MASSES,
    QE,
    REST,
    ROTATION,
    C,
    cg_squared,
    construct,
    converged_j,
    energy_grid,
    populations,
    threshold,
)


def analytic_bundle(target="N2", maximum=10, count=768):
    maximum_j = converged_j(target, 1000)
    edges = np.linspace(-1, 1, 17)
    integrals = np.diff(edges) / 2 + 0.15 * np.diff(edges**2)

    def reduced(energies):
        return np.broadcast_to(1e-20 * integrals, (len(energies), len(integrals)))

    def inclusive(energies):
        return np.broadcast_to(
            2e-20 * np.diff(edges) / 2, (len(energies), len(integrals))
        )

    energies = energy_grid(target, maximum_j, [2], maximum, count=count)
    return construct(
        target,
        energies,
        edges,
        {2: reduced},
        inclusive,
        maximum_j,
        model="analytic_test",
    )


def write_sampler_reference(bundle, folder):
    path = folder / f"{bundle.target}_analytic.rot"
    bundle.write(path)
    edges = bundle.edges
    mu1 = (edges[1:] + edges[:-1]) / 2
    mu2 = (edges[1:] ** 2 + edges[1:] * edges[:-1] + edges[:-1] ** 2) / 3
    mu4 = np.diff(edges**5) / (5 * np.diff(edges))
    loss = bundle.losses[None, :]
    with (folder / "reference.txt").open("a") as output:
        for temperature in [0, 100, 300, 1000]:
            rates = bundle.at_temperature(temperature)
            for energy in [0, 1e-5, 0.003, 0.05, 0.5, 5]:
                index = max(
                    0,
                    min(
                        np.searchsorted(bundle.energies, energy) - 1,
                        len(bundle.energies) - 2,
                    ),
                )
                fraction = (energy - bundle.energies[index]) / np.diff(bundle.energies)[
                    index
                ]
                joint = (1 - fraction) * rates[index] + fraction * rates[index + 1]
                rate = joint.sum()
                if rate == 0:
                    continue
                distribution = joint / rate
                values = [
                    loss,
                    loss**2,
                    mu1[:, None],
                    mu1[:, None] * loss,
                    loss == 0,
                    loss > 0,
                    loss < 0,
                    mu2[:, None],
                ]
                squares = [
                    loss**2,
                    loss**4,
                    mu2[:, None],
                    mu2[:, None] * loss**2,
                    loss == 0,
                    loss > 0,
                    loss < 0,
                    mu4[:, None],
                ]
                moments = [float((distribution * v).sum()) for v in values + squares]
                output.write(
                    f"{path.resolve()} {temperature} {energy:.17g} {rate:.17g} "
                    + " ".join(f"{v:.17g}" for v in moments)
                    + "\n"
                )


def main():
    if len(sys.argv) > 1:
        folder = Path(sys.argv[1])
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "reference.txt").write_text("")
    # The integral of three Legendre polynomials is twice the squared 3j;
    # this quadrature is independent of the implementation's factorial formula.
    x, w = np.polynomial.legendre.leggauss(128)
    for initial in [0, 1, 7, 23, 60]:
        for rank in [2, 4, 6]:
            for final in range(abs(initial - rank), initial + rank + 1, 2):
                independent = (
                    (2 * final + 1)
                    / 2
                    * np.dot(
                        w,
                        eval_legendre(initial, x)
                        * eval_legendre(rank, x)
                        * eval_legendre(final, x),
                    )
                )
                np.testing.assert_allclose(
                    cg_squared(initial, rank, final), independent, rtol=1e-9, atol=1e-12
                )
    for target in ["N2", "O2"]:
        for temperature in [0, 100, 300, 1000]:
            maximum_j = converged_j(target, temperature)
            p, tail = populations(target, temperature, maximum_j)
            reference, _ = populations(target, temperature, 1024)
            assert np.abs(p - reference[: len(p)]).sum() < 1e-10
            assert reference[len(p) :].sum() < tail + 1e-300
            assert abs(p.sum() - 1) < 1e-14
            if target == "O2":
                assert np.all(p[::2] == 0)
        bundle = analytic_bundle(target)
        bundle.validate()
        # Deterministic Maxwellian quadrature avoids Monte Carlo cancellation
        # noise in the zero-field power balance. Recoil and moving neutrals
        # are tested separately in the portable kinematics and MCC drivers.
        grid = np.r_[0, np.geomspace(1e-10, 10, 15000)]
        for temperature in [100, 300, 1000]:
            thermal = bundle.at_temperature(temperature).sum(axis=1)
            rates = np.array(
                [np.interp(grid, bundle.energies, row) for row in thermal.T]
            ).T
            loss = bundle.losses
            distribution = np.sqrt(grid) * np.exp(
                -grid / (8.617333262145e-5 * temperature)
            )
            cooling = np.trapezoid(
                distribution * (rates * np.maximum(loss, 0)).sum(axis=1), grid
            )
            heating = np.trapezoid(
                distribution * (rates * np.maximum(-loss, 0)).sum(axis=1), grid
            )
            assert abs(cooling - heating) / (cooling + heating) < 0.001
            for electron_temperature, sign in [
                (temperature / 2, -1),
                (temperature * 2, 1),
            ]:
                distribution = np.sqrt(grid) * np.exp(
                    -grid / (8.617333262145e-5 * electron_temperature)
                )
                power = np.trapezoid(distribution * (rates * loss).sum(axis=1), grid)
                assert sign * power > 0
        for channel in range(0, len(bundle.transitions), 2):
            initial, final = bundle.transitions[channel]
            loss = ROTATION[target] * (final * (final + 1) - initial * (initial + 1))
            mass = MASSES[target] * C**2 / QE + ROTATION[target] * initial * (
                initial + 1
            )
            down_e = 0.07
            up_e = (mass + loss) / mass * down_e + threshold(target, initial, final)
            # Direct Kallen momenta at one invariant s, using factored terms.
            invariant = (REST + mass) ** 2 + 2 * mass * up_e
            p_up2 = mass**2 * up_e * (up_e + 2 * REST) / invariant
            p_down2 = (mass + loss) ** 2 * down_e * (down_e + 2 * REST) / invariant
            velocity_up = C * np.sqrt(up_e * (up_e + 2 * REST)) / (up_e + REST)
            velocity_down = C * np.sqrt(down_e * (down_e + 2 * REST)) / (down_e + REST)
            rate_up = np.interp(
                up_e, bundle.energies, bundle.rates[:, :, channel + 1].sum(axis=1)
            )
            rate_down = np.interp(
                down_e, bundle.energies, bundle.rates[:, :, channel + 2].sum(axis=1)
            )
            balance = (
                (2 * initial + 1)
                * p_up2
                * rate_up
                / velocity_up
                / ((2 * final + 1) * p_down2 * rate_down / velocity_down)
            )
            assert abs(balance - 1) < 0.002, (target, initial, final, balance)
        if len(sys.argv) > 1:
            folder = Path(sys.argv[1])
            folder.mkdir(parents=True, exist_ok=True)
            write_sampler_reference(bundle, folder)
        print(
            target,
            "population, angular identities, thresholds and detailed balance passed",
        )


if __name__ == "__main__":
    main()
