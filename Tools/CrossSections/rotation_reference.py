# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Host reference kernels for a thermal rigid rotor.

Rates are integral v*sigma values in m3/s. Positive losses excite
the molecule. Reverse rates are constructed before thermal averaging. The
binary bundle stores state-resolved rates; WarpX applies the requested
Boltzmann populations and builds its sampling tables during initialization.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import gammaln

KB = 8.617333262145e-5
REST = 510998.95069
C = 299792458.0
A0 = 5.29177210544e-11
HARTREE = 27.211386245988
MASSES = {"N2": 28.0134 * 1.66053906660e-27, "O2": 31.9988 * 1.66053906660e-27}
QE = 1.602176634e-19
ROTATION = {"N2": 0.0002477204284695341, "O2": 0.00017828927734694198}


def weights(target, j):
    j = np.asarray(j)
    spin = np.where(j % 2 == 0, 6, 3) if target == "N2" else j % 2
    return (2 * j + 1) * spin


def populations(target, temperature, maximum_j):
    if temperature < 0 or not np.isfinite(temperature):
        raise ValueError("Rotational temperature must be finite and nonnegative")
    j = np.arange(maximum_j + 1)
    ground = 0 if target == "N2" else 1
    p = np.zeros(j.size)
    if temperature == 0:
        p[ground] = 1
        return p, 0.0
    a = ROTATION[target] / (KB * temperature)
    unnormalized = weights(target, j) * np.exp(
        -a * (j * (j + 1) - ground * (ground + 1))
    )
    partition = unnormalized.sum()
    # Upper bound from the first omitted term plus the decreasing continuum
    # integral. The factor six also bounds the O2 spin weight.
    first = maximum_j + 1
    tail = (
        6
        * (2 * first + 1 + 1 / a)
        * np.exp(-a * (first * (first + 1) - ground * (ground + 1)))
    )
    if a * (2 * first + 1) ** 2 < 2:
        tail = np.inf
    return unnormalized / partition, tail / partition


def converged_j(target, temperature, tolerance=1e-10):
    for maximum in range(8, 4097, 8):
        if populations(target, temperature, maximum)[1] < tolerance:
            return maximum
    raise ValueError("Rotational population did not converge below J=4096")


def cg_squared(initial, rank, final):
    """Squared C(initial,0;rank,0|final,0), from the factorial 3j identity."""
    if abs(initial - rank) > final or initial + rank < final:
        return 0.0
    total = initial + rank + final
    if total % 2:
        return 0.0
    half = total // 2
    differences = np.array([half - initial, half - rank, half - final])
    return (2 * final + 1) * np.exp(
        2 * (gammaln(half + 1) - gammaln(differences + 1).sum())
        + gammaln(2 * differences + 1).sum()
        - gammaln(total + 2)
    )


def speed(energy):
    energy = np.asarray(energy)
    return C * np.sqrt(energy * (energy + 2 * REST)) / (energy + REST)


def threshold(target, initial, final):
    loss = ROTATION[target] * (final * (final + 1) - initial * (initial + 1))
    mass = MASSES[target] * C**2 / QE + ROTATION[target] * initial * (initial + 1)
    return loss * (1 + REST / mass) + loss**2 / (2 * mass)


def phase_rate(energy, loss, mass):
    """v(E)*p_out_COM/p_in_COM, including a finite exothermic E=0 limit."""
    cutoff = loss * (1 + REST / mass) + loss**2 / (2 * mass)
    first = np.maximum(energy - cutoff, 0)
    second = np.maximum(
        energy + 2 * REST - loss * (1 - REST / mass) - loss**2 / (2 * mass), 0
    )
    return C * np.sqrt(first * second) / (energy + REST)


def thermal_rates(target, temperature, maximum_j, transitions, rates):
    p, tail = populations(target, temperature, maximum_j)
    if tail > 1e-10:
        raise ValueError(f"Unresolved rotational population: {tail}")
    factors = np.array([1.0, *(p[initial] for initial, _ in transitions)])
    return rates * factors[None, None, :]


@dataclass
class Bundle:
    target: str
    model: str
    maximum_j: int
    reference_temperature: float
    energies: np.ndarray
    edges: np.ndarray
    transitions: list
    rates: np.ndarray

    @property
    def losses(self):
        return np.array(
            [
                0,
                *(
                    ROTATION[self.target] * (f * (f + 1) - i * (i + 1))
                    for i, f in self.transitions
                ),
            ]
        )

    def at_temperature(self, temperature):
        return thermal_rates(
            self.target, temperature, self.maximum_j, self.transitions, self.rates
        )

    def validate(self):
        if self.model not in ("elastic_dcs", "analytic_test"):
            raise ValueError("Unknown rotational source model")
        if self.energies[0] != 0 or np.any(
            np.diff(self.energies.astype(np.float32)) <= 0
        ):
            raise ValueError(
                "Energy knots must start at zero and remain distinct in float32"
            )
        if len(self.edges) != 2:
            raise ValueError("Elastic-DCS rotation bundles contain only integral rates")
        if (
            self.edges[0] != -1
            or self.edges[-1] != 1
            or np.any(np.diff(self.edges) <= 0)
        ):
            raise ValueError("Cosine bins must span [-1,1]")
        if self.rates.shape != (
            len(self.energies),
            len(self.edges) - 1,
            len(self.transitions) + 1,
        ):
            raise ValueError("Invalid rate-array dimensions")
        if not np.isfinite(self.rates).all() or np.any(self.rates < 0):
            index = np.unravel_index(np.argmin(self.rates), self.rates.shape)
            raise ValueError(
                f"Negative or nonfinite inclusive decomposition at {index}: {self.rates[index]}"
            )
        for channel, loss in enumerate(self.losses):
            if loss > 0:
                physical_threshold = threshold(
                    self.target, *self.transitions[channel - 1]
                )
                if np.any(
                    self.rates[self.energies <= physical_threshold, :, channel] != 0
                ):
                    raise ValueError(
                        "Excitation rate is nonzero below its canonical threshold"
                    )
                # Including every threshold prevents mixtures from selecting a
                # discrete loss larger than the actual incident energy.
                if physical_threshold < self.energies[-1] and not np.any(
                    self.energies.astype(np.float32) == np.float32(physical_threshold)
                ):
                    raise ValueError("Missing excitation threshold knot")

    def write(self, path):
        self.validate()
        # Little-endian IEEE binary64 payload; text header remains inspectable.
        with Path(path).open("wb") as output:
            output.write(
                (
                    "WARPX_THERMAL_ROTATION_V2\n"
                    f"{self.target} {self.model} {self.maximum_j} {self.reference_temperature:.17g}\n"
                    f"{len(self.energies)} {len(self.transitions)}\n"
                ).encode("ascii")
            )
            for array, dtype in [
                (self.energies, "<f8"),
                (self.transitions, "<i4"),
                (self.rates, "<f8"),
            ]:
                output.write(np.asarray(array, dtype=dtype).tobytes())


def construct(
    target,
    energies,
    edges,
    reduced_elementary,
    inclusive,
    maximum_j,
    reference_temperature=0.0,
    model="elastic_dcs",
):
    """Construct integral Eq. 11.24 rates and their detailed-balance reverse rates.

    reduced_elementary maps rank to A(E)=sigma_0,rank*p_in/p_out.
    inclusive(E) returns the integral cross section at reference_temperature.
    Both have a length-one trailing axis for the reference quadrature helpers.
    Detailed balance uses finite-mass COM momenta at the same invariant energy.
    No rotational angular distribution is generated.
    """
    states = np.arange(maximum_j + 1)
    states = states[weights(target, states) > 0]
    transitions = []
    channels = []
    for initial in states:
        for final in states[states > initial]:
            coefficients = {
                rank: cg_squared(int(initial), rank, int(final))
                for rank in reduced_elementary
            }
            if not any(coefficients.values()):
                continue
            loss = ROTATION[target] * (final * (final + 1) - initial * (initial + 1))
            amplitude = sum(
                coef * reduced_elementary[rank](energies)
                for rank, coef in coefficients.items()
                if coef
            )
            mass = MASSES[target] * C**2 / QE + ROTATION[target] * initial * (
                initial + 1
            )
            upward = phase_rate(energies, loss, mass)[:, None] * amplitude
            shifted = (mass + loss) / mass * energies + threshold(
                target, initial, final
            )
            reverse_amplitude = sum(
                coef * reduced_elementary[rank](shifted)
                for rank, coef in coefficients.items()
                if coef
            )
            downward = (
                weights(target, initial)
                / weights(target, final)
                * mass
                / (mass + loss)
                * (C * np.sqrt(shifted * (shifted + 2 * REST)) / (energies + REST))[
                    :, None
                ]
                * reverse_amplitude
            )
            transitions.extend([(int(initial), int(final)), (int(final), int(initial))])
            channels.extend([upward, downward])
    rates = np.stack([np.zeros_like(channels[0]), *channels], axis=-1)
    weighted = thermal_rates(
        target, reference_temperature, maximum_j, transitions, rates
    )
    rates[:, :, 0] = speed(energies)[:, None] * inclusive(energies) - weighted.sum(
        axis=-1
    )
    return Bundle(
        target,
        model,
        maximum_j,
        reference_temperature,
        energies,
        edges,
        transitions,
        rates,
    )


def energy_grid(target, maximum_j, ranks, maximum, count=512):
    thresholds = [
        threshold(target, i, f)
        for i in range(maximum_j + 1)
        if weights(target, i) > 0
        for f in range(i + 2, min(maximum_j, i + max(ranks)) + 1, 2)
    ]
    values = np.unique(np.r_[0, np.geomspace(1e-9, maximum, count), thresholds])
    values = values[values <= maximum]
    # Different initial-state rest masses can split otherwise equal rotational
    # thresholds by much less than one float32 ulp. Keep the upper threshold
    # in each representable group, never an averaged energy-loss label.
    _, group = np.unique(values.astype(np.float32), return_inverse=True)
    merged = np.zeros(group.max() + 1)
    np.maximum.at(merged, group, values)
    return merged


def write_audit(path, report):
    Path(path).write_text(json.dumps(report, indent=2) + "\n")


def refine_grid(energies, moments, tolerance=0.001, relative_floor=1e-4):
    """Refine rates and absolute first/second transfer moments for linear mixtures.

    ``moments(E)`` returns columns for all validation temperatures. The floor
    bounds insignificant threshold tails relative to the largest moment in
    each column. Intervals narrower than four float32 ulps cannot be refined
    portably and belong to the separately checked threshold uncertainty band.
    """
    grid = np.asarray(energies)
    scale = np.maximum(moments(grid).max(axis=0), np.finfo(float).tiny)
    fractions = np.array([0.25, 0.5, 0.75])
    for _ in range(24):
        values = moments(grid)
        probes = grid[:-1, None] + np.diff(grid)[:, None] * fractions
        exact = moments(probes.ravel()).reshape(len(grid) - 1, 3, -1)
        linear = (
            values[:-1, None, :] * (1 - fractions[None, :, None])
            + values[1:, None, :] * fractions[None, :, None]
        )
        error = np.abs(linear - exact) > tolerance * np.maximum(
            np.abs(exact), relative_floor * scale
        )
        split = error.any(axis=(1, 2)) & (
            np.diff(grid) > 4 * np.finfo(np.float32).eps * np.maximum(grid[:-1], 1e-30)
        )
        # The first interval uses a sqrt(E) mixture to retain the threshold
        # law of finite elastic cross sections as E tends to zero.
        split[0] = False
        if not split.any():
            return grid
        grid = np.sort(np.r_[grid, probes[split, 1]])
    raise ValueError("Rotational rate/moment grid did not converge")
