# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""IAA integral models and adaptive tabulation for WarpX's linear reader.

Equations follow Schmalzried Chapter 11 and elmolcs Excitation.genA/genF.
All fitted parameters and original knots are supplied by physical source records.
"""

import numpy as np

REST = 9.1093837015e-31 * 299792458**2 / 1.602176634e-19
RYD = 13.605693122944
BOHR = 5.29177210903e-11
AREA = 4 * np.pi * BOHR**2
MAX_ENERGY = 1e9
TOLERANCE = 5e-4
PEAK_FLOOR = 1e-6
SR_JOIN = (9.955581, 10.213826)


def analytic(energy, fit):
    energy = np.asarray(energy, dtype=float)
    c = list(fit["coefficients"])
    threshold = fit["threshold_eV"]
    order = fit["asymptotic_order"]
    eq = fit["equation"]
    if order == 0 or "a" in eq:
        if len(c) == 5:
            threshold = c[4]
        a, b, d, power = c[:4]
        prefactor = AREA * a * RYD / (threshold * power)
        if len(c) == 6:
            x = np.maximum(energy / threshold - 1, 0)
            value = (
                prefactor
                * np.expm1(c[4] * np.log1p(x)) ** power
                / (x + b) ** (c[4] * power + 1)
                * (
                    np.log(4 * threshold / RYD * d * (x + 1))
                    + c[5] * RYD / threshold / (x + b)
                )
            )
        else:
            tau = energy / REST
            beta2 = tau * (tau + 2) / (tau + 1) ** 2
            x = np.maximum(0.5 * REST * beta2 / threshold - 1, 0)
            argument = x * (4 * threshold / RYD * d) * (1 + tau) ** 2
            value = prefactor / (x + b) * np.log1p(argument**power)
    else:
        normalization = AREA * (RYD / threshold) ** order
        if eq == "f0m":
            threshold = c.pop(2)
        x = np.maximum((energy - threshold) / threshold, 0)
        if len(c) in (4, 6) and c[2] == 1:
            c.pop(2)
        if len(c) == 3:
            a, b, p = c
            value = normalization * a * x**p / (x + b) ** (p + order)
        elif len(c) == 4:
            a, b, d, p = c
            numerator = np.expm1(d * np.log1p(x)) ** p
            denominator = (
                (x / b) ** (p * d + order) + 1
                if "l" in eq
                else (x + b) ** (p * d + order)
            )
            value = normalization * a * numerator / denominator
        elif len(c) == 5:
            a, b, d, p, f = c
            value = (
                normalization
                * a
                * f**d
                * x**p
                * (1 + (x / f) ** d)
                / (x + b) ** (p + d + order)
            )
        elif len(c) == 6:
            a, b, d, p, f, q = c
            value = (
                normalization
                * a
                * f**q
                * np.expm1(d * np.log1p(x)) ** p
                * (1 + (x / f) ** q)
                / (x + b) ** (q + p * d + order)
            )
        elif len(c) == 7:
            strength, phase, center, width, a, b, p = c
            z = 2 * (energy - center) / width
            resonance = 1 + strength * (np.cos(2 * phase) - z * np.sin(2 * phase)) / (
                1 + z * z
            )
            value = normalization * a * x**p / (x + b) ** (p + order) * resonance
        else:
            raise ValueError(f"Unsupported excitation equation: {fit}")
    return np.where(energy > threshold, value, 0.0)


class Rate:
    def __init__(self, record):
        self.record = record
        self.energy, self.sigma = np.asarray(record["table"], dtype=float).T
        self.maximum = self.energy[-1]
        self.minimum = min(self.energy[0], record["loss_eV"])
        self.knots = np.unique(np.r_[self.minimum, self.energy])
        if "partition" in record:
            self.minimum = min(record["loss_eV"], record["fit"]["threshold_eV"])
            self.knots = np.unique(np.r_[self.knots, self.minimum, SR_JOIN])
            # Both pieces describe one total; the repaired onset extends to
            # the end of the transition rather than stopping at an artificial zero.
            self.maximum = max(self.maximum, SR_JOIN[1])
        self.tail = None
        if record["model"] == "analytic":
            self.continuation = {"model": "source_analytic", "maximum_eV": MAX_ENERGY}
        elif self.sigma[-1] == 0:
            self.continuation = {"model": "zero", "join_eV": self.maximum}
        else:
            start = len(self.energy) - 1
            while start > 0 and self.sigma[start - 1] > self.sigma[start] > 0:
                start -= 1
            start = max(start, len(self.energy) - 5)
            if len(self.energy) - start < 3:
                raise ValueError(f"No supported declining tail for {record['id']}")
            self.tail = float(
                np.polyfit(np.log(self.energy[start:]), np.log(self.sigma[start:]), 1)[
                    0
                ]
            )
            if record["id"] == "complement_14":
                self.tail = -1.0
            if self.tail >= 0:
                raise ValueError(
                    f"Increasing production continuation for {record['id']}"
                )
            self.continuation = {
                "model": "endpoint_power_law",
                "exponent": self.tail,
                "join_eV": self.maximum,
                "fit_range_eV": [float(self.energy[start]), self.maximum],
                "assumption": "dipole-forbidden asymptote"
                if record["id"] == "complement_14"
                else "terminal decreasing source segment",
            }
        probes = np.unique(
            np.r_[self.knots, np.geomspace(max(self.minimum, 1e-5), MAX_ENERGY, 2000)]
        )
        self.peak = float(np.max(self(probes)))

    def __call__(self, energy):
        energy = np.asarray(energy, dtype=float)
        r = self.record
        if r["model"] == "analytic":
            value = analytic(energy, r["fit"])
            if "partition" in r:
                weight = np.clip(
                    (energy - SR_JOIN[0]) / (SR_JOIN[1] - SR_JOIN[0]), 0, 1
                )
                value *= 1 - weight if r["partition"] == "onset" else weight
        else:
            i = np.clip(
                np.searchsorted(self.energy, energy, side="right") - 1,
                0,
                len(self.energy) - 2,
            )
            e0, e1 = self.energy[i], self.energy[i + 1]
            s0, s1 = self.sigma[i], self.sigma[i + 1]
            f = np.clip((energy - e0) / (e1 - e0), 0, 1)
            value = (1 - f) * s0 + f * s1
            positive = (e0 > 0) & (s0 > 0) & (s1 > 0)
            with np.errstate(divide="ignore", invalid="ignore"):
                logf = np.clip(
                    np.log(np.maximum(energy, 1e-300) / e0) / np.log(e1 / e0), 0, 1
                )
                logarithmic = np.exp((1 - logf) * np.log(s0) + logf * np.log(s1))
            value = np.where(positive, logarithmic, value)
            value = np.where(energy == e0, s0, np.where(energy == e1, s1, value))
            if self.tail is not None:
                value = np.where(
                    energy > self.maximum,
                    self.sigma[-1]
                    * (np.maximum(energy, self.maximum) / self.maximum) ** self.tail,
                    value,
                )
            else:
                value = np.where(energy > self.maximum, 0.0, value)
            value = np.where(energy < self.energy[0], 0.0, value)
        value = np.where(energy < r["loss_eV"], 0.0, value)
        if not np.all(np.isfinite(value)) or np.any(value < 0):
            raise ValueError(
                f"Nonphysical source evaluation for {r['target']}/{r['id']}"
            )
        return value


def master_grid(rates, extra_knots=()):
    """Refine a common grid against every channel and transfer moments."""
    knots = np.unique(
        np.concatenate(
            [
                np.array([0.0, MAX_ENERGY]),
                np.geomspace(1e-5, MAX_ENERGY, 400),
                *(r.knots for r in rates),
                np.asarray(extra_knots),
            ]
        )
    )
    fractions = np.array([0.2113248654, 0.5, 0.7886751346])
    for _ in range(64):
        probes = knots[:-1, None] + np.diff(knots)[:, None] * fractions
        actual_fractions = (probes - knots[:-1, None]) / np.diff(knots)[:, None]
        split = np.zeros(len(knots) - 1, dtype=bool)
        for rate in rates:
            values = rate(knots)
            exact = rate(probes)
            interpolated = (1 - actual_fractions) * values[
                :-1, None
            ] + actual_fractions * values[1:, None]
            split |= np.any(
                np.abs(interpolated - exact)
                > 0.7 * TOLERANCE * np.maximum(exact, PEAK_FLOOR * rate.peak),
                axis=1,
            )
        if not split.any():
            return knots
        if np.any(split & ((probes[:, 1] == knots[:-1]) | (probes[:, 1] == knots[1:]))):
            bad = np.flatnonzero(
                split & ((probes[:, 1] == knots[:-1]) | (probes[:, 1] == knots[1:]))
            )[0]
            interval = knots[bad : bad + 2]
            details = [
                (r.record["id"], r(interval).tolist(), r(np.mean(interval)).tolist())
                for r in rates
                if np.ptp(r(interval)) > TOLERANCE * PEAK_FLOOR * r.peak
            ]
            raise ValueError(
                f"Unresolved source discontinuity at {interval}: {details}"
            )
        knots = np.sort(np.r_[knots, probes[split, 1]])
        if len(knots) > 100000:
            raise ValueError("Integral grid exceeds the production size budget")
    raise ValueError("Integral interpolation did not converge")


def write_rate(path, rate, grid, source=False):
    r = rate.record
    maximum = rate.maximum if source else MAX_ENERGY
    energy = np.unique(
        np.r_[rate.minimum, grid[(grid > rate.minimum) & (grid < maximum)], maximum]
    )
    # Each file uses a subset of the common grid. Preserve its source knots,
    # but omit knots introduced solely to resolve another collision channel.
    values = rate(energy)
    anchors = np.unique(
        np.r_[
            0,
            len(energy) - 1,
            np.searchsorted(
                energy,
                rate.knots[(rate.knots >= energy[0]) & (rate.knots <= energy[-1])],
            ),
        ]
    )
    keep = set(anchors.tolist())
    pending = list(zip(anchors[:-1], anchors[1:]))
    while pending:
        lo, hi = pending.pop()
        if hi - lo < 2:
            continue
        f = (energy[lo + 1 : hi] - energy[lo]) / (energy[hi] - energy[lo])
        approximate = (1 - f) * values[lo] + f * values[hi]
        errors = np.abs(approximate - values[lo + 1 : hi]) / np.maximum(
            values[lo + 1 : hi], PEAK_FLOOR * rate.peak
        )
        worst = int(np.argmax(errors)) + lo + 1
        probes = energy[lo] + np.array([0.2113248654, 0.5, 0.7886751346]) * (
            energy[hi] - energy[lo]
        )
        exact = rate(probes)
        estimate = np.interp(probes, energy[[lo, hi]], values[[lo, hi]])
        if errors[worst - lo - 1] > 0.25 * TOLERANCE or np.any(
            np.abs(estimate - exact)
            > 0.6 * TOLERANCE * np.maximum(exact, PEAK_FLOOR * rate.peak)
        ):
            keep.add(worst)
            pending.extend([(lo, worst), (worst, hi)])
    energy = energy[sorted(keep)]
    header = (
        f"IAA {r['target']}: {r['original_process']}\n"
        f"source_file = {r['source_file']}\nsource_sha256 = {r['source_sha256']}\n"
        f"configured_loss_eV = {r['loss_eV']:.17g}\nunits = {r['units']}\n"
        f"source_range_eV = {r['source_range_eV']}\n"
        f"continuation = {'none (source-bounded)' if source else rate.continuation}\n"
        f"native_energy_grid = resample\nenergy_min_eV = 0\nenergy_max_eV = {maximum:.17g}\noutside_energy_range = error\n"
        "Energy (eV), cross section; piecewise-linear runtime interpolation"
    )
    np.savetxt(
        path, np.column_stack((energy, rate(energy))), fmt="%.17e", header=header
    )
    return len(energy)
