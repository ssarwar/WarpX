# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Offline source reference: source-constrained reciprocal rigid rotors.

No production files are emitted. Rates have units m3/s/sr. The source primitive
X_b(E,theta)=p(E)*A_b(E,theta) is shared between forward and reverse kernels.
The cold unchanged background is explicit. State sums, special functions,
normalization, interpolation checks, and sampling work belong offline.
"""

import json

import numpy as np
from reference_paths import OUTPUT
from scipy.sparse import csr_matrix
from scipy.special import gammaln, spherical_jn

HERE = OUTPUT
REFERENCE_VERSION = 2
from elmolcs import reader
from elmolcs_source import Source
from rotation_reference import (
    A0,
    HARTREE,
    REST,
    ROTATION,
    C,
    converged_j,
    populations,
    weights,
)
from source_constraints import GOTE, JUNG


def smooth(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


def cg_array(initial, final, ranks):
    i = np.asarray(initial)[:, None]
    f = np.asarray(final)[:, None]
    L = np.asarray(ranks)[None, :]
    total = i + f + L
    half = total // 2
    ok = (abs(i - L) <= f) & (f <= i + L) & (total % 2 == 0)
    di = half - i
    df = half - f
    dl = half - L
    result = np.zeros_like(total, dtype=float)
    with np.errstate(invalid="ignore", over="ignore"):
        value = (2 * f + 1) * np.exp(
            2
            * (gammaln(half + 1) - gammaln(di + 1) - gammaln(df + 1) - gammaln(dl + 1))
            + gammaln(2 * di + 1)
            + gammaln(2 * df + 1)
            + gammaln(2 * dl + 1)
            - gammaln(total + 2)
        )
    result[ok] = value[ok]
    return result


class Elastic:
    def __init__(self, target):
        self.s = Source(target)
        self.target = target
        t = reader.readDCS(target, format="grid")[0]["data"]
        self.energy = np.asarray(t.columns, float)
        self.d = np.asarray(t).T
        self.y = np.sin(np.asarray(t.index, float) * np.pi / 360)
        dy = np.diff(self.y)
        self.omega = np.zeros(len(self.y))
        self.omega[:-1] += 8 * np.pi * dy * (2 * self.y[:-1] + self.y[1:]) / 6
        self.omega[1:] += 8 * np.pi * dy * (self.y[:-1] + 2 * self.y[1:]) / 6

    def shape(self, E, y):
        if E >= 8000:
            k2 = E * (E + 2 * REST) / (7.2973525693e-3 * REST) ** 2
            eta = 1 / (4 * {"N2": 0.6052, "O2": 0.5677}[self.target] ** 2 * k2)
            sr = eta * (1 + eta) / (4 * np.pi * (eta + y * y) ** 2)
            if E >= 1e4:
                return sr
        E = max(E, self.energy[0])
        i = np.clip(np.searchsorted(self.energy, E) - 1, 0, len(self.energy) - 2)
        t = np.clip(
            np.log(E / self.energy[i]) / np.log(self.energy[i + 1] / self.energy[i]),
            0,
            1,
        )
        d = (1 - t) * self.d[i] + t * self.d[i + 1]
        native = np.interp(y, self.y, d) / (d @ self.omega)
        if E > 8000:
            fraction = smooth(np.log(E / 8000) / np.log(10000 / 8000))
            return (1 - fraction) * native + fraction * sr
        return native

    def cdf(self, E, y):
        y = np.asarray(y)
        if E >= 8000:
            k2 = E * (E + 2 * REST) / (7.2973525693e-3 * REST) ** 2
            eta = 1 / (4 * {"N2": 0.6052, "O2": 0.5677}[self.target] ** 2 * k2)
            sr = (1 + eta) * y * y / (eta + y * y)
            if E >= 1e4:
                return sr
        E = max(E, self.energy[0])
        i = np.clip(np.searchsorted(self.energy, E) - 1, 0, len(self.energy) - 2)
        t = np.clip(
            np.log(E / self.energy[i]) / np.log(self.energy[i + 1] / self.energy[i]),
            0,
            1,
        )
        d = (1 - t) * self.d[i] + t * self.d[i + 1]
        dy = np.diff(self.y)
        slope = np.diff(d) / dy
        y0 = self.y[:-1]
        mass = (
            8
            * np.pi
            * (y0 * d[:-1] * dy + (d[:-1] + y0 * slope) * dy**2 / 2 + slope * dy**3 / 3)
        )
        prefix = np.r_[0, np.cumsum(mass)]
        j = np.clip(np.searchsorted(self.y, y) - 1, 0, len(dy) - 1)
        h = y - self.y[j]
        u = (
            prefix[j]
            + 8
            * np.pi
            * (
                self.y[j] * d[j] * h
                + (d[j] + self.y[j] * slope[j]) * h * h / 2
                + slope[j] * h**3 / 3
            )
        ) / prefix[-1]
        if E > 8000:
            w = smooth(np.log(E / 8000) / np.log(10000 / 8000))
            u = (1 - w) * u + w * sr
        return np.clip(u, 0, 1)

    def dcs(self, E, y):
        return self.s.elastic(E) * self.shape(E, y)


class Priors:
    def __init__(self, target, ranks, initial, final, y, measure):
        self.target = target
        self.elastic = Elastic(target)
        self.s = self.elastic.s
        self.ranks = ranks
        self.initial = initial
        self.final = final
        self.y = y
        self.measure = measure
        self.mu = 1 - 2 * y * y
        self.theta = np.rad2deg(2 * np.arcsin(y))
        self.R = {"N2": 2.068, "O2": 2.281}[target]
        self.nr = len(ranks)
        self.special_initial = (
            np.arange(1, int(max(initial)) + 1, 2)
            if target == "O2"
            else np.array([], dtype=int)
        )
        self.special_delta = ROTATION["O2"] * (4 * self.special_initial + 6)
        self.ns = len(self.special_initial)
        if target == "N2":
            data = JUNG
            self.jung = []
            for energy in [2.22, 2.47]:
                P = np.array(data[str(energy)])
                a = np.arange(15, 106, 15)
                m = np.cos(np.deg2rad(a))
                h_a = np.array(
                    [
                        5 * (3 * m * m - 1) ** 2 / (16 * np.pi),
                        5 * (9 * m**4 - 9 * m * m + 4) / (56 * np.pi),
                        5 * (m * m + 3) ** 2 / (224 * np.pi),
                    ]
                )
                m = self.mu
                h = np.array(
                    [
                        5 * (3 * m * m - 1) ** 2 / (16 * np.pi),
                        5 * (9 * m**4 - 9 * m * m + 4) / (56 * np.pi),
                        5 * (m * m + 3) ** 2 / (224 * np.pi),
                    ]
                )
                row = []
                for r in range(3):
                    v = np.interp(self.theta, a, P[r])
                    v[self.theta < 15] *= h[r, self.theta < 15] / h_a[r, 0]
                    v[self.theta > 105] *= h[r, self.theta > 105] / h_a[r, -1]
                    row.append(v)
                row = np.array(row)
                row /= row.sum(axis=0)
                mask = smooth((self.theta - 105) / 15) + smooth((15 - self.theta) / 15)
                for r, scale in [(1, 0.52823198), (2, 0.84013224)]:
                    removed = row[r] * (1 - scale) * mask
                    row[r] -= removed
                    row[0] += removed
                amp = np.r_[
                    self.s.unchanged(energy),
                    self.s.elementary[2].reduced(energy),
                    self.s.elementary[4].reduced(energy),
                ]
                self.jung.append(row * self.elastic.dcs(energy, y) / amp[:, None])
            self.gote_values = []
            self.gote_valid = []
            for e, rows in GOTE.items():
                rows = np.array(rows[:-1])
                lo = np.maximum(rows, 0).sum(axis=0)
                hi = lo + (rows < 0).sum(axis=0)
                valid = (lo <= 100.3) & (hi >= 99.7)
                self.gote_valid.append(valid)
                # An uncertainty fit: discard internally inconsistent rows, fill censored
                # terms with their interval midpoint, and interpolate from valid neighbors.
                # Printed values and rejected rows remain in the separate source audit.
                values = np.where(rows < 0, 0.5, rows) / 100
                if len(values) == 5:
                    values = np.vstack((values, np.zeros(16)))
                values /= values.sum(axis=0)
                for rank_index in range(6):
                    values[rank_index] = np.interp(
                        np.arange(16), np.flatnonzero(valid), values[rank_index, valid]
                    )
                self.gote_values.append(values)
            self.gote_values = np.array(self.gote_values)
            self.ge = np.array(list(GOTE))
        else:
            self.o2_cache = {}

    def bare(self, E, y=None):
        if y is None:
            y = self.y
        p = np.sqrt(E * (E + 2 * REST))
        z = self.R * p / (7.2973525693e-3 * REST) * y
        v = (2 * self.ranks[:, None] + 1) * spherical_jn(
            self.ranks[:, None], z[None, :]
        ) ** 2
        den = v.sum(axis=0)
        assert (den > 0).all()
        return v / den

    def gote(self, E, y=None):
        if y is None:
            y = self.y
        theta = np.rad2deg(2 * np.arcsin(y))
        i = np.clip(np.searchsorted(self.ge, E) - 1, 0, len(self.ge) - 2)
        t = np.clip(np.log(E / self.ge[i]) / np.log(self.ge[i + 1] / self.ge[i]), 0, 1)
        v = (1 - t) * self.gote_values[i] + t * self.gote_values[i + 1]
        p = np.zeros((self.nr, len(y)))
        p[:6] = np.array([np.interp(theta, np.arange(10, 161, 10), r) for r in v])
        bare = self.bare(E, y)
        low = theta < 10
        high = theta > 160
        a = smooth(theta[low] / 10)
        p[:, low] = (1 - a) * bare[:, low] + a * p[:, low]
        a = smooth((theta[high] - 160) / 20)
        p[:, high] = (1 - a) * p[:, high] + a * bare[:, high]
        return p / p.sum(axis=0)

    def o2_ranks(self, E):
        if not hasattr(self, "o2_completion"):
            from o2_completion import Completion

            self.o2_completion = Completion(
                self.ranks, self.y, self.measure, self.elastic
            )
        return self.o2_completion(E)

    def amplitude(self, E):
        arr = np.zeros((self.nr + self.ns, len(self.y)))
        inclusive = self.elastic.dcs(E, self.y)
        if E >= 1000:
            arr[: self.nr] = self.bare(E) * inclusive
            return arr
        if self.target == "N2":
            if E <= 10:
                low = np.zeros((self.nr, len(self.y)))
                low[0] = inclusive
                for idx, L in enumerate([2, 4, 6], 1):
                    low[idx] = self.s.elementary[L].reduced(E) / (4 * np.pi)
                t = np.clip((E - 2.22) / 0.25, 0, 1)
                h = (1 - t) * self.jung[0] + t * self.jung[1]
                res = np.zeros_like(low)
                amp = np.r_[
                    self.s.unchanged(E),
                    self.s.elementary[2].reduced(E),
                    self.s.elementary[4].reduced(E),
                ]
                res[:3] = amp[:, None] * h
                res[3] = low[3]
                w = smooth((E - 1) / 0.25)
                arr[: self.nr] = (1 - w) * low + w * res
                if E > 4:
                    w = smooth(np.log(E / 4) / np.log(10 / 4))
                    g = self.gote(E) * inclusive
                    arr[: self.nr] = (1 - w) * arr[: self.nr] + w * g
            elif E <= 200:
                arr[: self.nr] = self.gote(E) * inclusive
            else:
                k = np.sqrt(E * (E + 2 * REST))
                k200 = np.sqrt(200 * (200 + 2 * REST))
                q_ratio = self.y * k / k200
                mapped = np.minimum(q_ratio, 1.0)
                data = self.gote(200, mapped)
                bare = self.bare(E)
                # Beyond the measured q range, use the two-centre continuation.
                edge = smooth((q_ratio - 1) / 0.05)
                data = (1 - edge) * data + edge * bare
                w = smooth(np.log(E / 200) / np.log(1000 / 200))
                arr[: self.nr] = ((1 - w) * data + w * bare) * inclusive
        else:
            if E < 20:
                N = self.special_initial
                delta = self.special_delta
                ki = np.sqrt(2 * E / HARTREE)
                kf = np.sqrt(2 * np.maximum(E - delta, 0) / HARTREE)
                q = np.sqrt(
                    np.maximum(
                        ki * ki + kf[:, None] ** 2 - 2 * ki * kf[:, None] * self.mu, 0
                    )
                )
                reduced = (
                    (6 / 5 * (N + 1) * (N + 2) / (2 * N + 1) / (2 * N + 3))[:, None]
                    * (-0.29 / 3 + np.pi * 4.93 * q / 32) ** 2
                    * A0**2
                )
                w = smooth(np.log(max(E, 1) / 1) / np.log(20))
                arr[self.nr :] = (1 - w) * reduced
                high = self.o2_ranks(max(E, 1.0))
                arr[: self.nr] = w * high
                arr[0] = inclusive
            elif E <= 200:
                arr[: self.nr] = self.o2_ranks(E)
            else:
                # Continue the positive source-constrained 200 eV distribution at equal q.
                k = np.sqrt(E * (E + 2 * REST))
                k200 = np.sqrt(200 * (200 + 2 * REST))
                q_ratio = self.y * k / k200
                mapped = np.minimum(q_ratio, 1.0)
                P200 = self.o2_anchor
                data = np.array([np.interp(mapped, self.y, row) for row in P200])
                bare = self.bare(E)
                edge = smooth((q_ratio - 1) / 0.05)
                data = (1 - edge) * data + edge * bare
                w = smooth(np.log(E / 200) / np.log(1000 / 200))
                arr[: self.nr] = ((1 - w) * data + w * bare) * inclusive
        return arr


class Hybrid:
    def __init__(
        self,
        target,
        T=300,
        max_energy=1000,
        resolution=1,
        rank_max=40,
        angle_order=3,
        angular_override=None,
        energy_override=None,
        angular_resolution=None,
    ):
        if target not in ("N2", "O2") or max_energy <= 0:
            raise ValueError("Unsupported target or energy range")
        self.target = target
        self.T = T
        self.B = ROTATION[target]
        self.ranks = np.arange(0, rank_max + 1, 2)
        self.nr = len(self.ranks)
        self.jmax = converged_j(target, T)
        self.bound = int(
            np.floor(
                (np.sqrt(1 + 4 * {"N2": 9.759, "O2": 5.116}[target] / self.B) - 1) / 2
            )
        )
        if target == "O2":
            self.bound -= self.bound % 2 == 0
        if self.jmax > self.bound:
            raise ValueError("Rotational bath exceeds the bounded-rotor state space")
        self.pop = populations(target, T, self.jmax)[0]
        states = np.flatnonzero(weights(target, np.arange(self.jmax + 1)))
        pairs = [
            (i, f)
            for i in states
            for f in range(i + 2, min(i + rank_max, self.bound) + 1, 2)
        ]
        self.initial, self.final = np.array(pairs).T
        self.loss = self.B * (
            self.final * (self.final + 1) - self.initial * (self.initial + 1)
        )
        coeff = cg_array(self.initial, self.final, self.ranks)
        self.special_initial = (
            np.arange(1, self.jmax + 1, 2)
            if target == "O2"
            else np.array([], dtype=int)
        )
        self.ns = len(self.special_initial)
        self.nb = self.nr + self.ns
        if self.ns:
            extra = np.zeros((len(pairs), self.ns))
            for rank_index, j in enumerate(self.special_initial):
                extra[(self.initial == j) & (self.final == j + 2), rank_index] = 1
            coeff = np.column_stack((coeff, extra))
        # Group equal canonical losses to make future-energy evaluation sparse.
        off = np.rint(self.loss / (2 * self.B)).astype(int)
        self.offset = np.unique(off)
        self.delta = self.offset * 2 * self.B
        group = np.searchsorted(self.offset, off)
        popi = self.pop[self.initial]
        popf = np.where(
            self.final <= self.jmax, self.pop[np.minimum(self.final, self.jmax)], 0
        )
        ratio = weights(target, self.initial) / weights(target, self.final)
        up = np.zeros((len(self.delta), self.nb))
        down = np.zeros_like(up)
        np.add.at(up, group, popi[:, None] * coeff)
        np.add.at(down, group, (popf * ratio)[:, None] * coeff)
        self.up = csr_matrix(up)
        keep = np.flatnonzero(down.sum(axis=1) > 1e-25)
        self.down = down[keep]
        self.d_delta = self.delta[keep]
        self.guard = 2.0 + (float(self.d_delta.max()) if len(self.d_delta) else 0.0)
        diag = cg_array(states, states, self.ranks)
        self.diag = np.r_[self.pop[states] @ diag, np.zeros(self.ns)]
        ar = resolution if angular_resolution is None else angular_resolution
        edges = np.unique(
            np.r_[0, np.geomspace(1e-8, 1, 26 * ar), np.linspace(0, 1, 31 * ar)]
        )
        x, w = np.polynomial.legendre.leggauss(angle_order)
        y = (edges[:-1, None] + np.diff(edges)[:, None] * (x + 1) / 2).ravel()
        meas = (
            np.diff(edges)[:, None] * w / 2 * 8 * np.pi * y.reshape(-1, angle_order)
        ).ravel()
        extra = np.sin(
            np.deg2rad(np.r_[0, np.arange(10, 181, 10), 15, 45, 75, 105]) / 2
        )
        self.y = np.r_[y, extra]
        self.measure = np.r_[meas, np.zeros(len(extra))]
        self.extra = extra
        if angular_override is not None:
            self.y, self.measure = angular_override
        order = np.argsort(self.y)
        self.y = self.y[order]
        self.measure = self.measure[order]
        self.prior = Priors(
            target, self.ranks, self.initial, self.final, self.y, self.measure
        )
        self.elastic = self.prior.elastic
        if target == "O2" and (energy_override is None or energy_override[0] < 1000):
            a = self.prior.o2_ranks(200)
            inclusive = self.elastic.dcs(200, self.y)
            # For the high-energy angular continuation use sudden conservation exactly.
            a = self.prior.o2_ranks(200)
            a[0] = inclusive - a[1:].sum(axis=0)
            if (a[0] < 0).any():
                raise ValueError("O2 200 eV rank background is negative")
            self.prior.o2_anchor = a / inclusive
        thresholds = np.unique(self.loss[(self.final - self.initial) <= 6])
        thresholds = thresholds[thresholds < 1]
        s = self.prior.s
        E = np.unique(
            np.r_[
                np.geomspace(1e-9, max_energy + self.guard, 450 * resolution),
                np.linspace(0.001, 0.1, 100 * resolution),
                np.linspace(1, 4, 240 * resolution),
                thresholds,
                s.knots[s.knots <= min(max_energy + self.guard, 1000)],
                self.elastic.energy[self.elastic.energy <= max_energy + self.guard],
                0.001,
                1,
                1.25,
                2.22,
                2.47,
                4,
                10,
                20,
                30,
                40,
                50,
                60,
                70,
                100,
                150,
                200,
                220,
                8000,
                9000,
                10000,
                max_energy,
                max_energy + 1,
                max_energy + self.guard,
            ]
        )
        if energy_override is not None:
            E = np.asarray(energy_override)
        self.energy = E[(E > 0) & (E <= max_energy + self.guard)]
        self.boundary = np.searchsorted(self.energy, max_energy + 1)
        self.p = np.sqrt(self.energy * (self.energy + 2 * REST))
        self.kin = C / (self.energy + REST)
        self.X = np.zeros((len(self.energy), self.nb, len(self.y)))
        self.c = np.ones((len(self.energy), len(self.y)))
        self.unchanged = np.zeros((len(self.energy), len(self.y)))
        self.raw_residual_min = 1.0
        self.cold_anchor = None
        self.fail = []

    def up_coeff(self, E):
        out = np.maximum(E - self.delta, 0)
        return np.asarray(np.sqrt(out * (out + 2 * REST)) @ self.up) / np.sqrt(
            E * (E + 2 * REST)
        )

    def future_weights(self, E, index, power=0):
        future = E + self.d_delta
        lo = np.searchsorted(self.energy, future, side="right") - 1
        hi = lo + 1
        t = (future - self.energy[lo]) / (self.energy[hi] - self.energy[lo])
        nodes = np.unique(np.r_[lo, hi])
        W = np.zeros((len(nodes), self.nb))
        weighted = self.down * self.d_delta[:, None] ** power
        np.add.at(W, np.searchsorted(nodes, lo), (1 - t[:, None]) * weighted)
        np.add.at(W, np.searchsorted(nodes, hi), t[:, None] * weighted)
        same = nodes == index
        implicit = W[same].sum(axis=0)
        return nodes[~same], W[~same], implicit

    def solve(self):
        # Externally supplied asymptotic closure, followed by a two-eV buffer.
        for i in range(self.boundary, len(self.energy)):
            E = self.energy[i]
            A = self.prior.amplitude(E)
            P = self.p[i] * A
            up = self.up_coeff(E)
            native = (
                np.sqrt((E + self.d_delta) * (E + self.d_delta + 2 * REST))
                / self.p[i]
                @ self.down
            )
            raw = self.kin[i] * (self.diag + up + native) @ P
            target = self.kin[i] * self.p[i] * self.elastic.dcs(E, self.y)
            self.c[i] = target / raw
            self.X[i] = P * self.c[i]
        for i in range(self.boundary - 1, -1, -1):
            E = self.energy[i]
            A = self.prior.amplitude(E)
            P = self.p[i] * A
            up = self.up_coeff(E)
            nodes, W, implicit = self.future_weights(E, i)
            future = self.kin[i] * np.einsum("nb,nba->a", W, self.X[nodes])
            selfpart = self.kin[i] * (implicit @ P)
            target = self.kin[i] * self.p[i] * self.elastic.dcs(E, self.y)
            cutoff = 1.0 if self.target == "N2" else 200.0
            if E <= cutoff:
                # Protect the changing reference kernels; solve only the unchanged term.
                P[0] = 0
                up_rate = self.kin[i] * (up @ P)
                diag_rate = self.kin[i] * (self.diag @ P)
                selfpart = self.kin[i] * (implicit @ P)
                residual = target - future - selfpart - up_rate - diag_rate
                if E >= 0.001:
                    self.raw_residual_min = min(
                        self.raw_residual_min, float(np.min(residual / target))
                    )
                    if np.any(residual < 0):
                        self.fail.append((float(E), "unchanged"))
                        raise ValueError(self.fail[-1])
                    P[0] = residual / (self.kin[i] * self.diag[0])
                    if E == 0.001:
                        self.cold_anchor = residual / (
                            self.kin[i] * self.p[i] * self.diag[0]
                        )
                else:
                    if self.cold_anchor is None:
                        raise ValueError("missing cold anchor")
                    b = smooth(np.sqrt(E / 0.001))
                    a = {"N2": 0.44, "O2": 0.3}[self.target]
                    # Rank-zero s-wave background; rotational kernels retain their own limits.
                    P[0] = self.p[i] * (
                        (1 - b) * a * a * A0 * A0 + b * self.cold_anchor
                    )
                self.X[i] = P
            else:
                release = 1.25 if self.target == "N2" else 220.0
                if E < release:
                    protected = P.copy()
                    protected[0] = 0
                    other_forward = self.kin[i] * ((self.diag + up) @ protected)
                    original = self.kin[i] * self.diag[0] * P[0]
                    available = target - future
                    w = smooth((E - cutoff) / (release - cutoff))
                    # The implicit reverse contribution is c*selfpart. Include
                    # that same corrected contribution in the protected
                    # background, so the model does not depend on grid spacing.
                    linear = (
                        (1 - w) * available + w * (other_forward + original) + selfpart
                    )
                    discriminant = ((1 - w) * available - selfpart) ** 2 + (
                        2
                        * w
                        * (other_forward + original)
                        * ((1 - w) * available + selfpart)
                        + (w * (other_forward + original)) ** 2
                    )
                    c = 2 * available / (linear + np.sqrt(discriminant))
                    residual = available - c * selfpart - other_forward
                    if np.any(residual < 0):
                        raise ValueError(("negative transition background", E))
                    P[0] = (1 - w) * residual / (self.kin[i] * self.diag[0]) + w * P[0]
                else:
                    raw = self.kin[i] * ((self.diag + up) @ P) + selfpart
                    c = (target - future) / raw
                if np.any(~np.isfinite(c)) or np.any(c < 0):
                    self.fail.append((float(E), "common factor"))
                    raise ValueError(self.fail[-1])
                self.c[i] = c
                self.X[i] = P * c
        return self

    def interp(self, values, E):
        E = np.atleast_1d(E)
        i = np.clip(np.searchsorted(self.energy, E) - 1, 0, len(self.energy) - 2)
        t = (E - self.energy[i]) / (self.energy[i + 1] - self.energy[i])
        t = t.reshape((-1,) + (1,) * (values.ndim - 1))
        result = (1 - t) * values[i] + t * values[i + 1]
        below = E < self.energy[0]
        if np.any(below):
            if self.energy[0] > 1e-6:
                raise ValueError("Primitive requested below local reference window")
            ratio = (
                np.sqrt(np.maximum(E[below], 0) * (np.maximum(E[below], 0) + 2 * REST))
                / self.p[0]
            )
            result[below] = values[0] * ratio.reshape((-1,) + (1,) * (values.ndim - 1))
        return result

    def angular_moments(self, E):
        P = self.interp(self.X, E)[0]
        p = np.sqrt(E * (E + 2 * REST))
        kin = C / (E + REST)
        if E == 0:
            down = []
            for power in [0, 1, 2]:
                nodes, W, _ = self.future_weights(E, -1, power)
                down.append(kin * np.einsum("nb,nba->a", W, self.X[nodes]))
            zero = np.zeros(len(self.y))
            return np.array([zero, zero, down[0], zero, down[1], down[2]])
        out = np.maximum(E - self.delta, 0)
        phase = np.sqrt(out * (out + 2 * REST)) / p
        up = []
        down = []
        for power in [0, 1, 2]:
            coeff = np.asarray((phase * self.delta**power) @ self.up)
            up.append(kin * (coeff @ P))
            nodes, W, _ = self.future_weights(E, -1, power)
            down.append(kin * np.einsum("nb,nba->a", W, self.X[nodes]))
        un = kin * (self.diag @ P)
        return np.array([un, up[0], down[0], up[1], down[1], up[2] + down[2]])

    def evaluate(self, E):
        value = self.angular_moments(E)
        m = value @ self.measure
        total = value[:3].sum(axis=0)
        p = np.sqrt(E * (E + 2 * REST))
        kin = C / (E + REST)
        inclusive = kin * p * self.elastic.dcs(E, self.y)
        return {
            "E_eV": E,
            "moments": m.tolist(),
            "total_rate": float(m[:3].sum()),
            "max_inclusive_error": float(np.max(abs(total / inclusive - 1)))
            if E >= 0.001
            else None,
            "integral_inclusive_error": float(
                total @ self.measure / (kin * p * self.prior.s.elastic(E)) - 1
            )
            if E >= 0.001
            else None,
            "unchanged_fraction": float(m[0] / m[:3].sum())
            if m[:3].sum() > 0
            else None,
        }


if __name__ == "__main__":
    import argparse
    import time

    p = argparse.ArgumentParser()
    p.add_argument("--target", default="N2")
    p.add_argument("--resolution", type=int, default=1)
    p.add_argument("--rank", type=int, default=40)
    p.add_argument("--temperature", type=float, default=300)
    p.add_argument("--max-energy", type=float, default=1000)
    p.add_argument("--output", default="hybrid")
    p.add_argument("--no-save", action="store_true")
    a = p.parse_args()
    t = time.perf_counter()
    m = Hybrid(a.target, a.temperature, a.max_energy, a.resolution, a.rank).solve()
    values = [
        m.evaluate(e)
        for e in [
            1e-8,
            0.0005,
            0.001,
            0.003,
            0.01,
            0.025,
            0.1,
            0.5,
            1,
            1.1,
            1.25,
            1.5,
            2.22,
            2.47,
            3,
            4,
            5,
            8,
            10,
            15,
            20,
            30,
            50,
            100,
            150,
            200,
            300,
            500,
            900,
            1000,
            5000,
            9999,
            10000,
            10001,
            1e5,
            1e6,
            2.5e6,
            1e7,
            1e8,
            5e8,
        ]
        if e < a.max_energy
    ]
    report = {
        "target": a.target,
        "temperature_K": a.temperature,
        "energy_nodes": len(m.energy),
        "angle_nodes": len(m.y),
        "rank_max": a.rank,
        "table_bytes": m.X.nbytes,
        "seconds": time.perf_counter() - t,
        "minimum_raw_unchanged_fraction": m.raw_residual_min,
        "common_factor_range": [
            float(m.c[: m.boundary].min()),
            float(m.c[: m.boundary].max()),
        ],
        "evaluations": values,
    }
    (HERE / (a.output + "-" + a.target + ".json")).write_text(
        json.dumps(report, indent=2) + "\n"
    )
    if not a.no_save:
        np.savez_compressed(
            HERE / (a.output + "-" + a.target + ".npz"),
            energy=m.energy,
            y=m.y,
            measure=m.measure,
            X=m.X,
            c=m.c,
        )
    print(json.dumps(report, indent=2), flush=True)
