# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Check the reference energy grid against newly solved intermediate rows.

A point inserted in the advanced-energy system depends only on existing
higher-energy rows and on itself. This gives an independent off-grid check
without mistaking interpolation of the current grid for the source model.
"""

import numpy as np
from hybrid_reference import A0, REST, C, smooth
from reference_paths import OUTPUT

ROOT = OUTPUT


def inserted_moments(m, E):
    p = np.sqrt(E * (E + 2 * REST))
    kin = C / (E + REST)
    P = p * m.prior.amplitude(E)
    up = m.up_coeff(E)
    future = E + m.d_delta
    lo = np.searchsorted(m.energy, future, side="right") - 1
    hi = lo + 1
    same = m.energy[lo] < E
    left = m.energy[lo].copy()
    left[same] = E
    t = (future - left) / (m.energy[hi] - left)
    moments = []
    implicit = []
    for power in [0, 1, 2]:
        W = m.down * m.d_delta[:, None] ** power
        nodes = np.unique(np.r_[hi, lo[~same]])
        grouped = np.zeros((len(nodes), m.nb))
        np.add.at(grouped, np.searchsorted(nodes, hi), W * t[:, None])
        np.add.at(
            grouped, np.searchsorted(nodes, lo[~same]), W[~same] * (1 - t[~same, None])
        )
        known = np.einsum("nb,nba->a", grouped, m.X[nodes])
        moments.append(kin * known)
        implicit.append(((1 - t[same, None]) * W[same]).sum(axis=0))
    inclusive = kin * p * m.elastic.dcs(E, m.y)
    cutoff = 1 if m.target == "N2" else 200
    if E <= cutoff:
        P[0] = 0
        residual = inclusive - moments[0] - kin * ((m.diag + up + implicit[0]) @ P)
        if E >= 0.001:
            assert np.min(residual) >= 0, (m.target, E)
            P[0] = residual / (kin * m.diag[0])
        else:
            w = smooth(np.sqrt(E / 0.001))
            a = {"N2": 0.44, "O2": 0.3}[m.target]
            P[0] = p * ((1 - w) * a * a * A0 * A0 + w * m.cold_anchor)
    else:
        release = 1.25 if m.target == "N2" else 220
        if E < release:
            tmp = P.copy()
            tmp[0] = 0
            other_forward = kin * ((m.diag + up) @ tmp)
            implicit_rate = kin * (implicit[0] @ tmp)
            original = kin * m.diag[0] * P[0]
            available = inclusive - moments[0]
            w = smooth((E - cutoff) / (release - cutoff))
            lo = np.zeros_like(available)
            hi = (
                2
                * available
                / ((1 - w) * available + w * (other_forward + original) + implicit_rate)
            )
            # Independent bisection of the physical normalization equation;
            # do not reuse the production reference's quadratic root.
            for _ in range(56):
                mid = (lo + hi) / 2
                background = (1 - w) * (
                    available - mid * implicit_rate - other_forward
                ) + w * original
                total = mid * (other_forward + background + implicit_rate)
                lo = np.where(total < available, mid, lo)
                hi = np.where(total >= available, mid, hi)
            c = (lo + hi) / 2
            residual = available - c * implicit_rate - other_forward
            assert np.min(residual) >= 0
            P[0] = (1 - w) * residual / (kin * m.diag[0]) + w * P[0]
        else:
            c = (inclusive - moments[0]) / (kin * ((m.diag + up + implicit[0]) @ P))
        assert np.min(c) >= 0
        P *= c
    phase = (
        np.sqrt(np.maximum(E - m.delta, 0) * (np.maximum(E - m.delta, 0) + 2 * REST))
        / p
    )
    U = [kin * (np.asarray((phase * m.delta**n) @ m.up) @ P) for n in [0, 1, 2]]
    D = [moments[n] + kin * (implicit[n] @ P) for n in [0, 1, 2]]
    val = np.array([kin * m.diag @ P, U[0], D[0], U[1], D[1], U[2] + D[2]])
    return val, P
