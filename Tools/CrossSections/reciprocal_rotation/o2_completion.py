# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Spectator-informed, positive completion of the O2 published moment data.

Bhattacharyya & Goswami constrain N=1->3 and the sum of higher transitions,
not a unique elementary-rank spectrum. Higher-rank proportions are supplied
by the two-centre prior. An entropy projection retains the IAA angular
marginal while matching the inferred positive rank integrals and first angular
moments. This is an explicit model closure, not extra published measurements.
"""

import json

import numpy as np
from hybrid_reference import REST, Elastic
from reference_paths import output_path
from rotation_reference import A0, ROTATION, populations
from scipy.optimize import minimize
from scipy.special import logsumexp, spherical_jn
from source_constraints import ENERGY, MOMENT, TOTAL


class Completion:
    def __init__(self, ranks, y, measure, elastic=None):
        self.ranks = ranks
        self.y = y
        self.mu = 1 - 2 * y * y
        self.w = measure
        self.elastic = Elastic("O2") if elastic is None else elastic
        self.cache = {}
        self.diagnostics = []

    def __call__(self, E):
        if E in self.cache:
            return self.cache[E]
        e = np.clip(E, 20, 200)
        k = np.clip(np.searchsorted(ENERGY, e) - 1, 0, len(ENERGY) - 2)
        t = np.log(e / ENERGY[k]) / np.log(ENERGY[k + 1] / ENERGY[k])
        s = np.exp((1 - t) * np.log(TOTAL[k]) + t * np.log(TOTAL[k + 1]))
        mt = np.exp((1 - t) * np.log(MOMENT[k]) + t * np.log(MOMENT[k + 1]))
        F = self.elastic.shape(E, self.y)
        F /= F @ self.w
        z = 2.281 * np.sqrt(e * (e + 2 * REST)) / (7.2973525693e-3 * REST) * self.y
        bare = (2 * self.ranks[:, None] + 1) * spherical_jn(
            self.ranks[:, None], z[None, :]
        ) ** 2
        bare /= bare.sum(axis=0)
        S = bare @ (F * self.w)
        M = bare @ (F * self.w * (1 - self.mu))
        factors = np.ones(len(S))
        factors[:2] = 0
        factors[2] = 5 / 9
        residue = s[2] - s[0] - s[1]
        residue_m = mt[2] - mt[0] - mt[1]
        a = S * residue / (factors @ S)
        am = M * residue_m / (factors @ M)
        a[:2] = 0
        am[:2] = 0
        a[1] = (s[1] - 4 / 9 * a[2]) / (3 / 5)
        am[1] = (mt[1] - 4 / 9 * am[2]) / (3 / 5)
        total = self.elastic.s.elastic(E) / A0**2
        keep = np.flatnonzero(a > total * 1e-13)
        target = a[keep] / total
        means = 1 - am[keep] / a[keep]
        assert np.all((means > -1) & (means < 1)) and target.sum() < 1
        prior = bare[keep].copy()
        prior[keep == 1] = 1.0
        logprior = np.log(np.maximum(prior, 1e-300))
        qw = F * self.w
        goal = np.r_[target, target * means]
        n = len(keep)

        def objective(x):
            logits = logprior + x[:n, None] + x[n:, None] * self.mu
            norm = logsumexp(np.vstack([np.zeros(len(self.y)), logits]), axis=0)
            prob = np.exp(logits - norm)
            value = qw @ norm - goal @ x
            grad = np.r_[prob @ qw, prob @ (qw * self.mu)] - goal
            return value, grad

        x = np.r_[np.log(target / (1 - target.sum())) - np.log(prior @ qw), np.zeros(n)]
        result = minimize(
            objective,
            x,
            jac=True,
            method="BFGS",
            options={"gtol": 1e-11, "maxiter": 600},
        )
        logits = logprior + result.x[:n, None] + result.x[n:, None] * self.mu
        norm = logsumexp(np.vstack([np.zeros(len(self.y)), logits]), axis=0)
        P = np.zeros_like(bare)
        P[keep] = np.exp(logits - norm)
        P[0] = np.exp(-norm)
        achieved = P[keep] @ qw
        achieved_m = P[keep] @ (qw * self.mu)
        error = max(
            np.max(abs(achieved - target)), np.max(abs(achieved_m - target * means))
        )
        assert error < 2e-7, (E, result.message, error)
        value = self.elastic.s.elastic(E) * F * P
        self.cache[E] = value
        self.diagnostics.append(
            {
                "E": E,
                "constraint_absolute_error_fraction": float(error),
                "minimum_unchanged_rank_probability": float(P[0].min()),
                "rank_integrals_a02": a.tolist(),
                "rank_MT_a02": am.tolist(),
            }
        )
        return value


if __name__ == "__main__":
    x, w = np.polynomial.legendre.leggauss(256)
    y = np.sqrt((1 - x) / 2)
    measure = 2 * np.pi * w
    ranks = np.arange(0, 42, 2)
    model = Completion(ranks, y, measure)
    out = []
    pop = populations("O2", 300, 96)[0]
    N = np.arange(len(pop))
    J2 = pop @ (N * (N + 1))
    B = ROTATION["O2"]
    for E in ENERGY:
        a = model(E)
        diag = model.diagnostics[-1]
        S = np.array(diag["rank_integrals_a02"])
        L2 = ranks * (ranks + 1)
        diag["first_rotational_moment_integral_eV_a02"] = float(B * S @ L2)
        diag["second_rotational_moment_integral_eV2_a02"] = float(
            B * B * S @ (L2 * L2 + 2 * J2 * L2)
        )
        out.append(diag)
        print(
            E,
            diag["constraint_absolute_error_fraction"],
            diag["minimum_unchanged_rank_probability"],
            flush=True,
        )
    output_path("o2_completion.json").write_text(json.dumps(out, indent=2) + "\n")
