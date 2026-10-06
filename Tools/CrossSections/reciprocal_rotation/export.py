# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Export fixed-temperature joint distributions; never invoked by WarpX."""

import argparse
import json
from pathlib import Path

import numpy as np
from build_reference import build_reference
from bundle import ALIAS, write_bundle
from high_distribution import high_tables
from hybrid_reference import REFERENCE_VERSION, REST, C, Hybrid, cg_array
from lookup_index import lookup_arrays
from numba import njit
from rotation_reference import MASSES, weights


@njit(cache=True)
def alias_columns(probability):
    count = len(probability)
    scaled = probability * count
    cut = np.ones(count)
    alternate = np.arange(count)
    small = [i for i in range(count) if scaled[i] < 1]
    large = [i for i in range(count) if scaled[i] >= 1]
    while len(small) and len(large):
        low, high = small.pop(), large.pop()
        cut[low] = scaled[low]
        alternate[low] = high
        scaled[high] -= 1 - scaled[low]
        if scaled[high] < 1:
            small.append(high)
        else:
            large.append(high)
    return cut, alternate


def angular_cdf(y, density):
    """Integrate the positive linear differential kernel in sin(theta/2)."""
    dy = np.diff(y)
    slope = np.diff(density) / dy
    y0 = y[:-1]
    mass = (
        8
        * np.pi
        * (
            y0 * density[:-1] * dy
            + (density[:-1] + y0 * slope) * dy**2 / 2
            + slope * dy**3 / 3
        )
    )
    if np.min(mass) < 0 or mass.sum() <= 0:
        raise ValueError("Invalid angular marginal")
    prefix = np.r_[0.0, np.cumsum(mass)]
    return prefix / prefix[-1], prefix[-1]


def compress(u, distribution, changing, features):
    """Refine each interval independently, avoiding repeated global fits."""
    exact = features @ distribution
    joint = exact * changing
    scale = np.maximum(abs(joint), np.max(abs(joint), axis=1)[:, None] * 1e-8)
    scale = np.maximum(scale, 1e-100)
    selected = {0, len(u) - 1}
    stack = [(0, len(u) - 1)]
    while stack:
        lo, hi = stack.pop()
        if hi - lo <= 1:
            continue
        t = (u[lo + 1 : hi] - u[lo]) / (u[hi] - u[lo])
        fit = (1 - t) * distribution[:, lo, None] + t * distribution[:, hi, None]
        tv = 0.5 * abs(fit - distribution[:, lo + 1 : hi]).sum(axis=0)
        error = abs(features @ fit - exact[:, lo + 1 : hi])
        score = np.maximum(
            tv * changing[lo + 1 : hi] / 1e-5,
            (error * changing[lo + 1 : hi] / scale[:, lo + 1 : hi]).max(axis=0) / 2e-4,
        )
        worst = int(score.argmax())
        if score[worst] > 1:
            mid = lo + 1 + worst
            selected.add(mid)
            stack.extend([(lo, mid), (mid, hi)])
    return np.asarray(sorted(selected))


class Exporter:
    def __init__(self, model, sampling):
        self.model = model
        self.sampling = sampling
        self.high = high_tables(model.target, model.T)
        self.palette = np.column_stack((self.high["loss"], self.high["initial_energy"]))
        self.pair_ids = {
            (int(i), int(f)): k
            for k, (i, f) in enumerate(zip(self.high["initial"], self.high["final"]))
        }
        self.coeff = cg_array(model.initial, model.final, model.ranks)
        if model.ns:
            extra = np.zeros((len(model.initial), model.ns))
            for k, j in enumerate(model.special_initial):
                extra[(model.initial == j) & (model.final == j + 2), k] = 1
            self.coeff = np.column_stack((self.coeff, extra))
        self.reverse = model.final <= model.jmax
        self.ii = model.initial[self.reverse]
        self.ff = model.final[self.reverse]
        self.dd = model.loss[self.reverse]
        self.cc = self.coeff[self.reverse]
        self.ids = np.array(
            [self.pair_ids[int(i), int(f)] for i, f in zip(model.initial, model.final)]
            + [self.pair_ids[int(f), int(i)] for i, f in zip(self.ii, self.ff)]
        )
        self.loss = np.r_[model.loss, -self.dd]
        self.features = np.vstack(
            (
                self.loss > 0,
                self.loss < 0,
                np.maximum(self.loss, 0),
                np.maximum(-self.loss, 0),
                self.loss**2,
            )
        )
        self.aliases, self.cdfs, self.outcome_ids = [], [], []
        self.cell_offsets = [0]
        self.angular_offsets, self.conditional_offsets = [0], [0]
        self.angular_u, self.deflection, self.changing = [], [], []
        self.conditional_u, self.conditional_cells = [], []
        self.rates = []
        self.storage_error = 0.0
        self.pruned_probability = 0.0

    def add_cell(self, probability, ids):
        keep = probability > 1e-15
        if not keep.any():
            raise ValueError("Empty conditional distribution")
        dropped = probability[~keep].sum()
        self.pruned_probability = max(self.pruned_probability, float(dropped))
        probability = probability[keep].copy()
        ids = ids[keep]
        probability /= probability.sum()
        if len(probability) > 65535 or np.max(ids) > 65535:
            raise ValueError("Alias indices exceed uint16 capacity")
        cut, alternate = alias_columns(probability)
        entries = np.empty(len(probability), ALIAS)
        entries["cut"], entries["alias"], entries["outcome"] = cut, alternate, ids
        restored = entries["cut"].astype(float) / len(entries)
        np.add.at(
            restored, alternate, (1 - entries["cut"].astype(float)) / len(entries)
        )
        loss = self.palette[ids, 0]
        moments = np.vstack(
            (np.ones(len(loss)), np.maximum(loss, 0), np.maximum(-loss, 0), loss**2)
        )
        expected = moments @ probability
        self.storage_error = max(
            self.storage_error,
            float(
                np.max(abs(moments @ restored - expected) / np.maximum(expected, 1e-30))
            ),
        )
        index = len(self.cell_offsets) - 1
        self.cell_offsets.append(self.cell_offsets[-1] + len(entries))
        if self.sampling == "alias":
            self.aliases.append(entries)
        else:
            cdf = np.cumsum(probability)
            cdf[-1] = 1
            self.cdfs.append(cdf)
            self.outcome_ids.append(ids.astype("<u2"))
        return index

    def row(self, energy):
        m = self.model
        y = m.y
        if energy > 1000:
            self.rates.append(
                float(
                    C
                    * np.sqrt(energy * (energy + 2 * REST))
                    / (energy + REST)
                    * m.elastic.s.elastic(energy)
                )
            )
            if energy <= 10000:
                u = m.elastic.cdf(energy, y)
                keep = np.r_[True, np.diff(u) > 1e-15]
                self.angular_u.extend(u[keep])
                self.deflection.extend(2 * y[keep] ** 2)
                self.changing.extend(np.zeros(np.count_nonzero(keep)))
            self.angular_offsets.append(len(self.angular_u))
            self.conditional_offsets.append(len(self.conditional_u))
            return
        kin = C / (energy + REST)
        X = m.interp(m.X, energy)[0]
        phase = np.zeros(len(m.loss))
        if energy > 0:
            outgoing = np.maximum(energy - m.loss, 0)
            phase = np.sqrt(
                outgoing * (outgoing + 2 * REST) / (energy * (energy + 2 * REST))
            )
        up = (kin * m.pop[m.initial] * phase)[:, None] * (self.coeff @ X)
        down = np.empty((len(self.ii), len(y)))
        for start in range(0, len(self.ii), 32):
            sl = slice(start, start + 32)
            future = m.interp(m.X, energy + self.dd[sl])
            down[sl] = (
                kin
                * m.pop[self.ff[sl]]
                * weights(m.target, self.ii[sl])
                / weights(m.target, self.ff[sl])
            )[:, None] * np.einsum("pb,pba->pa", self.cc[sl], future)
        rates = np.vstack((up, down))
        changing = rates.sum(axis=0)
        unchanged = kin * (m.diag @ X)
        density = unchanged + changing
        if not np.any(density):
            self.rates.append(0.0)
            self.angular_u.extend([0.0, 1.0])
            self.deflection.extend([0.0, 2.0])
            self.changing.extend([0.0, 0.0])
            self.angular_offsets.append(len(self.angular_u))
            self.conditional_offsets.append(len(self.conditional_u))
            return
        u, total = angular_cdf(y, density)
        if energy >= 0.001:
            # The inclusive marginal is prescribed independently of the
            # rotational quadrature mesh. Integrate the actual IAA DCS instead
            # of approximating its angular knots on that mesh a second time.
            u = m.elastic.cdf(energy, y)
            total = (
                C
                * np.sqrt(energy * (energy + 2 * REST))
                / (energy + REST)
                * m.elastic.s.elastic(energy)
            )
        if energy == 0:
            # At rest there is no incident axis; the finite superelastic limit
            # has an isotropic angular marginal.
            u = y**2
        rho = np.divide(
            changing, density, out=np.zeros_like(changing), where=density > 0
        )
        distribution = np.divide(
            rates, changing, out=np.zeros_like(rates), where=changing > 0
        )
        keep = np.r_[True, np.diff(u) > 1e-15]
        u, rho, distribution = u[keep], rho[keep], distribution[:, keep]
        deflection = 2 * y[keep] ** 2
        if u[-1] != 1:
            raise ValueError("Missing angular endpoint")
        valid = np.flatnonzero(rho > 1e-14)
        if len(valid):
            for j in np.flatnonzero(rho <= 1e-14):
                distribution[:, j] = distribution[:, valid[np.argmin(abs(valid - j))]]
            selected = compress(u, distribution, rho, self.features)
            previous = None
            previous_cell = -1
            for j in selected:
                probability = distribution[:, j]
                if previous is None or not np.array_equal(probability, previous):
                    previous_cell = self.add_cell(probability, self.ids)
                    previous = probability.copy()
                self.conditional_u.append(u[j])
                self.conditional_cells.append(previous_cell)
        self.rates.append(total)
        self.angular_u.extend(u)
        self.deflection.extend(deflection)
        self.changing.extend(rho)
        self.angular_offsets.append(len(self.angular_u))
        self.conditional_offsets.append(len(self.conditional_u))

    def write(self, directory, energy, thresholds):
        high_cells = []
        high = self.high
        for offset, count in high["cells"]:
            offset, count = int(offset), int(count)
            cdf = high["cdf"][offset : offset + count]
            high_cells.append(
                self.add_cell(
                    np.diff(np.r_[0, cdf]),
                    high["alias"]["outcome"][offset : offset + count],
                )
            )
        flags = np.isin(energy, thresholds) | (energy == 0)
        arrays = [
            ("energies", "f64", energy),
            ("rates", "f64", self.rates),
            ("coordinates", "u32", flags),
            ("angular_offsets", "u32", self.angular_offsets),
            ("angular_u", "f64", self.angular_u),
            ("deflection", "f64", self.deflection),
            ("changing", "f32", self.changing),
            ("conditional_offsets", "u32", self.conditional_offsets),
            ("conditional_u", "f64", self.conditional_u),
            ("conditional_cells", "u32", self.conditional_cells),
            ("cell_offsets", "u32", self.cell_offsets),
            ("outcomes", "f64", self.palette),
            ("high_edges", "f64", high["edges"]),
            ("high_cells", "u32", high_cells),
        ]
        if self.sampling == "alias":
            arrays.append(("aliases", "alias", np.concatenate(self.aliases)))
        else:
            arrays.extend(
                [
                    ("cdf", "f64", np.concatenate(self.cdfs)),
                    ("outcome_ids", "u16", np.concatenate(self.outcome_ids)),
                ]
            )
        arrays.extend(lookup_arrays({name: values for name, _, values in arrays}))
        target = self.model.target
        metadata = [
            f"{target} reciprocal_hybrid {self.sampling}",
            f"{self.model.T:.17g} {MASSES[target]:.17g} {energy[-1]:.17g} 1000 10000 "
            f"{dict(N2=2.068, O2=2.281)[target]} {dict(N2=0.6052, O2=0.5677)[target]}",
        ]
        if self.storage_error > 2e-5 or self.pruned_probability > 1e-10:
            raise ValueError("Packed outcome accuracy budget exceeded")
        size = write_bundle(directory, metadata, arrays)
        return dict(
            target=target,
            bytes=size,
            energy_rows=len(energy),
            angular_nodes=len(self.angular_u),
            conditional_nodes=len(self.conditional_u),
            cells=len(self.cell_offsets) - 1,
            entries=self.cell_offsets[-1],
            storage_error=self.storage_error,
            pruned_probability=self.pruned_probability,
            high_partial_error=high["partial_error"],
        )


def load_reference(target, temperature, directory):
    checkpoint = directory / f"refined-reference-{target}.npz"
    maximum = (
        float(np.load(checkpoint)["maximum_energy"]) if checkpoint.exists() else 1000.0
    )
    model = Hybrid(target, temperature, maximum, 2, 48, angular_resolution=4)
    if not checkpoint.exists():
        return build_reference(target, temperature)[0]
    saved = np.load(checkpoint)
    if int(saved.get("reference_version", 0)) != REFERENCE_VERSION:
        raise ValueError(
            "Outdated normalization reference; run build_reference.py again"
        )
    if float(saved["temperature"]) != temperature or not np.array_equal(
        saved["y"], model.y
    ):
        raise ValueError("Reference checkpoint does not match the requested bath/grid")
    model.energy, model.X, model.c = saved["energy"], saved["X"], saved["c"]
    model.p = np.sqrt(model.energy * (model.energy + 2 * REST))
    model.kin = C / (model.energy + REST)
    model.boundary = np.searchsorted(model.energy, maximum + 1)
    anchor = np.flatnonzero(model.energy == 0.001)[0]
    model.cold_anchor = model.X[anchor, 0] / model.p[anchor]
    return model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=["N2", "O2"], required=True)
    parser.add_argument("--temperature", type=float, default=300)
    parser.add_argument("--maximum-energy", type=float, default=1e9)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sampling", choices=["alias", "cumulative"], default="alias")
    parser.add_argument(
        "--row-limit", type=int, help="Development-only prefix; never publish"
    )
    args = parser.parse_args()
    model = load_reference(args.target, args.temperature, args.reference_dir)
    low = np.load(args.reference_dir / f"adaptive-low-grid-{args.target}.npz")
    if int(low.get("reference_version", 0)) != REFERENCE_VERSION:
        raise ValueError("Outdated sampling grid; run adaptive_low_grid.py again")
    # The scalar tail is cheap; include its source joins and independently
    # refined rate knots. Angular joins receive additional logarithmic knots.
    energy = np.unique(
        np.r_[
            low["energy"],
            args.maximum_energy,
            np.geomspace(1000, 1e9, 301),
            np.geomspace(1000, 10000, 161),
            6000,
            8000,
            9000,
            10000,
        ]
    )
    energy = energy[(energy >= 0) & (energy <= args.maximum_energy)]
    if float(low["temperature"]) != args.temperature or low["energy"][-1] < min(
        args.maximum_energy, 1000
    ):
        raise ValueError("Sampling grid does not match the requested bath/range")
    for _ in range(20):
        hi = energy[energy >= 1000]
        probe = (hi[:-1] + hi[1:]) / 2

        def speed(e):
            return C * np.sqrt(e * (e + 2 * REST)) / (e + REST)

        exact = speed(probe) * model.elastic.s.elastic(probe)
        rows = speed(hi) * model.elastic.s.elastic(hi)
        bad = abs((rows[:-1] + rows[1:]) / 2 / exact - 1) > 5e-4
        if not bad.any():
            break
        energy = np.unique(np.r_[energy, probe[bad]])
    else:
        raise ValueError("High-energy rate refinement did not converge")
    if args.row_limit:
        energy = energy[: args.row_limit]
    exporter = Exporter(model, args.sampling)
    for i, value in enumerate(energy):
        exporter.row(value)
        if i % 25 == 0:
            print(
                args.target,
                i,
                "/",
                len(energy),
                "E",
                value,
                "entries",
                exporter.cell_offsets[-1],
                flush=True,
            )
    report = exporter.write(args.output, energy, low["thresholds"])
    (args.reference_dir / f"export-{args.target}-{args.sampling}.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(report, flush=True)


if __name__ == "__main__":
    main()
