# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Energy moments of the inclusive effective-pair PJG reference.

The electron kinetic-energy moment follows directly from the SDCS. Adding
the conditional PJG binding threshold defines the effective-pair ionization
loss, not a complete stopping model or a unique reconstruction of projectile
loss from inclusive electron yield. Neither discrete excitation nor nuclear
stopping is included. The production projectile remains prescribed/rigid.
"""

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pjg_model import (
    SpectrumGrid,
    endpoint_momentum_broadening,
    kinematics,
    molecular_endpoint,
)
from reference import (
    BETHE_CONSTANT,
    ELECTRON_REST_ENERGY,
    ORBITALS,
    PROTON_REST_ENERGY,
    RUDD_PARAMETERS,
    RYDBERG,
)
from scipy.integrate import simpson
from scipy.special import expit
from target_parameters import TARGETS

AVOGADRO = 6.02214076e23
MOLAR_MASS = {"N2": 28.0134, "O2": 31.9988}
PSTAR_DENSITY = {"N2": 1.16528e-3, "O2": 1.33151e-3}
PSTAR_MEAN_EXCITATION = {"N2": 82.0, "O2": 95.0}
PSTAR_FORM = "https://physics.nist.gov/PhysRefData/Star/Text/PSTAR-t.html"
PSTAR_ENDPOINT = "https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl"


def mass_stopping(target, loss_cross_section):
    """Convert eV cm^2 per molecule to MeV cm^2/g, with no density factor."""
    return np.asarray(loss_cross_section) * AVOGADRO / MOLAR_MASS[target] * 1e-6


@dataclass(frozen=True)
class PstarTable:
    target: str
    # MeV; electronic/nuclear/total MeV cm^2/g; ranges g/cm^2; detour factor.
    rows: np.ndarray
    sha256: str

    @property
    def energies(self):
        return self.rows[:, 0] * 1e6

    @property
    def electronic(self):
        return self.rows[:, 1]


def load_pstar(target, path=None):
    """Read the preserved NIST rows or an external text-only NIST HTML output."""
    path = Path(path) if path else Path(__file__).with_name("pstar_reference.json")
    raw = path.read_bytes()
    if path.suffix == ".json":
        data = json.loads(raw)
        if data["columns"][:2] != ["energy_MeV", "electronic_MeV_cm2_per_g"]:
            raise ValueError("The PSTAR energy or stopping-power units do not match")
        if data["materials"][target] != {"N2": "007", "O2": "008"}[target]:
            raise ValueError("The PSTAR material does not match")
        rows = np.asarray(data["rows"][target], dtype=float)
    else:
        html = raw.decode("latin1")
        material = {"N2": "NITROGEN", "O2": "OXYGEN"}[target]
        if (
            material not in html
            or "STOP(e) = electronic stopping power, MeV cm2/g" not in html
        ):
            raise ValueError("The PSTAR material or stopping-power units do not match")
        text = re.sub(r"<br\s*/?>", "\n", html, flags=re.IGNORECASE)
        pattern = r"^\s*((?:\d+\.\d+E[+-]\d+\s+){6}\d+\.\d+)\s*$"
        rows = np.array(
            [
                list(map(float, row.split()))
                for row in re.findall(pattern, text, re.MULTILINE)
            ]
        )
    if rows.ndim != 2 or rows.shape[1] != 7 or len(rows) < 100:
        raise ValueError("PSTAR did not return the expected default seven-column table")
    if np.any(np.diff(rows[:, 0]) <= 0) or np.any(rows <= 0):
        raise ValueError(
            "PSTAR energies must increase and the returned quantities must be positive"
        )
    if rows[0, 0] != 0.001 or rows[-1, 0] != 10000:
        raise ValueError("The default PSTAR energy range is incomplete")
    # Independent columns are rounded to four significant figures.
    if np.max(np.abs((rows[:, 1] + rows[:, 2]) / rows[:, 3] - 1)) > 1e-3:
        raise ValueError(
            "PSTAR electronic plus nuclear stopping does not match its total"
        )
    return PstarTable(target, rows, hashlib.sha256(raw).hexdigest())


@dataclass(frozen=True)
class LossMoments:
    total: np.ndarray
    kinetic: np.ndarray
    binding: np.ndarray
    kinetic_second: np.ndarray
    tail_total: np.ndarray
    tail_kinetic: np.ndarray
    tail_kinetic_second: np.ndarray

    @property
    def ionization(self):
        return self.kinetic + self.binding


class LossGrid:
    """Resolve the soft peak, both sides of the binary edge, and molecular tail.

    The first/fourth segments are logarithmic in 1+T; the middle segments
    are linear in T. All four retain their own exact Jacobians. Splitting
    at free Tmax makes the above-free moment independent of bin alignment.
    The outer tail is integrated through the full molecular endpoint.
    """

    def __init__(self, target, energies, points=1025, continuation="exponential"):
        if points < 3 or points % 2 != 1:
            raise ValueError("Use an odd number of quadrature points, at least three")
        self.target = target
        self.energies = np.atleast_1d(np.asarray(energies, dtype=float))
        _, _, equivalent, free = kinematics(self.energies)
        p = TARGETS[target]
        endpoint = molecular_endpoint(target, self.energies, min(p.thresholds))
        if np.any(free >= endpoint):
            raise ValueError(
                "This stopping audit requires molecular support above free Tmax"
            )
        width = endpoint_momentum_broadening(self.energies) * np.sqrt(
            equivalent * max(p.thresholds)
        )
        # Always resolve the soft peak logarithmically. At moderate energies
        # 32 widths can exceed Tmax: a linear interval starting at zero then
        # under-resolves the peak even if it resolves the binary edge well.
        left = np.maximum(free / 2, free - 32 * width)
        right = np.minimum(endpoint, free + 32 * width)
        u = np.linspace(0, 1, points)[None, :]
        segments, jacobians = [], []
        for lower, upper, logarithmic in (
            (np.zeros_like(free), left, True),
            (left, free, False),
            (free, right, False),
            (right, endpoint, True),
        ):
            if logarithmic:
                length = (np.log1p(upper) - np.log1p(lower))[:, None]
                x = np.log1p(lower)[:, None] + length * u
                t = np.expm1(x)
                jacobian = length * np.exp(x)
            else:
                length = (upper - lower)[:, None]
                t = lower[:, None] + length * u
                jacobian = np.broadcast_to(length, t.shape)
            segments.append(t)
            jacobians.append(jacobian)
        self.secondary = np.stack(segments, axis=1)
        self.jacobian = np.stack(jacobians, axis=1)
        self.dx = 1 / (points - 1)
        self.grid = SpectrumGrid(
            target, self.energies[:, None, None], self.secondary, continuation
        )

    def __call__(self, parameters, center_s=None):
        spectrum = self.grid(parameters, center_s=center_s)
        binding = self.grid.effective_binding(parameters)
        weighted = spectrum * self.jacobian
        integrals = []
        for factor in (
            np.ones_like(binding),
            self.secondary,
            binding,
            self.secondary**2,
        ):
            integrals.append(simpson(weighted * factor, dx=self.dx, axis=-1))
        return LossMoments(
            *(value.sum(axis=1) for value in integrals),
            integrals[0][:, 2:].sum(axis=1),
            integrals[1][:, 2:].sum(axis=1),
            integrals[3][:, 2:].sum(axis=1),
        )


def free_second_moment(target, energy):
    """Analytic integral of T^2 times the full free Bhabha SDCS on [0,Tmax]."""
    _, beta2, equivalent, maximum = kinematics(energy)
    integral = maximum * (1 - beta2 / 2)
    integral += maximum**3 / (6 * (np.asarray(energy) + PROTON_REST_ENERGY) ** 2)
    return TARGETS[target].electrons * BETHE_CONSTANT / equivalent * integral


def optical_strength(target, parameters, center_s=None):
    """Analytic effective continuum integral of (T+mean(I))*a(T), not an exact TRK sum."""
    p = TARGETS[target]
    peak = p.center_s if center_s is None else center_s
    center = peak - parameters.center_scale * p.center_numerator / (
        ELECTRON_REST_ENERGY / 2 + p.center_denominator
    )
    width = parameters.width * (
        1
        + (parameters.low_width_ratio - 1)
        * p.gamma_denominator
        / (ELECTRON_REST_ENERGY / 2 + p.gamma_denominator)
    )
    broad = np.hypot(width, parameters.broad_excess)
    area = (np.pi / 2 + np.arctan(center / width)) / width
    area -= (np.pi / 2 + np.arctan(center / broad)) / broad
    integral = 0.5 * np.log1p(parameters.broad_excess**2 / (center**2 + width**2))
    integral += (center + np.dot(p.fractions, p.thresholds)) * area
    return parameters.amplitude * p.k * width**2 / BETHE_CONSTANT * integral


def rudd_orbital_sdcs(target, energy, secondary):
    """Separate the original nonrelativistic Rudd SDCS into its nominal orbitals.

    This independently vectorizes Eqs. (41)--(48), using the published
    orbital thresholds/occupancies and empirical coefficients in reference.py.
    It is not a relativistic Rudd extension or an independent experiment.
    """
    energy, secondary = np.broadcast_arrays(energy, secondary)
    if np.any((energy < 5000) | (energy > 4e6)) or np.any(secondary < 0):
        raise ValueError(
            "This nonrelativistic comparison requires 5--4000 keV and T >= 0"
        )
    threshold, number = map(np.asarray, ORBITALS[target])
    coefficients = np.array(
        [
            RUDD_PARAMETERS["inner" if value > 2 * threshold[0] else target]
            for value in threshold
        ]
    ).T
    a1, b1, c1, d1, e1, a2, b2, c2, d2, alpha = coefficients
    v = np.sqrt(
        energy[..., None] * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY / threshold
    )
    w = secondary[..., None] / threshold
    first = a1 * np.log1p(v**2) / (v**2 + b1 / v**2)
    first += c1 * v**d1 / (1 + e1 * v ** (d1 + 4))
    high, low = a2 / v**2 + b2 / v**4, c2 * v**d2
    second = high * low / (high + low)
    cutoff = 4 * v**2 - 2 * v - RYDBERG / (4 * threshold)
    return (
        BETHE_CONSTANT
        * number
        / threshold**3
        * (first + second * w)
        / (1 + w) ** 3
        * expit(alpha * (cutoff - w) / v)
    )


def rudd_nominal_loss(target, energies, points=8193):
    """Nominal Rudd orbital-binding loss: a correlated, nonrelativistic diagnostic.

    Integrate the recommended curve through T=E; its empirical tail is
    negligible there in this energy range. This does not reconstruct unique
    projectile energy loss from inclusive yield. O2 SDCS values above
    300 keV and N2 values above 1700 keV are reference-curve extrapolations.
    """
    energy = np.atleast_1d(energies)
    x = np.linspace(0, 1, points)[None, :] * np.log1p(energy[:, None])
    t = np.expm1(x)
    channels = rudd_orbital_sdcs(target, energy[:, None], t)
    loss = np.sum(channels * (t[..., None] + np.array(ORBITALS[target][0])), axis=-1)
    return simpson(loss * np.exp(x), x=x, axis=-1)
