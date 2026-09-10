# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Calibrated PJG-type proton SDCS and parameterized offline fitting evaluator.

Independent Python reference for the production C++ implementation. Energies are in eV,
cross sections in cm^2, and the incident neutral is at rest. The free Bhabha
kernel is exact for a pointlike spin-1/2 proton at tree level. Molecular
distortion, optical shape, and binary-edge broadening remain empirical.
"""

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from reference import (
    BETHE_CONSTANT,
    ELECTRON_REST_ENERGY,
    PROTON_REST_ENERGY,
    RYDBERG,
    free_maximum_transfer,
)
from scipy.special import expit
from target_parameters import TARGETS, lorentzian_pair

MOLECULAR_REST_ENERGIES = {
    "N2": 28.0134 * 931494103.72,
    "O2": 31.9988 * 931494103.72,
}


@dataclass(frozen=True)
class MatchedParameters:
    j: float
    power: float
    amplitude: float
    width: float
    broad_excess: float
    bethe_scale: float = 1
    center_scale: float = 1
    low_width_ratio: float = 1
    edge_scale: float = 1

    def __post_init__(self):
        values = (
            self.j,
            self.power,
            self.amplitude,
            self.width,
            self.broad_excess,
            self.bethe_scale,
            self.low_width_ratio,
            self.edge_scale,
        )
        if any(not np.isfinite(value) or value <= 0 for value in values):
            raise ValueError(
                "Scales, widths, and distortion parameters must be positive"
            )
        if not np.isfinite(self.center_scale) or self.center_scale < 0:
            raise ValueError("The center numerator scale must be nonnegative")


def initial_parameters(target):
    p = TARGETS[target]
    return MatchedParameters(p.j, p.power, 1, p.gamma_s, p.broad_width)


def binding_log_scale(target):
    """Fix the logarithm to ln(4 E_e gamma^2 / mean(I_j) + e) - beta^2."""
    p = TARGETS[target]
    printed = np.exp(
        np.dot(p.fractions, np.log(np.array(p.bethe_constants) / p.thresholds))
    )
    return 1 / (printed * np.dot(p.fractions, p.thresholds))


PARAMETERS = MappingProxyType(
    {
        "N2": MatchedParameters(
            j=13.8235317299573,
            power=0.807,
            amplitude=0.8812026712177866,
            width=12.503227296612915,
            broad_excess=88.54007881169514,
            bethe_scale=binding_log_scale("N2"),
            center_scale=0.17909741025355236,
            low_width_ratio=1,
            edge_scale=1.1911400606628986,
        ),
        "O2": MatchedParameters(
            j=8.76635921598519,
            power=1.314,
            amplitude=0.7815344574102883,
            width=16.757498368297895,
            broad_excess=159.2697988004284,
            bethe_scale=binding_log_scale("O2"),
            center_scale=0.37449064801445814,
            low_width_ratio=1,
            edge_scale=1,
        ),
    }
)


def kinematics(energy):
    energy = np.asarray(energy, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy <= 0):
        raise ValueError("The incident kinetic energy must be finite and positive")
    mass = PROTON_REST_ENERGY
    gamma = 1 + energy / mass
    beta2 = (energy / (energy + mass)) * ((energy + 2 * mass) / (energy + mass))
    equivalent = ELECTRON_REST_ENERGY * beta2 / 2
    maximum = free_maximum_transfer(energy)
    return gamma, beta2, equivalent, maximum


def endpoint_momentum_broadening(energy):
    """Relativistic/NR ratio of free-endpoint sensitivity to target momentum.

    Differentiate exact, collinear two-body scattering at zero initial
    electron momentum. Normalize the derivative by its NR value at the
    same projectile speed. This supplies an energy-dependent broadening
    multiplier, not a bound-state momentum distribution or a new fit parameter.
    """
    gamma = 1 + np.asarray(energy) / PROTON_REST_ENERGY
    ratio = ELECTRON_REST_ENERGY / PROTON_REST_ENERGY
    denominator = 1 + 2 * gamma * ratio + ratio**2
    return gamma * (gamma + ratio) * (1 + ratio) ** 3 / denominator**2


def hard_factor(energy, secondary, continuation="exponential"):
    """Complete Bhabha bracket inside the free domain, positive continuation outside.

    Beyond free Tmax the exponential continuation matches the value and first
    derivative at the endpoint. It is explicitly an empirical bound-tail
    continuation, not free scattering at a forbidden transfer. The constant
    alternative is exposed only for a tail-model sensitivity check.
    """
    energy, t = np.broadcast_arrays(energy, secondary)
    gamma, beta2, _, maximum = kinematics(energy)
    if np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError("The secondary energy must be finite and nonnegative")
    mass = PROTON_REST_ENERGY
    electron = ELECTRON_REST_ENERGY
    total = energy + mass
    # Do not evaluate the free polynomial at an out-of-domain argument.
    bounded = np.minimum(t, maximum)
    x = bounded / maximum
    value = (1 - x) + x / gamma**2 + 0.5 * (bounded / total) ** 2
    if continuation == "constant":
        return value
    if continuation != "exponential":
        raise ValueError("Unknown bound-tail continuation")
    momentum = np.sqrt(energy * (energy + 2 * mass))
    invariant = mass**2 + electron**2 + 2 * electron * total
    # Factored 1-(Tmax/p)^2 keeps the terminal slope nonpositive at large gamma.
    minus = mass**2 + electron**2 + 2 * electron * mass**2 / (total + momentum)
    slope = -beta2 / maximum * (minus / invariant)
    slope *= 1 + 2 * electron * momentum / invariant
    terminal = 1 / gamma**2 + 0.5 * (maximum / total) ** 2
    return value * np.exp(slope / terminal * np.maximum(t - maximum, 0))


def molecular_endpoint(
    target, energy, binding, projectile_rest_energy=PROTON_REST_ENERGY
):
    """Exact upper electron-energy endpoint for a specified three-body channel."""
    energy = np.asarray(energy, dtype=float)
    mass = projectile_rest_energy
    electron = ELECTRON_REST_ENERGY
    neutral = MOLECULAR_REST_ENERGIES[target]
    threshold = binding * (1 + mass / neutral) + binding**2 / (2 * neutral)
    root_s = np.sqrt((mass + neutral) ** 2 + 2 * neutral * energy)
    excess = 2 * neutral * np.maximum(energy - threshold, 0)
    excess /= root_s + mass + neutral + binding
    kinetic = excess * (1 - (electron + excess / 2) / root_s)
    p2 = energy * (energy + 2 * mass)
    gm1 = p2 / (root_s * (energy + mass + neutral + root_s))
    endpoint = kinetic + gm1 * (kinetic + electron)
    endpoint += np.sqrt(p2) / root_s * np.sqrt(kinetic * (kinetic + 2 * electron))
    return np.where(energy >= threshold, endpoint, 0)


def molecular_angular_fraction(
    target, energy, secondary, binding, projectile_rest_energy=PROTON_REST_ENERGY
):
    """Fraction of isotropic directions allowed by exact three-body support.

    Used only as a parameter-free molecular-endpoint taper. It is not the
    ionization angular distribution or a near-threshold dynamical theory.
    It is identically one over the spectra used in the refit to numerical
    significance; the empirical binary gate suppresses the distant endpoint.
    """
    energy, t = np.broadcast_arrays(energy, secondary)
    mass = projectile_rest_energy
    electron = ELECTRON_REST_ENERGY
    neutral = MOLECULAR_REST_ENERGIES[target]
    threshold = binding * (1 + mass / neutral) + binding**2 / (2 * neutral)
    constant = (neutral - electron) * (energy - threshold)
    constant -= electron * binding * (mass + binding / 2) / neutral
    numerator = (energy + mass + neutral) * t - constant
    denominator = np.sqrt(energy * (energy + 2 * mass) * t * (t + 2 * electron))
    fraction = np.divide(
        denominator - numerator,
        2 * denominator,
        out=np.where(numerator <= 0, 1.0, 0.0),
        where=denominator > 0,
    )
    return np.where(
        (energy >= threshold)
        & (t < molecular_endpoint(target, energy, binding, projectile_rest_energy)),
        np.clip(fraction, 0, 1),
        0,
    )


class SpectrumGrid:
    """Cache parameter-independent physics for repeated host-side fitting."""

    def __init__(self, target, energy, secondary, continuation="exponential"):
        self.target = target
        self.energy, self.secondary = np.broadcast_arrays(energy, secondary)
        if np.any(~np.isfinite(self.secondary)) or np.any(self.secondary < 0):
            raise ValueError("The secondary energy must be finite and nonnegative")
        gamma, self.beta2, self.equivalent, maximum = kinematics(self.energy)
        self.hard_factor = hard_factor(self.energy, self.secondary, continuation)
        self.distortion_energy = self.energy * ELECTRON_REST_ENERGY / PROTON_REST_ENERGY
        p = TARGETS[target]
        # The printed C_j/I_j differ only by table rounding. Replace their
        # independently named values by one weighted geometric mean. The
        # selected model fixes its multiplier from the mean binding energy.
        c = np.exp(
            np.dot(p.fractions, np.log(np.array(p.bethe_constants) / p.thresholds))
        )
        self.log_argument = 4 * self.equivalent * c * gamma**2
        broadening = endpoint_momentum_broadening(self.energy)
        alpha = {"N2": 0.70, "O2": 0.59}[target]
        self.gate = np.zeros_like(self.secondary, dtype=float)
        self.gate_terms = []
        for fraction, binding in zip(p.fractions, p.thresholds, strict=True):
            width = broadening * np.sqrt(self.equivalent * binding)
            center = maximum - 2 * width - RYDBERG / 4
            argument = alpha * (center - self.secondary) / width
            weight = fraction * molecular_angular_fraction(
                target, self.energy, self.secondary, binding
            )
            self.gate_terms.append((argument, weight, binding))
            self.gate += weight * expit(argument)

    def parts(self, parameters, hard_distortion="fading", center_s=None):
        p = TARGETS[self.target]
        # The optional override studies refitting the existing PJG T_s in
        # place of its T_a multiplier, not adding a seventh fit parameter.
        peak = p.center_s if center_s is None else center_s
        center = peak - parameters.center_scale * p.center_numerator / (
            self.equivalent + p.center_denominator
        )
        width = parameters.width * (
            1
            + (parameters.low_width_ratio - 1)
            * p.gamma_denominator
            / (self.equivalent + p.gamma_denominator)
        )
        broad, difference = lorentzian_pair(
            self.secondary, width, center, parameters.broad_excess
        )
        distortion = expit(
            parameters.power * np.log(self.distortion_energy / parameters.j)
        )
        if hard_distortion == "fading":
            hard_d = 1 - (1 - distortion) / (
                1 + (self.secondary / parameters.broad_excess) ** 2
            )
        elif hard_distortion == "common":
            hard_d = distortion
        elif hard_distortion == "none":
            hard_d = 1
        else:
            raise ValueError("Unknown hard-distortion model")
        logarithm = (
            np.log(parameters.bethe_scale * self.log_argument + np.e) - self.beta2
        )
        gate = (
            self.gate
            if parameters.edge_scale == 1
            else sum(
                weight * expit(parameters.edge_scale * argument)
                for argument, weight, _ in self.gate_terms
            )
        )
        scale = gate / self.equivalent
        soft = scale * distortion * parameters.amplitude * p.k * width**2
        soft *= logarithm * difference
        hard = scale * hard_d * p.electrons * BETHE_CONSTANT * broad * self.hard_factor
        return soft, hard

    def __call__(self, parameters, hard_distortion="fading", center_s=None):
        soft, hard = self.parts(parameters, hard_distortion, center_s)
        return soft + hard

    def effective_binding(self, parameters):
        """Conditional mean of the retained PJG continuum thresholds at (E,T)."""
        gates = [
            weight * expit(parameters.edge_scale * argument)
            for argument, weight, _ in self.gate_terms
        ]
        denominator = sum(gates)
        numerator = sum(
            gate * binding
            for gate, (_, _, binding) in zip(gates, self.gate_terms, strict=True)
        )
        return np.divide(
            numerator,
            denominator,
            out=np.full_like(denominator, min(TARGETS[self.target].thresholds)),
            where=denominator > 0,
        )


def sdcs(target, energy, secondary, parameters, **kwargs):
    energy = np.asarray(energy, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0):
        raise ValueError("The incident energy must be finite and nonnegative")
    value = SpectrumGrid(target, np.where(energy > 0, energy, 1), secondary, **kwargs)(
        parameters
    )
    return np.where(energy > 0, value, 0)


def optical_coefficient(
    target, secondary, parameters, relativistic=True, center_s=None
):
    p = TARGETS[target]
    equivalent = ELECTRON_REST_ENERGY / 2 if relativistic else np.inf
    peak = p.center_s if center_s is None else center_s
    center = peak - parameters.center_scale * p.center_numerator / (
        equivalent + p.center_denominator
    )
    width = parameters.width * (
        1
        + (parameters.low_width_ratio - 1)
        * p.gamma_denominator
        / (equivalent + p.gamma_denominator)
    )
    _, difference = lorentzian_pair(secondary, width, center, parameters.broad_excess)
    return parameters.amplitude * p.k * width**2 / BETHE_CONSTANT * difference
