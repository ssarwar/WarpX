#!/usr/bin/env python3
"""Compare actual PIC products with the independent final Python SDCS.

Python uses vectorized physics and Simpson integration, not the C++ host
quadrature or device inverse table. Experimental-fit tests live separately.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    str(
        Path(__file__).resolve().parents[3] / "Tools/Algorithms/ProtonImpactIonization"
    ),
)
from calibrated_pjg import PARAMETERS, total_cross_section
from pjg_model import SpectrumGrid, molecular_endpoint
from pjg_moments import LossGrid
from pjg_properties import mapped_cdf
from reference import free_maximum_transfer
from target_parameters import TARGETS

C = 299_792_458.0
K_B = 1.380_649e-23
Q_E = 1.602_176_634e-19


def kinetic_energy(ux, uy, uz, mass):
    proper_speed_squared = ux**2 + uy**2 + uz**2
    gamma = np.sqrt(1 + proper_speed_squared / C**2)
    return mass * proper_speed_squared / ((gamma + 1) * Q_E)


data = np.load("proton_impact_ionization_results.npz")
projectile_mass = float(data["projectile_rest_energy"]) * Q_E / C**2
steps = int(data["steps"])
for target in ("N2", "O2"):
    for component in ("x", "y", "z", "ux", "uy", "uz", "w", "id"):
        np.testing.assert_array_equal(
            data[f"{target}_beam_initial_{component}"],
            data[f"{target}_beam_{component}"],
        )
    electron_weight = data[f"{target}_electrons_w"]
    ion_weight = data[f"{target}_ions_w"]
    assert 10000 < len(electron_weight) <= int(data["max_products_per_cell"]) * steps
    np.testing.assert_array_equal(electron_weight, ion_weight)
    assert np.all(electron_weight == float(data["fixed_product_weight"]))
    for position in ("x", "y", "z"):
        np.testing.assert_array_equal(
            data[f"{target}_electrons_{position}"], data[f"{target}_ions_{position}"]
        )
    for species in ("electrons", "ions"):
        ids = data[f"{target}_{species}_id"]
        assert np.unique(ids).size == ids.size

    # Quantize only roundoff in prescribed input energies, not a spectrum.
    beam_energy = np.round(
        kinetic_energy(
            *(data[f"{target}_beam_{c}"] for c in ("ux", "uy", "uz")), projectile_mass
        ),
        4,
    )
    energies, inverse = np.unique(beam_energy, return_inverse=True)
    sigma = total_cross_section(target, energies)
    rest = float(data["projectile_rest_energy"])
    speed = C * np.sqrt(energies * (energies + 2 * rest)) / (energies + rest)
    beam_weights = np.bincount(inverse, weights=data[f"{target}_beam_w"])
    score = beam_weights * sigma * 1e-4 * speed
    expected_weight = (
        score.sum()
        * float(data["background_density"])
        * float(data["time_step"])
        * steps
    )
    assert np.isclose(electron_weight.sum(), expected_weight, rtol=2e-3)

    energy = kinetic_energy(
        *(data[f"{target}_electrons_{c}"] for c in ("ux", "uy", "uz")),
        9.109_383_7139e-31,
    )
    assert np.all(energy >= 0)
    endpoint = max(
        molecular_endpoint(target, energies, min(TARGETS[target].thresholds))
    )
    assert energy.max() <= endpoint * (1 + 2e-12)
    grid = LossGrid(target, energies, 8193)
    density = grid.grid(PARAMETERS[target])
    cdfs = [
        mapped_cdf(grid.secondary[i], density[i] * grid.jacobian[i], grid.dx)[:2]
        for i in range(len(energies))
    ]
    mixture = score / score.sum()
    sorted_energy = np.sort(energy)
    probability = sum(
        f * np.interp(sorted_energy, x, cdf)
        for f, (x, cdf) in zip(mixture, cdfs, strict=True)
    )
    empirical = (np.arange(len(energy)) + 0.5) / len(energy)
    ks = np.max(np.abs(probability - empirical))
    assert ks < 4e-3
    x = np.unique(np.concatenate([values[0] for values in cdfs]))
    cdf = sum(f * np.interp(x, t, p) for f, (t, p) in zip(mixture, cdfs, strict=True))
    quantiles = np.array([0.01, 0.1, 0.5, 0.9, 0.99, 0.999])
    expected_quantiles = np.interp(quantiles, cdf, x)
    observed_quantiles = np.quantile(energy, quantiles)
    np.testing.assert_allclose(
        observed_quantiles[:-1], expected_quantiles[:-1], rtol=2e-2, atol=0.03
    )
    assert np.isclose(observed_quantiles[-1], expected_quantiles[-1], rtol=8e-2)
    moments = grid(PARAMETERS[target])
    expected_mean = np.dot(mixture, moments.kinetic / moments.total)
    assert np.isclose(np.mean(energy), expected_mean, rtol=1e-2)
    if len(energies) == 1:
        tail = np.mean(energy > free_maximum_transfer(energies[0]))
        expected_tail = float(moments.tail_total[0] / moments.total[0])
        assert abs(tail - expected_tail) < 2e-3

    # Check the retained closure's geometry, not agreement with measured DDCS.
    electron_u = np.column_stack(
        [data[f"{target}_electrons_{c}"] for c in ("ux", "uy", "uz")]
    )
    direction = data[f"{target}_direction"]
    unit = electron_u / np.linalg.norm(electron_u, axis=1)[:, None]
    cosine = unit @ direction
    assert np.all(abs(cosine) <= 1 + 2e-14)
    if len(energies) == 1:
        binding = SpectrumGrid(target, energies[0], energy).effective_binding(
            PARAMETERS[target]
        )
        maximum = free_maximum_transfer(energies[0])
        free_energy = np.minimum(energy, maximum)
        free_cosine = np.sqrt(
            free_energy
            * (maximum + 2 * 510998.95069)
            / (maximum * (free_energy + 2 * 510998.95069))
        )
        center = free_cosine * (energy + binding / 2) / (energy + binding)
        width = binding / (energy + binding)
        # Exact moments/CDF of the conditioned uniform closure, not a DDCS fit.
        lower = np.maximum(-1, center - width)
        upper = np.minimum(1, center + width)
        expected_cosine = (lower + upper) / 2
        assert abs(cosine.mean() - expected_cosine.mean()) < 1.5e-2
        expected_second = (lower**2 + lower * upper + upper**2) / 3
        assert abs(np.mean(cosine**2) - expected_second.mean()) < 1.5e-2
        conditional_quantile = (cosine - lower) / (upper - lower)
        assert np.all(
            (conditional_quantile > -2e-5) & (conditional_quantile < 1 + 2e-5)
        )
        angular_ks = np.max(np.abs(np.sort(conditional_quantile) - empirical))
        assert angular_ks < 4e-3
        assert np.count_nonzero(cosine == 1) == 0
        assert cosine[energy >= np.quantile(energy, 0.9)].mean() > (
            cosine[energy <= np.quantile(energy, 0.5)].mean()
        )
    transverse = unit.mean(axis=0) - np.dot(unit.mean(axis=0), direction) * direction
    assert np.linalg.norm(transverse) < 1.5e-2
    thermal_speed = math.sqrt(
        K_B
        * float(data[f"{target}_temperature"])
        / float(data[f"{target}_neutral_mass"])
    )
    ion_u = np.column_stack([data[f"{target}_ions_{c}"] for c in ("ux", "uy", "uz")])
    assert np.all(abs(ion_u.mean(axis=0)) < 0.03 * thermal_speed)
    np.testing.assert_allclose(ion_u.std(axis=0), thermal_speed, rtol=0.035)
    print(
        f"{target}: pairs={len(energy)}, energies={energies}, KS={ks:.3e}, mean={energy.mean():.6f}/{expected_mean:.6f} eV"
    )

assert float(data["initialization_elapsed"]) < 30
assert float(data["step_elapsed"]) < 30
