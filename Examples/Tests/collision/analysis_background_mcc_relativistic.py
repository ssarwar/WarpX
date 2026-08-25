#!/usr/bin/env python3

import math
import re
from pathlib import Path

import numpy as np


HISTOGRAM_PATH = Path("hist_uz.txt")
PARTICLE_COUNT = 131072
INITIAL_UZ = 0.5810524930349901
BACKGROUND_DENSITY = 1.0e20
BACKGROUND_MASS_RATIO = 100.0
TIME_STEP = 1.0062975508741132e-9

C = 299792458.0
E_CHARGE = 1.602176634e-19
ELECTRON_MASS = 9.1093837139e-31

ENERGY_GRID = np.array([0.0, 1000.0, 5000.0, 100000.0])
CROSS_SECTION_GRID = np.array([0.0, 1.0e-20, 1.0e-20, 1.0e-21])


def load_histogram(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(encoding="utf-8") as stream:
        header = stream.readline()
    centers = np.array(
        [float(value) for value in re.findall(r"bin\d+=([-+0-9.eE]+)", header)]
    )
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return centers, data[-1, 2:]


def expected_collision_fraction() -> float:
    scan_energy = np.arange(0.0, 100000.0, 0.2)
    scan_sigma = np.interp(scan_energy, ENERGY_GRID, CROSS_SECTION_GRID)
    scan_speed = np.sqrt(2.0 * E_CHARGE * scan_energy / ELECTRON_MASS)
    nu_max = BACKGROUND_DENSITY * np.max(scan_sigma * scan_speed)

    gamma = math.sqrt(1.0 + INITIAL_UZ**2)
    electron_speed = C * INITIAL_UZ / gamma
    kinetic_energy = (gamma - 1.0) * ELECTRON_MASS * C**2 / E_CHARGE
    sigma = np.interp(kinetic_energy, ENERGY_GRID, CROSS_SECTION_GRID)

    candidate_probability = -math.expm1(-nu_max * TIME_STEP)
    return candidate_probability * BACKGROUND_DENSITY * sigma * electron_speed / nu_max


def expected_backward_uz() -> float:
    gamma = math.sqrt(1.0 + INITIAL_UZ**2)
    beta_com = INITIAL_UZ / (gamma + BACKGROUND_MASS_RATIO)
    gamma_com = 1.0 / math.sqrt(1.0 - beta_com**2)
    p_star = gamma_com * (INITIAL_UZ - beta_com * gamma)
    energy_star = gamma_com * (gamma - beta_com * INITIAL_UZ)
    return gamma_com * (-p_star + beta_com * energy_star)


bin_centers, histogram = load_histogram(HISTOGRAM_PATH)
assert np.isclose(np.sum(histogram), 1.0, rtol=0.0, atol=1.0e-12)

negative = bin_centers < 0.0
collision_fraction = np.sum(histogram[negative]) / np.sum(histogram)
mean_collided_uz = np.sum(histogram[negative] * bin_centers[negative]) / np.sum(
    histogram[negative]
)

expected_fraction = expected_collision_fraction()
standard_error = math.sqrt(
    expected_fraction * (1.0 - expected_fraction) / PARTICLE_COUNT
)
assert abs(collision_fraction - expected_fraction) < 6.0 * standard_error, (
    collision_fraction,
    expected_fraction,
    standard_error,
)

expected_uz = expected_backward_uz()
assert abs(mean_collided_uz - expected_uz) < 7.5e-4, (
    mean_collided_uz,
    expected_uz,
)

print(
    f"collision fraction: observed={collision_fraction:.8f}, "
    f"expected={expected_fraction:.8f}"
)
print(f"backward uz: observed={mean_collided_uz:.8f}, expected={expected_uz:.8f}")
