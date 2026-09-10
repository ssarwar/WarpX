# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Independent energy-moment checks, not acceptance of an empirical SDCS fit."""

import unittest

import numpy as np
from original_pjg import bhabha
from pjg_model import (
    PARAMETERS,
    SpectrumGrid,
    endpoint_momentum_broadening,
    kinematics,
    molecular_endpoint,
    optical_coefficient,
)
from pjg_moments import (
    AVOGADRO,
    LossGrid,
    free_second_moment,
    load_pstar,
    mass_stopping,
    optical_strength,
    rudd_orbital_sdcs,
)
from reference import rudd_sdcs
from scipy.integrate import quad
from scipy.special import expit
from target_parameters import TARGETS


class StoppingChecks(unittest.TestCase):
    def test_separate_rudd_orbitals_reproduce_reference_sdcs(self):
        energy = np.geomspace(5000, 4e6, 23)[:, None]
        t = np.expm1(np.linspace(0, np.log1p(1e4), 97))[None, :]
        for target in TARGETS:
            np.testing.assert_allclose(
                rudd_orbital_sdcs(target, energy, t).sum(axis=-1),
                rudd_sdcs(target, energy, t),
                rtol=2e-14,
                atol=0,
            )

    def test_reference_units_and_materials(self):
        for target, expected, atomic_mass in (
            ("N2", 225.9, 14.0067),
            ("O2", 216.1, 15.9994),
        ):
            table = load_pstar(target)
            self.assertEqual(
                table.electronic[np.flatnonzero(table.energies == 1e6)[0]], expected
            )
            # Two atoms per target molecule: using atomic mass would double
            # the stopping inferred from this per-molecule cross section.
            self.assertAlmostEqual(
                mass_stopping(target, 1e-15), 1e-21 * AVOGADRO / (2 * atomic_mass)
            )

    def test_composite_quadrature_convergence(self):
        energies = [5e3, 1e5, 1.75e6, 8e8, 1e10]
        for target, parameters in PARAMETERS.items():
            reference = LossGrid(target, energies, 4097)(parameters)
            result = LossGrid(target, energies, 1025)(parameters)
            for field in reference.__dataclass_fields__:
                np.testing.assert_allclose(
                    getattr(result, field), getattr(reference, field), rtol=1e-7, atol=0
                )

    def test_moments_against_adaptive_integrals(self):
        # Integrate in log(1+T), independently of the composite-grid mapping.
        # Form the binding sum explicitly from its channel weights instead
        # of using effective_binding in the expected result.
        for target, parameters in PARAMETERS.items():
            for energy in (5e3, 1.75e6, 1e10):
                _, _, equivalent, free = kinematics(energy)
                endpoint = molecular_endpoint(
                    target, energy, min(TARGETS[target].thresholds)
                )
                width = endpoint_momentum_broadening(energy) * np.sqrt(
                    equivalent * max(TARGETS[target].thresholds)
                )
                boundaries = np.unique(
                    np.clip([free - 10 * width, free, free + 10 * width], 0, endpoint)
                )
                expected = []
                for moment in range(4):

                    def integrand(x):
                        t = np.expm1(x)
                        grid = SpectrumGrid(target, energy, t)
                        value = float(grid(parameters)) * np.exp(x) * 1e16
                        if moment == 2:
                            terms = [
                                (
                                    float(
                                        weight * expit(parameters.edge_scale * argument)
                                    ),
                                    binding,
                                )
                                for argument, weight, binding in grid.gate_terms
                            ]
                            gate = sum(weight for weight, _ in terms)
                            factor = (
                                sum(weight * binding for weight, binding in terms)
                                / gate
                                if gate
                                else 0
                            )
                        else:
                            factor = t ** (2 if moment == 3 else moment)
                        return value * factor

                    value, _ = quad(
                        integrand,
                        0,
                        np.log1p(endpoint),
                        points=np.log1p(boundaries),
                        epsabs=1e-10,
                        epsrel=1e-9,
                        limit=250,
                    )
                    expected.append(value * 1e-16)
                result = LossGrid(target, [energy])(parameters)
                actual = [
                    result.total[0],
                    result.kinetic[0],
                    result.binding[0],
                    result.kinetic_second[0],
                ]
                np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=0)

    def test_free_second_moment_against_bhabha_quadrature(self):
        for target in TARGETS:
            for energy in (5e3, 1e6, 8e8, 1e10):
                _, _, _, free = kinematics(energy)
                result, _ = quad(
                    lambda x: (x * free) ** 2 * bhabha(target, energy, x * free),
                    0,
                    1,
                    epsabs=0,
                    epsrel=1e-12,
                )
                self.assertAlmostEqual(
                    result * free / free_second_moment(target, energy), 1, delta=1e-12
                )

    def test_optical_strength_against_quadrature(self):
        for target, parameters in PARAMETERS.items():
            binding = np.dot(TARGETS[target].fractions, TARGETS[target].thresholds)
            for peak in (None, -5, 10):
                result, _ = quad(
                    lambda x: (
                        (np.expm1(x) + binding)
                        * optical_coefficient(
                            target, np.expm1(x), parameters, center_s=peak
                        )
                        * np.exp(x)
                    ),
                    0,
                    np.log1p(1e9),
                    epsabs=1e-9,
                    epsrel=1e-10,
                )
                self.assertAlmostEqual(
                    result / optical_strength(target, parameters, peak), 1, delta=1e-9
                )
            self.assertLess(
                optical_strength(target, parameters), TARGETS[target].electrons
            )

    def test_positive_moment_and_available_energy_budgets(self):
        energy = np.geomspace(5e3, 1e10, 31)
        for target, parameters in PARAMETERS.items():
            result = LossGrid(target, energy)(parameters)
            self.assertTrue(np.all(result.ionization > result.kinetic))
            self.assertTrue(np.all(result.ionization < energy * result.total))
            self.assertTrue(
                np.all(result.kinetic_second * result.total >= result.kinetic**2)
            )
            self.assertTrue(np.all(result.kinetic_second < energy * result.kinetic))
            low, high = min(TARGETS[target].thresholds), max(TARGETS[target].thresholds)
            self.assertTrue(np.all(result.binding >= low * result.total))
            self.assertTrue(np.all(result.binding <= high * result.total))
            yield_fraction = result.tail_total / result.total
            energy_fraction = result.tail_kinetic / result.kinetic
            second_fraction = result.tail_kinetic_second / result.kinetic_second
            self.assertTrue(
                np.all((yield_fraction >= 0) & (yield_fraction <= energy_fraction))
            )
            self.assertTrue(
                np.all((energy_fraction <= second_fraction) & (second_fraction <= 1))
            )


if __name__ == "__main__":
    unittest.main()
