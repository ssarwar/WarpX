#!/usr/bin/env python3
"""Check N2/O2 source budgets, sub-particle remainders, caps and bare-ion scaling."""

import sys
from pathlib import Path

import numpy as np

from pywarpx import picmi

sys.path.insert(
    0,
    str(
        Path(__file__).resolve().parents[3] / "Tools/Algorithms/ProtonImpactIonization"
    ),
)
from calibrated_pjg import total_cross_section

DT = 1.0e-9
ENERGY = 50.0e3
C = picmi.constants.c
ME = picmi.constants.m_e
MP = picmi.constants.m_p
QE = picmi.constants.q_e
MU = 1.66053906660e-27
REST = MP * C**2 / QE
PROPER_SPEED = C * np.sqrt(ENERGY * (ENERGY + 2 * REST)) / REST
SPEED = PROPER_SPEED / (1 + ENERGY / REST)
BEAM_DENSITY = 1.0e8
STEPS = 7

grid = picmi.Cartesian3DGrid(
    number_of_cells=[1] * 3,
    lower_bound=[0.0] * 3,
    upper_bound=[1.0] * 3,
    lower_boundary_conditions=["periodic"] * 3,
    upper_boundary_conditions=["periodic"] * 3,
    lower_boundary_conditions_particles=["periodic"] * 3,
    upper_boundary_conditions_particles=["periodic"] * 3,
    warpx_max_grid_size=1,
    warpx_blocking_factor=1,
)
solver = picmi.ElectromagneticSolver(grid=grid, method="Yee")
cases = {}
collisions = []
options = dict(
    warpx_do_not_push=True, warpx_do_not_deposit=True, warpx_do_not_gather=True
)
for target, mass_number in (("N2", 28.0134), ("O2", 31.9988)):
    sigma = float(total_cross_section(target, ENERGY)) * 1e-4
    density = 0.37 / (BEAM_DENSITY * sigma * SPEED * DT)
    for mode in ("coarse", "fine", "pulse", "alpha"):
        name = f"{target}_{mode}"
        alpha = mode == "alpha"
        beam = picmi.Species(
            name=f"beam_{name}",
            charge=(2 if alpha else 1) * QE,
            mass=(4 if alpha else 1) * MP,
            initial_distribution=picmi.UniformDistribution(
                density=BEAM_DENSITY, directed_velocity=[0.0, 0.0, PROPER_SPEED]
            ),
            **options,
        )
        electrons = picmi.Species(
            name=f"electrons_{name}", particle_type="electron", **options
        )
        ions = picmi.Species(
            name=f"ions_{name}", charge=QE, mass=mass_number * MU - ME, **options
        )
        weight = 0.125 if mode == "fine" else 1.0
        cap = 1 if alpha else 8
        collision = picmi.ProtonImpactIonizationCollisions(
            name=name,
            species=beam,
            product_species=[electrons, ions],
            ionization_target=target,
            background_density=f"{density:.17g}*(t < 1.5e-9)"
            if mode == "pulse"
            else density,
            background_temperature=0.0,
            fixed_product_weight=weight,
            max_products_per_cell=cap,
        )
        collisions.append(collision)
        cases[name] = dict(
            beam=beam, electrons=electrons, ions=ions, mode=mode, weight=weight, cap=cap
        )

sim = picmi.Simulation(
    solver=solver,
    time_step_size=DT,
    max_steps=STEPS,
    warpx_collisions=collisions,
    warpx_random_seed=42,
    verbose=0,
)
for case in cases.values():
    for key in ("beam", "electrons", "ions"):
        sim.add_species(
            case[key],
            layout=picmi.GriddedLayout(
                n_macroparticle_per_cell=[2 if key == "beam" else 0] * 3, grid=grid
            ),
        )
sim.initialize_inputs()
sim.initialize_warpx()


def host(values):
    return values.get() if hasattr(values, "get") else np.asarray(values)


def weights(name):
    arrays = [host(tile["w"]) for tile in sim.particles.get(name).iterator(level=0)]
    return np.concatenate(arrays) if arrays else np.empty(0)


previous = {name: (0.0, 0.0, 0) for name in cases}
for step in range(1, STEPS + 1):
    sim.step(1)
    for name, case in cases.items():
        ew = weights(f"electrons_{name}")
        iw = weights(f"ions_{name}")
        np.testing.assert_array_equal(ew, iw)
        remainder = float(
            host(sim.fields.get(f"{name}_product_weight_remainder", level=0)[...]).sum()
        )
        budget = float(ew.sum()) + remainder
        active_steps = min(step, 2) if case["mode"] == "pulse" else step
        expected = 0.37 * active_steps * (4 if case["mode"] == "alpha" else 1)
        # Alpha/proton rates at fixed speed differ slightly beyond Z^2 through
        # the exact finite-projectile-mass endpoints and Bhabha factor.
        assert np.isclose(budget, expected, rtol=2e-3), (name, step, budget, expected)
        assert 0 <= remainder < case["weight"]
        assert len(ew) - previous[name][2] <= case["cap"]
        if case["mode"] == "alpha":
            assert len(ew) == step and np.all(ew > case["weight"]) and remainder == 0
        else:
            assert np.all(ew == case["weight"])
        if case["mode"] == "pulse" and step > 2:
            assert (float(ew.sum()), remainder, len(ew)) == previous[name]
        previous[name] = (float(ew.sum()), remainder, len(ew))
        # The background is cold and the beam remains rigid.
        for tile in sim.particles.get(f"ions_{name}").iterator(level=0):
            for component in ("ux", "uy", "uz"):
                assert np.all(host(tile[component]) == 0)
        for tile in sim.particles.get(f"beam_{name}").iterator(level=0):
            assert np.all(host(tile["uz"]) == PROPER_SPEED)

for target in ("N2", "O2"):
    coarse = sum(previous[f"{target}_coarse"][:2])
    fine = sum(previous[f"{target}_fine"][:2])
    assert np.isclose(coarse, fine, rtol=2e-14)
print(
    "PASS: both-gas source budgets, fractional carry, zero-density pause, weight/cap independence and Z^2 scaling"
)
sim.finalize()
