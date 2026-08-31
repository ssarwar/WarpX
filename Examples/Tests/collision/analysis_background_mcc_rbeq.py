#!/usr/bin/env python3
"""Compare sampled N2/O2 ionization energies with independent RBEQ moments."""

import math

import numpy as np

PARTICLE_COUNT = 32768
MC2_EV = 510998.95069
LN2 = math.log(2.0)
MAJORANT_MARGIN = 1.001

CASES = [
    ("n2_near", "N2", 15.78),
    ("n2_partial_guard", "N2", 19.50),
    ("n2_mid", "N2", 100.0),
    ("n2_high", "N2", 1000.0),
    ("n2_gev", "N2", 1.0e9),
    ("o2_near", "O2", 12.27),
    ("o2_sdcs_guard", "O2", 18.702),
    ("o2_partial_guard", "O2", 21.50),
    ("o2_mid", "O2", 100.0),
    ("o2_high", "O2", 1000.0),
    ("o2_gev", "O2", 1.0e9),
]

SHELLS = {
    "N2": np.array(
        [
            [409.50, 603.30, 4.0, 1.000],
            [37.30, 71.13, 2.0, 0.760],
            [18.72, 63.18, 2.0, 1.000],
            [16.74, 44.30, 4.0, 0.938],
            [15.58, 54.91, 2.0, 0.792],
        ]
    ),
    "O2": np.array(
        [
            [543.80, 796.20, 4.0, 1.0000],
            [40.33, 79.73, 2.0, 0.9600],
            [27.05, 90.92, 2.0, 1.0000],
            [20.30, 71.84, 2.0, 1.0000],
            [17.08, 59.89, 4.0, 1.0000],
            [12.07, 84.88, 2.0, 0.9314],
        ]
    ),
}


def rbeq_terms(energy, shell):
    binding, kinetic, occupation, q = shell
    if energy <= binding:
        return None
    t = energy / MC2_EV
    b = binding / MC2_EV
    u = kinetic / MC2_EV
    ratio = energy / binding
    gamma_tilde = 1.0 + t + u + b
    beta_tilde_sq = 1.0 - 1.0 / gamma_tilde**2
    correction = -(1.0 + q - (5.0 - 3.0 * q) * LN2) / q
    prefactor = occupation / (2.0 * b * beta_tilde_sq)
    binding_sq = (b / gamma_tilde) ** 2
    exchange = (2.0 * gamma_tilde - 1.0) / ((1.0 + ratio) * gamma_tilde**2)
    beta_sq = t * (t + 2.0) / (1.0 + t) ** 2
    bethe_log = math.log(t * (t + 2.0)) - beta_sq - math.log(2.0 * b)
    bethe_log += correction
    total = prefactor * (
        0.5 * q * bethe_log * (1.0 - 1.0 / ratio**2)
        + (2.0 - q)
        * (
            1.0
            - 1.0 / ratio
            - math.log(ratio) * exchange
            + 0.5 * binding_sq * (ratio - 1.0)
        )
    )
    return prefactor, ratio, binding_sq, exchange, bethe_log, total


def differential_shape(secondary_energy, shell, terms):
    binding, _, _, q = shell
    _, ratio, binding_sq, exchange_coefficient, bethe_log, _ = terms
    w = np.asarray(secondary_energy) / binding
    forward = w + 1.0
    exchange = ratio - w
    return q * bethe_log * (forward**-3 + exchange**-3) + (2.0 - q) * (
        forward**-2
        + exchange**-2
        + binding_sq
        - exchange_coefficient * (forward**-1 + exchange**-1)
    )


def has_nonnegative_differential(energy, shell):
    terms = rbeq_terms(energy, shell)
    if terms is None or terms[-1] <= 0.0:
        return False
    maximum = 0.5 * (energy - shell[0])
    grid = np.linspace(0.0, maximum, 16385)
    return np.min(differential_shape(grid, shell, terms)) >= 0.0


def find_threshold(shell, predicate):
    lower = shell[0]
    upper = np.nextafter(lower, math.inf)
    while not predicate(upper, shell):
        lower = upper
        upper *= 1.1
        assert upper <= 1.0e9
    for _ in range(64):
        midpoint = 0.5 * (lower + upper)
        if predicate(midpoint, shell):
            upper = midpoint
        else:
            lower = midpoint
    return upper


def positive_partial(energy, shell):
    terms = rbeq_terms(energy, shell)
    return terms is not None and terms[-1] > 0.0


THRESHOLDS = {}
for target, shells in SHELLS.items():
    target_thresholds = []
    for shell in shells:
        positive = find_threshold(shell, positive_partial)
        nonnegative = find_threshold(shell, has_nonnegative_differential)
        uniform = max(nonnegative, shell[0] * (1.0 + 1.0e-3))
        target_thresholds.append((positive, uniform))
    THRESHOLDS[target] = target_thresholds


def conditional_cdf(
    secondary_energy, incident_energy, shell, terms, uniform_threshold
):
    binding, _, _, q = shell
    prefactor, ratio, binding_sq, exchange, bethe_log, total = terms
    if total <= 0.0:
        return 0.0
    maximum = 0.5 * (incident_energy - binding)
    values = np.asarray(secondary_energy)
    if incident_energy <= uniform_threshold:
        return np.clip(values / maximum, 0.0, 1.0)
    w = np.maximum(values, 0.0) / binding
    c1 = (
        0.5
        * bethe_log
        * q
        * (1.0 - 1.0 / (w + 1.0) ** 2 + 1.0 / (ratio - w) ** 2 - 1.0 / ratio**2)
    )
    c2 = 1.0 - 1.0 / (w + 1.0) + 1.0 / (ratio - w) - 1.0 / ratio + w * binding_sq
    c3 = np.log((w + 1.0) * ratio / (ratio - w)) * exchange
    cdf = np.clip(prefactor * (c1 + (2.0 - q) * (c2 - c3)) / total, 0.0, 1.0)
    return np.where(values >= maximum, 1.0, cdf)


def expected_statistics(target, incident_energy, cdf_probes):
    shells = SHELLS[target]
    terms = [rbeq_terms(incident_energy, shell) for shell in shells]
    partials = np.array(
        [max(term[-1], 0.0) if term is not None else 0.0 for term in terms]
    )
    probabilities = partials / np.sum(partials)

    binding_mean = np.sum(probabilities * shells[:, 0])
    binding_variance = np.sum(probabilities * (shells[:, 0] - binding_mean) ** 2)
    secondary_mean = 0.0
    secondary_second_moment = 0.0
    cdf_at_probes = np.zeros(len(cdf_probes))
    primary_cosine_mean = 0.0
    primary_cosine_second_moment = 0.0
    secondary_cosine_mean = 0.0
    secondary_cosine_second_moment = 0.0
    for probability, shell, term, (_, uniform_threshold) in zip(
        probabilities, shells, terms, THRESHOLDS[target]
    ):
        if probability == 0.0:
            continue
        maximum = 0.5 * (incident_energy - shell[0])
        # A logarithmic grid resolves the O(binding-energy) peak and the long
        # relativistic tail simultaneously, including for a 1 GeV projectile.
        positive_grid = np.geomspace(
            max(maximum * 1.0e-14, np.finfo(float).tiny), maximum, 40001
        )
        grid = np.concatenate(([0.0], positive_grid))
        cdf = conditional_cdf(
            grid, incident_energy, shell, term, uniform_threshold
        )
        secondary_mean += probability * np.trapezoid(1.0 - cdf, grid)
        secondary_second_moment += probability * np.trapezoid(
            2.0 * grid * (1.0 - cdf), grid
        )
        cdf_at_probes += probability * conditional_cdf(
            cdf_probes, incident_energy, shell, term, uniform_threshold
        )
        pdf = np.maximum(np.gradient(cdf, grid, edge_order=2), 0.0)
        pdf /= np.trapezoid(pdf, grid)
        available = incident_energy - shell[0]
        primary_energy = available - grid
        primary_cosine = np.sqrt(
            primary_energy
            * (available + 2.0 * MC2_EV)
            / (available * (primary_energy + 2.0 * MC2_EV))
        )
        denominator = grid + shell[0]
        free_component = np.sqrt(grid / incident_energy) * (
            grid + 0.5 * shell[0]
        ) / denominator
        bound_weight = shell[0] / denominator

        # Eq. (11.132) is uniform on [a-b, a+b], followed by the same physical
        # cosine clamp as the device implementation. Integrate the clipped
        # interval analytically so the near-threshold test is sensitive to the
        # exact IAA expression rather than only to its unbounded first moment.
        lower = free_component - bound_weight
        upper = free_component + bound_weight
        clipped = upper > 1.0
        secondary_conditional_mean = free_component.copy()
        secondary_conditional_second_moment = (
            free_component**2 + bound_weight**2 / 3.0
        )
        secondary_conditional_mean[clipped] = (
            0.5 * (1.0 - lower[clipped] ** 2) + upper[clipped] - 1.0
        ) / (2.0 * bound_weight[clipped])
        secondary_conditional_second_moment[clipped] = (
            (1.0 - lower[clipped] ** 3) / 3.0 + upper[clipped] - 1.0
        ) / (2.0 * bound_weight[clipped])
        primary_cosine_mean += probability * np.trapezoid(primary_cosine * pdf, grid)
        primary_cosine_second_moment += probability * np.trapezoid(
            primary_cosine**2 * pdf, grid
        )
        secondary_cosine_mean += probability * np.trapezoid(
            secondary_conditional_mean * pdf, grid
        )
        secondary_cosine_second_moment += probability * np.trapezoid(
            secondary_conditional_second_moment * pdf, grid
        )
    secondary_variance = secondary_second_moment - secondary_mean**2
    primary_cosine_variance = primary_cosine_second_moment - primary_cosine_mean**2
    secondary_cosine_variance = (
        secondary_cosine_second_moment - secondary_cosine_mean**2
    )
    return (
        probabilities,
        binding_mean,
        binding_variance,
        secondary_mean,
        secondary_variance,
        cdf_at_probes,
        primary_cosine_mean,
        primary_cosine_variance,
        secondary_cosine_mean,
        secondary_cosine_variance,
    )


data = np.load("background_mcc_rbeq_results.npz")
particle_real_bytes = int(data["particle_real_bytes"])
for name, target, incident_energy in CASES:
    particle_count = int(data[f"{name}_particle_count"])
    event_count = int(data[f"{name}_event_count"])
    secondary = data[f"{name}_secondary_energies"]
    binding_losses = data[f"{name}_binding_losses"]
    primary_cosines = data[f"{name}_primary_cosines"]
    secondary_cosines = data[f"{name}_secondary_cosines"]
    momentum_residual = data[f"{name}_momentum_residual"]
    observed_binding_mean = float(data[f"{name}_binding_mean"])
    outer_binding = np.min(SHELLS[target][:, 0])
    cdf_probes = [min(0.1 * (incident_energy - outer_binding), 10.0 * outer_binding)]
    if incident_energy >= 1.0e8:
        cdf_probes += [1.0e3, 1.0e6]
    cdf_probes = np.array(cdf_probes)
    (
        expected_shell_probabilities,
        expected_binding_mean,
        expected_binding_variance,
        expected_secondary_mean,
        expected_secondary_variance,
        expected_cdfs,
        expected_primary_cosine,
        expected_primary_cosine_variance,
        expected_secondary_cosine,
        expected_secondary_cosine_variance,
    ) = expected_statistics(target, incident_energy, cdf_probes)

    expected_event_fraction = -math.expm1(-5.0) / MAJORANT_MARGIN
    event_standard_error = math.sqrt(
        expected_event_fraction * (1.0 - expected_event_fraction) / particle_count
    )
    assert abs(event_count / particle_count - expected_event_fraction) < (
        7.0 * event_standard_error
    )
    assert int(data[f"{name}_ion_count"]) == event_count
    assert primary_cosines.size == event_count
    assert secondary_cosines.size == event_count
    assert np.isfinite(secondary).all()
    assert np.all(secondary >= 0.0)
    assert np.max(secondary) <= 0.5 * (incident_energy - outer_binding)

    if incident_energy < 1.0e8:
        shell_distances = np.abs(binding_losses[:, None] - SHELLS[target][:, 0])
        sampled_shells = np.argmin(shell_distances, axis=1)
        # Event-by-event energy closure must recover one of the discrete target
        # binding energies, including in single-particle-precision builds.
        assert np.max(np.min(shell_distances, axis=1)) < 5.0e-2
        observed_shell_probabilities = np.bincount(
            sampled_shells, minlength=len(SHELLS[target])
        ) / event_count
        shell_standard_errors = np.sqrt(
            expected_shell_probabilities
            * (1.0 - expected_shell_probabilities)
            / event_count
        )
        assert np.all(
            np.abs(observed_shell_probabilities - expected_shell_probabilities)
            < 7.0 * shell_standard_errors + 1.0e-2
        )

        # The guarded shell has an analytically negative partial at these
        # energies and therefore must never be selected by grid interpolation.
        forbidden_shell = {"n2_partial_guard": 2, "o2_partial_guard": 3}.get(name)
        if forbidden_shell is not None:
            assert np.count_nonzero(sampled_shells == forbidden_shell) == 0

        # Isolate a shell whose total partial is positive while its published
        # SDCS is locally negative. Its conservative continuation must be
        # uniform rather than an inversion of a non-monotone CDF.
        guarded_shell = {"o2_sdcs_guard": 4}.get(name)
        if guarded_shell is not None:
            guarded = secondary[sampled_shells == guarded_shell]
            assert guarded.size >= 20
            maximum = 0.5 * (
                incident_energy - SHELLS[target][guarded_shell, 0]
            )
            uniform_samples = np.sort(guarded / maximum)
            empirical_upper = np.arange(1, guarded.size + 1) / guarded.size
            empirical_lower = np.arange(guarded.size) / guarded.size
            ks_distance = max(
                np.max(empirical_upper - uniform_samples),
                np.max(uniform_samples - empirical_lower),
            )
            assert ks_distance < 2.0 / math.sqrt(guarded.size) + 2.0e-2

    if particle_real_bytes > 4 or incident_energy < 1.0e8:
        binding_standard_error = math.sqrt(expected_binding_variance / event_count)
        assert abs(observed_binding_mean - expected_binding_mean) < (
            7.0 * binding_standard_error + 5.0e-3 * max(1.0, expected_binding_mean)
        )
    observed_secondary_mean = float(np.mean(secondary))
    if incident_energy < 1.0e8:
        secondary_standard_error = math.sqrt(expected_secondary_variance / event_count)
        secondary_tolerance = 7.0 * secondary_standard_error + 5.0e-3 * max(
            1.0, expected_secondary_mean
        )
        assert (
            abs(observed_secondary_mean - expected_secondary_mean) < secondary_tolerance
        ), (
            f"{name} secondary mean {observed_secondary_mean} differs from "
            f"{expected_secondary_mean} by more than {secondary_tolerance}"
        )

    observed_cdfs = np.array([np.mean(secondary <= probe) for probe in cdf_probes])
    for probe_index, (observed_cdf, expected_cdf) in enumerate(
        zip(observed_cdfs, expected_cdfs)
    ):
        cdf_standard_error = math.sqrt(
            expected_cdf * (1.0 - expected_cdf) / event_count
        )
        interpolation_tolerance = 4.0e-3 if probe_index == 0 else 1.0e-3
        assert abs(observed_cdf - expected_cdf) < (
            7.0 * cdf_standard_error + interpolation_tolerance
        )

    assert np.isfinite(primary_cosines).all()
    assert np.isfinite(secondary_cosines).all()
    assert np.all(np.abs(primary_cosines) <= 1.0)
    assert np.all(np.abs(secondary_cosines) <= 1.0)
    primary_cosine_standard_error = math.sqrt(
        expected_primary_cosine_variance / event_count
    )
    secondary_cosine_standard_error = math.sqrt(
        expected_secondary_cosine_variance / event_count
    )
    assert abs(np.mean(primary_cosines) - expected_primary_cosine) < (
        7.0 * primary_cosine_standard_error + 5.0e-3
    )
    assert abs(np.mean(secondary_cosines) - expected_secondary_cosine) < (
        7.0 * secondary_cosine_standard_error + 2.0e-3
    )
    assert np.linalg.norm(momentum_residual) < 2.0e-5

    print(
        f"{name}: events={event_count}, binding mean="
        f"{observed_binding_mean:.5f}/{expected_binding_mean:.5f} eV, "
        f"secondary mean={observed_secondary_mean:.5f}/{expected_secondary_mean:.5f} eV, "
        f"CDF={np.array2string(observed_cdfs, precision=5)}/"
        f"{np.array2string(expected_cdfs, precision=5)}, "
        f"<cos1>={np.mean(primary_cosines):.5f}/{expected_primary_cosine:.5f}, "
        f"<cos2>={np.mean(secondary_cosines):.5f}/{expected_secondary_cosine:.5f}, "
        f"|dp|/p0={np.linalg.norm(momentum_residual):.3e}"
    )
