/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_RBEQ_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_RBEQ_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

/** Host-only source definitions shared by the sampler and table exporter. */
namespace BackgroundMCCRBEQ
{
enum class Model { IAAThesis2023, ElmolcsB8643810 };

inline std::string
name (Model const model)
{
    return model == Model::IAAThesis2023 ? "iaa_thesis_2023" : "elmolcs_b8643810";
}

inline Model
parse (std::string const& value)
{
    if (value == "iaa_thesis_2023") {
        return Model::IAAThesis2023;
    }
    if (value == "elmolcs_b8643810") {
        return Model::ElmolcsB8643810;
    }
    throw std::invalid_argument("Unknown RBEQ model: " + value);
}

struct RBEQShell
{
    double binding_energy;
    double kinetic_energy;
    double occupation;
    double oscillator_strength_ratio;
    double dipole_correction = 0.0;
};

inline std::vector<RBEQShell>
shells (bool const nitrogen, Model const model)
{
    std::vector<RBEQShell> result;
    if (nitrogen && model == Model::IAAThesis2023) {
        result = {{409.50, 603.30, 4, 1},
                  {37.30, 71.13, 2, .760},
                  {18.72, 63.18, 2, 1},
                  {16.74, 44.30, 4, .938},
                  {15.58, 54.91, 2, .792}};
    } else if (nitrogen) {
        result = {{409.50, 603.30, 4, 1}, {37.8, 71.13, 2, .760},   {23.6, 63.18, 1, 1},
                  {18.746, 63.18, 1, 1},  {16.716, 44.30, 4, .938}, {15.58, 54.91, 2, .792}};
    } else {
        result = {{543.80, 796.20, 4, 1}, {40.33, 79.73, 2, .9600}, {27.05, 90.92, 2, 1},
                  {20.30, 71.84, 2, 1},   {17.08, 59.89, 4, 1},     {12.07, 84.88, 2, .9314}};
        if (model == Model::ElmolcsB8643810) {
            result[0].binding_energy = 531;
        }
    }
    for (auto& shell : result) {
        auto const q = shell.oscillator_strength_ratio;
        shell.dipole_correction = -(1 + q - (5 - 3 * q) * std::log(2.0)) / q;
        if (model == Model::ElmolcsB8643810) {
            // Literal Cd(1 MeV,B,Q) from the archived lnBm expression.
            // Algebraically cancel log(B/Ryd) to improve numerical stability.
            auto const t = 1.0e6 / shell.binding_energy;
            auto const norm = 1 + .5 / t * (1 - 1 / t) - 2 / (t + 1);
            auto const moment = (1 + .5 / t - 3 / (t + 1) - 2 * std::log((t + 1) / 2) / (t + 1) +
                                 std::log(t) / (4 * (t + 1))) /
                                norm;
            shell.dipole_correction = -2 * moment + 1 - 1 / q + (5 / q - 3) * std::log(2.0);
        }
    }
    return result;
}

// CODATA constants in SI; normalization is 4*pi*a0^2*alpha^4 in m^2.
constexpr double normalization = 4 * 3.14159265358979323846 * 5.29177210544e-11 *
                                 5.29177210544e-11 * 7.2973525643e-3 * 7.2973525643e-3 *
                                 7.2973525643e-3 * 7.2973525643e-3;

struct RBEQTerms
{
    double prefactor = 0.0;
    double incident_to_binding = 0.0;
    double binding_sq = 0.0;
    double exchange = 0.0;
    double bethe_log = 0.0;
    double total = 0.0;
};

/** Evaluate the dimensionless terms in the RBEQ* model.
 *
 * This is Schmalzried, Eqs. (11.119)--(11.121), with the corrected
 * high-energy dipole coefficient. The common dimensional normalization
 * cancels in shell selection and conditional energy sampling.
 */
inline RBEQTerms
rbeqTerms (double const incident_energy, RBEQShell const& shell)
{
    if (incident_energy <= shell.binding_energy) {
        return {};
    }

    constexpr double electron_rest_energy = 510998.95069; // eV

    auto const t = incident_energy / electron_rest_energy;
    auto const b = shell.binding_energy / electron_rest_energy;
    auto const u = shell.kinetic_energy / electron_rest_energy;
    auto const energy_ratio = incident_energy / shell.binding_energy;
    auto const gamma_tilde = 1.0 + t + u + b;
    auto const beta_tilde_sq =
        (gamma_tilde - 1.0) * (gamma_tilde + 1.0) / (gamma_tilde * gamma_tilde);
    auto const q = shell.oscillator_strength_ratio;
    auto const dipole_correction = shell.dipole_correction;

    RBEQTerms result;
    result.prefactor = shell.occupation / (2.0 * b * beta_tilde_sq);
    result.incident_to_binding = energy_ratio;
    result.binding_sq = (b / gamma_tilde) * (b / gamma_tilde);
    result.exchange =
        (2.0 * gamma_tilde - 1.0) / ((1.0 + energy_ratio) * gamma_tilde * gamma_tilde);
    result.bethe_log = std::log(t * (t + 2.0)) - t * (t + 2.0) / ((1.0 + t) * (1.0 + t)) -
                       std::log(2.0 * b) + dipole_correction;

    result.total =
        result.prefactor *
        (0.5 * q * result.bethe_log * (1.0 - 1.0 / (energy_ratio * energy_ratio)) +
         (2.0 - q) * (1.0 - 1.0 / energy_ratio - std::log(energy_ratio) * result.exchange +
                      0.5 * result.binding_sq * (energy_ratio - 1.0)));
    return result;
}

inline double
rbeqCumulative (double const secondary_energy, RBEQShell const& shell, RBEQTerms const& terms)
{
    auto const w = secondary_energy / shell.binding_energy;
    auto const t = terms.incident_to_binding;
    auto const q = shell.oscillator_strength_ratio;
    auto const c1 =
        0.5 * terms.bethe_log * q *
        (1.0 - 1.0 / ((w + 1.0) * (w + 1.0)) + 1.0 / ((t - w) * (t - w)) - 1.0 / (t * t));
    auto const c2 = 1.0 - 1.0 / (w + 1.0) + 1.0 / (t - w) - 1.0 / t + w * terms.binding_sq;
    auto const c3 = std::log((w + 1.0) * t / (t - w)) * terms.exchange;
    return terms.prefactor * (c1 + (2.0 - q) * (c2 - c3));
}

inline bool
hasNonnegativeDifferential (double const incident_energy, RBEQShell const& shell)
{
    auto const terms = rbeqTerms(incident_energy, shell);
    if (terms.total <= 0.0) {
        return false;
    }

    // Set y = 1/(w+1) + 1/(t-w). The differential RBEQ shape is the cubic
    // A*y^3+B*y^2+C*y+D, so its exact minimum on the symmetric half-domain is
    // attained at an endpoint or one of at most two stationary points.
    auto const t = terms.incident_to_binding;
    auto const sum = t + 1.0;
    auto const q = shell.oscillator_strength_ratio;
    auto const a = q * terms.bethe_log;
    auto const b = (2.0 - q) - 3.0 * q * terms.bethe_log / sum;
    auto const c = -(2.0 - q) * (2.0 / sum + terms.exchange);
    auto const d = (2.0 - q) * terms.binding_sq;
    auto const evaluate = [=] (double const y) { return ((a * y + b) * y + c) * y + d; };

    auto const lower = 4.0 / sum;
    auto const upper = sum / t;
    auto minimum = std::min(evaluate(lower), evaluate(upper));
    auto const evaluate_stationary_point = [&] (double const y) {
        if (y > lower && y < upper) {
            minimum = std::min(minimum, evaluate(y));
        }
    };

    auto const derivative_a = 3.0 * a;
    auto const derivative_b = 2.0 * b;
    auto const derivative_c = c;
    auto const coefficient_scale =
        std::max({std::abs(derivative_a), std::abs(derivative_b), std::abs(derivative_c)});
    if (std::abs(derivative_a) <= std::numeric_limits<double>::epsilon() * coefficient_scale) {
        if (std::abs(derivative_b) > std::numeric_limits<double>::epsilon() * coefficient_scale) {
            evaluate_stationary_point(-derivative_c / derivative_b);
        }
    } else {
        auto const discriminant = derivative_b * derivative_b - 4.0 * derivative_a * derivative_c;
        if (discriminant >= 0.0) {
            auto const root = std::sqrt(discriminant);
            evaluate_stationary_point((-derivative_b - root) / (2.0 * derivative_a));
            evaluate_stationary_point((-derivative_b + root) / (2.0 * derivative_a));
        }
    }
    return minimum >= 0.0;
}

template <typename Predicate>
double
findPositiveThreshold (RBEQShell const& shell, Predicate&& predicate)
{
    auto lower = shell.binding_energy;
    auto upper = std::nextafter(lower, std::numeric_limits<double>::infinity());
    while (!predicate(upper)) {
        lower = upper;
        upper *= 1.1;
        if (!std::isfinite(upper) || upper > 1.0e9) {
            throw std::runtime_error("No physical RBEQ continuation below 1 GeV.");
        }
    }

    constexpr int bisection_iterations = 64;
    for (int iteration = 0; iteration < bisection_iterations; ++iteration) {
        auto const midpoint = 0.5 * (lower + upper);
        if (predicate(midpoint)) {
            upper = midpoint;
        } else {
            lower = midpoint;
        }
    }
    return upper;
}

inline double
crossSection (double const energy, std::vector<RBEQShell> const& parameters,
              bool const positive_part = true)
{
    double total = 0;
    for (auto const& shell : parameters) {
        auto const partial = rbeqTerms(energy, shell).total;
        total += positive_part ? std::max(partial, 0.0) : partial;
    }
    return normalization * total;
}
} // namespace BackgroundMCCRBEQ
#endif
