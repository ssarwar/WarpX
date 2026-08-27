/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCIonization.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Gpu.H>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>

namespace
{
struct RBEQShell
{
    double binding_energy;
    double kinetic_energy;
    double occupation;
    double oscillator_strength_ratio;
};

constexpr std::array<RBEQShell, 5> n2_shells{{{409.50, 603.30, 4.0, 1.000},
                                              {37.30, 71.13, 2.0, 0.760},
                                              {18.72, 63.18, 2.0, 1.000},
                                              {16.74, 44.30, 4.0, 0.938},
                                              {15.58, 54.91, 2.0, 0.792}}};

constexpr std::array<RBEQShell, 6> o2_shells{{{543.80, 796.20, 4.0, 1.0000},
                                              {40.33, 79.73, 2.0, 0.9600},
                                              {27.05, 90.92, 2.0, 1.0000},
                                              {20.30, 71.84, 2.0, 1.0000},
                                              {17.08, 59.89, 4.0, 1.0000},
                                              {12.07, 84.88, 2.0, 0.9314}}};

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
RBEQTerms rbeqTerms (double const incident_energy, RBEQShell const& shell)
{
    if (incident_energy <= shell.binding_energy)
    {
        return {};
    }

    constexpr double electron_rest_energy = 510998.95069; // eV
    constexpr double log_two = 0.69314718055994530942;

    auto const t = incident_energy / electron_rest_energy;
    auto const b = shell.binding_energy / electron_rest_energy;
    auto const u = shell.kinetic_energy / electron_rest_energy;
    auto const energy_ratio = incident_energy / shell.binding_energy;
    auto const gamma_tilde = 1.0 + t + u + b;
    auto const beta_tilde_sq =
        (gamma_tilde - 1.0) * (gamma_tilde + 1.0) / (gamma_tilde * gamma_tilde);
    auto const q = shell.oscillator_strength_ratio;
    auto const dipole_correction = -(1.0 + q - (5.0 - 3.0 * q) * log_two) / q;

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

double rbeqCumulative (double const secondary_energy, RBEQShell const& shell,
                       RBEQTerms const& terms)
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

template <std::size_t N>
void initializeTable (std::array<RBEQShell, N> const& shells, double const energy_min,
                      double const log_energy_step,
                      amrex::Gpu::HostVector<amrex::ParticleReal>& shell_probabilities,
                      amrex::Gpu::HostVector<amrex::ParticleReal>& inverse_cdf)
{
    constexpr int inverse_iterations = 36;
    constexpr double uniform_threshold_ratio = 1.0e-3;
    constexpr int max_shells = BackgroundMCCIonizationModel::max_shell_count;
    constexpr int energy_count = BackgroundMCCIonizationModel::energy_grid_size;
    constexpr int quantile_count = BackgroundMCCIonizationModel::quantile_grid_size;

    for (int energy_index = 0; energy_index < energy_count; ++energy_index)
    {
        auto const incident_energy =
            energy_min * std::exp(log_energy_step * static_cast<double>(energy_index));

        std::array<RBEQTerms, N> terms;
        double total_cross_section = 0.0;
        for (int shell_index = 0; shell_index < static_cast<int>(N); ++shell_index)
        {
            terms[shell_index] = rbeqTerms(incident_energy, shells[shell_index]);
            // Some Q=1 partial cross sections are negative immediately above
            // their thresholds because the RBEQ dipole correction is too
            // strong. A negative partial cannot be sampled; it enters only once
            // positive.
            auto const partial = std::max(terms[shell_index].total, 0.0);
            shell_probabilities[energy_index * max_shells + shell_index] =
                static_cast<amrex::ParticleReal>(partial);
            total_cross_section += partial;
        }

        if (total_cross_section > 0.0)
        {
            for (int shell_index = 0; shell_index < static_cast<int>(N); ++shell_index)
            {
                shell_probabilities[energy_index * max_shells + shell_index] /=
                    static_cast<amrex::ParticleReal>(total_cross_section);
            }
        }

        for (int shell_index = 0; shell_index < static_cast<int>(N); ++shell_index)
        {
            auto const& shell = shells[shell_index];
            auto const available_energy = std::max(incident_energy - shell.binding_energy, 0.0);
            auto const maximum_secondary_energy = 0.5 * available_energy;
            auto const use_uniform =
                maximum_secondary_energy <= 0.0 || terms[shell_index].total <= 0.0 ||
                available_energy / shell.binding_energy < uniform_threshold_ratio;

            for (int quantile_index = 0; quantile_index < quantile_count; ++quantile_index)
            {
                auto const coordinate =
                    static_cast<double>(quantile_index) / static_cast<double>(quantile_count - 1);
                auto const coordinate_fourth =
                    coordinate * coordinate * coordinate * coordinate;
                auto const complement = 1.0 - coordinate;
                auto const complement_fourth =
                    complement * complement * complement * complement;
                auto const quantile =
                    coordinate_fourth / (coordinate_fourth + complement_fourth);
                double fraction = quantile;

                if (!use_uniform && quantile > 0.0 && quantile < 1.0)
                {
                    auto const target = quantile * terms[shell_index].total;
                    double lower = 0.0;
                    double upper = maximum_secondary_energy;
                    for (int iteration = 0; iteration < inverse_iterations; ++iteration)
                    {
                        auto const midpoint = 0.5 * (lower + upper);
                        if (rbeqCumulative(midpoint, shell, terms[shell_index]) < target)
                        {
                            lower = midpoint;
                        }
                        else
                        {
                            upper = midpoint;
                        }
                    }
                    fraction = 0.5 * (lower + upper) / maximum_secondary_energy;
                }

                auto const index =
                    (energy_index * max_shells + shell_index) * quantile_count + quantile_index;
                inverse_cdf[index] = static_cast<amrex::ParticleReal>(fraction);
            }
        }
    }
}
} // namespace

BackgroundMCCIonizationModel::BackgroundMCCIonizationModel (
    BackgroundMCCIonizationTarget const target, amrex::ParticleReal const maximum_energy)
{
    using namespace amrex::literals;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(target == BackgroundMCCIonizationTarget::N2 ||
                                         target == BackgroundMCCIonizationTarget::O2,
                                     "RBEQ energy sharing supports only N2 and O2 targets.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(static_cast<double>(maximum_energy)) &&
                                         maximum_energy > 0.0_prt,
                                     "The RBEQ maximum energy must be finite and positive.");

    auto const energy_min = static_cast<double>(outerBindingEnergy(target)) * (1.0 + 1.0e-6);
    // RBEQ is analytic at high energy. Extending the initialization grid to 1
    // GeV avoids a device-side extrapolation in the intended subthermal-to-MeV
    // range.
    auto const energy_max = std::max(static_cast<double>(maximum_energy), 1.0e9);
    auto const log_energy_step =
        std::log(energy_max / energy_min) / static_cast<double>(energy_grid_size - 1);

    m_shell_probabilities_h.assign(energy_grid_size * max_shell_count, 0.0_prt);
    m_inverse_cdf_h.assign(energy_grid_size * max_shell_count * quantile_grid_size, 0.0_prt);

    m_executor_h.m_model = IonizationEnergySharingModel::RBEQ;
    m_executor_h.m_energy_grid_size = energy_grid_size;
    m_executor_h.m_quantile_grid_size = quantile_grid_size;
    m_executor_h.m_energy_min = static_cast<amrex::ParticleReal>(energy_min);
    m_executor_h.m_log_energy_min = static_cast<amrex::ParticleReal>(std::log(energy_min));
    m_executor_h.m_inverse_log_energy_step =
        static_cast<amrex::ParticleReal>(1.0 / log_energy_step);

    if (target == BackgroundMCCIonizationTarget::N2)
    {
        m_executor_h.m_shell_count = static_cast<int>(n2_shells.size());
        for (int i = 0; i < m_executor_h.m_shell_count; ++i)
        {
            m_executor_h.m_binding_energies[i] =
                static_cast<amrex::ParticleReal>(n2_shells[i].binding_energy);
        }
        initializeTable(n2_shells, energy_min, log_energy_step, m_shell_probabilities_h,
                        m_inverse_cdf_h);
    }
    else
    {
        m_executor_h.m_shell_count = static_cast<int>(o2_shells.size());
        for (int i = 0; i < m_executor_h.m_shell_count; ++i)
        {
            m_executor_h.m_binding_energies[i] =
                static_cast<amrex::ParticleReal>(o2_shells[i].binding_energy);
        }
        initializeTable(o2_shells, energy_min, log_energy_step, m_shell_probabilities_h,
                        m_inverse_cdf_h);
    }

    m_executor_h.m_shell_probabilities = m_shell_probabilities_h.data();
    m_executor_h.m_inverse_cdf = m_inverse_cdf_h.data();

#ifdef AMREX_USE_GPU
    m_shell_probabilities_d.resize(m_shell_probabilities_h.size());
    m_inverse_cdf_d.resize(m_inverse_cdf_h.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_shell_probabilities_h.begin(),
                          m_shell_probabilities_h.end(), m_shell_probabilities_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_inverse_cdf_h.begin(), m_inverse_cdf_h.end(),
                          m_inverse_cdf_d.begin());

    m_executor_d = m_executor_h;
    m_executor_d.m_shell_probabilities = m_shell_probabilities_d.data();
    m_executor_d.m_inverse_cdf = m_inverse_cdf_d.data();
    amrex::Gpu::streamSynchronize();
#endif
}

BackgroundMCCIonizationTarget BackgroundMCCIonizationModel::parseTarget (std::string const& target)
{
    auto normalized = target;
    std::transform(
        normalized.begin(), normalized.end(), normalized.begin(),
        [] (unsigned char const character) { return static_cast<char>(std::toupper(character)); });

    if (normalized == "N2")
    {
        return BackgroundMCCIonizationTarget::N2;
    }
    if (normalized == "O2")
    {
        return BackgroundMCCIonizationTarget::O2;
    }
    return BackgroundMCCIonizationTarget::None;
}

amrex::ParticleReal
BackgroundMCCIonizationModel::outerBindingEnergy (BackgroundMCCIonizationTarget const target)
{
    if (target == BackgroundMCCIonizationTarget::N2)
    {
        return static_cast<amrex::ParticleReal>(n2_shells.back().binding_energy);
    }
    if (target == BackgroundMCCIonizationTarget::O2)
    {
        return static_cast<amrex::ParticleReal>(o2_shells.back().binding_energy);
    }
    return 0.0;
}
