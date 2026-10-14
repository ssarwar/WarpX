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
using namespace BackgroundMCCRBEQ;

void
findShellThresholds (std::vector<RBEQShell> const& shells, std::vector<double>& positive_thresholds,
                     std::vector<double>& uniform_thresholds)
{
    constexpr double uniform_threshold_ratio = 1.0e-3;
    for (int shell_index = 0; shell_index < static_cast<int>(shells.size()); ++shell_index) {
        auto const& shell = shells[shell_index];
        positive_thresholds[shell_index] = findPositiveThreshold(
            shell, [&shell] (double const energy) { return rbeqTerms(energy, shell).total > 0.0; });
        auto const differential_threshold = findPositiveThreshold(
            shell, [&shell] (double const energy) {
                return hasNonnegativeDifferential(energy, shell);
            });
        uniform_thresholds[shell_index] = std::max(
            differential_threshold, shell.binding_energy * (1.0 + uniform_threshold_ratio));
    }
}

void
initializeShellCrossSections (std::vector<RBEQShell> const& shells, double const energy_min,
                              double const shell_log_energy_step,
                              amrex::Gpu::HostVector<amrex::ParticleReal>& shell_cross_sections)
{
    constexpr int max_shells = BackgroundMCCIonizationModel::max_shell_count;
    constexpr int shell_energy_count =
        BackgroundMCCIonizationModel::shell_energy_grid_size;
    for (int energy_index = 0; energy_index < shell_energy_count; ++energy_index)
    {
        auto const incident_energy = energy_min *
            std::exp(shell_log_energy_step * static_cast<double>(energy_index));
        for (int shell_index = 0; shell_index < static_cast<int>(shells.size()); ++shell_index) {
            // Some Q=1 partial cross sections are negative immediately above
            // their thresholds because the RBEQ dipole correction is too
            // strong. A negative partial cannot be sampled; it enters only once
            // positive.
            auto const partial =
                std::max(rbeqTerms(incident_energy, shells[shell_index]).total, 0.0);
            shell_cross_sections[energy_index * max_shells + shell_index] =
                static_cast<amrex::ParticleReal>(partial);
        }
    }
}

void
initializeInverseCdf (std::vector<RBEQShell> const& shells, double const energy_min,
                      double const log_energy_step, std::vector<double> const& uniform_thresholds,
                      amrex::Gpu::HostVector<amrex::ParticleReal>& inverse_cdf)
{
    constexpr int inverse_iterations = 36;
    constexpr int max_shells = BackgroundMCCIonizationModel::max_shell_count;
    constexpr int energy_count = BackgroundMCCIonizationModel::energy_grid_size;
    constexpr int quantile_count = BackgroundMCCIonizationModel::quantile_grid_size;

    for (int energy_index = 0; energy_index < energy_count; ++energy_index)
    {
        auto const incident_energy =
            energy_min * std::exp(log_energy_step * static_cast<double>(energy_index));

        std::vector<RBEQTerms> terms(shells.size());
        for (int shell_index = 0; shell_index < static_cast<int>(shells.size()); ++shell_index) {
            terms[shell_index] = rbeqTerms(incident_energy, shells[shell_index]);
        }

        for (int shell_index = 0; shell_index < static_cast<int>(shells.size()); ++shell_index) {
            auto const& shell = shells[shell_index];
            auto const available_energy = std::max(incident_energy - shell.binding_energy, 0.0);
            auto const maximum_secondary_energy = 0.5 * available_energy;
            auto const use_uniform =
                maximum_secondary_energy <= 0.0 || terms[shell_index].total <= 0.0 ||
                incident_energy <= uniform_thresholds[shell_index];

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
    BackgroundMCCIonizationTarget const target, amrex::ParticleReal const maximum_energy,
    BackgroundMCCRBEQ::Model const model)
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
    auto const shell_log_energy_step =
        std::log(energy_max / energy_min) /
        static_cast<double>(shell_energy_grid_size - 1);

    m_shell_cross_sections_h.assign(shell_energy_grid_size * max_shell_count, 0.0_prt);
    m_inverse_cdf_h.assign(energy_grid_size * max_shell_count * quantile_grid_size, 0.0_prt);

    m_executor_h.m_model = IonizationEnergySharingModel::RBEQ;
    m_executor_h.m_energy_grid_size = energy_grid_size;
    m_executor_h.m_shell_energy_grid_size = shell_energy_grid_size;
    m_executor_h.m_quantile_grid_size = quantile_grid_size;
    m_executor_h.m_energy_min = static_cast<amrex::ParticleReal>(energy_min);
    m_executor_h.m_log_energy_min = static_cast<amrex::ParticleReal>(std::log(energy_min));
    m_executor_h.m_inverse_log_energy_step =
        static_cast<amrex::ParticleReal>(1.0 / log_energy_step);

    auto const parameters =
        BackgroundMCCRBEQ::shells(target == BackgroundMCCIonizationTarget::N2, model);
    std::vector<double> positive_thresholds(parameters.size());
    std::vector<double> uniform_thresholds(parameters.size());
    findShellThresholds(parameters, positive_thresholds, uniform_thresholds);
    m_executor_h.m_shell_count = static_cast<int>(parameters.size());
    for (int i = 0; i < m_executor_h.m_shell_count; ++i) {
        m_executor_h.m_binding_energies[i] =
            static_cast<amrex::ParticleReal>(parameters[i].binding_energy);
        m_executor_h.m_positive_threshold_coordinates[i] = static_cast<amrex::ParticleReal>(
            std::log(positive_thresholds[i] / energy_min) / shell_log_energy_step);
        m_executor_h.m_uniform_threshold_coordinates[i] = static_cast<amrex::ParticleReal>(
            std::log(uniform_thresholds[i] / energy_min) / log_energy_step);
    }
    initializeShellCrossSections(parameters, energy_min, shell_log_energy_step,
                                 m_shell_cross_sections_h);
    initializeInverseCdf(parameters, energy_min, log_energy_step, uniform_thresholds,
                         m_inverse_cdf_h);

    m_executor_h.m_shell_cross_sections = m_shell_cross_sections_h.data();
    m_executor_h.m_inverse_cdf = m_inverse_cdf_h.data();

#ifdef AMREX_USE_GPU
    m_shell_cross_sections_d.resize(m_shell_cross_sections_h.size());
    m_inverse_cdf_d.resize(m_inverse_cdf_h.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_shell_cross_sections_h.begin(),
                          m_shell_cross_sections_h.end(), m_shell_cross_sections_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_inverse_cdf_h.begin(), m_inverse_cdf_h.end(),
                          m_inverse_cdf_d.begin());

    m_executor_d = m_executor_h;
    m_executor_d.m_shell_cross_sections = m_shell_cross_sections_d.data();
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
        return static_cast<amrex::ParticleReal>(15.58);
    }
    if (target == BackgroundMCCIonizationTarget::O2)
    {
        return static_cast<amrex::ParticleReal>(12.07);
    }
    return 0.0;
}
