/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCElasticScattering.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{
    void readDifferentialCrossSection (
        std::string const& file_name, amrex::Gpu::HostVector<amrex::ParticleReal>& energies,
        std::vector<std::vector<double>>& differential_cross_sections,
        double& screening_radius)
    {
        std::ifstream input(file_name);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            input.is_open(), "Could not open the elastic differential-cross-section file '" +
                                 file_name + "'.");

        std::string line;
        std::size_t angle_count = 0;
        while (std::getline(input, line))
        {
            std::istringstream row_stream(line);
            double energy;
            if (!(row_stream >> energy))
            {
                std::istringstream metadata_stream(line);
                std::string label;
                std::string projectile;
                std::string separator;
                std::string target;
                metadata_stream >> label;
                if (label == "#") { metadata_stream >> label; }
                if (metadata_stream >> projectile >> separator >> target &&
                    label == "SPECIES:" && projectile == "e" && separator == "/")
                {
                    if (target == "N2") {
                        screening_radius = 0.6052;
                    } else if (target == "O2") {
                        screening_radius = 0.5677;
                    }
                }
                continue;
            }

            std::vector<double> row;
            double value;
            while (row_stream >> value) { row.push_back(value); }

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                row.size() >= 3u,
                "Each elastic differential-cross-section row must contain one energy "
                "and at least three angular values.");
            if (angle_count == 0u) { angle_count = row.size(); }
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                row.size() == angle_count,
                "Elastic differential-cross-section rows must have the same number "
                "of angular values.");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(energy) && energy > 0.0,
                "Elastic differential-cross-section energies must be finite and positive.");
            auto const particle_energy = static_cast<amrex::ParticleReal>(energy);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(static_cast<double>(particle_energy)),
                "Elastic differential-cross-section energies exceed particle precision.");
            if (!energies.empty())
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    particle_energy > energies.back(),
                    "Elastic differential-cross-section energies must remain strictly "
                    "increasing in particle precision.");
            }
            for (auto const dcs : row)
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(dcs) && dcs >= 0.0,
                    "Elastic differential cross sections must be finite and non-negative.");
            }

            energies.push_back(particle_energy);
            differential_cross_sections.push_back(std::move(row));
            if (screening_radius > 0.0 &&
                particle_energy >=
                    BackgroundMCCElasticScatteringModel::Executor::
                        m_screened_rutherford_energy)
            {
                // The first row at or above the analytic cutoff completes the
                // final interpolation interval. Later rows are never sampled.
                break;
            }
        }

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            energies.size() >= 2u,
            "Elastic differential-cross-section tables must contain at least two energy rows.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            energies.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
            "Elastic differential-cross-section table contains too many energy rows.");
    }

    void initializeInverseCdf (
        std::vector<std::vector<double>> const& differential_cross_sections,
        amrex::Gpu::HostVector<amrex::ParticleReal>& inverse_cdf_deflection,
        amrex::Gpu::HostVector<amrex::ParticleReal>& row_integrals)
    {
        constexpr int quantile_count =
            BackgroundMCCElasticScatteringModel::quantile_grid_size;
        auto const angle_count = differential_cross_sections.front().size();
        auto const angle_step = static_cast<double>(MathConst::pi) /
                                static_cast<double>(angle_count - 1u);
        inverse_cdf_deflection.resize(
            differential_cross_sections.size() * quantile_count);
        row_integrals.resize(differential_cross_sections.size());

        std::vector<double> half_sines(angle_count);
        for (std::size_t angle_index = 0; angle_index < angle_count; ++angle_index)
        {
            auto const theta = static_cast<double>(angle_index) * angle_step;
            half_sines[angle_index] = std::sin(0.5 * theta);
        }
        std::vector<double> cumulative(angle_count);
        for (std::size_t energy_index = 0;
             energy_index < differential_cross_sections.size(); ++energy_index)
        {
            auto const& dcs = differential_cross_sections[energy_index];
            cumulative[0] = 0.0;
            for (std::size_t angle_index = 1; angle_index < angle_count; ++angle_index)
            {
                auto const y_lo = half_sines[angle_index - 1u];
                auto const y_hi = half_sines[angle_index];
                auto const width = y_hi - y_lo;
                auto const dcs_lo = dcs[angle_index - 1u];
                auto const slope = (dcs[angle_index] - dcs_lo) / width;
                auto const segment_integral = y_lo * dcs_lo * width +
                    0.5 * (y_lo * slope + dcs_lo) * width * width +
                    slope * width * width * width / 3.0;
                cumulative[angle_index] =
                    cumulative[angle_index - 1u] + segment_integral;
            }

            auto const total = cumulative.back();
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(total) && total > 0.0,
                "Every elastic differential-cross-section energy row must have a "
                "finite, positive angular integral.");
            auto const particle_total = static_cast<amrex::ParticleReal>(total);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(static_cast<double>(particle_total)) && particle_total > 0.0,
                "Elastic differential-cross-section row integrals exceed particle "
                "precision.");
            row_integrals[energy_index] = particle_total;

            auto const row_offset = energy_index * quantile_count;
            inverse_cdf_deflection[row_offset] = 0.0;
            std::size_t angle_index = 1u;
            for (int quantile = 1; quantile < quantile_count - 1; ++quantile)
            {
                auto const coordinate = static_cast<double>(quantile) /
                                        static_cast<double>(quantile_count - 1);
                auto const probability =
                    1.0 - (1.0 - coordinate) * (1.0 - coordinate);
                auto const target = total * probability;
                while (angle_index + 1u < angle_count && cumulative[angle_index] < target)
                {
                    ++angle_index;
                }

                auto const y_lo = half_sines[angle_index - 1u];
                auto const y_hi = half_sines[angle_index];
                auto const dcs_lo = dcs[angle_index - 1u];
                auto const slope =
                    (dcs[angle_index] - dcs_lo) / (y_hi - y_lo);
                auto const partialIntegral = [=] (double const width)
                {
                    return y_lo * dcs_lo * width +
                        0.5 * (y_lo * slope + dcs_lo) * width * width +
                        slope * width * width * width / 3.0;
                };
                auto lower = 0.0;
                auto upper = y_hi - y_lo;
                auto const local_target = target - cumulative[angle_index - 1u];
                for (int iteration = 0; iteration < 60; ++iteration)
                {
                    auto const midpoint = 0.5 * (lower + upper);
                    if (partialIntegral(midpoint) < local_target) {
                        lower = midpoint;
                    } else {
                        upper = midpoint;
                    }
                }
                auto const half_sine = y_lo + 0.5 * (lower + upper);
                inverse_cdf_deflection[row_offset + quantile] =
                    static_cast<amrex::ParticleReal>(2.0 * half_sine * half_sine);
            }
            inverse_cdf_deflection[row_offset + quantile_count - 1] = 2.0;
        }
    }
} // namespace

BackgroundMCCElasticScatteringModel::BackgroundMCCElasticScatteringModel (
    std::string const& differential_cross_section)
{
    std::vector<std::vector<double>> differential_cross_sections;
    amrex::Gpu::HostVector<amrex::ParticleReal> energies;
    double screening_radius = 0.0;
    readDifferentialCrossSection(
        differential_cross_section, energies, differential_cross_sections,
        screening_radius);
    initializeInverseCdf(
        differential_cross_sections, m_inverse_cdf_deflection_h, m_row_integrals_h);

    m_log_energies_h.resize(energies.size());
    for (std::size_t i = 0; i < energies.size(); ++i) {
        m_log_energies_h[i] = std::log(energies[i]);
    }

    m_executor_h.m_log_energies = m_log_energies_h.data();
    m_executor_h.m_row_integrals = m_row_integrals_h.data();
    m_executor_h.m_inverse_cdf_deflection = m_inverse_cdf_deflection_h.data();
    m_executor_h.m_energy_grid_size = static_cast<int>(energies.size());
    m_executor_h.m_energy_lo = energies.front();
    m_executor_h.m_energy_hi = energies.back();
    m_executor_h.m_screening_radius = screening_radius;

#ifdef AMREX_USE_GPU
    m_log_energies_d.resize(m_log_energies_h.size());
    m_row_integrals_d.resize(m_row_integrals_h.size());
    m_inverse_cdf_deflection_d.resize(m_inverse_cdf_deflection_h.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_log_energies_h.begin(), m_log_energies_h.end(),
        m_log_energies_d.begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_row_integrals_h.begin(), m_row_integrals_h.end(),
        m_row_integrals_d.begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_inverse_cdf_deflection_h.begin(),
        m_inverse_cdf_deflection_h.end(), m_inverse_cdf_deflection_d.begin());

    m_executor_d = m_executor_h;
    m_executor_d.m_log_energies = m_log_energies_d.data();
    m_executor_d.m_row_integrals = m_row_integrals_d.data();
    m_executor_d.m_inverse_cdf_deflection = m_inverse_cdf_deflection_d.data();
    amrex::Gpu::streamSynchronize();
#endif
}
