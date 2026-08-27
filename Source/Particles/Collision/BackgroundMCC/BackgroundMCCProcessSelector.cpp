/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCProcessSelector.H"

#include "Utils/TextMsg.H"

#include <AMReX_Gpu.H>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

BackgroundMCCProcessSelector::BackgroundMCCProcessSelector (
    amrex::Vector<ScatteringProcess> const& processes)
{
    if (processes.empty()) { return; }

    std::size_t knot_count = 0u;
    for (auto const& process : processes)
    {
        auto const grid_size = static_cast<std::size_t>(process.getEnergyGrid().size());
        if (grid_size > std::numeric_limits<std::size_t>::max() - knot_count) { return; }
        knot_count += grid_size;
    }

    std::vector<amrex::ParticleReal> union_energies;
    union_energies.reserve(knot_count);
    for (auto const& process : processes)
    {
        auto const& energies = process.getEnergyGrid();
        union_energies.insert(union_energies.end(), energies.begin(), energies.end());
    }
    std::sort(union_energies.begin(), union_energies.end());
    union_energies.erase(
        std::unique(union_energies.begin(), union_energies.end()), union_energies.end());

    auto const energy_count = union_energies.size();
    auto const process_count = static_cast<std::size_t>(processes.size());
    constexpr std::size_t max_table_bytes = 64u * 1024u * 1024u;
    constexpr std::size_t max_table_entries =
        max_table_bytes / sizeof(amrex::ParticleReal);
    if (energy_count < 2u ||
        energy_count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        process_count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        process_count > max_table_entries / energy_count)
    {
        // Very heterogeneous tables can make the exact prefix representation
        // quadratic in input size. The collision kernel retains its linear-scan
        // fallback for this uncommon case.
        return;
    }

    m_energies_h.resize(energy_count);
    std::copy(union_energies.begin(), union_energies.end(), m_energies_h.begin());
    m_cumulative_cross_sections_h.resize(energy_count * process_count);
    std::vector<long double> cumulative(energy_count, 0.0L);

    for (std::size_t process_index = 0; process_index < process_count; ++process_index)
    {
        auto const& energies = processes[process_index].getEnergyGrid();
        auto const& cross_sections = processes[process_index].getCrossSectionGrid();
        auto const process_energy_count = static_cast<std::size_t>(energies.size());
        std::size_t interval = 0u;
        for (std::size_t energy_index = 0; energy_index < energy_count; ++energy_index)
        {
            auto const energy = union_energies[energy_index];
            long double cross_section;
            if (energy <= energies.front())
            {
                cross_section = static_cast<long double>(cross_sections.front());
            }
            else if (energy >= energies.back())
            {
                cross_section = static_cast<long double>(cross_sections.back());
            }
            else
            {
                while (interval + 1u < process_energy_count &&
                       energy > energies[interval + 1u]) {
                    ++interval;
                }
                auto const fraction =
                    static_cast<long double>(energy - energies[interval]) /
                    static_cast<long double>(energies[interval + 1u] - energies[interval]);
                cross_section = static_cast<long double>(cross_sections[interval]) +
                    fraction * static_cast<long double>(
                        cross_sections[interval + 1u] - cross_sections[interval]);
            }

            cumulative[energy_index] += cross_section;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(cumulative[energy_index]) &&
                    cumulative[energy_index] <= static_cast<long double>(
                        std::numeric_limits<amrex::ParticleReal>::max()),
                "The total Background MCC cross section is not representable as ParticleReal."
            );
            m_cumulative_cross_sections_h[process_index * energy_count + energy_index] =
                static_cast<amrex::ParticleReal>(cumulative[energy_index]);
        }
    }

    m_executor_h.m_energies = m_energies_h.data();
    m_executor_h.m_cumulative_cross_sections = m_cumulative_cross_sections_h.data();
    m_executor_h.m_energy_count = static_cast<int>(energy_count);
    m_executor_h.m_process_count = static_cast<int>(process_count);
    m_executor_h.m_energy_lo = m_energies_h.front();
    m_executor_h.m_energy_hi = m_energies_h.back();

#ifdef AMREX_USE_GPU
    m_energies_d.resize(m_energies_h.size());
    m_cumulative_cross_sections_d.resize(m_cumulative_cross_sections_h.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_energies_h.begin(), m_energies_h.end(),
                          m_energies_d.begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_cumulative_cross_sections_h.begin(),
        m_cumulative_cross_sections_h.end(), m_cumulative_cross_sections_d.begin());

    m_executor_d = m_executor_h;
    m_executor_d.m_energies = m_energies_d.data();
    m_executor_d.m_cumulative_cross_sections = m_cumulative_cross_sections_d.data();
    amrex::Gpu::streamSynchronize();
#endif
}
