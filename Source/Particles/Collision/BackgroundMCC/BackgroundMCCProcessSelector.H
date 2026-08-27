/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_PROCESS_SELECTOR_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_PROCESS_SELECTOR_H_

#include "Particles/Collision/ScatteringProcess.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/** Exact cumulative cross-section table for Background MCC process selection.
 *
 * The host constructor merges all process energy knots. A sum of piecewise
 * linear cross sections is itself linear between consecutive knots in this
 * union grid, so runtime interpolation remains exact. Device selection uses
 * one energy bisection and one process-prefix bisection instead of one energy
 * bisection per process.
 */
class BackgroundMCCProcessSelector
{
public:
    struct Executor
    {
        /** Select the process whose cumulative cross section contains a draw.
         *
         * @param collision_energy electron collision energy in eV
         * @param cross_section_draw draw expressed in cross-section units
         * @param selected_process selected process index, or -1 for a null event
         * @param total_cross_section interpolated sum of all process cross sections
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void select (
            amrex::ParticleReal const collision_energy,
            amrex::ParticleReal const cross_section_draw,
            int& selected_process,
            amrex::ParticleReal& total_cross_section) const noexcept
        {
            using namespace amrex::literals;

            selected_process = -1;
            total_cross_section = 0.0_prt;
            if (!enabled()) { return; }

            int energy_index = 0;
            amrex::ParticleReal energy_fraction = 0.0_prt;
            if (collision_energy >= m_energy_hi)
            {
                energy_index = m_energy_count - 2;
                energy_fraction = 1.0_prt;
            }
            else if (collision_energy > m_energy_lo)
            {
                energy_index = amrex::bisect(
                    m_energies, 0, m_energy_count - 1, collision_energy);
                energy_fraction =
                    (collision_energy - m_energies[energy_index]) /
                    (m_energies[energy_index + 1] - m_energies[energy_index]);
            }

            total_cross_section = cumulativeCrossSection(
                m_process_count - 1, energy_index, energy_fraction);
            if (!(cross_section_draw < total_cross_section)) { return; }

            int process_lo = 0;
            int process_hi = m_process_count - 1;
            while (process_lo < process_hi)
            {
                auto const process_mid = process_lo + (process_hi - process_lo) / 2;
                auto const cumulative = cumulativeCrossSection(
                    process_mid, energy_index, energy_fraction);
                if (cross_section_draw < cumulative) {
                    process_hi = process_mid;
                } else {
                    process_lo = process_mid + 1;
                }
            }
            selected_process = process_lo;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
        enabled () const noexcept
        {
            return m_energy_count >= 2 && m_process_count > 0;
        }

        amrex::ParticleReal const* m_energies = nullptr;
        amrex::ParticleReal const* m_cumulative_cross_sections = nullptr;
        int m_energy_count = 0;
        int m_process_count = 0;
        amrex::ParticleReal m_energy_lo = 0.0;
        amrex::ParticleReal m_energy_hi = 0.0;

    private:
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
        cumulativeCrossSection (
            int const process_index, int const energy_index,
            amrex::ParticleReal const energy_fraction) const noexcept
        {
            auto const index = process_index * m_energy_count + energy_index;
            auto const value_lo = m_cumulative_cross_sections[index];
            return value_lo + energy_fraction *
                (m_cumulative_cross_sections[index + 1] - value_lo);
        }
    };

    explicit BackgroundMCCProcessSelector (
        amrex::Vector<ScatteringProcess> const& processes);

    [[nodiscard]] Executor const& executor () const
    {
#ifdef AMREX_USE_GPU
        return m_executor_d;
#else
        return m_executor_h;
#endif
    }

private:
    amrex::Gpu::HostVector<amrex::ParticleReal> m_energies_h;
    amrex::Gpu::HostVector<amrex::ParticleReal> m_cumulative_cross_sections_h;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_energies_d;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_cumulative_cross_sections_d;
    Executor m_executor_d;
#endif
    Executor m_executor_h;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_PROCESS_SELECTOR_H_
