/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_THERMAL_ROTATION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_THERMAL_ROTATION_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>

class ScatteringProcess;

/** One thermal vibrationally elastic family, including signed rotational
 * changes.
 *
 * State populations, angular marginals and sparse conditional alias tables are
 * built on the host. No rotational-state loop or special function is evaluated
 * on the device. The reference data contain angular-bin integrals of v*sigma.
 */
class BackgroundMCCThermalRotation
{
  public:
    struct Alias
    {
        amrex::ParticleReal m_probability = 1;
        int m_alternate = 0;
        int m_outcome = 0;
    };

    struct Cell
    {
        int m_offset = 0;
        int m_count = 0;
    };

    struct Outcome
    {
        amrex::ParticleReal m_loss = 0;
        amrex::ParticleReal m_initial_energy = 0;
    };

    struct Executor
    {
        struct Interpolation
        {
            int m_index = 0;
            amrex::ParticleReal m_fraction = 0;
            amrex::ParticleReal m_rate = 0;
        };

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
        enabled () const noexcept
        {
            return m_energy_count != 0;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
        inRange (amrex::ParticleReal energy) const noexcept
        {
            return energy >= 0 && energy <= m_maximum_energy;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Interpolation
        interpolate (amrex::ParticleReal energy) const noexcept
        {
            if (!enabled()) {
                return {};
            }
            int const index = energy >= m_maximum_energy
                                  ? m_energy_count - 2
                                  : amrex::bisect(m_energies, 0, m_energy_count - 1, energy);
            auto fraction = amrex::max(amrex::ParticleReal(0),
                                       amrex::min((energy - m_energies[index]) /
                                                      (m_energies[index + 1] - m_energies[index]),
                                                  amrex::ParticleReal(1)));
            if (index == 0) {
                fraction = std::sqrt(fraction);
            }
            return {index, fraction,
                    (amrex::ParticleReal(1) - fraction) * m_rates[index] +
                        fraction * m_rates[index + 1]};
        }

        /** Sample one alias and recycle its conditional uniform within the
         * outcome. This supplies the within-bin angular coordinate without
         * another draw.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static int
        aliasSample (Alias const* table, Cell cell, double& draw) noexcept
        {
            auto const coordinate =
                amrex::min(draw, 1.0 - std::numeric_limits<double>::epsilon()) * cell.m_count;
            int const index = static_cast<int>(coordinate);
            auto const fraction = coordinate - index;
            auto const& entry = table[cell.m_offset + index];
            if (fraction < entry.m_probability) {
                draw = fraction / entry.m_probability;
                return entry.m_outcome;
            }
            draw = (fraction - entry.m_probability) / (1.0 - entry.m_probability);
            return table[cell.m_offset + entry.m_alternate].m_outcome;
        }

        /** Explicit cumulative reference for physics-equivalent performance
         * comparisons. */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
        select (Alias const* table, Cell cell, double& draw) const noexcept
        {
            draw = amrex::max(0.0, amrex::min(draw, 1.0 - std::numeric_limits<double>::epsilon()));
            if (!m_cumulative) {
                return aliasSample(table, cell, draw);
            }
            int lo = 0, hi = cell.m_count - 1;
            while (lo < hi) {
                int const mid = lo + (hi - lo) / 2;
                if (draw < table[cell.m_offset + mid].m_probability) {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            double const left = lo > 0 ? table[cell.m_offset + lo - 1].m_probability : 0;
            double const right = table[cell.m_offset + lo].m_probability;
            draw = right > left ? (draw - left) / (right - left) : 0;
            return table[cell.m_offset + lo].m_outcome;
        }

        /** Mixtures interpolate probabilities; discrete loss labels remain
         * exact. */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Outcome
        sample (Interpolation state, double row_draw, double angle_draw, double loss_draw,
                double& cosine) const noexcept
        {
            auto const lower_rate =
                (amrex::ParticleReal(1) - state.m_fraction) * m_rates[state.m_index];
            int const row = state.m_index + (row_draw * state.m_rate >= lower_rate ? 1 : 0);
            int const angle =
                select(m_angle_alias, {row * m_angle_count, m_angle_count}, angle_draw);
            cosine = (1 - angle_draw) * m_edges[angle] + angle_draw * m_edges[angle + 1];
            int const outcome =
                select(m_loss_alias, m_cells[row * m_angle_count + angle], loss_draw);
            return m_outcomes[outcome];
        }

        amrex::ParticleReal const* m_energies = nullptr;
        amrex::ParticleReal const* m_rates = nullptr;
        amrex::ParticleReal const* m_edges = nullptr;
        Alias const* m_angle_alias = nullptr;
        Alias const* m_loss_alias = nullptr;
        Cell const* m_cells = nullptr;
        Outcome const* m_outcomes = nullptr;
        int m_energy_count = 0;
        int m_angle_count = 0;
        amrex::ParticleReal m_maximum_energy = 0;
        bool m_cumulative = false;
    };

    static std::shared_ptr<BackgroundMCCThermalRotation> get (std::string const& file,
                                                              std::string const& model,
                                                              double temperature,
                                                              bool cumulative = false);

    BackgroundMCCThermalRotation (std::string const& file, std::string const& model,
                                  double temperature, bool cumulative = false);

    /** Check the source-temperature inclusive rate before replacing its elastic
     * channel. */
    void checkInclusiveRate (ScatteringProcess const& process) const;

    [[nodiscard]] Executor const&
    executor () const
    {
#ifdef AMREX_USE_GPU
        return m_executor_d;
#else
        return m_executor_h;
#endif
    }
    [[nodiscard]] Executor const&
    hostExecutor () const
    {
        return m_executor_h;
    }
    [[nodiscard]] double
    neutralMass () const
    {
        return m_neutral_mass;
    }
    [[nodiscard]] double
    maximumRate () const
    {
        return m_maximum_rate;
    }
    [[nodiscard]] std::size_t
    tableBytes () const
    {
        return m_table_bytes;
    }
    [[nodiscard]] amrex::Gpu::HostVector<amrex::ParticleReal> const&
    energies () const
    {
        return m_energies_h;
    }

  private:
    amrex::Gpu::HostVector<amrex::ParticleReal> m_energies_h, m_rates_h, m_edges_h;
    amrex::Gpu::HostVector<Alias> m_angle_alias_h, m_loss_alias_h;
    amrex::Gpu::HostVector<Cell> m_cells_h;
    amrex::Gpu::HostVector<Outcome> m_outcomes_h;
    amrex::Gpu::HostVector<double> m_reference_rates;
    Executor m_executor_h;
#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_energies_d, m_rates_d, m_edges_d;
    amrex::Gpu::DeviceVector<Alias> m_angle_alias_d, m_loss_alias_d;
    amrex::Gpu::DeviceVector<Cell> m_cells_d;
    amrex::Gpu::DeviceVector<Outcome> m_outcomes_d;
    Executor m_executor_d;
#endif
    double m_neutral_mass = 0;
    double m_maximum_rate = 0;
    std::size_t m_table_bytes = 0;
};

#endif
