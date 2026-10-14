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
 * State populations and sparse energy-change alias tables are built on the host.
 * The elastic DCS supplies an independent angle. Only the internal-energy
 * threshold restricts rotational outcomes; there is no device rotational-state loop.
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
        // Keep canonical thresholds in double precision in both particle builds.
        double m_loss = 0;
        double m_initial_energy = 0;
    };

    struct Executor
    {
        struct Interpolation
        {
            int m_index = 0;
            amrex::ParticleReal m_fraction = 0;
            amrex::ParticleReal m_rate = 0;
            double m_energy = 0;
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
                        fraction * m_rates[index + 1],
                    double(energy)};
        }

        /** Sample one alias without a transition search or rejection. */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static int
        aliasSample (Alias const *table, Cell cell, double draw) noexcept
        {
            auto const coordinate =
                amrex::min(draw, 1.0 - std::numeric_limits<double>::epsilon()) * cell.m_count;
            int const index = static_cast<int>(coordinate);
            auto const fraction = coordinate - index;
            auto const& entry = table[cell.m_offset + index];
            if (fraction < entry.m_probability) {
                return entry.m_outcome;
            }
            return table[cell.m_offset + entry.m_alternate].m_outcome;
        }

        /** Mix energy rows, then sample an angle-independent internal change.
         *
         * A subthreshold excitation caused by energy-grid roundoff becomes an
         * unchanged event. This retains the inclusive rate without resampling
         * or increasing any other rotational channel's rate.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Outcome
        sample (Interpolation state, double row_draw, double loss_draw) const noexcept
        {
            double const lower_rate = (1 - double(state.m_fraction)) * m_rates[state.m_index];
            int const row = state.m_index + (row_draw * state.m_rate >= lower_rate ? 1 : 0);
            auto const cell = m_cells[row];
            int outcome;
            if (!m_cumulative) {
                outcome = aliasSample(m_loss_alias, cell, loss_draw);
            } else {
                int lo = 0, hi = cell.m_count - 1;
                while (lo < hi) {
                    int const mid = lo + (hi - lo) / 2;
                    if (loss_draw < m_cdf[cell.m_offset + mid]) {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                outcome = m_loss_alias[cell.m_offset + lo].m_outcome;
            }
            auto const result = m_outcomes[outcome];
            return result.m_loss > state.m_energy ? Outcome{} : result;
        }

        amrex::ParticleReal const* m_energies = nullptr;
        amrex::ParticleReal const* m_rates = nullptr;
        double const* m_cdf = nullptr;
        Alias const* m_loss_alias = nullptr;
        Cell const* m_cells = nullptr;
        Outcome const* m_outcomes = nullptr;
        int m_energy_count = 0;
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
    amrex::Gpu::HostVector<amrex::ParticleReal> m_energies_h, m_rates_h;
    amrex::Gpu::HostVector<double> m_cdf_h;
    amrex::Gpu::HostVector<Alias> m_loss_alias_h;
    amrex::Gpu::HostVector<Cell> m_cells_h;
    amrex::Gpu::HostVector<Outcome> m_outcomes_h;
    amrex::Gpu::HostVector<double> m_reference_rates;
    Executor m_executor_h;
#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_energies_d, m_rates_d;
    amrex::Gpu::DeviceVector<double> m_cdf_d;
    amrex::Gpu::DeviceVector<Alias> m_loss_alias_d;
    amrex::Gpu::DeviceVector<Cell> m_cells_d;
    amrex::Gpu::DeviceVector<Outcome> m_outcomes_d;
    Executor m_executor_d;
#endif
    double m_neutral_mass = 0;
    double m_maximum_rate = 0;
    std::size_t m_table_bytes = 0;
};

#endif
