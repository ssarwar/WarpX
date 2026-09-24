/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_THERMAL_ROTATION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_THERMAL_ROTATION_H_

#include "BackgroundMCCElasticKinematics.H"

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
 * The elastic DCS supplies the angle. Only exact kinematic accessibility couples
 * that angle to the rotational draw; no rotational DCS or device state loop is used.
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

        /** Sample one alias and retain its conditional uniform in draw. */
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

        struct Accessible
        {
            int m_count = 0;
            double m_probability = 0;
        };

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
        accessible (int outcome, double energy, double cosine) const noexcept
        {
            auto const value = m_outcomes[outcome];
            return value.m_loss <= 0 || energy >= value.m_loss * m_all_angle_factor ||
                   value.m_loss <= BackgroundMCCElasticKinematics::maximumInternalEnergyChange(
                                       energy, m_electron_rest_energy,
                                       m_neutral_rest_energy + value.m_initial_energy, cosine);
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Accessible
        allowed (int row, double energy, double cosine) const noexcept
        {
            auto const cell = m_cells[row];
            if (accessible(m_loss_alias[cell.m_offset + cell.m_count - 1].m_outcome, energy,
                           cosine)) {
                return {cell.m_count, 1};
            }
            // Outcomes are ordered by increasing loss, then decreasing initial
            // mass for equal losses. The host verifies the heavy-rotor ordering.
            int lo = 0, hi = cell.m_count;
            while (lo < hi) {
                int const mid = lo + (hi - lo) / 2;
                if (accessible(m_loss_alias[cell.m_offset + mid].m_outcome, energy, cosine)) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            return {lo, lo > 0 ? double(m_cdf[cell.m_offset + lo - 1]) : 0};
        }

        /** Condition the mixture on accessibility, retaining the elastic angle.
         *
         * The usual path uses one alias. Restricted threshold draws use a prefix
         * CDF with bounded bisections; there is no rejection or rotational loop.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Outcome
        sample (Interpolation state, double cosine, double row_draw,
                double loss_draw) const noexcept
        {
            auto const lower = allowed(state.m_index, state.m_energy, cosine);
            auto const upper = allowed(state.m_index + 1, state.m_energy, cosine);
            double const lo_rate =
                (1 - double(state.m_fraction)) * m_rates[state.m_index] * lower.m_probability;
            double const hi_rate =
                double(state.m_fraction) * m_rates[state.m_index + 1] * upper.m_probability;
            int const offset = row_draw * (lo_rate + hi_rate) >= lo_rate ? 1 : 0;
            auto const valid = offset == 0 ? lower : upper;
            auto const cell = m_cells[state.m_index + offset];
            if (valid.m_count == cell.m_count && !m_cumulative) {
                return m_outcomes[aliasSample(m_loss_alias, cell, loss_draw)];
            }
            double const target = loss_draw * valid.m_probability;
            int lo = 0, hi = valid.m_count - 1;
            while (lo < hi) {
                int const mid = lo + (hi - lo) / 2;
                if (target < m_cdf[cell.m_offset + mid]) {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            return m_outcomes[m_loss_alias[cell.m_offset + lo].m_outcome];
        }

        amrex::ParticleReal const* m_energies = nullptr;
        amrex::ParticleReal const* m_rates = nullptr;
        double const* m_cdf = nullptr;
        Alias const* m_loss_alias = nullptr;
        Cell const* m_cells = nullptr;
        Outcome const* m_outcomes = nullptr;
        int m_energy_count = 0;
        amrex::ParticleReal m_maximum_energy = 0;
        double m_neutral_rest_energy = 0;
        double m_electron_rest_energy =
            PhysConst::m_e_v<double> * PhysConst::c2_v<double> / PhysConst::q_e_v<double>;
        double m_all_angle_factor = std::numeric_limits<double>::infinity();
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
    maximumLoss () const
    {
        return m_maximum_loss;
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
    double m_maximum_loss = 0;
    std::size_t m_table_bytes = 0;
};

#endif
