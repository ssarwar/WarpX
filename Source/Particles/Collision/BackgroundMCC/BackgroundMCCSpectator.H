/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_SPECTATOR_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_SPECTATOR_H_

#include "BackgroundMCCThermalRotation.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <cstddef>
#include <memory>
#include <string>

/** Discrete IAA spectator outcomes conditional on the inclusive elastic angle.
 *
 * All differential weights are tabulated offline. A Boltzmann alias selects a
 * virtual initial level, followed by seven possible final levels. The amount
 * of device work is independent of the number of thermally populated levels;
 * no neutral state is stored or evolved, and no special function is evaluated.
 */
class BackgroundMCCSpectator
{
  public:
    using Alias = BackgroundMCCThermalRotation::Alias;
    using Outcome = BackgroundMCCThermalRotation::Outcome;

    struct Executor
    {
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool enabled () const noexcept
        {
            return m_energy_count > 0;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
        inRange (double energy) const noexcept
        {
            return energy >= 0 && energy <= m_maximum_energy;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Outcome
        sample (double energy, double cosine, double row_draw, double initial_draw,
                double outcome_draw) const noexcept
        {
            int const lower = energy >= m_maximum_energy
                                  ? m_energy_count - 2
                                  : amrex::bisect(m_energies, 0, m_energy_count - 1,
                                                  static_cast<amrex::ParticleReal>(energy));
            double fraction =
                (energy - m_energies[lower]) / (m_energies[lower + 1] - m_energies[lower]);
            fraction = amrex::max(0.0, amrex::min(fraction, 1.0));
            if (lower == 0) {
                fraction = std::sqrt(fraction);
            }
            int const row = lower + (row_draw < fraction ? 1 : 0);
            int initial;
            if (!m_cumulative) {
                initial = BackgroundMCCThermalRotation::Executor::aliasSample(
                    m_initial_alias, {0, m_population_count}, initial_draw);
            } else {
                int lo = 0, hi = m_population_count - 1;
                while (lo < hi) {
                    int const mid = lo + (hi - lo) / 2;
                    if (initial_draw < m_initial_cdf[mid]) {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                initial = m_initial_alias[lo].m_outcome;
            }

            double const y = std::sqrt(amrex::max(0.0, (1 - cosine) * 0.5));
            int const first = m_angle_offsets[row];
            int const last = m_angle_offsets[row + 1] - 1;
            int const angle = y >= 1 ? last - 1 : amrex::bisect(m_angles, first, last, float(y));
            double const angular_fraction = amrex::max(
                0.0,
                amrex::min((y - m_angles[angle]) / (m_angles[angle + 1] - m_angles[angle]), 1.0));
            double basis[4];
            for (int rank = 0; rank < 4; ++rank) {
                basis[rank] = (1 - angular_fraction) * m_angular_basis[4 * angle + rank] +
                              angular_fraction * m_angular_basis[4 * (angle + 1) + rank];
            }
            auto const* weights = m_weights + (row * m_state_count + initial) * 28;
            double scores[7];
            double total = 0;
            for (int choice = 0; choice < 7; ++choice) {
                auto const* w = weights + 4 * choice;
                scores[choice] =
                    w[0] * basis[0] + w[1] * basis[1] + w[2] * basis[2] + w[3] * basis[3];
                total += scores[choice];
            }
            if (total <= 0) {
                return {};
            }
            double const draw = outcome_draw * total;
            double prefix = 0;
            int choice = 3;
            for (int i = 0; i < 7; ++i) {
                prefix += scores[i];
                if (draw < prefix) {
                    choice = i;
                    break;
                }
            }
            int const final = initial + 2 * (choice - 3);
            constexpr double b = 0.0002477204284695341; // N2, eV
            double const loss = b * (final * (final + 1) - initial * (initial + 1));
            // Energy-grid roundoff never permits a subthreshold excitation.
            if (loss > energy) {
                return {};
            }
            return {loss, b * initial * (initial + 1)};
        }

        amrex::ParticleReal const* m_energies = nullptr;
        int const* m_angle_offsets = nullptr;
        float const* m_angles = nullptr;
        float const* m_angular_basis = nullptr;
        float const* m_weights = nullptr;
        Alias const* m_initial_alias = nullptr;
        double const* m_initial_cdf = nullptr;
        int m_energy_count = 0;
        int m_state_count = 0;
        int m_population_count = 0;
        double m_maximum_energy = 0;
        bool m_cumulative = false;
    };

    BackgroundMCCSpectator (std::string const& file, double temperature, bool cumulative);

    static std::shared_ptr<BackgroundMCCSpectator> get (std::string const& file, double temperature,
                                                        bool cumulative);

    [[nodiscard]] Executor const& executor () const
    {
#ifdef AMREX_USE_GPU
        return m_executor_d;
#else
        return m_executor_h;
#endif
    }

    [[nodiscard]] Executor const& hostExecutor () const
    {
        return m_executor_h;
    }
    [[nodiscard]] std::size_t tableBytes () const
    {
        return m_table_bytes;
    }
    [[nodiscard]] static double neutralMass ()
    {
        return 28.0134 * 1.66053906660e-27;
    }

  private:
    struct Data;
    std::shared_ptr<Data const> m_data;
    amrex::Gpu::HostVector<Alias> m_alias_h;
    amrex::Gpu::HostVector<double> m_cdf_h;
    Executor m_executor_h;
#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<Alias> m_alias_d;
    amrex::Gpu::DeviceVector<double> m_cdf_d;
    Executor m_executor_d;
#endif
    std::size_t m_table_bytes = 0;
};

#endif
