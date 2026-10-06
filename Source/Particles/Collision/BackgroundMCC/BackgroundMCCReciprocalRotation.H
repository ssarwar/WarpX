/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_RECIPROCAL_ROTATION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_RECIPROCAL_ROTATION_H_

#include "BackgroundMCCThermalRotation.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

class ScatteringProcess;

/** Read-only, offline thermal rotational distributions and their elastic
 * marginal. */
class BackgroundMCCReciprocalRotation
{
public:
    struct Alias
    {
        float m_probability;
        std::uint16_t m_alternate;
        std::uint16_t m_outcome;
    };
    static_assert(sizeof(Alias) == 8, "The V6 alias layout requires eight bytes.");
    using Outcome = BackgroundMCCThermalRotation::Outcome;

    struct Executor
    {
        static constexpr int lookup_bins = 256;
        struct Interpolation
        {
            int m_index = 0;
            double m_fraction = 0;
            double m_rate = 0;
            double m_energy = 0;
        };
        struct Sample
        {
            Outcome m_outcome;
            double m_deflection = 0;
            bool m_valid = true;
        };

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool enabled () const noexcept
        {
            return m_energy_count > 0;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
        inRange (double energy) const noexcept
        {
            constexpr double tolerance = 8 * std::numeric_limits<amrex::ParticleReal>::epsilon();
            return energy >= 0 && energy <= m_maximum_energy * (1 + tolerance);
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Interpolation
        interpolate (double energy) const noexcept
        {
            if (!enabled()) {
                return {};
            }
            energy = amrex::min(energy, m_maximum_energy);
            int const index = energy >= m_maximum_energy
                                  ? m_energy_count - 2
                                  : amrex::bisect(m_energies, 0, m_energy_count - 1, energy);
            double fraction =
                (energy - m_energies[index]) / (m_energies[index + 1] - m_energies[index]);
            if (m_coordinates[index] != 0) {
                fraction = std::sqrt(fraction);
            }
            return {index, fraction,
                    (1 - fraction) * m_rates[index] + fraction * m_rates[index + 1], energy};
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Outcome
        sampleCell (std::uint32_t cell, double draw) const noexcept
        {
            auto const first = m_cell_offsets[cell];
            auto const count = m_cell_offsets[cell + 1] - first;
            draw = amrex::min(draw, 1.0 - 0.5 * std::numeric_limits<double>::epsilon());
            std::uint32_t outcome;
            if (!m_cumulative) {
                double const coordinate = draw * count;
                auto const index = static_cast<std::uint32_t>(coordinate);
                auto const entry = m_aliases[first + index];
                outcome = coordinate - index < entry.m_probability
                              ? entry.m_outcome
                              : m_aliases[first + entry.m_alternate].m_outcome;
            } else {
                auto lo = first;
                auto hi = first + count - 1;
                while (lo < hi) {
                    auto const mid = lo + (hi - lo) / 2;
                    if (draw < m_cdf[mid]) {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                outcome = m_outcome_ids[lo];
            }
            return m_outcomes[outcome];
        }

        /** Use optional offline bounds to shorten an exact CDF-grid search. */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE std::uint32_t
        quantileInterval (double const* coordinates, std::uint32_t const* lookup, int row,
                          std::uint32_t first, std::uint32_t last, double quantile) const noexcept
        {
            if (quantile >= 1) {
                return last - 1;
            }
            if (lookup != nullptr) {
                auto const bin = static_cast<int>(quantile * lookup_bins);
                auto const offset = row * (lookup_bins + 1) + bin;
                first = lookup[offset];
                last = lookup[offset + 1] + 1;
            }
            return amrex::bisect(coordinates, first, last, quantile);
        }

        /** The row mixture and its conditional uniform share a single draw.
         * The unchanged outcome still produces the sampled deflection.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Sample
        sample (Interpolation const &state, double angle_draw, double outcome_draw) const noexcept
        {
            double deflection;
            int row = state.m_index;
            double quantile = angle_draw;
            double changing = 0;
            if (state.m_energy >= m_rutherford_energy) {
                constexpr double rest = 510998.95069;
                constexpr double alpha = 7.2973525693e-3;
                double const tau = state.m_energy / rest;
                double const k2 = tau * (tau + 2) / (alpha * alpha);
                double const eta = 1 / (4 * m_screening_radius * m_screening_radius * k2);
                deflection = 2 * eta * angle_draw / (1 - angle_draw + eta);
            } else {
                double const lower = (1 - state.m_fraction) * m_rates[row];
                double const upper = state.m_fraction * m_rates[row + 1];
                double const draw = angle_draw * state.m_rate;
                if (upper == 0) {
                    quantile = angle_draw;
                } else if (lower == 0) {
                    ++row;
                    quantile = angle_draw;
                } else if (draw < lower) {
                    quantile = draw / lower;
                } else {
                    ++row;
                    quantile = (draw - lower) / upper;
                }
                quantile = amrex::max(0.0, amrex::min(quantile, 1.0));
                auto const first = m_angular_offsets[row];
                auto const last = m_angular_offsets[row + 1] - 1;
                auto const index =
                    quantileInterval(m_angular_u, m_angular_lookup, row, first, last, quantile);
                double const fraction =
                    (quantile - m_angular_u[index]) / (m_angular_u[index + 1] - m_angular_u[index]);
                deflection =
                    (1 - fraction) * m_deflection[index] + fraction * m_deflection[index + 1];
                changing = (1 - fraction) * m_changing[index] + fraction * m_changing[index + 1];
            }
            Outcome outcome;
            if (state.m_energy >= m_high_energy) {
                constexpr double rest = 510998.95069;
                constexpr double alpha = 7.2973525693e-3;
                double const q =
                    m_separation *
                    std::sqrt(state.m_energy * (state.m_energy + 2 * rest) * (0.5 * deflection)) /
                    (alpha * rest);
                int const bin = q >= m_high_edges[m_high_count]
                                    ? m_high_count - 1
                                    : amrex::bisect(m_high_edges, 0, m_high_count, q);
                outcome = sampleCell(m_high_cells[bin], outcome_draw);
            } else if (outcome_draw < changing) {
                auto const first = m_conditional_offsets[row];
                auto const last = m_conditional_offsets[row + 1] - 1;
                auto const index = quantileInterval(m_conditional_u, m_conditional_lookup, row,
                                                    first, last, quantile);
                double const fraction = (quantile - m_conditional_u[index]) /
                                        (m_conditional_u[index + 1] - m_conditional_u[index]);
                double draw = outcome_draw / changing;
                auto selected = index;
                if (draw < 1 - fraction) {
                    draw /= 1 - fraction;
                } else {
                    ++selected;
                    draw = fraction > 0 ? (draw - (1 - fraction)) / fraction : 0;
                }
                outcome = sampleCell(m_conditional_cells[selected], draw);
            }
            return {outcome, deflection, outcome.m_loss <= state.m_energy};
        }

        double const *m_energies = nullptr;
        double const *m_rates = nullptr;
        std::uint32_t const *m_coordinates = nullptr;
        std::uint32_t const *m_angular_offsets = nullptr;
        double const *m_angular_u = nullptr;
        double const *m_deflection = nullptr;
        float const *m_changing = nullptr;
        std::uint32_t const *m_conditional_offsets = nullptr;
        double const *m_conditional_u = nullptr;
        std::uint32_t const *m_conditional_cells = nullptr;
        std::uint32_t const* m_angular_lookup = nullptr;
        std::uint32_t const* m_conditional_lookup = nullptr;
        std::uint32_t const *m_cell_offsets = nullptr;
        Alias const *m_aliases = nullptr;
        double const *m_cdf = nullptr;
        std::uint16_t const *m_outcome_ids = nullptr;
        Outcome const *m_outcomes = nullptr;
        double const *m_high_edges = nullptr;
        std::uint32_t const *m_high_cells = nullptr;
        int m_energy_count = 0;
        int m_high_count = 0;
        double m_maximum_energy = 0;
        double m_high_energy = 1000;
        double m_rutherford_energy = 10000;
        double m_separation = 0;
        double m_screening_radius = 0;
        bool m_cumulative = false;
    };

    static std::shared_ptr<BackgroundMCCReciprocalRotation>
    get (std::string const &file, double temperature, bool cumulative);
    BackgroundMCCReciprocalRotation (std::string const &file, double temperature, bool cumulative);
    ~BackgroundMCCReciprocalRotation ();
    [[nodiscard]] Executor const &executor () const
    {
        return m_executor;
    }
    [[nodiscard]] Executor const &hostExecutor () const
    {
        return m_host_executor;
    }
    [[nodiscard]] std::vector<double> const &energies () const
    {
        return m_energies;
    }
    [[nodiscard]] double neutralMass () const
    {
        return m_neutral_mass;
    }
    [[nodiscard]] std::size_t tableBytes () const
    {
        return m_table_bytes;
    }
    void checkInclusiveRate (ScatteringProcess const &process) const;

private:
    struct Data;
    std::unique_ptr<Data> m_data;
    Executor m_executor;
    Executor m_host_executor;
    std::vector<double> m_energies, m_rates;
    std::vector<std::uint32_t> m_coordinates;
    double m_neutral_mass = 0;
    std::size_t m_table_bytes = 0;
};

#endif
