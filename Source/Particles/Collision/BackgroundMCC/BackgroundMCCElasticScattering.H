/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_SCATTERING_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_SCATTERING_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>

/** Precomputed inverse angular CDF for an IAA/elmolcs elastic DCS table.
 *
 * Below 10 keV, differential cross sections are linear in sin(theta/2) and
 * log(E), following the IAA interpolation scheme. The inverse CDF stores
 * 1-cos(theta), which preserves small forward angles in single precision.
 * Expensive integration and CDF inversion are performed once on the host. For
 * recognized N2 and O2 elmolcs tables, the unresolved high-energy forward lobe
 * is sampled from the analytic screened-Rutherford continuation. Runtime
 * sampling has no allocation or iterative solve.
 */
class BackgroundMCCElasticScatteringModel
{
public:
    static constexpr int quantile_grid_size = 1025;

    struct Executor
    {
        /** Sample the cosine of the neutral-rest-frame electron scattering angle.
         *
         * @param collision_energy electron collision energy in eV
         * @param angle_draw uniform random number in [0, 1)
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
        sampleCosine (amrex::ParticleReal const collision_energy,
                      double const angle_draw) const noexcept
        {
            using namespace amrex::literals;

            if (m_screening_radius > 0.0 &&
                collision_energy >= m_screened_rutherford_energy)
            {
                // Schmalzried Eq. (3.11), with the molecular screening radii
                // fitted in Table 11.4. Here k is in inverse Bohr radii.
                constexpr double fine_structure = 7.2973525693e-3;
                constexpr double electron_rest_energy = 510998.95069; // eV
                auto const tau = static_cast<double>(collision_energy) /
                    electron_rest_energy;
                auto const k_sq = tau * (tau + 2.0) /
                    (fine_structure * fine_structure);
                auto const eta = 1.0 /
                    (4.0 * m_screening_radius * m_screening_radius * k_sq);
                auto const deflection = 2.0 * eta * angle_draw /
                    (1.0 - angle_draw + eta);
                return amrex::max(-1.0, amrex::min(1.0 - deflection, 1.0));
            }

            if (m_energy_grid_size < 2)
            {
                return 1.0 - 2.0 * angle_draw;
            }

            int row_index = 0;
            double row_draw = angle_draw;
            if (collision_energy >= m_energy_hi)
            {
                row_index = m_energy_grid_size - 1;
            }
            else if (collision_energy > m_energy_lo)
            {
                auto const log_energy = static_cast<amrex::ParticleReal>(
                    std::log(static_cast<double>(collision_energy)));
                auto const energy_index = amrex::bisect(
                    m_log_energies, 0, m_energy_grid_size - 1, log_energy);
                auto const energy_fraction = static_cast<double>(
                    (log_energy - m_log_energies[energy_index]) /
                    (m_log_energies[energy_index + 1] -
                     m_log_energies[energy_index]));

                // A DCS interpolated linearly in log(E) is an unnormalized
                // mixture of the two adjacent angular rows. Select one row
                // with its integrated weight, then reuse the same variate as
                // a conditional draw. This is exact and loads only one row.
                auto const low_weight = (1.0 - energy_fraction) *
                    static_cast<double>(m_row_integrals[energy_index]);
                auto const high_weight = energy_fraction *
                    static_cast<double>(m_row_integrals[energy_index + 1]);
                auto const weighted_draw = angle_draw * (low_weight + high_weight);
                if (weighted_draw < low_weight)
                {
                    row_index = energy_index;
                    row_draw = low_weight > 0.0 ? weighted_draw / low_weight : 0.0;
                }
                else
                {
                    row_index = energy_index + 1;
                    row_draw = high_weight > 0.0
                        ? (weighted_draw - low_weight) / high_weight
                        : 0.0;
                }
            }

            // The transformed grid q = 1-(1-x)^2 resolves rare large-angle
            // events without enlarging the table. Invert it here to locate the
            // two stored quantiles bracketing the uniform draw q.
            auto const transformed_draw = 1.0 -
                std::sqrt(amrex::max(1.0 - row_draw, 0.0));
            auto quantile_coordinate = transformed_draw *
                static_cast<double>(quantile_grid_size - 1);
            quantile_coordinate = amrex::max(
                0.0, amrex::min(quantile_coordinate,
                                static_cast<double>(quantile_grid_size - 1)));
            int quantile_index =
                static_cast<int>(amrex::Math::floor(quantile_coordinate));
            auto quantile_fraction = quantile_coordinate - quantile_index;
            if (quantile_index >= quantile_grid_size - 1)
            {
                quantile_index = quantile_grid_size - 2;
                quantile_fraction = 1.0;
            }

            auto const index = row_index * quantile_grid_size + quantile_index;
            auto const deflection = static_cast<double>(m_inverse_cdf_deflection[index]) +
                static_cast<double>(quantile_fraction) *
                    (static_cast<double>(m_inverse_cdf_deflection[index + 1]) -
                     static_cast<double>(m_inverse_cdf_deflection[index]));
            return amrex::max(-1.0, amrex::min(1.0 - deflection, 1.0));
        }

        amrex::ParticleReal const* m_log_energies = nullptr;
        amrex::ParticleReal const* m_row_integrals = nullptr;
        amrex::ParticleReal const* m_inverse_cdf_deflection = nullptr;
        int m_energy_grid_size = 0;
        amrex::ParticleReal m_energy_lo = 0.0;
        amrex::ParticleReal m_energy_hi = 0.0;
        double m_screening_radius = 0.0;
        static constexpr amrex::ParticleReal m_screened_rutherford_energy = 1.0e4;
    };

    explicit BackgroundMCCElasticScatteringModel (std::string const& differential_cross_section);

    [[nodiscard]] Executor const& executor () const
    {
#ifdef AMREX_USE_GPU
        return m_executor_d;
#else
        return m_executor_h;
#endif
    }

private:
    amrex::Gpu::HostVector<amrex::ParticleReal> m_log_energies_h;
    amrex::Gpu::HostVector<amrex::ParticleReal> m_row_integrals_h;
    amrex::Gpu::HostVector<amrex::ParticleReal> m_inverse_cdf_deflection_h;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_log_energies_d;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_row_integrals_d;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_inverse_cdf_deflection_d;
    Executor m_executor_d;
#endif
    Executor m_executor_h;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_SCATTERING_H_
