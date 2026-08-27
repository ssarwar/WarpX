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
 * Differential cross sections are integrated in polar-angle space with the
 * sin(theta) solid-angle Jacobian. The inverse CDF stores 1-cos(theta), which
 * preserves small forward angles in single precision. Expensive integration
 * and CDF inversion are performed once on the host. Each accepted device event
 * uses a bisection in the energy grid followed by fixed-size interpolation,
 * without allocation or an iterative solve.
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
                      amrex::ParticleReal const angle_draw) const noexcept
        {
            using namespace amrex::literals;

            if (m_energy_grid_size < 2)
            {
                return 1.0 - 2.0 * static_cast<double>(angle_draw);
            }

            int energy_index = 0;
            amrex::ParticleReal energy_fraction = 0.0_prt;
            if (collision_energy >= m_energy_hi)
            {
                energy_index = m_energy_grid_size - 2;
                energy_fraction = 1.0_prt;
            }
            else if (collision_energy > m_energy_lo)
            {
                energy_index = amrex::bisect(
                    m_energies, 0, m_energy_grid_size - 1, collision_energy);
                energy_fraction =
                    (collision_energy - m_energies[energy_index]) /
                    (m_energies[energy_index + 1] - m_energies[energy_index]);
            }

            // The transformed grid q = 1-(1-x)^2 resolves rare large-angle
            // events without enlarging the table. Invert it here to locate the
            // two stored quantiles bracketing the uniform draw q.
            auto const transformed_draw = 1.0_prt -
                std::sqrt(amrex::max(1.0_prt - angle_draw, 0.0_prt));
            auto quantile_coordinate = transformed_draw *
                static_cast<amrex::ParticleReal>(quantile_grid_size - 1);
            quantile_coordinate = amrex::max(
                0.0_prt,
                amrex::min(quantile_coordinate,
                           static_cast<amrex::ParticleReal>(quantile_grid_size - 1)));
            int quantile_index =
                static_cast<int>(amrex::Math::floor(quantile_coordinate));
            auto quantile_fraction = quantile_coordinate - quantile_index;
            if (quantile_index >= quantile_grid_size - 1)
            {
                quantile_index = quantile_grid_size - 2;
                quantile_fraction = 1.0_prt;
            }

            auto const index = energy_index * quantile_grid_size + quantile_index;
            auto const deflection_lo = static_cast<double>(m_inverse_cdf_deflection[index]) +
                static_cast<double>(quantile_fraction) *
                    (static_cast<double>(m_inverse_cdf_deflection[index + 1]) -
                     static_cast<double>(m_inverse_cdf_deflection[index]));
            auto const deflection_hi =
                static_cast<double>(m_inverse_cdf_deflection[index + quantile_grid_size]) +
                static_cast<double>(quantile_fraction) *
                    (static_cast<double>(
                         m_inverse_cdf_deflection[index + quantile_grid_size + 1]) -
                     static_cast<double>(
                         m_inverse_cdf_deflection[index + quantile_grid_size]));
            auto const deflection = deflection_lo + static_cast<double>(energy_fraction) *
                (deflection_hi - deflection_lo);
            return amrex::max(-1.0, amrex::min(1.0 - deflection, 1.0));
        }

        amrex::ParticleReal const* m_energies = nullptr;
        amrex::ParticleReal const* m_inverse_cdf_deflection = nullptr;
        int m_energy_grid_size = 0;
        amrex::ParticleReal m_energy_lo = 0.0;
        amrex::ParticleReal m_energy_hi = 0.0;
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
    amrex::Gpu::HostVector<amrex::ParticleReal> m_energies_h;
    amrex::Gpu::HostVector<amrex::ParticleReal> m_inverse_cdf_deflection_h;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_energies_d;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_inverse_cdf_deflection_d;
    Executor m_executor_d;
#endif
    Executor m_executor_h;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_SCATTERING_H_
