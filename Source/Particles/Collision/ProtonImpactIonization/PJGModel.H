/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_PJGMODEL_H_
#define WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_PJGMODEL_H_

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>

namespace ProtonImpactIonization
{
    enum class PJGTarget
    {
        N2,
        O2
    };

    /** Corrected Porter-Jackman-Green proton-impact ionization model.
     *
     * The model uses the 1977 erratum, exact relativistic two-body maximum
     * energy transfer, and the Bhabha hard-collision remainder. Its distortion
     * and normalization are refitted to recommended proton data. A table of
     * total cross sections and inverse SDCS quantiles is built on the host, so
     * device calls require only constant-time interpolation.
     */
    class PJGModel
    {
    public:
        static constexpr int table_energy_points = 256;
        static constexpr int table_quantile_points = 513;

        struct Executor
        {
            amrex::ParticleReal const* m_cross_section = nullptr;
            amrex::ParticleReal const* m_log_secondary_energy = nullptr;
            amrex::ParticleReal const* m_binding_energy = nullptr;
            amrex::ParticleReal m_log_projectile_energy_min = 0.0;
            amrex::ParticleReal m_inv_log_projectile_energy_step = 0.0;
            amrex::ParticleReal m_projectile_energy_min = 0.0;
            amrex::ParticleReal m_projectile_energy_max = 0.0;
            amrex::ParticleReal m_projectile_rest_energy = 0.0;

            /** Return the unit-charge total ionization cross section in m^2. */
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
            crossSection (amrex::ParticleReal projectile_energy) const noexcept;

            /** Sample ejected-electron kinetic and effective binding energies in eV. */
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
            sample (amrex::ParticleReal projectile_energy, amrex::ParticleReal quantile,
                    amrex::ParticleReal& secondary_energy,
                    amrex::ParticleReal& binding_energy) const noexcept;

            /** Return the exact free-electron maximum energy transfer in eV. */
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
            maximumEnergyTransfer (amrex::ParticleReal projectile_energy) const noexcept;

        private:
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
            energyCoordinate (amrex::ParticleReal projectile_energy, int& index,
                              amrex::ParticleReal& fraction) const noexcept;
        };

        PJGModel (PJGTarget target, amrex::ParticleReal projectile_rest_energy,
                  amrex::ParticleReal projectile_energy_min,
                  amrex::ParticleReal projectile_energy_max);

        [[nodiscard]] Executor executor () const noexcept;

        [[nodiscard]] static PJGTarget parseTarget (std::string const& target);
        [[nodiscard]] static std::string targetName (PJGTarget target);

        /** Host reference SDCS for a unit-charge projectile, in m^2/eV. */
        [[nodiscard]] static double differentialCrossSection (PJGTarget target,
                                                              double projectile_energy,
                                                              double secondary_energy,
                                                              double projectile_rest_energy);

        /** Host reference total cross section for a unit-charge projectile, in m^2. */
        [[nodiscard]] static double integratedCrossSection (PJGTarget target,
                                                            double projectile_energy,
                                                            double projectile_rest_energy);

        /** Exact relativistic maximum transfer to a free stationary electron, in eV. */
        [[nodiscard]] static double maximumEnergyTransfer (double projectile_energy,
                                                           double projectile_rest_energy);

    private:
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_cross_section;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_log_secondary_energy;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_binding_energy;
        amrex::ParticleReal m_log_projectile_energy_min = 0.0;
        amrex::ParticleReal m_inv_log_projectile_energy_step = 0.0;
        amrex::ParticleReal m_projectile_energy_min = 0.0;
        amrex::ParticleReal m_projectile_energy_max = 0.0;
        amrex::ParticleReal m_projectile_rest_energy = 0.0;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    PJGModel::Executor::energyCoordinate (amrex::ParticleReal const projectile_energy, int& index,
                                          amrex::ParticleReal& fraction) const noexcept
    {
        using namespace amrex::literals;
        using std::log;

        if (projectile_energy < m_projectile_energy_min ||
            projectile_energy > m_projectile_energy_max) {
            return false;
        }

        auto const coordinate = (log(projectile_energy) - m_log_projectile_energy_min) *
                                m_inv_log_projectile_energy_step;
        index = static_cast<int>(coordinate);
        index = amrex::max(0, amrex::min(index, table_energy_points - 2));
        fraction = amrex::max(
            0.0_prt, amrex::min(coordinate - static_cast<amrex::ParticleReal>(index), 1.0_prt));
        return true;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
    PJGModel::Executor::crossSection (amrex::ParticleReal const projectile_energy) const noexcept
    {
        using namespace amrex::literals;
        int energy_index;
        amrex::ParticleReal energy_fraction;
        if (!energyCoordinate(projectile_energy, energy_index, energy_fraction)) {
            return 0.0_prt;
        }

        return (1.0_prt - energy_fraction) * m_cross_section[energy_index] +
               energy_fraction * m_cross_section[energy_index + 1];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
    PJGModel::Executor::maximumEnergyTransfer (
        amrex::ParticleReal const projectile_energy) const noexcept
    {
        using namespace amrex::literals;

        constexpr auto electron_rest_energy = 510998.95069_prt;
        auto const gamma = 1.0_prt + projectile_energy / m_projectile_rest_energy;
        auto const beta_squared = 1.0_prt - 1.0_prt / (gamma * gamma);
        auto const mass_ratio = electron_rest_energy / m_projectile_rest_energy;
        return 2.0_prt * electron_rest_energy * beta_squared * gamma * gamma /
               (1.0_prt + 2.0_prt * gamma * mass_ratio + mass_ratio * mass_ratio);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    PJGModel::Executor::sample (amrex::ParticleReal const projectile_energy,
                                amrex::ParticleReal const quantile,
                                amrex::ParticleReal& secondary_energy,
                                amrex::ParticleReal& binding_energy) const noexcept
    {
        using namespace amrex::literals;
        using std::exp;

        int energy_index;
        amrex::ParticleReal energy_fraction;
        if (!energyCoordinate(projectile_energy, energy_index, energy_fraction)) {
            secondary_energy = 0.0_prt;
            binding_energy = 0.0_prt;
            return;
        }

        auto const bounded_quantile = amrex::max(0.0_prt, amrex::min(quantile, 1.0_prt));
        // Invert q = x^4/(x^4 + (1-x)^4). Uniform x resolves the long
        // Bhabha tail without increasing the device table or per-event search.
        auto const probability_root = std::sqrt(std::sqrt(bounded_quantile));
        auto const complement_root =
            std::sqrt(std::sqrt(amrex::max(1.0_prt - bounded_quantile, 0.0_prt)));
        auto const root_sum = probability_root + complement_root;
        auto const transformed_quantile =
            root_sum > 0.0_prt ? probability_root / root_sum : 0.0_prt;
        auto const quantile_coordinate =
            transformed_quantile * static_cast<amrex::ParticleReal>(table_quantile_points - 1);
        auto quantile_index = static_cast<int>(quantile_coordinate);
        quantile_index = amrex::max(0, amrex::min(quantile_index, table_quantile_points - 2));
        auto const quantile_fraction =
            amrex::max(0.0_prt, amrex::min(quantile_coordinate -
                                               static_cast<amrex::ParticleReal>(quantile_index),
                                           1.0_prt));

        auto const lower_offset = energy_index * table_quantile_points + quantile_index;
        auto const upper_offset = lower_offset + table_quantile_points;
        auto const lower_log_energy =
            (1.0_prt - quantile_fraction) * m_log_secondary_energy[lower_offset] +
            quantile_fraction * m_log_secondary_energy[lower_offset + 1];
        auto const upper_log_energy =
            (1.0_prt - quantile_fraction) * m_log_secondary_energy[upper_offset] +
            quantile_fraction * m_log_secondary_energy[upper_offset + 1];
        auto const log_energy =
            (1.0_prt - energy_fraction) * lower_log_energy + energy_fraction * upper_log_energy;
        secondary_energy =
            amrex::min(exp(log_energy) - 1.0_prt, maximumEnergyTransfer(projectile_energy));

        auto const lower_binding = (1.0_prt - quantile_fraction) * m_binding_energy[lower_offset] +
                                   quantile_fraction * m_binding_energy[lower_offset + 1];
        auto const upper_binding = (1.0_prt - quantile_fraction) * m_binding_energy[upper_offset] +
                                   quantile_fraction * m_binding_energy[upper_offset + 1];
        binding_energy =
            (1.0_prt - energy_fraction) * lower_binding + energy_fraction * upper_binding;
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_PJGMODEL_H_
