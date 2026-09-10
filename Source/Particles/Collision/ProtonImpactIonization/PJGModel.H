/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_PJGMODEL_H_
#define WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_PJGMODEL_H_

#include "RelativisticIonization.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>
#include <string>

namespace ProtonImpactIonization
{
    enum class PJGTarget
    {
        N2,
        O2
    };

    /** Calibrated PJG-type inclusive electron source for N2 and O2.
     *
     * The positive soft Lorentzian difference and fixed complete Bhabha
     * hard kernel share a smooth molecular gate, including the bound tail.
     * Initialization integrates the same SDCS for its total and inverse CDF.
     * Device events use fixed-size interpolation, not quadrature or root solves.
     * See the proton-impact-ionization theory documentation for the calibrated
     * range, original-model comparison, sources, and effective-pair convention.
     */
    class PJGModel
    {
      public:
        static constexpr int table_energy_points = 256;
        static constexpr int table_quantile_points = 1025;

        template <typename Real> struct ExecutorT
        {
            Real const* m_cross_section = nullptr;
            Real const* m_log_secondary_energy = nullptr;
            Real const* m_binding_energy = nullptr;
            Real m_log_projectile_energy_min = 0.0;
            Real m_inv_log_projectile_energy_step = 0.0;
            Real m_projectile_energy_min = 0.0;
            Real m_projectile_energy_max = 0.0;
            Real m_projectile_rest_energy = 0.0;
            Real m_neutral_rest_energy = 0.0;
            Real m_minimum_binding_energy = 0.0;

            /** Return the unit-charge total ionization cross section in m^2. */
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
            crossSection (Real projectile_energy) const noexcept;

            /** Sample ejected-electron kinetic and effective binding energies in eV. */
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
            sample (Real projectile_energy, double quantile, Real& secondary_energy,
                    Real& binding_energy) const noexcept;

            /** Return the exact free-electron maximum energy transfer in eV. */
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
            maximumEnergyTransfer (Real projectile_energy) const noexcept;

          private:
            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
            energyCoordinate (Real projectile_energy, int& index, Real& fraction) const noexcept;
        };

        using Executor = ExecutorT<amrex::ParticleReal>;

        struct Moments
        {
            double m_total = 0.0;          //!< Integral S dT, m^2.
            double m_kinetic = 0.0;        //!< Integral T S dT, eV m^2.
            double m_kinetic_second = 0.0; //!< Integral T^2 S dT, eV^2 m^2.
            double m_binding = 0.0;        //!< Integral I_eff S dT, eV m^2.
            double m_above_free = 0.0;     //!< Integral above free Tmax, m^2.
        };

        PJGModel (PJGTarget target, amrex::ParticleReal projectile_rest_energy,
                  amrex::ParticleReal projectile_energy_min,
                  amrex::ParticleReal projectile_energy_max);

        [[nodiscard]] Executor executor () const noexcept;

        [[nodiscard]] static PJGTarget parseTarget (std::string const& target);
        [[nodiscard]] static std::string targetName (PJGTarget target);

        /** Host reference SDCS in m^2/eV; zero outside the calibrated incident range. */
        [[nodiscard]] static double differentialCrossSection (PJGTarget target,
                                                              double projectile_energy,
                                                              double secondary_energy,
                                                              double projectile_rest_energy);

        /** Host reference total cross section for a unit-charge projectile, in m^2. */
        [[nodiscard]] static double integratedCrossSection (PJGTarget target,
                                                            double projectile_energy,
                                                            double projectile_rest_energy);

        /** Host energy moments from positive quadrature of the calibrated SDCS. */
        [[nodiscard]] static Moments integratedMoments (PJGTarget target, double projectile_energy,
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
        amrex::ParticleReal m_neutral_rest_energy = 0.0;
        amrex::ParticleReal m_minimum_binding_energy = 0.0;
    };

    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    PJGModel::ExecutorT<Real>::energyCoordinate (Real const projectile_energy, int& index,
                                                 Real& fraction) const noexcept
    {
        using std::log;

        constexpr auto tolerance = Real(8) * std::numeric_limits<Real>::epsilon();
        if (!(projectile_energy >= m_projectile_energy_min * (Real(1) - tolerance) &&
              projectile_energy <= m_projectile_energy_max * (Real(1) + tolerance))) {
            return false;
        }

        auto const coordinate = (log(projectile_energy) - m_log_projectile_energy_min) *
                                m_inv_log_projectile_energy_step;
        index = static_cast<int>(coordinate);
        index = amrex::max(0, amrex::min(index, table_energy_points - 2));
        fraction = amrex::max(Real(0), amrex::min(coordinate - static_cast<Real>(index), Real(1)));
        return true;
    }

    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    PJGModel::ExecutorT<Real>::crossSection (Real const projectile_energy) const noexcept
    {
        int energy_index;
        Real energy_fraction;
        if (!energyCoordinate(projectile_energy, energy_index, energy_fraction)) {
            return Real(0);
        }

        return (Real(1) - energy_fraction) * m_cross_section[energy_index] +
               energy_fraction * m_cross_section[energy_index + 1];
    }

    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    PJGModel::ExecutorT<Real>::maximumEnergyTransfer (Real const projectile_energy) const noexcept
    {
        return relativisticMaximumTransfer(projectile_energy, m_projectile_rest_energy);
    }

    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    PJGModel::ExecutorT<Real>::sample (Real const projectile_energy, double const quantile,
                                       Real& secondary_energy, Real& binding_energy) const noexcept
    {
        using std::expm1;

        int energy_index;
        Real energy_fraction;
        if (!(quantile >= 0.0 && quantile <= 1.0) ||
            !energyCoordinate(projectile_energy, energy_index, energy_fraction)) {
            secondary_energy = Real(0);
            binding_energy = Real(0);
            return;
        }

        // Invert q = x^4/(x^4 + (1-x)^4). Endpoint clustering resolves
        // hard-electron moments; there is no per-event CDF search. Form
        // 1-q in double even with float particles to retain rare hard tails.
        auto const probability_root = std::sqrt(std::sqrt(static_cast<Real>(quantile)));
        auto const complement_root = std::sqrt(std::sqrt(static_cast<Real>(1.0 - quantile)));
        auto const root_sum = probability_root + complement_root;
        auto const transformed_quantile =
            root_sum > Real(0) ? probability_root / root_sum : Real(0);
        auto const quantile_coordinate =
            transformed_quantile * static_cast<Real>(table_quantile_points - 1);
        auto quantile_index = static_cast<int>(quantile_coordinate);
        quantile_index = amrex::max(0, amrex::min(quantile_index, table_quantile_points - 2));
        auto const quantile_fraction = amrex::max(
            Real(0), amrex::min(quantile_coordinate - static_cast<Real>(quantile_index), Real(1)));

        auto const lower_offset = energy_index * table_quantile_points + quantile_index;
        auto const upper_offset = lower_offset + table_quantile_points;
        auto const lower_log_energy =
            (Real(1) - quantile_fraction) * m_log_secondary_energy[lower_offset] +
            quantile_fraction * m_log_secondary_energy[lower_offset + 1];
        auto const upper_log_energy =
            (Real(1) - quantile_fraction) * m_log_secondary_energy[upper_offset] +
            quantile_fraction * m_log_secondary_energy[upper_offset + 1];
        auto const log_energy =
            (Real(1) - energy_fraction) * lower_log_energy + energy_fraction * upper_log_energy;
        // Guard the molecular endpoint against inter-row roundoff, never
        // clip the physical bound-electron tail at the free binary limit.
        secondary_energy = amrex::min(
            expm1(log_energy),
            molecularMaximumSecondaryEnergy(projectile_energy, m_projectile_rest_energy,
                                            m_neutral_rest_energy, m_minimum_binding_energy));

        auto const lower_binding = (Real(1) - quantile_fraction) * m_binding_energy[lower_offset] +
                                   quantile_fraction * m_binding_energy[lower_offset + 1];
        auto const upper_binding = (Real(1) - quantile_fraction) * m_binding_energy[upper_offset] +
                                   quantile_fraction * m_binding_energy[upper_offset + 1];
        binding_energy =
            (Real(1) - energy_fraction) * lower_binding + energy_fraction * upper_binding;
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_PJGMODEL_H_
