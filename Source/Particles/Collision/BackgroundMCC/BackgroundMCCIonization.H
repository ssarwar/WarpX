/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_IONIZATION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_IONIZATION_H_

#include "Utils/WarpXAlgorithmSelection.H"

#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>

/** Molecular target supported by the RBEQ ionization model. */
enum class BackgroundMCCIonizationTarget
{
    None,
    N2,
    O2
};

/** Precomputed RBEQ shell-selection and inverse-CDF data.
 *
 * The relatively expensive RBEQ logarithms and CDF inversions are evaluated
 * once on the host. An accepted ionization event then needs one logarithm,
 * fixed-size interpolation and no device allocation or iterative root solve.
 */
class BackgroundMCCIonizationModel
{
public:
    static constexpr int max_shell_count = 6;
    static constexpr int energy_grid_size = 257;
    // A symmetric fourth-power probability transform resolves both inverse-CDF
    // tails. 513 coordinates retain the former central probability spacing.
    static constexpr int quantile_grid_size = 513;

    struct Executor
    {
        /** Sample the binding-energy loss and secondary-electron energy.
         *
         * @param incident_energy incident electron energy in eV
         * @param fixed_binding_energy fixed loss for the equal-sharing model,
         * in eV
         * @param shell_draw uniform random number for subshell selection
         * @param energy_draw uniform random number for inverse-CDF sampling
         * @param binding_energy sampled binding-energy loss in eV
         * @param secondary_energy sampled secondary-electron energy in eV
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
        sampleEnergy (amrex::ParticleReal const incident_energy,
                      amrex::ParticleReal const fixed_binding_energy,
                      amrex::ParticleReal const shell_draw, amrex::ParticleReal const energy_draw,
                      amrex::ParticleReal& binding_energy,
                      amrex::ParticleReal& secondary_energy) const noexcept
        {
            using namespace amrex::literals;

            if (m_model != IonizationEnergySharingModel::RBEQ)
            {
                binding_energy = fixed_binding_energy;
                secondary_energy = 0.5_prt * amrex::max(incident_energy - binding_energy, 0.0_prt);
                return;
            }

            amrex::ParticleReal energy_coordinate = 0.0_prt;
            if (incident_energy > m_energy_min)
            {
                energy_coordinate =
                    (std::log(incident_energy) - m_log_energy_min) * m_inverse_log_energy_step;
            }
            energy_coordinate = amrex::max(
                0.0_prt, amrex::min(energy_coordinate,
                                    static_cast<amrex::ParticleReal>(m_energy_grid_size - 1)));

            int energy_index = static_cast<int>(amrex::Math::floor(energy_coordinate));
            amrex::ParticleReal energy_fraction = energy_coordinate - energy_index;
            if (energy_index >= m_energy_grid_size - 1)
            {
                energy_index = m_energy_grid_size - 2;
                energy_fraction = 1.0_prt;
            }

            amrex::ParticleReal shell_weights[max_shell_count];
            amrex::ParticleReal total_weight = 0.0_prt;
            for (int shell = 0; shell < m_shell_count; ++shell)
            {
                auto const index_lo = energy_index * max_shell_count + shell;
                auto const index_hi = index_lo + max_shell_count;
                auto const interpolated = m_shell_probabilities[index_lo] +
                                          energy_fraction * (m_shell_probabilities[index_hi] -
                                                             m_shell_probabilities[index_lo]);
                shell_weights[shell] = incident_energy > m_binding_energies[shell]
                                           ? amrex::max(interpolated, 0.0_prt)
                                           : 0.0_prt;
                total_weight += shell_weights[shell];
            }

            int selected_shell = m_shell_count - 1;
            if (total_weight > 0.0_prt)
            {
                auto const target_weight = shell_draw * total_weight;
                amrex::ParticleReal cumulative_weight = 0.0_prt;
                for (int shell = 0; shell < m_shell_count; ++shell)
                {
                    cumulative_weight += shell_weights[shell];
                    if (target_weight <= cumulative_weight)
                    {
                        selected_shell = shell;
                        break;
                    }
                }
            }
            else
            {
                // This path is only reachable for a cross-section table that
                // accepts an event below every positive RBEQ partial. Choose
                // the outer shell so the outcome remains finite and respects
                // the lowest threshold.
                for (int shell = m_shell_count - 1; shell >= 0; --shell)
                {
                    if (incident_energy > m_binding_energies[shell])
                    {
                        selected_shell = shell;
                        break;
                    }
                }
            }

            binding_energy = m_binding_energies[selected_shell];
            auto const available_energy = amrex::max(incident_energy - binding_energy, 0.0_prt);

            // Invert q = x^4/(x^4 + (1-x)^4). The transformed coordinate
            // prevents the physical W_max endpoint from contaminating the
            // high-energy distribution while retaining fixed device work.
            auto const probability = amrex::max(
                0.0_prt, amrex::min(energy_draw, 1.0_prt));
            auto const probability_root = std::sqrt(std::sqrt(probability));
            auto const complement_root =
                std::sqrt(std::sqrt(amrex::max(1.0_prt - probability, 0.0_prt)));
            auto const root_sum = probability_root + complement_root;
            auto const transformed_draw =
                root_sum > 0.0_prt ? probability_root / root_sum : 0.0_prt;
            auto quantile_coordinate = transformed_draw *
                static_cast<amrex::ParticleReal>(m_quantile_grid_size - 1);
            quantile_coordinate = amrex::max(
                0.0_prt, amrex::min(quantile_coordinate,
                                    static_cast<amrex::ParticleReal>(m_quantile_grid_size - 1)));
            int quantile_index = static_cast<int>(amrex::Math::floor(quantile_coordinate));
            amrex::ParticleReal quantile_fraction = quantile_coordinate - quantile_index;
            if (quantile_index >= m_quantile_grid_size - 1)
            {
                quantile_index = m_quantile_grid_size - 2;
                quantile_fraction = 1.0_prt;
            }

            auto const row_size = max_shell_count * m_quantile_grid_size;
            auto const index =
                energy_index * row_size + selected_shell * m_quantile_grid_size + quantile_index;
            auto const fraction_lo =
                m_inverse_cdf[index] +
                quantile_fraction * (m_inverse_cdf[index + 1] - m_inverse_cdf[index]);
            auto const fraction_hi = m_inverse_cdf[index + row_size] +
                                     quantile_fraction * (m_inverse_cdf[index + row_size + 1] -
                                                          m_inverse_cdf[index + row_size]);
            auto const secondary_fraction = amrex::max(
                0.0_prt,
                amrex::min(fraction_lo + energy_fraction * (fraction_hi - fraction_lo), 1.0_prt));

            secondary_energy = 0.5_prt * available_energy * secondary_fraction;
        }

        IonizationEnergySharingModel m_model = IonizationEnergySharingModel::Equal;
        int m_shell_count = 0;
        int m_energy_grid_size = 0;
        int m_quantile_grid_size = 0;
        amrex::GpuArray<amrex::ParticleReal, max_shell_count> m_binding_energies{};
        amrex::ParticleReal const* m_shell_probabilities = nullptr;
        amrex::ParticleReal const* m_inverse_cdf = nullptr;
        amrex::ParticleReal m_energy_min = 0.0;
        amrex::ParticleReal m_log_energy_min = 0.0;
        amrex::ParticleReal m_inverse_log_energy_step = 0.0;
    };

    BackgroundMCCIonizationModel (BackgroundMCCIonizationTarget target,
                                  amrex::ParticleReal maximum_energy);

    [[nodiscard]] Executor const& executor () const
    {
#ifdef AMREX_USE_GPU
        return m_executor_d;
#else
        return m_executor_h;
#endif
    }

    [[nodiscard]] static BackgroundMCCIonizationTarget parseTarget (std::string const& target);

    [[nodiscard]] static amrex::ParticleReal
    outerBindingEnergy (BackgroundMCCIonizationTarget target);

private:
    amrex::Gpu::HostVector<amrex::ParticleReal> m_shell_probabilities_h;
    amrex::Gpu::HostVector<amrex::ParticleReal> m_inverse_cdf_h;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_shell_probabilities_d;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_inverse_cdf_d;
    Executor m_executor_d;
#endif
    Executor m_executor_h;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_IONIZATION_H_
