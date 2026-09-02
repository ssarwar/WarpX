/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PJGModel.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace ProtonImpactIonization
{
    namespace
    {
        constexpr int max_continua = 7;
        constexpr double square_centimeter_to_square_meter = 1.0e-4;
        constexpr double electron_rest_energy = 510998.95069;
        constexpr double euler_number = 2.71828182845904523536;

        struct Parameters
        {
            int num_continua;
            int num_electrons;
            double b0;
            double b1;
            double e0;
            double t1;
            double gamma1_fixed;
            double k;
            double gamma_s;
            double gamma_numerator;
            double gamma_denominator;
            double t_s;
            double t_numerator;
            double t_denominator;
            double j;
            double nu;
            std::array<double, max_continua> thresholds;
            std::array<double, max_continua> fractions;
            std::array<double, max_continua> bethe_constants;
        };

        // J is refitted to the Rudd recommended total over 5--4000 keV.
        constexpr Parameters n2_parameters{6,
                                           14,
                                           0.029,
                                           1.035,
                                           8239.0,
                                           53.3,
                                           115.0,
                                           7.58e-16,
                                           11.1,
                                           1.27e4,
                                           1.81e3,
                                           4.0,
                                           2.03e4,
                                           1.97e3,
                                           59.65,
                                           -1.93e-1,
                                           {15.58, 16.73, 18.75, 22.0, 23.6, 40.0, 0.0},
                                           {0.456, 0.2, 0.104, 0.07, 0.07, 0.1, 0.0},
                                           {2.48, 2.66, 2.99, 3.50, 3.76, 6.37, 0.0}};

        // O2 additionally requires K because J cannot change the asymptote.
        constexpr Parameters o2_parameters{7,
                                           16,
                                           0.030,
                                           1.035,
                                           8239.0,
                                           68.3,
                                           189.1,
                                           4.581e-16,
                                           13.1,
                                           5.0e5,
                                           7.60e4,
                                           6.34,
                                           2.52e3,
                                           1.28e2,
                                           19.28,
                                           3.14e-1,
                                           {12.1, 16.1, 16.9, 18.2, 20.3, 23.0, 37.0},
                                           {0.08, 0.19, 0.19, 0.17, 0.11, 0.16, 0.1},
                                           {1.93, 2.56, 2.69, 2.90, 3.23, 3.66, 5.89}};

        Parameters const&
        parameters (PJGTarget const target)
        {
            return target == PJGTarget::N2 ? n2_parameters : o2_parameters;
        }

        double
        piElectronChargeFourth ()
        {
            auto const classical_radius_cm = PhysConst::r_e_v<double> * 100.0;
            return static_cast<double>(MathConst::pi) *
                   std::pow(classical_radius_cm * electron_rest_energy, 2);
        }

        struct ProjectileState
        {
            double gamma;
            double beta_squared;
            double equivalent_electron_energy;
            double maximum_transfer;
            double total_energy;
        };

        ProjectileState
        projectileState (double const projectile_energy, double const projectile_rest_energy)
        {
            auto const gamma = 1.0 + projectile_energy / projectile_rest_energy;
            auto const beta_squared = 1.0 - 1.0 / (gamma * gamma);
            auto const mass_ratio = electron_rest_energy / projectile_rest_energy;
            auto const maximum_transfer =
                2.0 * electron_rest_energy * beta_squared * gamma * gamma /
                (1.0 + 2.0 * gamma * mass_ratio + mass_ratio * mass_ratio);
            return {gamma, beta_squared, 0.5 * electron_rest_energy * beta_squared,
                    maximum_transfer, projectile_energy + projectile_rest_energy};
        }

        struct EnergyFunctions
        {
            double gamma_width;
            double peak_energy;
            double relativistic_reduction;
            double distortion;
        };

        EnergyFunctions
        energyFunctions (Parameters const& p, ProjectileState const& state)
        {
            auto const equivalent_energy = state.equivalent_electron_energy;
            auto const gamma_width =
                p.gamma_numerator / (equivalent_energy + p.gamma_denominator) + p.gamma_s;
            auto const peak_energy = p.t_s - p.t_numerator / (equivalent_energy + p.t_denominator);
            auto const log_energy = std::log(equivalent_energy / p.e0);
            auto const relativistic_reduction = p.b0 * (log_energy * log_energy + p.b1);
            auto const power = p.nu + 1.0;
            auto const energy_power = std::pow(equivalent_energy, power);
            auto const distortion = energy_power / (std::pow(p.j, power) + energy_power);
            return {gamma_width, peak_energy, relativistic_reduction, distortion};
        }

        double
        continuumDifferentialCrossSectionCm (Parameters const& p, int const continuum,
                                             ProjectileState const& state,
                                             EnergyFunctions const& functions,
                                             double const secondary_energy)
        {
            auto const threshold = p.thresholds[continuum];
            auto const amplitude =
                p.fractions[continuum] * functions.distortion / state.equivalent_electron_energy;
            auto const gamma_squared = functions.gamma_width * functions.gamma_width;
            auto const line_shape =
                1.0 / (std::pow(secondary_energy - functions.peak_energy, 2) + gamma_squared) -
                functions.relativistic_reduction /
                    (std::pow(secondary_energy - p.t1, 2) + p.gamma1_fixed * p.gamma1_fixed);
            auto const bethe_factor =
                std::log(4.0 * state.equivalent_electron_energy * p.bethe_constants[continuum] /
                             (threshold * (1.0 - state.beta_squared)) +
                         euler_number) -
                state.beta_squared;
            auto const soft_term = p.k * gamma_squared * bethe_factor * line_shape;

            // Bhabha's exact hard-collision remainder after the common PJG
            // 1/(m beta^2 c^2/2) factor is extracted. The threshold shift is
            // the remaining PJG bound-electron continuation.
            auto const hard_term =
                static_cast<double>(p.num_electrons) * piElectronChargeFourth() *
                (1.0 / (2.0 * state.total_energy * state.total_energy) -
                 state.beta_squared / ((state.maximum_transfer + threshold) *
                                       (secondary_energy + threshold)));
            return amplitude * (soft_term + hard_term);
        }

        double
        integratedContinuumCrossSectionCm (Parameters const& p, int const continuum,
                                           ProjectileState const& state,
                                           EnergyFunctions const& functions,
                                           double const upper_energy)
        {
            auto const threshold = p.thresholds[continuum];
            auto const amplitude =
                p.fractions[continuum] * functions.distortion / state.equivalent_electron_energy;
            auto const gamma_width = functions.gamma_width;
            auto const gamma_squared = gamma_width * gamma_width;
            auto const line_integral =
                (std::atan((upper_energy - functions.peak_energy) / gamma_width) -
                 std::atan(-functions.peak_energy / gamma_width)) /
                    gamma_width -
                functions.relativistic_reduction *
                    (std::atan((upper_energy - p.t1) / p.gamma1_fixed) -
                     std::atan(-p.t1 / p.gamma1_fixed)) /
                    p.gamma1_fixed;
            auto const bethe_factor =
                std::log(4.0 * state.equivalent_electron_energy * p.bethe_constants[continuum] /
                             (threshold * (1.0 - state.beta_squared)) +
                         euler_number) -
                state.beta_squared;
            auto const soft_integral = p.k * gamma_squared * bethe_factor * line_integral;
            auto const hard_integral =
                static_cast<double>(p.num_electrons) * piElectronChargeFourth() *
                (upper_energy / (2.0 * state.total_energy * state.total_energy) -
                 state.beta_squared / (state.maximum_transfer + threshold) *
                     std::log1p(upper_energy / threshold));
            return amplitude * (soft_integral + hard_integral);
        }

        double
        cumulativeCrossSectionCm (PJGTarget const target, double const projectile_energy,
                                  double const secondary_energy,
                                  double const projectile_rest_energy)
        {
            if (projectile_energy <= 0.0 || secondary_energy <= 0.0) {
                return 0.0;
            }
            auto const& p = parameters(target);
            auto const state = projectileState(projectile_energy, projectile_rest_energy);
            auto const functions = energyFunctions(p, state);
            auto const upper_energy = std::min(secondary_energy, state.maximum_transfer);
            double result = 0.0;
            for (int continuum = 0; continuum < p.num_continua; ++continuum) {
                result +=
                    integratedContinuumCrossSectionCm(p, continuum, state, functions, upper_energy);
            }
            return result;
        }

        double
        effectiveBindingEnergy (PJGTarget const target, double const projectile_energy,
                                double const secondary_energy, double const projectile_rest_energy)
        {
            auto const& p = parameters(target);
            auto const state = projectileState(projectile_energy, projectile_rest_energy);
            auto const functions = energyFunctions(p, state);
            double weighted_binding = 0.0;
            double positive_sdcs = 0.0;
            for (int continuum = 0; continuum < p.num_continua; ++continuum) {
                auto const contribution =
                    std::max(continuumDifferentialCrossSectionCm(p, continuum, state, functions,
                                                                 secondary_energy),
                             0.0);
                positive_sdcs += contribution;
                weighted_binding += contribution * p.thresholds[continuum];
            }
            return positive_sdcs > 0.0 ? weighted_binding / positive_sdcs : p.thresholds[0];
        }
    } // namespace

    PJGModel::PJGModel (PJGTarget const target, amrex::ParticleReal const projectile_rest_energy,
                        amrex::ParticleReal const projectile_energy_min,
                        amrex::ParticleReal const projectile_energy_max)
        : m_projectile_energy_min(projectile_energy_min),
          m_projectile_energy_max(projectile_energy_max),
          m_projectile_rest_energy(projectile_rest_energy)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(projectile_rest_energy)) &&
                projectile_rest_energy > 0.0,
            "PJG projectile rest energy must be finite and positive.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(projectile_energy_min)) &&
                projectile_energy_min > 0.0,
            "PJG minimum projectile energy must be finite and positive.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(projectile_energy_max)) &&
                projectile_energy_max > projectile_energy_min,
            "PJG maximum projectile energy must exceed the minimum energy.");

        auto const log_energy_min = std::log(static_cast<double>(projectile_energy_min));
        auto const log_energy_max = std::log(static_cast<double>(projectile_energy_max));
        auto const log_energy_step =
            (log_energy_max - log_energy_min) / static_cast<double>(table_energy_points - 1);
        m_log_projectile_energy_min = static_cast<amrex::ParticleReal>(log_energy_min);
        m_inv_log_projectile_energy_step = static_cast<amrex::ParticleReal>(1.0 / log_energy_step);

        amrex::Vector<amrex::ParticleReal> host_cross_section(table_energy_points);
        amrex::Vector<amrex::ParticleReal> host_log_secondary_energy(table_energy_points *
                                                                     table_quantile_points);
        amrex::Vector<amrex::ParticleReal> host_binding_energy(table_energy_points *
                                                               table_quantile_points);

        for (int energy_index = 0; energy_index < table_energy_points; ++energy_index) {
            auto const projectile_energy =
                std::exp(log_energy_min + static_cast<double>(energy_index) * log_energy_step);
            auto const maximum_transfer = maximumEnergyTransfer(
                projectile_energy, static_cast<double>(projectile_rest_energy));
            auto const total_cross_section = integratedCrossSection(
                target, projectile_energy, static_cast<double>(projectile_rest_energy));
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(total_cross_section) &&
                                                 total_cross_section > 0.0,
                                             "The corrected PJG model produced a "
                                             "non-positive total cross section.");
            host_cross_section[energy_index] =
                static_cast<amrex::ParticleReal>(total_cross_section);

            double lower_energy = 0.0;
            for (int quantile_index = 0; quantile_index < table_quantile_points; ++quantile_index) {
                auto const offset = energy_index * table_quantile_points + quantile_index;
                auto const transformed_quantile = static_cast<double>(quantile_index) /
                                                  static_cast<double>(table_quantile_points - 1);
                auto const transformed_fourth_power = std::pow(transformed_quantile, 4);
                auto const complement_fourth_power = std::pow(1.0 - transformed_quantile, 4);
                auto const quantile =
                    transformed_fourth_power / (transformed_fourth_power + complement_fourth_power);
                double secondary_energy;
                if (quantile_index == 0) {
                    secondary_energy = 0.0;
                } else if (quantile_index == table_quantile_points - 1) {
                    secondary_energy = maximum_transfer;
                } else {
                    auto const target_cross_section = quantile * total_cross_section;
                    double upper_energy = maximum_transfer;
                    // A fixed iteration count makes table generation reproducible.
                    for (int iteration = 0; iteration < 48; ++iteration) {
                        auto const midpoint = 0.5 * (lower_energy + upper_energy);
                        auto const cumulative =
                            square_centimeter_to_square_meter *
                            cumulativeCrossSectionCm(target, projectile_energy, midpoint,
                                                     static_cast<double>(projectile_rest_energy));
                        if (cumulative < target_cross_section) {
                            lower_energy = midpoint;
                        } else {
                            upper_energy = midpoint;
                        }
                    }
                    secondary_energy = 0.5 * (lower_energy + upper_energy);
                    lower_energy = secondary_energy;
                }

                host_log_secondary_energy[offset] =
                    static_cast<amrex::ParticleReal>(std::log1p(secondary_energy));
                host_binding_energy[offset] = static_cast<amrex::ParticleReal>(
                    effectiveBindingEnergy(target, projectile_energy, secondary_energy,
                                           static_cast<double>(projectile_rest_energy)));
            }
        }

        m_cross_section.resize(host_cross_section.size());
        m_log_secondary_energy.resize(host_log_secondary_energy.size());
        m_binding_energy.resize(host_binding_energy.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_cross_section.begin(),
                         host_cross_section.end(), m_cross_section.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_log_secondary_energy.begin(),
                         host_log_secondary_energy.end(), m_log_secondary_energy.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, host_binding_energy.begin(),
                         host_binding_energy.end(), m_binding_energy.begin());
    }

    PJGModel::Executor
    PJGModel::executor () const noexcept
    {
        return {m_cross_section.dataPtr(),        m_log_secondary_energy.dataPtr(),
                m_binding_energy.dataPtr(),       m_log_projectile_energy_min,
                m_inv_log_projectile_energy_step, m_projectile_energy_min,
                m_projectile_energy_max,          m_projectile_rest_energy};
    }

    PJGTarget
    PJGModel::parseTarget (std::string const& target)
    {
        if (target == "N2" || target == "n2") {
            return PJGTarget::N2;
        }
        if (target == "O2" || target == "o2") {
            return PJGTarget::O2;
        }
        WARPX_ABORT_WITH_MESSAGE("PJG ionization target must be N2 or O2.");
        return PJGTarget::N2;
    }

    std::string
    PJGModel::targetName (PJGTarget const target)
    {
        return target == PJGTarget::N2 ? "N2" : "O2";
    }

    double
    PJGModel::differentialCrossSection (PJGTarget const target, double const projectile_energy,
                                        double const secondary_energy,
                                        double const projectile_rest_energy)
    {
        if (projectile_energy <= 0.0 || secondary_energy < 0.0) {
            return 0.0;
        }
        auto const& p = parameters(target);
        auto const state = projectileState(projectile_energy, projectile_rest_energy);
        if (secondary_energy > state.maximum_transfer) {
            return 0.0;
        }
        auto const functions = energyFunctions(p, state);
        double result = 0.0;
        for (int continuum = 0; continuum < p.num_continua; ++continuum) {
            result += continuumDifferentialCrossSectionCm(p, continuum, state, functions,
                                                          secondary_energy);
        }
        return square_centimeter_to_square_meter * std::max(result, 0.0);
    }

    double
    PJGModel::integratedCrossSection (PJGTarget const target, double const projectile_energy,
                                      double const projectile_rest_energy)
    {
        if (projectile_energy <= 0.0) {
            return 0.0;
        }
        auto const maximum_transfer =
            maximumEnergyTransfer(projectile_energy, projectile_rest_energy);
        return square_centimeter_to_square_meter *
               cumulativeCrossSectionCm(target, projectile_energy, maximum_transfer,
                                        projectile_rest_energy);
    }

    double
    PJGModel::maximumEnergyTransfer (double const projectile_energy,
                                     double const projectile_rest_energy)
    {
        if (projectile_energy <= 0.0 || projectile_rest_energy <= 0.0) {
            return 0.0;
        }
        return projectileState(projectile_energy, projectile_rest_energy).maximum_transfer;
    }
} // namespace ProtonImpactIonization
