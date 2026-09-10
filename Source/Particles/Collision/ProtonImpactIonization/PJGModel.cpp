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
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace ProtonImpactIonization
{
    namespace
    {
        constexpr int max_continua = 7;
        constexpr int intervals_per_segment = 1024;
        constexpr double electron_rest_energy = 510998.95069;
        constexpr double proton_rest_energy = 938272088.16;
        constexpr double atomic_rest_energy = 931494103.72;
        constexpr double rydberg = 13.605693122994;
        constexpr double bohr_radius = 5.29177210544e-11;
        constexpr double bethe_constant =
            4.0 * MathConst::pi * bohr_radius * bohr_radius * rydberg * rydberg;

        bool
        inCalibrationRange (double const energy, double const mass)
        {
            if (!std::isfinite(energy) || !std::isfinite(mass) ||
                mass < proton_rest_energy / 2.0) {
                return false;
            }
            auto const equivalent_proton_energy = energy / mass * proton_rest_energy;
            return equivalent_proton_energy >= 5.0e3 * (1.0 - 1.e-6) &&
                   equivalent_proton_energy <= 1.0e10 * (1.0 + 1.e-6);
        }

        struct Parameters
        {
            int m_num_continua;
            int m_num_electrons;
            double m_neutral_mass_number;
            double m_k;
            double m_width;
            double m_broad_excess;
            double m_t_s;
            double m_t_numerator;
            double m_t_denominator;
            double m_j;
            double m_power;
            double m_edge;
            std::array<double, max_continua> m_thresholds;
            std::array<double, max_continua> m_fractions;
        };

        // Frozen joint calibration. K is converted from cm^2 to m^2; the
        // peak numerator has units eV^2. No empirical hard normalization.
        constexpr Parameters n2_parameters{6,
                                           14,
                                           28.0134,
                                           6.679516247830822e-20,
                                           12.503227296612915,
                                           88.54007881169514,
                                           4.0,
                                           3635.677428147113,
                                           1970.0,
                                           13.8235317299573,
                                           0.807,
                                           0.70 * 1.1911400606628986,
                                           {15.58, 16.73, 18.75, 22.0, 23.6, 40.0, 0.0},
                                           {0.456, 0.2, 0.104, 0.07, 0.07, 0.1, 0.0}};
        constexpr Parameters o2_parameters{7,
                                           16,
                                           31.9988,
                                           5.119050696037388e-20,
                                           16.757498368297895,
                                           159.2697988004284,
                                           6.34,
                                           943.7164329964345,
                                           128.0,
                                           8.76635921598519,
                                           1.314,
                                           0.59,
                                           {12.1, 16.1, 16.9, 18.2, 20.3, 23.0, 37.0},
                                           {0.08, 0.19, 0.19, 0.17, 0.11, 0.16, 0.1}};

        Parameters const&
        parameters (PJGTarget const target)
        {
            return target == PJGTarget::N2 ? n2_parameters : o2_parameters;
        }

        double
        logistic (double const value)
        {
            auto const exponential = std::exp(-std::abs(value));
            return value >= 0.0 ? 1.0 / (1.0 + exponential) : exponential / (1.0 + exponential);
        }

        /** Cache projectile-dependent factors once per host integration. */
        struct Spectrum
        {
            Parameters const& m_p;
            double m_energy;
            double m_mass;
            double m_neutral;
            double m_equivalent;
            double m_maximum;
            double m_endpoint;
            double m_peak;
            double m_width_squared;
            double m_broad_squared;
            double m_distortion;
            double m_soft_scale;
            double m_terminal_factor;
            double m_terminal_log_slope;
            std::array<double, max_continua> m_edge_center{};
            std::array<double, max_continua> m_edge_inverse_width{};
            std::array<double, max_continua> m_channel_endpoint{};
            double m_largest_width = 0.0;

            Spectrum (PJGTarget const target, double const e, double const m)
                : m_p(parameters(target)), m_energy(e), m_mass(m),
                  m_neutral(m_p.m_neutral_mass_number * atomic_rest_energy)
            {
                auto const gamma = 1.0 + m_energy / m_mass;
                auto const beta_squared = relativisticBetaSquared(m_energy, m_mass);
                m_equivalent = 0.5 * electron_rest_energy * beta_squared;
                m_maximum = relativisticMaximumTransfer(m_energy, m_mass);
                m_endpoint = molecularMaximumSecondaryEnergy(m_energy, m_mass, m_neutral,
                                                             m_p.m_thresholds[0]);
                m_peak = m_p.m_t_s - m_p.m_t_numerator / (m_equivalent + m_p.m_t_denominator);
                m_width_squared = m_p.m_width * m_p.m_width;
                m_broad_squared = m_p.m_broad_excess * m_p.m_broad_excess;
                m_distortion = logistic(
                    m_p.m_power * std::log(m_energy * electron_rest_energy / (m_mass * m_p.m_j)));
                double mean_binding = 0.0;
                auto const ratio = electron_rest_energy / m_mass;
                auto const denominator = 1.0 + 2.0 * gamma * ratio + ratio * ratio;
                auto const broadening = gamma * (gamma + ratio) * std::pow(1.0 + ratio, 3) /
                                        (denominator * denominator);
                for (int j = 0; j < m_p.m_num_continua; ++j) {
                    mean_binding += m_p.m_fractions[j] * m_p.m_thresholds[j];
                    auto const width = broadening * std::sqrt(m_equivalent * m_p.m_thresholds[j]);
                    m_largest_width = std::max(m_largest_width, width);
                    m_edge_center[j] = m_maximum - 2.0 * width - rydberg / 4.0;
                    m_edge_inverse_width[j] = m_p.m_edge / width;
                    m_channel_endpoint[j] = molecularMaximumSecondaryEnergy(
                        m_energy, m_mass, m_neutral, m_p.m_thresholds[j]);
                }
                auto const logarithm =
                    std::log(4.0 * m_equivalent * gamma * gamma / mean_binding + std::exp(1.0)) -
                    beta_squared;
                m_soft_scale = m_distortion * m_p.m_k * m_width_squared * logarithm;

                m_terminal_factor = bhabhaSpinHalfFactor(m_energy, m_mass, m_maximum);
                auto const total = m_energy + m_mass;
                auto const momentum = std::sqrt(m_energy * (m_energy + 2.0 * m_mass));
                auto const invariant =
                    m_mass * m_mass + electron_rest_energy * (electron_rest_energy + 2.0 * total);
                // Factor 1-(Tmax/m_p)^2, keeping the terminal derivative
                // nonpositive without cancellation at relativistic energies.
                auto const minus =
                    m_mass * m_mass + electron_rest_energy * electron_rest_energy +
                    2.0 * electron_rest_energy * m_mass * m_mass / (total + momentum);
                m_terminal_log_slope = -beta_squared / m_maximum * (minus / invariant) *
                                       (1.0 + 2.0 * electron_rest_energy * momentum / invariant) /
                                       m_terminal_factor;
            }

            double
            phaseFraction (double const t, int const j) const
            {
                if (t >= m_channel_endpoint[j]) {
                    return 0.0;
                }
                auto const binding = m_p.m_thresholds[j];
                auto const threshold =
                    binding * (1.0 + m_mass / m_neutral) + binding * binding / (2.0 * m_neutral);
                auto const constant =
                    (m_neutral - electron_rest_energy) * (m_energy - threshold) -
                    electron_rest_energy * binding * (m_mass + binding / 2.0) / m_neutral;
                auto const numerator = (m_energy + m_mass + m_neutral) * t - constant;
                auto const denominator = std::sqrt(m_energy * (m_energy + 2.0 * m_mass) * t *
                                                   (t + 2.0 * electron_rest_energy));
                return denominator > 0.0
                           ? std::clamp((denominator - numerator) / (2.0 * denominator), 0.0, 1.0)
                           : (numerator <= 0.0 ? 1.0 : 0.0);
            }

            /** Return the SDCS and its binding-weighted value, in SI units. */
            std::array<double, 2>
            evaluate (double const t) const
            {
                if (!(t >= 0.0 && t < m_endpoint)) {
                    return {0.0, 0.0};
                }
                double gate = 0.0;
                double binding_gate = 0.0;
                for (int j = 0; j < m_p.m_num_continua; ++j) {
                    auto const weight = m_p.m_fractions[j] * phaseFraction(t, j) *
                                        logistic((m_edge_center[j] - t) * m_edge_inverse_width[j]);
                    gate += weight;
                    binding_gate += weight * m_p.m_thresholds[j];
                }
                auto const narrow = (t - m_peak) * (t - m_peak) + m_width_squared;
                auto const broad = 1.0 / (narrow + m_broad_squared);
                auto const difference = m_broad_squared * broad / narrow;
                auto const hard_distortion =
                    1.0 - (1.0 - m_distortion) / (1.0 + t * t / m_broad_squared);
                auto const factor =
                    t <= m_maximum
                        ? bhabhaSpinHalfFactor(m_energy, m_mass, t)
                        : m_terminal_factor * std::exp(m_terminal_log_slope * (t - m_maximum));
                auto const value =
                    (m_soft_scale * difference +
                     hard_distortion * m_p.m_num_electrons * bethe_constant * broad * factor) /
                    m_equivalent;
                return {gate * value, binding_gate * value};
            }
        };

        /** Positive quadrature with both CDF and survival integrals.
         *
         * The four segments resolve the peak, either side of free Tmax,
         * and the molecular tail. Independent reverse accumulation avoids
         * subtracting a nearly unit CDF to resolve rare hard electrons.
         */
        struct IntegratedSpectrum
        {
            Spectrum const& m_spectrum;
            std::vector<double> m_x;
            std::vector<double> m_density;
            std::vector<double> m_area;
            std::vector<double> m_cumulative;
            std::vector<double> m_survival;
            PJGModel::Moments m_moments{};

            explicit IntegratedSpectrum (Spectrum const& state) : m_spectrum(state)
            {
                auto const left =
                    std::max(state.m_maximum / 2.0, state.m_maximum - 32.0 * state.m_largest_width);
                auto const right =
                    std::min(state.m_endpoint, state.m_maximum + 32.0 * state.m_largest_width);
                std::array<double, 5> const edges{0.0, left, state.m_maximum, right,
                                                  state.m_endpoint};
                m_x.reserve(4 * intervals_per_segment + 1);
                m_x.push_back(0.0);
                for (int segment = 0; segment < 4; ++segment) {
                    auto const lower = edges[segment];
                    auto const upper = edges[segment + 1];
                    if (upper <= lower) {
                        continue;
                    }
                    for (int i = 1; i <= intervals_per_segment; ++i) {
                        auto const fraction = static_cast<double>(i) / intervals_per_segment;
                        auto const coordinate =
                            segment == 0 || segment == 3
                                ? std::log1p(lower) +
                                      fraction * (std::log1p(upper) - std::log1p(lower))
                                : std::log1p(lower + fraction * (upper - lower));
                        m_x.push_back(coordinate);
                    }
                }
                m_density.reserve(m_x.size());
                for (auto const coordinate : m_x) {
                    m_density.push_back(state.evaluate(std::expm1(coordinate))[0] *
                                        std::exp(coordinate));
                }
                m_area.resize(m_x.size() - 1);
                m_cumulative.resize(m_x.size(), 0.0);
                m_survival.resize(m_x.size(), 0.0);
                constexpr std::array<double, 4> nodes{-0.8611363115940526, -0.3399810435848563,
                                                      0.3399810435848563, 0.8611363115940526};
                constexpr std::array<double, 4> weights{0.3478548451374538, 0.6521451548625461,
                                                        0.6521451548625461, 0.3478548451374538};
                for (std::size_t i = 0; i < m_area.size(); ++i) {
                    auto const half = (m_x[i + 1] - m_x[i]) / 2.0;
                    for (int node = 0; node < 4; ++node) {
                        auto const coordinate = m_x[i] + half * (nodes[node] + 1.0);
                        auto const t = std::expm1(coordinate);
                        auto const values = state.evaluate(t);
                        auto const jacobian = half * weights[node] * std::exp(coordinate);
                        auto const contribution = jacobian * values[0];
                        m_area[i] += contribution;
                        m_moments.m_kinetic += contribution * t;
                        m_moments.m_kinetic_second += contribution * t * t;
                        m_moments.m_binding += jacobian * values[1];
                        if (t > state.m_maximum) {
                            m_moments.m_above_free += contribution;
                        }
                    }
                    m_cumulative[i + 1] = m_cumulative[i] + m_area[i];
                }
                for (std::size_t i = m_area.size(); i-- > 0;) {
                    m_survival[i] = m_survival[i + 1] + m_area[i];
                }
                m_moments.m_total = m_cumulative.back();
            }

            double
            quantile (double const coordinate) const
            {
                if (coordinate <= 0.0) {
                    return 0.0;
                }
                if (coordinate >= 1.0) {
                    return m_spectrum.m_endpoint;
                }
                auto const forward = coordinate <= 0.5;
                auto const a = std::pow(coordinate, 4);
                auto const b = std::pow(1.0 - coordinate, 4);
                auto const probability = (forward ? a : b) / (a + b);
                auto const target = probability * m_moments.m_total;
                auto const count = static_cast<std::ptrdiff_t>(m_area.size());
                auto const index =
                    forward
                        ? std::upper_bound(m_cumulative.begin(), m_cumulative.end(), target) -
                              m_cumulative.begin() - 1
                        : count -
                              (std::upper_bound(m_survival.rbegin(), m_survival.rend(), target) -
                               m_survival.rbegin());
                auto const i =
                    static_cast<std::size_t>(std::clamp(index, std::ptrdiff_t{0}, count - 1));
                if (m_area[i] <= 0.0) {
                    return std::expm1(m_x[i]);
                }
                auto const fraction =
                    std::clamp(forward ? (target - m_cumulative[i]) / m_area[i]
                                       : 1.0 - (target - m_survival[i + 1]) / m_area[i],
                               0.0, 1.0);
                auto const width = m_x[i + 1] - m_x[i];
                auto slope_a = m_density[i] * width / m_area[i];
                auto slope_b = m_density[i + 1] * width / m_area[i];
                auto const norm = std::hypot(slope_a, slope_b);
                if (norm > 3.0) {
                    // Monotone Hermite CDF: never create a negative m_density
                    // while interpolating an exponentially small tail cell.
                    slope_a *= 3.0 / norm;
                    slope_b *= 3.0 / norm;
                }
                double lower = 0.0;
                double upper = 1.0;
                for (int iteration = 0; iteration < 36; ++iteration) {
                    auto const s = (lower + upper) / 2.0;
                    auto const value = s * (slope_a + s * ((3.0 - 2.0 * slope_a - slope_b) +
                                                           s * (slope_a + slope_b - 2.0)));
                    if (value < fraction) {
                        lower = s;
                    } else {
                        upper = s;
                    }
                }
                return std::expm1(m_x[i] + (lower + upper) / 2.0 * width);
            }
        };
    } // namespace

    PJGModel::PJGModel (PJGTarget const target, amrex::ParticleReal const projectile_rest_energy,
                        amrex::ParticleReal const projectile_energy_min,
                        amrex::ParticleReal const projectile_energy_max)
        : m_projectile_energy_min(projectile_energy_min),
          m_projectile_energy_max(projectile_energy_max),
          m_projectile_rest_energy(projectile_rest_energy)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(projectile_rest_energy) &&
                                             projectile_rest_energy >= proton_rest_energy / 2,
                                         "PJG requires a finite heavy-projectile rest energy.");
        auto const mass_scale = static_cast<double>(projectile_rest_energy) / proton_rest_energy;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(projectile_energy_min) && std::isfinite(projectile_energy_max) &&
                projectile_energy_min >= 5.0e3 * mass_scale * (1.0 - 1.e-6) &&
                projectile_energy_max <= 1.0e10 * mass_scale * (1.0 + 1.e-6) &&
                projectile_energy_max > projectile_energy_min,
            "PJG table bounds require 5 keV <= (m_p/M) E <= 10 GeV and E_max > E_min.");
        auto const& p = parameters(target);
        m_neutral_rest_energy =
            static_cast<amrex::ParticleReal>(p.m_neutral_mass_number * atomic_rest_energy);
        m_minimum_binding_energy = static_cast<amrex::ParticleReal>(p.m_thresholds[0]);

        auto const log_min = std::log(static_cast<double>(projectile_energy_min));
        auto const step = std::log(static_cast<double>(projectile_energy_max) /
                                   static_cast<double>(projectile_energy_min)) /
                          (table_energy_points - 1);
        m_log_projectile_energy_min = static_cast<amrex::ParticleReal>(log_min);
        m_inv_log_projectile_energy_step = static_cast<amrex::ParticleReal>(1.0 / step);
        amrex::Vector<amrex::ParticleReal> cross_section(table_energy_points);
        amrex::Vector<amrex::ParticleReal> log_secondary(table_energy_points *
                                                         table_quantile_points);
        amrex::Vector<amrex::ParticleReal> binding(table_energy_points * table_quantile_points);
        for (int i = 0; i < table_energy_points; ++i) {
            auto const energy = std::exp(log_min + i * step);
            Spectrum const spectrum(target, energy, projectile_rest_energy);
            IntegratedSpectrum const integrated(spectrum);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(integrated.m_moments.m_total) && integrated.m_moments.m_total > 0.0,
                "The calibrated PJG integral must be finite and positive.");
            cross_section[i] = static_cast<amrex::ParticleReal>(integrated.m_moments.m_total);
            for (int j = 0; j < table_quantile_points; ++j) {
                auto const t =
                    integrated.quantile(static_cast<double>(j) / (table_quantile_points - 1));
                auto const offset = i * table_quantile_points + j;
                log_secondary[offset] = static_cast<amrex::ParticleReal>(std::log1p(t));
                auto const values = spectrum.evaluate(t);
                binding[offset] = static_cast<amrex::ParticleReal>(
                    values[0] > 0.0 ? values[1] / values[0] : p.m_thresholds[0]);
            }
        }
        m_cross_section.resize(cross_section.size());
        m_log_secondary_energy.resize(log_secondary.size());
        m_binding_energy.resize(binding.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, cross_section.begin(), cross_section.end(),
                         m_cross_section.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, log_secondary.begin(), log_secondary.end(),
                         m_log_secondary_energy.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, binding.begin(), binding.end(),
                         m_binding_energy.begin());
    }

    PJGModel::Executor
    PJGModel::executor () const noexcept
    {
        return {m_cross_section.dataPtr(),
                m_log_secondary_energy.dataPtr(),
                m_binding_energy.dataPtr(),
                m_log_projectile_energy_min,
                m_inv_log_projectile_energy_step,
                m_projectile_energy_min,
                m_projectile_energy_max,
                m_projectile_rest_energy,
                m_neutral_rest_energy,
                m_minimum_binding_energy};
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
        if (!inCalibrationRange(projectile_energy, projectile_rest_energy)) {
            return 0.0;
        }
        return Spectrum(target, projectile_energy, projectile_rest_energy)
            .evaluate(secondary_energy)[0];
    }

    PJGModel::Moments
    PJGModel::integratedMoments (PJGTarget const target, double const projectile_energy,
                                 double const projectile_rest_energy)
    {
        if (!inCalibrationRange(projectile_energy, projectile_rest_energy)) {
            return {};
        }
        Spectrum const spectrum(target, projectile_energy, projectile_rest_energy);
        if (spectrum.m_endpoint <= spectrum.m_maximum) {
            return {};
        }
        return IntegratedSpectrum(spectrum).m_moments;
    }

    double
    PJGModel::integratedCrossSection (PJGTarget const target, double const projectile_energy,
                                      double const projectile_rest_energy)
    {
        return integratedMoments(target, projectile_energy, projectile_rest_energy).m_total;
    }

    double
    PJGModel::maximumEnergyTransfer (double const projectile_energy,
                                     double const projectile_rest_energy)
    {
        return projectile_energy > 0.0 && projectile_rest_energy > 0.0 &&
                       std::isfinite(projectile_energy) && std::isfinite(projectile_rest_energy)
                   ? relativisticMaximumTransfer(projectile_energy, projectile_rest_energy)
                   : 0.0;
    }
} // namespace ProtonImpactIonization
