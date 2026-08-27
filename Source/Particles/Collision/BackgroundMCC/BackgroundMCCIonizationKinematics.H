/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_IONIZATION_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_IONIZATION_KINEMATICS_H_

#include "BackgroundMCCKinematics.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Math.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace BackgroundMCCIonizationKinematics
{
    using BackgroundMCCKinematics::Vector3;
    using BackgroundMCCKinematics::add;
    using BackgroundMCCKinematics::boostToLab;
    using BackgroundMCCKinematics::directionFromPolarAngle;
    using BackgroundMCCKinematics::dot;
    using BackgroundMCCKinematics::kineticEnergyFromProperVelocity;
    using BackgroundMCCKinematics::properSpeedFromKineticEnergy;
    using BackgroundMCCKinematics::scale;
    using BackgroundMCCKinematics::subtract;
    using BackgroundMCCKinematics::transverseDirections;

    /** Compute electron-ion impact-ionization kinematics in the neutral rest frame.
     *
     * For the IAA model, the primary polar angle follows Schmalzried Eq. (2.60),
     * and the secondary uses the bound/free interpolation of Eqs. (2.66) and
     * (11.132). The two electron azimuths differ by pi, which gives the natural
     * coplanar binary-encounter geometry while the ion absorbs the residual
     * momentum. Two fixed recoil corrections conserve kinetic plus binding
     * energy to far below particle precision for molecular targets.
     *
     * @param incident_ux,incident_uy,incident_uz incident lab-frame proper velocity
     * @param neutral_vx,neutral_vy,neutral_vz sampled neutral lab-frame velocity
     * @param electron_mass electron rest mass in kg
     * @param ion_mass product-ion rest mass in kg
     * @param binding_energy binding-energy loss in eV
     * @param secondary_fraction fraction of available electron energy given to
     *        the lower-energy electron, in [0, 1/2]
     * @param angle_model selected scattering-angle model
     * @param engine random-number engine
     * @param primary_ux,primary_uy,primary_uz outgoing primary proper velocity
     * @param secondary_ux,secondary_uy,secondary_uz outgoing secondary proper velocity
     * @param ion_ux,ion_uy,ion_uz outgoing ion proper velocity
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute (amrex::ParticleReal const incident_ux, amrex::ParticleReal const incident_uy,
             amrex::ParticleReal const incident_uz, amrex::ParticleReal const neutral_vx,
             amrex::ParticleReal const neutral_vy, amrex::ParticleReal const neutral_vz,
             double const electron_mass, double const ion_mass,
             amrex::ParticleReal const binding_energy, amrex::ParticleReal const secondary_fraction,
             ScatteringAngleModel const angle_model, amrex::RandomEngine const& engine,
             amrex::ParticleReal& primary_ux, amrex::ParticleReal& primary_uy,
             amrex::ParticleReal& primary_uz, amrex::ParticleReal& secondary_ux,
             amrex::ParticleReal& secondary_uy, amrex::ParticleReal& secondary_uz,
             amrex::ParticleReal& ion_ux, amrex::ParticleReal& ion_uy,
             amrex::ParticleReal& ion_uz) noexcept
    {
        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double elementary_charge = PhysConst::q_e_v<double>;
        constexpr double two_pi = 2.0 * static_cast<double>(MathConst::pi);

        Vector3 const neutral_velocity{static_cast<double>(neutral_vx),
                                       static_cast<double>(neutral_vy),
                                       static_cast<double>(neutral_vz)};
        Vector3 const incident_lab{static_cast<double>(incident_ux),
                                   static_cast<double>(incident_uy),
                                   static_cast<double>(incident_uz)};

        auto const neutral_speed_sq = dot(neutral_velocity, neutral_velocity);
        auto const neutral_gamma =
            1.0 / std::sqrt(amrex::max(1.0 - neutral_speed_sq * inv_c2,
                                       std::numeric_limits<double>::epsilon()));
        // Stable form of (gamma_n - 1)/V_n^2, including V_n = 0.
        auto const boost_spatial_factor =
            neutral_gamma * neutral_gamma / ((neutral_gamma + 1.0) * c2);

        auto const incident_rest = BackgroundMCCKinematics::boostToRest(
            incident_lab, neutral_velocity, neutral_gamma, boost_spatial_factor);

        auto const incident_speed = std::sqrt(dot(incident_rest, incident_rest));
        Vector3 incident_direction{0.0, 0.0, 1.0};
        if (incident_speed > std::numeric_limits<double>::min())
        {
            incident_direction = scale(incident_rest, 1.0 / incident_speed);
        }

        Vector3 transverse_1;
        Vector3 transverse_2;
        transverseDirections(incident_direction, transverse_1, transverse_2);

        auto const incident_energy = kineticEnergyFromProperVelocity(incident_rest, electron_mass);
        auto const binding = static_cast<double>(binding_energy);
        auto const available_energy = amrex::max(incident_energy - binding, 0.0);
        auto const sharing =
            amrex::max(0.0, amrex::min(static_cast<double>(secondary_fraction), 0.5));
        auto const unrecoiled_secondary_energy = sharing * available_energy;
        auto const unrecoiled_primary_energy = available_energy - unrecoiled_secondary_energy;

        Vector3 primary_direction;
        Vector3 secondary_direction;
        if (angle_model == ScatteringAngleModel::IAA && available_energy > 0.0)
        {
            auto const electron_rest_energy = electron_mass * c2 / elementary_charge;
            auto const primary_cosine = std::sqrt(
                amrex::max(0.0, amrex::min(unrecoiled_primary_energy *
                                               (available_energy + 2.0 * electron_rest_energy) /
                                               (available_energy * (unrecoiled_primary_energy +
                                                                    2.0 * electron_rest_energy)),
                                           1.0)));
            auto const secondary_kinematic_cosine = std::sqrt(
                amrex::max(0.0, amrex::min(unrecoiled_secondary_energy *
                                               (available_energy + 2.0 * electron_rest_energy) /
                                               (available_energy * (unrecoiled_secondary_energy +
                                                                    2.0 * electron_rest_energy)),
                                           1.0)));
            auto const denominator = unrecoiled_secondary_energy + binding;
            auto const free_weight =
                denominator > 0.0 ? unrecoiled_secondary_energy / denominator : 0.0;
            auto const bound_cosine = 1.0 - 2.0 * static_cast<double>(amrex::Random(engine));
            auto const secondary_cosine =
                amrex::max(-1.0, amrex::min(free_weight * secondary_kinematic_cosine +
                                                (1.0 - free_weight) * bound_cosine,
                                            1.0));
            auto const azimuth = two_pi * static_cast<double>(amrex::Random(engine));

            primary_direction = directionFromPolarAngle(incident_direction, transverse_1,
                                                        transverse_2, primary_cosine, azimuth);
            secondary_direction = directionFromPolarAngle(
                incident_direction, transverse_1, transverse_2, secondary_cosine,
                azimuth + static_cast<double>(MathConst::pi));
        }
        else if (angle_model == ScatteringAngleModel::Forward)
        {
            primary_direction = incident_direction;
            secondary_direction = incident_direction;
        }
        else if (angle_model == ScatteringAngleModel::Backward)
        {
            primary_direction = scale(incident_direction, -1.0);
            secondary_direction = scale(incident_direction, -1.0);
        }
        else
        {
            auto const primary_cosine = 1.0 - 2.0 * static_cast<double>(amrex::Random(engine));
            auto const secondary_cosine = 1.0 - 2.0 * static_cast<double>(amrex::Random(engine));
            primary_direction = directionFromPolarAngle(
                incident_direction, transverse_1, transverse_2, primary_cosine,
                two_pi * static_cast<double>(amrex::Random(engine)));
            secondary_direction = directionFromPolarAngle(
                incident_direction, transverse_1, transverse_2, secondary_cosine,
                two_pi * static_cast<double>(amrex::Random(engine)));
        }

        Vector3 primary_rest{};
        Vector3 secondary_rest{};
        Vector3 ion_rest{};
        auto electron_energy = available_energy;

        // Ion recoil is a small O(m_e/M) correction. Two fixed iterations leave
        // an O((m_e/M)^3) residual and avoid a divergent device-side loop.
        for (int iteration = 0; iteration < 2; ++iteration)
        {
            auto const secondary_energy = sharing * electron_energy;
            auto const primary_energy = electron_energy - secondary_energy;
            primary_rest = scale(primary_direction,
                                 properSpeedFromKineticEnergy(primary_energy, electron_mass));
            secondary_rest = scale(secondary_direction,
                                   properSpeedFromKineticEnergy(secondary_energy, electron_mass));
            ion_rest = scale(subtract(subtract(incident_rest, primary_rest), secondary_rest),
                             electron_mass / ion_mass);
            auto const ion_energy = kineticEnergyFromProperVelocity(ion_rest, ion_mass);
            electron_energy = amrex::max(available_energy - ion_energy, 0.0);
        }

        auto const final_secondary_energy = sharing * electron_energy;
        auto const final_primary_energy = electron_energy - final_secondary_energy;
        primary_rest = scale(primary_direction,
                             properSpeedFromKineticEnergy(final_primary_energy, electron_mass));
        secondary_rest = scale(secondary_direction,
                               properSpeedFromKineticEnergy(final_secondary_energy, electron_mass));
        ion_rest = scale(subtract(subtract(incident_rest, primary_rest), secondary_rest),
                         electron_mass / ion_mass);

        auto const primary_lab =
            boostToLab(primary_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        auto const secondary_lab =
            boostToLab(secondary_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        auto const ion_lab =
            boostToLab(ion_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);

        primary_ux = static_cast<amrex::ParticleReal>(primary_lab.x);
        primary_uy = static_cast<amrex::ParticleReal>(primary_lab.y);
        primary_uz = static_cast<amrex::ParticleReal>(primary_lab.z);
        secondary_ux = static_cast<amrex::ParticleReal>(secondary_lab.x);
        secondary_uy = static_cast<amrex::ParticleReal>(secondary_lab.y);
        secondary_uz = static_cast<amrex::ParticleReal>(secondary_lab.z);
        ion_ux = static_cast<amrex::ParticleReal>(ion_lab.x);
        ion_uy = static_cast<amrex::ParticleReal>(ion_lab.y);
        ion_uz = static_cast<amrex::ParticleReal>(ion_lab.z);
    }
} // namespace BackgroundMCCIonizationKinematics

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_IONIZATION_KINEMATICS_H_
