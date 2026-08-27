/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_KINEMATICS_H_

#include "BackgroundMCCKinematics.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace BackgroundMCCElasticKinematics
{
    /** Apply an elastic electron-neutral collision in the neutral rest frame.
     *
     * The sampled polar angle is the electron angle represented by the
     * laboratory-frame IAA/elmolcs DCS for a stationary target. Relativistic
     * two-body kinematics set the outgoing electron energy, and the neutral
     * receives the exact recoil momentum. Both products are then boosted back
     * to the simulation frame.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute (amrex::ParticleReal const incident_ux, amrex::ParticleReal const incident_uy,
             amrex::ParticleReal const incident_uz, amrex::ParticleReal const neutral_vx,
             amrex::ParticleReal const neutral_vy, amrex::ParticleReal const neutral_vz,
             double const electron_mass, double const neutral_mass,
             amrex::ParticleReal const scattering_cosine, amrex::RandomEngine const& engine,
             amrex::ParticleReal& electron_ux, amrex::ParticleReal& electron_uy,
             amrex::ParticleReal& electron_uz, amrex::ParticleReal& recoil_ux,
             amrex::ParticleReal& recoil_uy, amrex::ParticleReal& recoil_uz) noexcept
    {
        using BackgroundMCCKinematics::Vector3;
        using BackgroundMCCKinematics::boostToLab;
        using BackgroundMCCKinematics::boostToRest;
        using BackgroundMCCKinematics::directionFromPolarAngle;
        using BackgroundMCCKinematics::dot;
        using BackgroundMCCKinematics::kineticEnergyFromProperVelocity;
        using BackgroundMCCKinematics::properSpeedFromKineticEnergy;
        using BackgroundMCCKinematics::scale;
        using BackgroundMCCKinematics::subtract;
        using BackgroundMCCKinematics::transverseDirections;

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
        auto const boost_spatial_factor =
            neutral_gamma * neutral_gamma / ((neutral_gamma + 1.0) * c2);
        auto const incident_rest = boostToRest(
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

        auto const cosine = amrex::max(
            -1.0, amrex::min(static_cast<double>(scattering_cosine), 1.0));
        auto const sine_sq = amrex::max(1.0 - cosine * cosine, 0.0);
        auto const incident_energy =
            kineticEnergyFromProperVelocity(incident_rest, electron_mass);
        auto const electron_rest_energy = electron_mass * c2 / elementary_charge;
        auto const neutral_rest_energy = neutral_mass * c2 / elementary_charge;

        // This stable form avoids subtracting two nearly equal neutral rest
        // energies for strongly forward scattering from a heavy molecule.
        auto const recoil_root = std::sqrt(amrex::max(
            neutral_rest_energy * neutral_rest_energy -
                electron_rest_energy * electron_rest_energy * sine_sq,
            0.0));
        auto const recoil_factor =
            (incident_energy + electron_rest_energy) * sine_sq +
            neutral_rest_energy * (1.0 - cosine) +
            cosine * electron_rest_energy * electron_rest_energy * sine_sq /
                (neutral_rest_energy + recoil_root);
        auto const incident_pc_sq =
            incident_energy * (incident_energy + 2.0 * electron_rest_energy);
        auto const total_energy =
            incident_energy + electron_rest_energy + neutral_rest_energy;
        auto const denominator =
            total_energy * total_energy - incident_pc_sq * cosine * cosine;
        auto const energy_loss = denominator > 0.0
            ? recoil_factor * incident_pc_sq / denominator
            : incident_energy;
        auto const outgoing_energy =
            amrex::max(incident_energy - energy_loss, 0.0);

        auto const azimuth = two_pi * static_cast<double>(amrex::Random(engine));
        auto const outgoing_direction = directionFromPolarAngle(
            incident_direction, transverse_1, transverse_2, cosine, azimuth);
        auto const electron_rest = scale(
            outgoing_direction,
            properSpeedFromKineticEnergy(outgoing_energy, electron_mass));
        auto const recoil_rest = scale(
            subtract(incident_rest, electron_rest), electron_mass / neutral_mass);

        auto const electron_lab =
            boostToLab(electron_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        auto const recoil_lab =
            boostToLab(recoil_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        electron_ux = static_cast<amrex::ParticleReal>(electron_lab.x);
        electron_uy = static_cast<amrex::ParticleReal>(electron_lab.y);
        electron_uz = static_cast<amrex::ParticleReal>(electron_lab.z);
        recoil_ux = static_cast<amrex::ParticleReal>(recoil_lab.x);
        recoil_uy = static_cast<amrex::ParticleReal>(recoil_lab.y);
        recoil_uz = static_cast<amrex::ParticleReal>(recoil_lab.z);
    }
} // namespace BackgroundMCCElasticKinematics

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_KINEMATICS_H_
