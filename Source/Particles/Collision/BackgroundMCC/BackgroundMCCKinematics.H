/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_KINEMATICS_H_

#include "Utils/WarpXConst.H"

#include <AMReX_Algorithm.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace BackgroundMCCKinematics
{
    struct Vector3
    {
        double x;
        double y;
        double z;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double dot (Vector3 const& a,
                                                         Vector3 const& b) noexcept
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vector3 scale (Vector3 const& vector,
                                                            double const factor) noexcept
    {
        return {factor * vector.x, factor * vector.y, factor * vector.z};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vector3 add (Vector3 const& a,
                                                          Vector3 const& b) noexcept
    {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vector3 subtract (Vector3 const& a,
                                                               Vector3 const& b) noexcept
    {
        return {a.x - b.x, a.y - b.y, a.z - b.z};
    }

    /** Boost a proper velocity from a moving frame to the laboratory frame. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vector3
    boostToLab (Vector3 const& rest_frame_velocity, Vector3 const& frame_velocity,
                double const frame_gamma, double const spatial_factor) noexcept
    {
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;

        auto const rest_gamma =
            std::sqrt(1.0 + dot(rest_frame_velocity, rest_frame_velocity) * inv_c2);
        auto const factor =
            spatial_factor * dot(frame_velocity, rest_frame_velocity) + frame_gamma * rest_gamma;
        return add(rest_frame_velocity, scale(frame_velocity, factor));
    }

    /** Boost a proper velocity from the laboratory to a moving frame. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vector3
    boostToRest (Vector3 const& lab_velocity, Vector3 const& frame_velocity,
                 double const frame_gamma, double const spatial_factor) noexcept
    {
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;

        auto const lab_gamma = std::sqrt(1.0 + dot(lab_velocity, lab_velocity) * inv_c2);
        auto const factor =
            spatial_factor * dot(frame_velocity, lab_velocity) - frame_gamma * lab_gamma;
        return add(lab_velocity, scale(frame_velocity, factor));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    properSpeedFromKineticEnergy (double const energy_ev, double const mass) noexcept
    {
        constexpr double c = PhysConst::c_v<double>;
        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double elementary_charge = PhysConst::q_e_v<double>;

        auto const normalized_energy = amrex::max(energy_ev * elementary_charge / (mass * c2), 0.0);
        return c * std::sqrt(normalized_energy * (normalized_energy + 2.0));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    kineticEnergyFromProperVelocity (Vector3 const& velocity, double const mass) noexcept
    {
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double elementary_charge = PhysConst::q_e_v<double>;

        auto const proper_speed_sq = dot(velocity, velocity);
        auto const gamma = std::sqrt(1.0 + proper_speed_sq * inv_c2);
        return mass * proper_speed_sq / ((gamma + 1.0) * elementary_charge);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Vector3 directionFromPolarAngle (
        Vector3 const& incident_direction, Vector3 const& transverse_1, Vector3 const& transverse_2,
        double const cosine, double const azimuth) noexcept
    {
        auto const sine = std::sqrt(amrex::max(1.0 - cosine * cosine, 0.0));
        auto const transverse =
            add(scale(transverse_1, std::cos(azimuth)), scale(transverse_2, std::sin(azimuth)));
        return add(scale(incident_direction, cosine), scale(transverse, sine));
    }

    /** Construct a right-handed transverse basis around a unit direction. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void transverseDirections (
        Vector3 const& direction, Vector3& transverse_1, Vector3& transverse_2) noexcept
    {
        auto const transverse_xy = std::sqrt(direction.x * direction.x + direction.y * direction.y);
        if (transverse_xy > std::numeric_limits<double>::epsilon())
        {
            transverse_1 = {-direction.y / transverse_xy, direction.x / transverse_xy, 0.0};
        }
        else
        {
            transverse_1 = {1.0, 0.0, 0.0};
        }
        transverse_2 = {direction.y * transverse_1.z - direction.z * transverse_1.y,
                        direction.z * transverse_1.x - direction.x * transverse_1.z,
                        direction.x * transverse_1.y - direction.y * transverse_1.x};
    }
} // namespace BackgroundMCCKinematics

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_KINEMATICS_H_
