/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_BACKGROUND_MCC_UTILS_H_
#define WARPX_BACKGROUND_MCC_UTILS_H_

#include "Utils/WarpXConst.H"

#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace BackgroundMCCUtils
{
    /**
     * \brief Return target-rest-frame energy and invariant flux speed for an
     * electron colliding with a non-relativistic neutral background particle.
     *
     * The projectile components are WarpX normalized momenta, i.e. proper
     * velocities u = gamma*v. The neutral components are ordinary lab-frame
     * velocities sampled from a classical Maxwellian. A scalar relative
     * Lorentz factor is sufficient for the cross-section energy and Moller
     * speed, so this avoids a full three-vector Lorentz transformation.
     *
     * The neutral Lorentz factor uses its non-relativistic expansion through
     * order (v_n/c)^2. This is substantially more accurate than the neutral-gas
     * model while avoiding an additional square root for every collision
     * candidate.
     *
     * @param[in] ux,uy,uz projectile proper-velocity components [m/s]
     * @param[in] neutral_vx,neutral_vy,neutral_vz neutral lab velocity [m/s]
     * @param[in] projectile_mass projectile rest mass [kg]
     * @param[out] collision_energy projectile energy in the neutral frame [eV]
     * @param[out] collision_speed invariant flux speed [m/s]
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void getElectronNeutralCollisionParameters (
        const amrex::ParticleReal ux, const amrex::ParticleReal uy,
        const amrex::ParticleReal uz, const amrex::ParticleReal neutral_vx,
        const amrex::ParticleReal neutral_vy, const amrex::ParticleReal neutral_vz,
        const amrex::ParticleReal projectile_mass,
        amrex::ParticleReal& collision_energy, amrex::ParticleReal& collision_speed)
    {
        constexpr double c = PhysConst::c_v<double>;
        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double q_e = PhysConst::q_e_v<double>;

        const double ux_d = static_cast<double>(ux);
        const double uy_d = static_cast<double>(uy);
        const double uz_d = static_cast<double>(uz);
        const double neutral_vx_d = static_cast<double>(neutral_vx);
        const double neutral_vy_d = static_cast<double>(neutral_vy);
        const double neutral_vz_d = static_cast<double>(neutral_vz);

        const double u_sq = ux_d*ux_d + uy_d*uy_d + uz_d*uz_d;
        const double gamma = std::sqrt(1.0 + u_sq*inv_c2);
        const double gamma_minus_one = u_sq*inv_c2/(gamma + 1.0);

        const double neutral_v_sq = neutral_vx_d*neutral_vx_d
            + neutral_vy_d*neutral_vy_d + neutral_vz_d*neutral_vz_d;
        const double neutral_gamma_minus_one = 0.5*neutral_v_sq*inv_c2;
        const double neutral_gamma = 1.0 + neutral_gamma_minus_one;

        const double u_dot_neutral_v = ux_d*neutral_vx_d
            + uy_d*neutral_vy_d + uz_d*neutral_vz_d;
        const double relative_gamma_minus_one = amrex::max(
            gamma_minus_one + neutral_gamma_minus_one*gamma
                - neutral_gamma*u_dot_neutral_v*inv_c2,
            0.0);

        collision_energy = static_cast<amrex::ParticleReal>(
            static_cast<double>(projectile_mass)*c2*relative_gamma_minus_one/q_e);
        collision_speed = static_cast<amrex::ParticleReal>(
            c*std::sqrt(relative_gamma_minus_one*(relative_gamma_minus_one + 2.0))
                /(gamma*neutral_gamma));
    }
}

#endif // WARPX_BACKGROUND_MCC_UTILS_H_
