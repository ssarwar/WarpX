/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCKINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCKINEMATICS_H_

#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

/**
 * Relativistic projectile-target kinematics used by background MCC.
 *
 * WarpX stores proper velocity u = gamma*v. The neutral velocity sampled from
 * a Maxwellian is an ordinary velocity because the neutral background is
 * non-relativistic in the lab frame. This helper converts the target velocity
 * to proper velocity, Lorentz-transforms the projectile to the target rest
 * frame, and returns the physical relative speed and collision energy.
 */
struct BackgroundMCCKinematics
{
    amrex::ParticleReal target_ux = 0.0;
    amrex::ParticleReal target_uy = 0.0;
    amrex::ParticleReal target_uz = 0.0;
    amrex::ParticleReal relative_speed = 0.0;
    amrex::ParticleReal collision_energy_eV = 0.0;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
BackgroundMCCKinematics getBackgroundMCCKinematics (
    amrex::ParticleReal projectile_ux,
    amrex::ParticleReal projectile_uy,
    amrex::ParticleReal projectile_uz,
    amrex::ParticleReal target_vx,
    amrex::ParticleReal target_vy,
    amrex::ParticleReal target_vz,
    double projectile_mass,
    double target_mass) noexcept
{
    using namespace amrex::literals;

    BackgroundMCCKinematics result;

    amrex::ParticleReal target_v2 =
        target_vx*target_vx + target_vy*target_vy + target_vz*target_vz;
    // The background model is intended for non-relativistic neutrals. Protect
    // the Lorentz factor from an exceedingly rare unphysical Maxwellian draw.
    const amrex::ParticleReal max_v2 =
        (1.0_prt - 32.0_prt*std::numeric_limits<amrex::ParticleReal>::epsilon())
        * PhysConst::c2;
    if (target_v2 > max_v2) {
        const amrex::ParticleReal scale = std::sqrt(max_v2/target_v2);
        target_vx *= scale;
        target_vy *= scale;
        target_vz *= scale;
        target_v2 = max_v2;
    }

    const amrex::ParticleReal target_gamma =
        1.0_prt/std::sqrt(1.0_prt - target_v2*PhysConst::inv_c2);
    result.target_ux = target_gamma*target_vx;
    result.target_uy = target_gamma*target_vy;
    result.target_uz = target_gamma*target_vz;

    amrex::ParticleReal relative_ux = projectile_ux;
    amrex::ParticleReal relative_uy = projectile_uy;
    amrex::ParticleReal relative_uz = projectile_uz;
    ParticleUtils::doLorentzTransformWithU(
        relative_ux, relative_uy, relative_uz,
        result.target_ux, result.target_uy, result.target_uz);

    const amrex::ParticleReal relative_u2 =
        relative_ux*relative_ux + relative_uy*relative_uy + relative_uz*relative_uz;
    const amrex::ParticleReal relative_gamma =
        std::sqrt(1.0_prt + relative_u2*PhysConst::inv_c2);
    result.relative_speed =
        (relative_gamma > 0.0_prt) ? std::sqrt(relative_u2)/relative_gamma : 0.0_prt;

    double gamma_unused;
    double collision_energy_eV;
    ParticleUtils::getCollisionEnergy(
        relative_u2, projectile_mass, target_mass,
        gamma_unused, collision_energy_eV);
    result.collision_energy_eV =
        static_cast<amrex::ParticleReal>(collision_energy_eV);

    return result;
}

/**
 * Return the projectile speed in the target rest frame corresponding to a
 * center-of-momentum collision energy in eV.
 *
 * For a target initially at rest,
 *
 *   s = m^2 c^4 + M^2 c^4 + 2 gamma m M c^4,
 *
 * while sqrt(s) = (m+M)c^2 + E_coll. The factored expression below avoids
 * cancellation when E_coll is small.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::ParticleReal relativeSpeedFromCollisionEnergy (
    amrex::ParticleReal collision_energy_eV,
    double projectile_mass,
    double target_mass) noexcept
{
    using namespace amrex::literals;

    if (collision_energy_eV <= 0.0_prt) {
        return 0.0_prt;
    }

    const double collision_energy =
        static_cast<double>(collision_energy_eV)*PhysConst::q_e;
    const double mass_sum_c2 =
        (projectile_mass + target_mass)*PhysConst::c2;
    const double gamma_minus_one =
        collision_energy*(2.0*mass_sum_c2 + collision_energy)
        /(2.0*projectile_mass*target_mass*PhysConst::c2*PhysConst::c2);
    const double gamma = 1.0 + gamma_minus_one;
    const double beta2 =
        gamma_minus_one*(gamma + 1.0)/(gamma*gamma);

    return static_cast<amrex::ParticleReal>(
        PhysConst::c*std::sqrt(amrex::max(beta2, 0.0)));
}

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCCKINEMATICS_H_
