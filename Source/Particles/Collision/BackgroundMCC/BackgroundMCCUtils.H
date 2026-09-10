/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_BACKGROUND_MCC_UTILS_H_
#define WARPX_BACKGROUND_MCC_UTILS_H_

#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>

#include <cmath>

namespace BackgroundMCCUtils
{
    /**
     * \brief Return cross-section energy and collision-rate speed for an
     * electron colliding with a non-relativistic neutral background particle.
     *
     * The projectile components are WarpX normalized momenta, i.e. proper
     * velocities u = gamma*v. The neutral components are ordinary lab-frame
     * velocities sampled from a classical Maxwellian.
     *
     * For speed, the neutral velocity is subtracted directly from the electron
     * proper velocity and the result is treated as an approximate relative
     * proper velocity. This is correct in the non-relativistic limit, exact for
     * a stationary neutral, and has a relative error of order V_n/c when the
     * electron is relativistic. The latter is negligible for a classical gas.
     *
     * The lookup energy is evaluated directly from the approximate relative
     * proper velocity. For an electron incident on a heavy neutral, the
     * target-rest kinetic energy and total center-of-momentum kinetic energy
     * differ only by order m_e/M until extreme relativistic energies. Calling
     * ParticleUtils::getCollisionEnergy would add an extra square root without
     * a useful increase in accuracy for this model.
     *
     * @param[in] ux,uy,uz projectile proper-velocity components [m/s]
     * @param[in] neutral_vx,neutral_vy,neutral_vz neutral lab velocity [m/s]
     * @param[in] projectile_mass projectile rest mass [kg]
     * @param[out] collision_energy cross-section lookup energy [eV]
     * @param[out] collision_speed physical collision-rate speed [m/s]
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void getElectronNeutralCollisionParameters (
        const amrex::ParticleReal ux, const amrex::ParticleReal uy,
        const amrex::ParticleReal uz, const amrex::ParticleReal neutral_vx,
        const amrex::ParticleReal neutral_vy, const amrex::ParticleReal neutral_vz,
        const amrex::ParticleReal projectile_mass,
        amrex::ParticleReal& collision_energy, amrex::ParticleReal& collision_speed)
    {
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double q_e = PhysConst::q_e_v<double>;

        const double relative_ux = static_cast<double>(ux)
            - static_cast<double>(neutral_vx);
        const double relative_uy = static_cast<double>(uy)
            - static_cast<double>(neutral_vy);
        const double relative_uz = static_cast<double>(uz)
            - static_cast<double>(neutral_vz);
        const double relative_u_sq = relative_ux*relative_ux
            + relative_uy*relative_uy + relative_uz*relative_uz;
        const double relative_gamma = std::sqrt(1.0 + relative_u_sq*inv_c2);
        const double relative_u = std::sqrt(relative_u_sq);

        collision_energy = static_cast<amrex::ParticleReal>(
            static_cast<double>(projectile_mass)*relative_u_sq
                / ((relative_gamma + 1.0)*q_e));
        collision_speed = static_cast<amrex::ParticleReal>(
            relative_u/relative_gamma);
    }
}

#endif // WARPX_BACKGROUND_MCC_UTILS_H_
