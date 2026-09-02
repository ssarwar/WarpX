/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_

#include <AMReX_REAL.H>

#include <cmath>

namespace ProtonImpactIonization
{
    /** Free-electron binary-encounter polar cosine for a heavy projectile. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
    freeElectronCosine (amrex::ParticleReal const secondary_energy,
                        amrex::ParticleReal const maximum_transfer) noexcept
    {
        using namespace amrex::literals;
        using std::sqrt;

        if (secondary_energy <= 0.0_prt || maximum_transfer <= 0.0_prt) {
            return 0.0_prt;
        }
        constexpr auto electron_rest_energy = 510998.95069_prt;
        return sqrt(
            amrex::max(secondary_energy * (maximum_transfer + 2.0_prt * electron_rest_energy) /
                           (maximum_transfer * (secondary_energy + 2.0_prt * electron_rest_energy)),
                       0.0_prt));
    }

    /** Practical IAA bound/free angular closure for heavy-projectile ionization.
     *
     * The free-electron factor is replaced by the exact heavy-projectile
     * binary-encounter cosine. The original IAA binding interpolation and its
     * isotropic bound component are retained.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
    polarCosine (amrex::ParticleReal const secondary_energy,
                 amrex::ParticleReal const binding_energy,
                 amrex::ParticleReal const maximum_transfer,
                 amrex::ParticleReal const bound_uniform) noexcept
    {
        using namespace amrex::literals;

        auto const denominator = secondary_energy + binding_energy;
        auto const bound_cosine = 2.0_prt * bound_uniform - 1.0_prt;
        if (denominator <= 0.0_prt) {
            return bound_cosine;
        }

        auto const free_cosine = freeElectronCosine(secondary_energy, maximum_transfer);
        return amrex::max(
            -1.0_prt,
            amrex::min(free_cosine * (secondary_energy + 0.5_prt * binding_energy) / denominator +
                           binding_energy / denominator * bound_cosine,
                       1.0_prt));
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_
