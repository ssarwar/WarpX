/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

namespace ProtonImpactIonization
{
    /** Free-electron binary-encounter polar cosine for a heavy projectile. */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    freeElectronCosine (Real const secondary_energy, Real const maximum_transfer) noexcept
    {
        using std::sqrt;

        if (secondary_energy <= Real(0) || maximum_transfer <= Real(0)) {
            return Real(0);
        }
        constexpr auto electron_rest_energy = Real(510998.95069);
        // The energy model has a molecular tail. Its free reference direction
        // is forward above Tmax; do not evaluate free kinematics out of domain.
        auto const free_energy = amrex::min(secondary_energy, maximum_transfer);
        return sqrt(amrex::min(
            Real(1), free_energy * (maximum_transfer + Real(2) * electron_rest_energy) /
                         (maximum_transfer * (free_energy + Real(2) * electron_rest_energy))));
    }

    /** Practical IAA bound/free angular closure for heavy-projectile ionization.
     *
     * The free-electron factor is replaced by the exact heavy-projectile
     * binary-encounter cosine. The original IAA binding interpolation and its
     * isotropic bound component are retained. Condition the resulting uniform
     * interval on [-1,1], rather than clipping draws to a delta peak at mu=1.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    polarCosine (Real const secondary_energy, Real const binding_energy,
                 Real const maximum_transfer, Real const bound_uniform) noexcept
    {
        auto const denominator = secondary_energy + binding_energy;
        if (denominator <= Real(0)) {
            return Real(2) * bound_uniform - Real(1);
        }

        auto const free_cosine = freeElectronCosine(secondary_energy, maximum_transfer);
        auto const center =
            free_cosine * (secondary_energy + Real(.5) * binding_energy) / denominator;
        auto const width = binding_energy / denominator;
        auto const lower = amrex::max(Real(-1), center - width);
        auto const upper = amrex::min(Real(1), center + width);
        return (Real(1) - bound_uniform) * lower + bound_uniform * upper;
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_
