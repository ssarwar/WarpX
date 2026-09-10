/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_IONIZATION_SAMPLING_H_
#define WARPX_PARTICLES_COLLISION_IONIZATION_SAMPLING_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <cstdint>
#include <limits>

namespace ProtonImpactIonization
{
    /** Shifted fixed-point Kronecker sequence on the open unit interval.
     *
     * Unsigned wraparound evaluates (shift + index*increment) modulo 2^32
     * exactly. An odd increment visits every phase before repeating, without
     * losing fractional bits at large indices in single-precision kernels.
     * Midpoint conversion avoids log(0) in the thermal sampler. Rounding the
     * largest float phases down by one ulp keeps the upper endpoint open.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    shiftedKronecker (std::uint32_t const index, std::uint32_t const shift,
                      std::uint32_t const increment) noexcept
    {
        auto const phase = shift + index * increment;
        auto const value = (static_cast<Real>(phase) + Real(.5)) * Real(0x1p-32);
        return value < Real(1) ? value : Real(1) - Real(.5) * std::numeric_limits<Real>::epsilon();
    }

    /** Base-2 radical inverse with an independent cell-wide uniform shift.
     *
     * Each index has a uniform marginal when shift is uniform on [0,1).
     * Consecutive blocks retain low discrepancy, including blocks assigned
     * to different parent energies by ordered source-weight selection.
     * The integer sequence has period 2^32. No per-product RNG is required.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    shiftedRadicalInverse (std::uint32_t index, Real const shift) noexcept
    {
        index = ((index & 0x55555555u) << 1) | ((index >> 1) & 0x55555555u);
        index = ((index & 0x33333333u) << 2) | ((index >> 2) & 0x33333333u);
        index = ((index & 0x0f0f0f0fu) << 4) | ((index >> 4) & 0x0f0f0f0fu);
        index = ((index & 0x00ff00ffu) << 8) | ((index >> 8) & 0x00ff00ffu);
        index = (index << 16) | (index >> 16);
        auto value = static_cast<Real>(index) * Real(0x1p-32) + shift;
        value -= value >= Real(1) ? Real(1) : Real(0);
        // A float conversion/addition can round the original sum to two.
        return value < Real(1) ? value : Real(0);
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_IONIZATION_SAMPLING_H_
