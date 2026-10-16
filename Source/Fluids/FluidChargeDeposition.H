/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FLUID_CHARGE_DEPOSITION_H_
#define WARPX_FLUID_CHARGE_DEPOSITION_H_

#include "Utils/WarpXConst.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

namespace warpx::fluids
{
    /** Radial factor removed by ApplyInverseVolumeScalingToChargeDensity. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real RadialVolumeFactor (
        [[maybe_unused]] amrex::Real r, [[maybe_unused]] amrex::Real dr,
        [[maybe_unused]] bool use_axis_correction)
    {
        using namespace amrex::literals;
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER)
        r = amrex::Math::abs(r);
        const amrex::Real axis_factor = use_axis_correction ? 1.0_rt/3.0_rt : 1.0_rt/4.0_rt;
        return (r == 0.0_rt) ? MathConst::pi*dr*axis_factor : 2.0_rt*MathConst::pi*r;
#elif defined(WARPX_DIM_RSPHERE)
        r = amrex::Math::abs(r);
        const amrex::Real axis_factor = use_axis_correction ? 1.0_rt/4.0_rt : 1.0_rt/8.0_rt;
        return (r == 0.0_rt) ? 4.0_rt/3.0_rt*MathConst::pi*dr*dr*axis_factor
                            : 4.0_rt*MathConst::pi*r*r;
#else
        return 1.0_rt;
#endif
    }
}

#endif
