/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_IMMOBILE_DENSITY_DEPOSIT_H_
#define WARPX_IMMOBILE_DENSITY_DEPOSIT_H_

#include "Particles/ShapeFactors.H"

#include <AMReX_Array4.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuQualifiers.H>

/** Deposit a signed number increment using the kinetic electron's charge shape.
 *
 * Call from amrex::For, never ParallelFor: events can touch the same grid node.
 * HostDevice atomics also protect tiles of a shared FAB on OpenMP builds. The
 * caller applies RZ volume scaling and synchronization once to the fresh field.
 */
struct ImmobileDensityDeposit
{
    amrex::Array4<amrex::Real> m_increment;
    amrex::Real m_rmin = 0.0;
    amrex::Real m_zmin = 0.0;
    amrex::Real m_dr_inv = 0.0;
    amrex::Real m_dz_inv = 0.0;
    int m_rlo = 0;
    int m_zlo = 0;
    int m_rnodal = 0;
    int m_znodal = 0;
    int m_order = 1;

    template <int Order>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    deposit (amrex::Real r, amrex::Real z, amrex::Real weight) const noexcept
    {
#ifdef WARPX_DIM_RZ
        Compute_shape_factor<Order> const compute;
        amrex::Real sr[Order+1], sz[Order+1];
        int const ir = m_rlo + compute(sr, (r-m_rmin)*m_dr_inv - 0.5*(1-m_rnodal));
        int const iz = m_zlo + compute(sz, (z-m_zmin)*m_dz_inv - 0.5*(1-m_znodal));
        auto const density = weight*m_dr_inv*m_dz_inv;
        for (int j = 0; j <= Order; ++j) {
            for (int i = 0; i <= Order; ++i) {
                amrex::HostDevice::Atomic::Add(&m_increment(ir+i, iz+j, 0),
                                               sr[i]*sz[j]*density);
            }
        }
#else
        amrex::ignore_unused(r, z, weight);
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator() (amrex::Real r, amrex::Real z, amrex::Real weight) const noexcept
    {
        switch (m_order) {
        case 1: deposit<1>(r, z, weight); break;
        case 2: deposit<2>(r, z, weight); break;
        case 3: deposit<3>(r, z, weight); break;
        case 4: deposit<4>(r, z, weight); break;
        }
    }
};

#endif
