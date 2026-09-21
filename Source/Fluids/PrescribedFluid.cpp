/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXFluidContainer.H"

#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <ablastr/utils/Communication.H>

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

void
WarpXFluidContainer::InitPrescribedDensity (ablastr::fields::MultiFabRegister& fields, int lev)
{
    auto& density = *fields.get(name_mf_N, lev);
    if (!h_inj_rho) {
        return;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !h_inj_rho->needPreparation() && !h_inj_rho->distributed(),
        "Immobile initial density supports constant and parser profiles.");
    auto const& geom = WarpX::GetInstance().Geom(lev);
    auto const dx = geom.CellSizeArray();
    auto const lo = geom.ProbLoArray();
    auto const type = density.ixType().toIntVect();
    auto const* injector = d_inj_rho;
    for (amrex::MFIter mfi(density, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const array = density.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef WARPX_DIM_RZ
            auto const r = lo[0] + (i + 0.5 * (1 - type[0])) * dx[0];
            auto const z = lo[1] + (j + 0.5 * (1 - type[1])) * dx[1];
            array(i, j, k) = injector->getDensity(r, 0.0, z);
#else
            amrex::ignore_unused(array, injector, type, dx, lo, i, j, k);
#endif
        });
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        density.min(0) >= 0.0 && !density.contains_nan() && !density.contains_inf(),
        "Immobile initial number density must be finite and non-negative.");
    density.FillBoundary(geom.periodicity());
}

void
WarpXFluidContainer::CommitDensityIncrement (ablastr::fields::MultiFabRegister& fields, int lev)
{
    auto& warpx = WarpX::GetInstance();
    auto& increment = *fields.get(DensityIncrementName(), lev);
    auto& density = *fields.get(name_mf_N, lev);
#ifdef WARPX_DIM_RZ
    warpx.ApplyInverseVolumeScalingToChargeDensity(&increment, lev);
    // The axis fold has consumed these raw deposits. They are not a physical
    // negative-radius population and must not enter persistent state.
    if (warpx.Geom(lev).ProbLo(0) == 0.0) {
        for (amrex::MFIter mfi(increment, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const values = increment.array(mfi);
            auto box = mfi.growntilebox();
            box.setBig(0, std::min(box.bigEnd(0), warpx.Geom(lev).Domain().smallEnd(0) - 1));
            amrex::ParallelFor(
                box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept { values(i, j, k) = 0.0; });
        }
    }
#endif
    // Only fresh deposits are summed. Summing an accumulated nodal population
    // would duplicate it at every grid interface on each collision call.
    ablastr::utils::communication::SumBoundary(
        increment, 0, 1, increment.nGrowVect(), increment.nGrowVect(),
        WarpX::do_single_precision_comms, warpx.Geom(lev).periodicity());
    // Retain the charge-shape support beyond physical walls as well. Boundary
    // reflection and filtering act on the transient total charge, just as for
    // frozen kinetic ions, rather than truncating each persistent increment.
    amrex::MultiFab::Add(density, increment, 0, 0, 1, density.nGrowVect());
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        density.min(0) >= 0.0 && !density.contains_nan() && !density.contains_inf(),
        "An immobile-fluid update produced a negative or non-finite number density.");
    density.FillBoundary(warpx.Geom(lev).periodicity());
    increment.setVal(0.0);
}
