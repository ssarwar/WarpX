/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXFluidContainer.H"

#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>

#include <cmath>

std::array<amrex::Real, 5>
WarpXFluidContainer::PhysicalTotals () const
{
    auto const& warpx = WarpX::GetInstance();
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
                     amrex::ReduceOpSum>
        reduce;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> data(reduce);
    using Tuple = decltype(data)::Type;
    auto const model = m_model;
    auto const species_mass = mass;
    double const velocity = m_rigid_beam ? m_rigid_beam->velocity() : 0.0;
    double const beam_energy =
        m_rigid_beam ? m_rigid_beam->kineticEnergyEV() * PhysConst::q_e : 0.0;
    double const beam_momentum =
        species_mass * velocity / std::sqrt(1 - velocity * velocity / PhysConst::c2);
    bool correction = true;
    amrex::ParmParse("boundary").query("verboncoeur_axis_correction", correction);
    double const axis_factor = correction ? 1.0 / 3.0 : 1.0 / 4.0;
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev) {
        auto const& number = *warpx.m_fields.get(name_mf_N, lev);
        auto const& geom = warpx.Geom(lev);
        auto const dx = geom.CellSizeArray();
        auto const lo = geom.ProbLoArray();
        auto const domain = geom.Domain();
        auto const nodal = number.ixType().toIntVect();
        auto const ng = model == FluidModel::Immobile ? number.nGrowVect() : amrex::IntVect(0);
        auto const owners = amrex::OwnerMask(number, geom.periodicity(), ng);
        amrex::iMultiFab fine_mask;
        if (lev < warpx.finestLevel()) {
            fine_mask =
                amrex::makeFineMask(number, warpx.boxArray(lev + 1), warpx.refRatio(lev), 1, 0);
        }
        for (amrex::MFIter mfi(number, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const n = number.const_array(mfi);
            auto const owner = owners->const_array(mfi);
            auto const mask =
                fine_mask.ok() ? fine_mask.const_array(mfi) : amrex::Array4<int const>{};
            amrex::Array4<amrex::Real const> momentum[3];
            if (model == FluidModel::ColdRelativistic) {
                for (int dir = 0; dir < 3; ++dir) {
                    momentum[dir] =
                        warpx.m_fields.get(name_mf_NU, ablastr::fields::Direction{dir}, lev)
                            ->const_array(mfi);
                }
            }
            auto const mx = momentum[0], my = momentum[1], mz = momentum[2];
            reduce.eval(mfi.growntilebox(ng), data,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> Tuple {
                            if (!owner(i, j, k) || (mask && !mask(i, j, k))) {
                                return {0, 0, 0, 0, 0};
                            }
                            double volume = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
#ifdef WARPX_DIM_RZ
                            double const r =
                                lo[0] + (i - domain.smallEnd(0) + 0.5 * (1 - nodal[0])) * dx[0];
                            if (r < 0.0) {
                                return {0, 0, 0, 0, 0};
                            }
                            volume *= r == 0.0 ? MathConst::pi * dx[0] * axis_factor
                                               : 2 * MathConst::pi * r;
#else
                    amrex::ignore_unused(lo, domain, nodal, axis_factor);
#endif
                            double const count = n(i, j, k) * volume;
                            if (model == FluidModel::RigidBeam) {
                                return {amrex::Real(count), amrex::Real(count * beam_energy), 0, 0,
                                        amrex::Real(count * beam_momentum)};
                            }
                            if (model == FluidModel::Immobile || count == 0.0) {
                                return {amrex::Real(count), 0, 0, 0, 0};
                            }
                            double const ux = mx(i, j, k) / n(i, j, k),
                                         uy = my(i, j, k) / n(i, j, k),
                                         uz = mz(i, j, k) / n(i, j, k);
                            double const u2 = ux * ux + uy * uy + uz * uz;
                            double const energy =
                                species_mass * u2 / (1 + std::sqrt(1 + u2 / PhysConst::c2));
#ifdef WARPX_DIM_RZ
                            // Azimuthal integration of an axisymmetric fluid gives zero
                            // net transverse Cartesian momentum.
                            return {amrex::Real(count), amrex::Real(count * energy), 0, 0,
                                    amrex::Real(count * species_mass * uz)};
#else
                    return {amrex::Real(count), amrex::Real(count*energy),
                            amrex::Real(count*species_mass*ux), amrex::Real(count*species_mass*uy),
                            amrex::Real(count*species_mass*uz)};
#endif
                        });
        }
    }
    auto const value = data.value();
    std::array<amrex::Real, 5> totals{amrex::get<0>(value), amrex::get<1>(value),
                                      amrex::get<2>(value), amrex::get<3>(value),
                                      amrex::get<4>(value)};
    amrex::ParallelDescriptor::ReduceRealSum(totals.data(), 5);
    return totals;
}
