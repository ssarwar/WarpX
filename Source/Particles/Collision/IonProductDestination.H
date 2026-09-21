/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ION_PRODUCT_DESTINATION_H_
#define WARPX_ION_PRODUCT_DESTINATION_H_

#include "Fluids/ImmobileDensityDeposit.H"
#include "Fluids/MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX_MFIter.H>

#include <string>

/** Collision products may be kinetic particles or an immobile number density. */
class IonProductDestination
{
public:
    IonProductDestination (std::string const& name, MultiParticleContainer const& particles)
    {
        auto& warpx = WarpX::GetInstance();
        if (warpx.DoFluidSpecies()) { m_fluid = warpx.GetFluidContainer().FindSpecies(name); }
        if (m_fluid) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_fluid->getModel() == FluidModel::Immobile,
                "Collision ion products require a kinetic or immobile species: '" + name + "'.");
        } else {
            m_particles = &particles.GetParticleContainerFromName(name);
        }
    }

    [[nodiscard]] amrex::Real getMass () const
        { return m_fluid ? m_fluid->getMass() : m_particles->getMass(); }
    [[nodiscard]] amrex::Real getCharge () const
        { return m_fluid ? m_fluid->getCharge() : m_particles->getCharge(); }
    [[nodiscard]] WarpXParticleContainer* particles () const { return m_particles; }
    [[nodiscard]] bool isFluid () const { return m_fluid != nullptr; }

    [[nodiscard]] ImmobileDensityDeposit deposit (int lev, amrex::MFIter const& mfi) const
    {
        if (!m_fluid) { return {}; }
        auto& warpx = WarpX::GetInstance();
        auto& increment = *warpx.m_fields.get(m_fluid->DensityIncrementName(), lev);
        auto const& geom = warpx.Geom(lev);
        ImmobileDensityDeposit result;
        result.m_increment = increment.array(mfi);
#ifdef WARPX_DIM_RZ
        auto const lo = amrex::lbound(increment[mfi].box());
        result.m_rmin = geom.ProbLo(0) +
            (lo.x-geom.Domain().smallEnd(0))*geom.CellSize(0);
        result.m_zmin = geom.ProbLo(1) +
            (lo.y-geom.Domain().smallEnd(1))*geom.CellSize(1);
        result.m_dr_inv = geom.InvCellSize(0);
        result.m_dz_inv = geom.InvCellSize(1);
        result.m_rlo = lo.x;
        result.m_zlo = lo.y;
        result.m_rnodal = increment.ixType()[0];
        result.m_znodal = increment.ixType()[1];
        result.m_order = WarpX::nox;
#else
        amrex::ignore_unused(geom);
#endif
        return result;
    }

    void commit (int lev) const
    {
        if (m_fluid) {
            m_fluid->CommitDensityIncrement(WarpX::GetInstance().m_fields, lev);
        }
    }

private:
    WarpXParticleContainer* m_particles = nullptr;
    WarpXFluidContainer* m_fluid = nullptr;
};

#endif
