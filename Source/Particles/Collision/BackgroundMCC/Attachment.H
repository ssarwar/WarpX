/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_ATTACHMENT_H_
#define WARPX_PARTICLES_COLLISION_ATTACHMENT_H_

#include "Particles/WarpXParticleContainer.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Particles.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

/**
 * \brief Transform an attached electron into one negative-ion macroparticle.
 *
 * The neutral and any neutral fragment or third body remain part of the fixed
 * background reservoir. The negative ion is initialized from the local gas
 * temperature and the source electron is invalidated. Position, weight and
 * common runtime attributes are copied before this transform is called.
 */
class AttachmentTransformFunc
{
public:
    AttachmentTransformFunc (
        amrex::ParticleReal sqrt_kb_m_product,
        amrex::ParserExecutor<4> const& temperature_func,
        amrex::Real t)
        : m_sqrt_kb_m_product(sqrt_kb_m_product),
          m_temperature_func(temperature_func),
          m_t(t)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (
        DstData& dst, SrcData& src, int const i_src, int const i_dst,
        amrex::RandomEngine const& engine) const noexcept
    {
        using namespace amrex::literals;

        const auto& p = src.getSuperParticle(i_src);
        amrex::ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const amrex::ParticleReal temperature =
            amrex::max(m_temperature_func(x, y, z, m_t), 0.0_prt);
        const amrex::ParticleReal velocity_std =
            m_sqrt_kb_m_product * std::sqrt(temperature);

        dst.m_rdata[PIdx::ux][i_dst] =
            velocity_std * amrex::RandomNormal(0.0_prt, 1.0_prt, engine);
        dst.m_rdata[PIdx::uy][i_dst] =
            velocity_std * amrex::RandomNormal(0.0_prt, 1.0_prt, engine);
        dst.m_rdata[PIdx::uz][i_dst] =
            velocity_std * amrex::RandomNormal(0.0_prt, 1.0_prt, engine);

        src.m_idcpu[i_src] = amrex::ParticleIdCpus::Invalid;
    }

private:
    amrex::ParticleReal m_sqrt_kb_m_product;
    amrex::ParserExecutor<4> m_temperature_func;
    amrex::Real m_t;
};

#endif // WARPX_PARTICLES_COLLISION_ATTACHMENT_H_
