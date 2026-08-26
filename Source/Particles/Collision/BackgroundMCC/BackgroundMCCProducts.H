/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_

#include "Particles/Algorithms/KineticEnergy.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Particle.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

/** Select particles whose previously chosen process belongs to one product group. */
class BackgroundMCCProductFilterFunc
{
public:
    BackgroundMCCProductFilterFunc (
        int const product_group,
        int const* selected_process,
        int const* process_product_group)
        : m_product_group(product_group),
          m_selected_process(selected_process),
          m_process_product_group(process_product_group)
    {}

    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (
        PData const&, int const i, amrex::RandomEngine const&) const noexcept
    {
        int const process = m_selected_process[i];
        return process >= 0 && m_process_product_group[process] == m_product_group;
    }

private:
    int m_product_group;
    int const* m_selected_process;
    int const* m_process_product_group;
};

namespace BackgroundMCCProducts
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void setProperVelocityFromPhysicalVelocity (
        amrex::ParticleReal const vx,
        amrex::ParticleReal const vy,
        amrex::ParticleReal const vz,
        amrex::ParticleReal& ux,
        amrex::ParticleReal& uy,
        amrex::ParticleReal& uz) noexcept
    {
        using namespace amrex::literals;
        using std::sqrt;

        amrex::ParticleReal const beta2 =
            (vx*vx + vy*vy + vz*vz)/PhysConst::c2;
        amrex::ParticleReal const gamma = 1.0_prt/sqrt(1.0_prt - beta2);
        ux = gamma*vx;
        uy = gamma*vy;
        uz = gamma*vz;
    }
}

/** Apply the existing equal-energy-sharing ionization outcome for a selected channel. */
class ImpactIonizationTransformFunc
{
public:
    ImpactIonizationTransformFunc (
        ScatteringProcess::Executor const* processes,
        int const* selected_process,
        amrex::ParticleReal const* neutral_vx,
        amrex::ParticleReal const* neutral_vy,
        amrex::ParticleReal const* neutral_vz,
        double const mass1)
        : m_processes(processes),
          m_selected_process(selected_process),
          m_neutral_vx(neutral_vx),
          m_neutral_vy(neutral_vy),
          m_neutral_vz(neutral_vz),
          m_mass1(mass1)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (
        DstData& dst1,
        DstData& dst2,
        SrcData& src,
        int const i_src,
        int const i_dst1,
        int const i_dst2,
        amrex::RandomEngine const& engine) const noexcept
    {
        using namespace amrex;
        using std::sqrt;

        auto& ux = src.m_rdata[PIdx::ux][i_src];
        auto& uy = src.m_rdata[PIdx::uy][i_src];
        auto& uz = src.m_rdata[PIdx::uz][i_src];

        auto& e_ux = dst1.m_rdata[PIdx::ux][i_dst1];
        auto& e_uy = dst1.m_rdata[PIdx::uy][i_dst1];
        auto& e_uz = dst1.m_rdata[PIdx::uz][i_dst1];
        auto& i_ux = dst2.m_rdata[PIdx::ux][i_dst2];
        auto& i_uy = dst2.m_rdata[PIdx::uy][i_dst2];
        auto& i_uz = dst2.m_rdata[PIdx::uz][i_dst2];

        int const process_index = m_selected_process[i_src];
        ParticleReal const energy_cost =
            m_processes[process_index].m_energy_penalty;

        constexpr auto eV = PhysConst::q_e;
        double const E_coll =
            Algorithms::KineticEnergy<double>(ux, uy, uz, m_mass1)/eV;
        double const E_available = E_coll > static_cast<double>(energy_cost)
            ? E_coll - static_cast<double>(energy_cost) : 0.0;
        ParticleReal const E_out = static_cast<ParticleReal>(0.5*E_available*eV);

        constexpr auto c2 = PhysConst::c2;
        auto const mc2 = m_mass1*c2;
        ParticleReal const up =
            sqrt(E_out*(E_out + 2.0_prt*mc2)/c2)/m_mass1;

        ParticleUtils::RandomizeVelocity(ux, uy, uz, up, engine);
        ParticleUtils::RandomizeVelocity(e_ux, e_uy, e_uz, up, engine);

        BackgroundMCCProducts::setProperVelocityFromPhysicalVelocity(
            m_neutral_vx[i_src],
            m_neutral_vy[i_src],
            m_neutral_vz[i_src],
            i_ux,
            i_uy,
            i_uz);
    }

private:
    ScatteringProcess::Executor const* m_processes;
    int const* m_selected_process;
    amrex::ParticleReal const* m_neutral_vx;
    amrex::ParticleReal const* m_neutral_vy;
    amrex::ParticleReal const* m_neutral_vz;
    double m_mass1;
};

/** Create a negative ion at the sampled neutral velocity and remove the electron. */
class AttachmentTransformFunc
{
public:
    AttachmentTransformFunc (
        amrex::ParticleReal const* neutral_vx,
        amrex::ParticleReal const* neutral_vy,
        amrex::ParticleReal const* neutral_vz)
        : m_neutral_vx(neutral_vx),
          m_neutral_vy(neutral_vy),
          m_neutral_vz(neutral_vz)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (
        DstData& dst,
        SrcData& src,
        int const i_src,
        int const i_dst,
        amrex::RandomEngine const&) const noexcept
    {
        BackgroundMCCProducts::setProperVelocityFromPhysicalVelocity(
            m_neutral_vx[i_src],
            m_neutral_vy[i_src],
            m_neutral_vz[i_src],
            dst.m_rdata[PIdx::ux][i_dst],
            dst.m_rdata[PIdx::uy][i_dst],
            dst.m_rdata[PIdx::uz][i_dst]);

        src.m_idcpu[i_src] = amrex::ParticleIdCpus::Invalid;
    }

private:
    amrex::ParticleReal const* m_neutral_vx;
    amrex::ParticleReal const* m_neutral_vy;
    amrex::ParticleReal const* m_neutral_vz;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_
