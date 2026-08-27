/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_

#include "Particles/Algorithms/KineticEnergy.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCIonization.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Particle.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>

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

/** Apply the configured ionization energy-sharing outcome for a selected channel. */
class ImpactIonizationTransformFunc
{
public:
    ImpactIonizationTransformFunc (
        ScatteringProcess::Executor const* processes,
        BackgroundMCCIonizationModel::Executor const* ionization_processes,
        int const* selected_process,
        amrex::ParticleReal const* neutral_vx,
        amrex::ParticleReal const* neutral_vy,
        amrex::ParticleReal const* neutral_vz,
        double const mass1)
        : m_processes(processes),
          m_ionization_processes(ionization_processes),
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
        using namespace amrex::literals;

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
        ParticleReal const relative_ux = ux - m_neutral_vx[i_src];
        ParticleReal const relative_uy = uy - m_neutral_vy[i_src];
        ParticleReal const relative_uz = uz - m_neutral_vz[i_src];
        double const E_coll = Algorithms::KineticEnergy<double>(
            relative_ux, relative_uy, relative_uz, m_mass1)/eV;

        ParticleReal binding_energy;
        ParticleReal secondary_energy;
        m_ionization_processes[process_index].sampleEnergy(
            static_cast<ParticleReal>(E_coll),
            energy_cost,
            amrex::Random(engine),
            amrex::Random(engine),
            binding_energy,
            secondary_energy);
        auto const primary_energy = amrex::max(
            static_cast<ParticleReal>(E_coll) - binding_energy - secondary_energy,
            0.0_prt);

        constexpr double c2 = PhysConst::c2_v<double>;
        auto const mc2 = m_mass1*c2;
        auto const primary_energy_joule =
            static_cast<double>(primary_energy)*eV;
        auto const secondary_energy_joule =
            static_cast<double>(secondary_energy)*eV;
        auto const primary_up = static_cast<ParticleReal>(
            std::sqrt(primary_energy_joule*(primary_energy_joule + 2.0*mc2)/c2)/
            m_mass1);
        auto const secondary_up = static_cast<ParticleReal>(
            std::sqrt(
                secondary_energy_joule*(secondary_energy_joule + 2.0*mc2)/c2)/
            m_mass1);

        ParticleUtils::RandomizeVelocity(ux, uy, uz, primary_up, engine);
        ParticleUtils::RandomizeVelocity(
            e_ux, e_uy, e_uz, secondary_up, engine);

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
    BackgroundMCCIonizationModel::Executor const* m_ionization_processes;
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
