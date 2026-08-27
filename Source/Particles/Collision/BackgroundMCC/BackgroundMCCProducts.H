/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCPRODUCTS_H_

#include "Particles/Collision/BackgroundMCC/BackgroundMCCIonization.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCIonizationKinematics.H"
#include "Particles/Collision/ScatteringProcess.H"
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
        double const mass1,
        double const ion_mass)
        : m_processes(processes),
          m_ionization_processes(ionization_processes),
          m_selected_process(selected_process),
          m_neutral_vx(neutral_vx),
          m_neutral_vy(neutral_vy),
          m_neutral_vz(neutral_vz),
          m_mass1(mass1),
          m_ion_mass(ion_mass)
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

        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double elementary_charge = PhysConst::q_e_v<double>;
        double const relative_ux = static_cast<double>(ux) - m_neutral_vx[i_src];
        double const relative_uy = static_cast<double>(uy) - m_neutral_vy[i_src];
        double const relative_uz = static_cast<double>(uz) - m_neutral_vz[i_src];
        double const relative_u_sq =
            relative_ux*relative_ux + relative_uy*relative_uy +
            relative_uz*relative_uz;
        double const relative_gamma = std::sqrt(
            1.0 + relative_u_sq*inv_c2);
        double const collision_energy =
            m_mass1*relative_u_sq/
            ((relative_gamma + 1.0)*elementary_charge);

        ParticleReal binding_energy;
        ParticleReal secondary_energy;
        m_ionization_processes[process_index].sampleEnergy(
            static_cast<ParticleReal>(collision_energy),
            energy_cost,
            amrex::Random(engine),
            amrex::Random(engine),
            binding_energy,
            secondary_energy);
        auto const available_energy = amrex::max(
            static_cast<ParticleReal>(collision_energy) - binding_energy,
            0.0_prt);
        auto const secondary_fraction = available_energy > 0.0_prt
            ? secondary_energy/available_energy : 0.0_prt;

        BackgroundMCCIonizationKinematics::compute(
            ux,
            uy,
            uz,
            m_neutral_vx[i_src],
            m_neutral_vy[i_src],
            m_neutral_vz[i_src],
            m_mass1,
            m_ion_mass,
            binding_energy,
            secondary_fraction,
            m_processes[process_index].m_scattering_angle_model,
            engine,
            ux,
            uy,
            uz,
            e_ux,
            e_uy,
            e_uz,
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
    double m_ion_mass;
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
