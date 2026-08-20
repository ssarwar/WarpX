/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCCollision.H"

#include "Attachment.H"
#include "BackgroundMCCKinematics.H"
#include "ImpactIonization.H"
#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Particles/Collision/BinaryCollision/TwoProductUtil.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/ParticleUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "WarpX.H"

#include <ablastr/profiler/ProfilerWrapper.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_INT.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

void
BackgroundMCCCollision::doCollisions (
    amrex::Real cur_time, amrex::Real dt,
    MultiParticleContainer* mypc)
{
    ABLASTR_PROFILE("BackgroundMCCCollision::doCollisions()");
    using namespace amrex::literals;

    if (!m_initialized) {
        initialize(mypc);
    }

    if (m_nu_max <= 0.0_prt) {
        return;
    }

    const amrex::ParticleReal coll_n = m_nu_max*dt;
    const amrex::ParticleReal total_collision_probability =
        1.0_prt - std::exp(-coll_n);

    if (coll_n > 0.1_prt && !m_warned_large_step) {
        ablastr::warn_manager::WMRecordWarning(
            "BackgroundMCC Collisions",
            "nu_max*dt = " + std::to_string(coll_n)
            + " is greater than 0.1. Use <collision>.ndt_subcycle so each "
              "collision substep has nu_max*dt_sub <= 0.1.");
        m_warned_large_step = true;
    }

    auto& source_species =
        mypc->GetParticleContainerFromName(m_species_names[0]);

    const int finest_level = source_species.finestLevel();
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto* cost = WarpX::getCosts(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(source_species, lev); pti.isValid(); ++pti) {
            if (cost &&
                WarpX::load_balance_costs_update_algo ==
                    LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            auto wt = static_cast<amrex::Real>(amrex::second());

            const long source_particle_count = pti.numParticles();
            amrex::Gpu::DeviceVector<int> event_process(source_particle_count);

            selectCollisionEventsWithinTile(
                pti, cur_time, total_collision_probability, event_process);
            createReactionProductsWithinTile(
                lev, pti, source_particle_count, event_process,
                source_species, mypc, cur_time);

            if (cost &&
                WarpX::load_balance_costs_update_algo ==
                    LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = static_cast<amrex::Real>(amrex::second()) - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }
    }
}

void
BackgroundMCCCollision::selectCollisionEventsWithinTile (
    WarpXParIter& pti,
    amrex::Real t,
    amrex::ParticleReal total_collision_probability,
    amrex::Gpu::DeviceVector<int>& event_process) const
{
    using namespace amrex::literals;

    const long np = pti.numParticles();
    if (np == 0) {
        return;
    }

    int* event_process_ptr = event_process.dataPtr();
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (long ip) noexcept
        {
            event_process_ptr[ip] = -1;
        });

    auto n_a_func = m_background_density_func;
    auto T_a_func = m_background_temperature_func;
    auto* processes = m_processes_exe.data();
    auto* rate_multipliers = m_rate_multipliers_d.data();
    const int process_count = static_cast<int>(m_processes_exe.size());
    const amrex::ParticleReal nu_max = m_nu_max;
    const double projectile_mass = m_projectile_mass;
    const double background_mass = m_background_mass;

    auto GetPosition = GetParticlePosition<PIdx>(pti);
    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux =
        attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy =
        attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz =
        attribs[PIdx::uz].dataPtr();
    auto* const AMREX_RESTRICT idcpu =
        pti.GetStructOfArrays().GetIdCPUData().data();

    amrex::ParallelForRNG(
        np,
        [=] AMREX_GPU_HOST_DEVICE (
            long ip, amrex::RandomEngine const& engine) noexcept
        {
            if (idcpu[ip] == amrex::ParticleIdCpus::Invalid) {
                return;
            }
            if (amrex::Random(engine) > total_collision_probability) {
                return;
            }

            amrex::ParticleReal x, y, z;
            GetPosition(ip, x, y, z);
            const amrex::ParticleReal background_density =
                amrex::max(n_a_func(x, y, z, t), 0.0_prt);
            const amrex::ParticleReal background_temperature =
                amrex::max(T_a_func(x, y, z, t), 0.0_prt);

            const amrex::ParticleReal velocity_std = std::sqrt(
                PhysConst::kb*background_temperature/background_mass);
            const amrex::ParticleReal target_vx =
                velocity_std*amrex::RandomNormal(0.0_prt, 1.0_prt, engine);
            const amrex::ParticleReal target_vy =
                velocity_std*amrex::RandomNormal(0.0_prt, 1.0_prt, engine);
            const amrex::ParticleReal target_vz =
                velocity_std*amrex::RandomNormal(0.0_prt, 1.0_prt, engine);

            const auto kinematics = getBackgroundMCCKinematics(
                ux[ip], uy[ip], uz[ip],
                target_vx, target_vy, target_vz,
                projectile_mass, background_mass);
            if (kinematics.relative_speed <= 0.0_prt) {
                return;
            }

            const amrex::ParticleReal channel_draw = amrex::Random(engine);
            amrex::ParticleReal cumulative_probability = 0.0_prt;

            for (int i = 0; i < process_count; ++i) {
                auto const& process = processes[i];
                const amrex::ParticleReal cross_section = amrex::max(
                    process.getCrossSection(kinematics.collision_energy_eV),
                    0.0_prt);
                cumulative_probability +=
                    background_density*rate_multipliers[i]*cross_section
                    *kinematics.relative_speed/nu_max;

                if (channel_draw > cumulative_probability) {
                    continue;
                }

                if (process.m_type == ScatteringProcessType::IONIZATION ||
                    process.m_type == ScatteringProcessType::ATTACHMENT)
                {
                    event_process_ptr[ip] = i;
                    return;
                }

                amrex::ParticleReal projectile_ux_out;
                amrex::ParticleReal projectile_uy_out;
                amrex::ParticleReal projectile_uz_out;
                amrex::ParticleReal target_ux_out;
                amrex::ParticleReal target_uy_out;
                amrex::ParticleReal target_uz_out;
                TwoProductComputeProductMomenta(
                    ux[ip], uy[ip], uz[ip], projectile_mass,
                    kinematics.target_ux,
                    kinematics.target_uy,
                    kinematics.target_uz,
                    background_mass,
                    projectile_ux_out,
                    projectile_uy_out,
                    projectile_uz_out,
                    projectile_mass,
                    target_ux_out,
                    target_uy_out,
                    target_uz_out,
                    background_mass,
                    -process.m_energy_penalty*PhysConst::q_e,
                    process.m_scattering_angle_model,
                    engine);

                ux[ip] = projectile_ux_out;
                uy[ip] = projectile_uy_out;
                uz[ip] = projectile_uz_out;
                return;
            }
        });
}
