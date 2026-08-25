/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCCollision.H"

#include "BackgroundMCCUtils.H"
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
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <string>
#include <vector>

namespace
{
    struct CrossSectionKnot
    {
        amrex::ParticleReal energy;
        int process_index;
        int knot_index;
    };

    struct CompareCrossSectionKnots
    {
        bool operator() (CrossSectionKnot const& lhs, CrossSectionKnot const& rhs) const
        {
            if (lhs.energy != rhs.energy) { return lhs.energy > rhs.energy; }
            if (lhs.process_index != rhs.process_index) {
                return lhs.process_index > rhs.process_index;
            }
            return lhs.knot_index > rhs.knot_index;
        }
    };
}

BackgroundMCCCollision::BackgroundMCCCollision (std::string const& collision_name)
    : CollisionBase(collision_name)
{
    using namespace amrex::literals;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_species_names.size() == 1,
        "Background MCC must have exactly one incident species."
    );

    const amrex::ParmParse pp_collision_name(collision_name);

    amrex::ParticleReal background_density = 0;
    if (utils::parser::queryWithParser(
            pp_collision_name, "background_density", background_density))
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            background_density > 0,
            "The background density must be greater than 0."
        );
        m_background_density_parser = utils::parser::makeParser(
            std::to_string(background_density), {"x", "y", "z", "t"});
    }
    else
    {
        std::string background_density_str;
        utils::parser::Store_parserString(
            pp_collision_name, "background_density(x,y,z,t)", background_density_str);
        m_background_density_parser = utils::parser::makeParser(
            background_density_str, {"x", "y", "z", "t"});
    }

    amrex::ParticleReal background_temperature;
    if (utils::parser::queryWithParser(
            pp_collision_name, "background_temperature", background_temperature))
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            background_temperature >= 0,
            "The background temperature must be non-negative."
        );
        m_background_temperature_parser = utils::parser::makeParser(
            std::to_string(background_temperature), {"x", "y", "z", "t"});
    }
    else
    {
        std::string background_temperature_str;
        utils::parser::Store_parserString(
            pp_collision_name,
            "background_temperature(x,y,z,t)",
            background_temperature_str);
        m_background_temperature_parser = utils::parser::makeParser(
            background_temperature_str, {"x", "y", "z", "t"});
    }

    m_background_density_func = m_background_density_parser.compile<4>();
    m_background_temperature_func = m_background_temperature_parser.compile<4>();

    m_user_nu_max = utils::parser::queryWithParser(
        pp_collision_name, "nu_max", m_nu_max);
    if (m_user_nu_max)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(m_nu_max)) && m_nu_max > 0.0_prt,
            "Background MCC nu_max must be finite and greater than 0."
        );
    }

    utils::parser::queryWithParser(
        pp_collision_name, "max_background_density", m_max_background_density);
    if (m_max_background_density == 0 && background_density != 0) {
        m_max_background_density = background_density;
    }
    if (!m_user_nu_max)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_max_background_density > 0,
            "The maximum background density must be greater than 0 when nu_max "
            "is calculated automatically."
        );
    }

    utils::parser::queryWithParser(
        pp_collision_name, "background_mass", m_background_mass);

    auto processes = BinaryCollisionUtils::parse_scattering_processes(collision_name);
    for (auto& process : processes)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            process.type() != ScatteringProcessType::INVALID,
            "Cannot add an unknown scattering process type."
        );

        if (process.type() == ScatteringProcessType::IONIZATION)
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                !ionization_flag,
                "Background MCC currently supports one ionization process per object."
            );
            ionization_flag = true;
            m_ionization_process_index = static_cast<int>(m_processes.size());

            std::string secondary_species;
            pp_collision_name.get("ionization_species", secondary_species);
            m_species_names.push_back(secondary_species);
        }
        m_processes.push_back(std::move(process));
    }

#ifdef AMREX_USE_GPU
    amrex::Gpu::HostVector<ScatteringProcess::Executor> host_processes;
    host_processes.reserve(m_processes.size());
    for (auto const& process : m_processes) {
        host_processes.push_back(process.executor());
    }
    m_processes_exe.resize(host_processes.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice,
        host_processes.begin(),
        host_processes.end(),
        m_processes_exe.begin());
    amrex::Gpu::streamSynchronize();
#else
    for (auto const& process : m_processes) {
        m_processes_exe.push_back(process.executor());
    }
#endif
}

amrex::ParticleReal
BackgroundMCCCollision::get_nu_max (
    amrex::Vector<ScatteringProcess> const& processes) const
{
    using Accumulator = long double;

    if (processes.empty()) { return 0; }

    std::priority_queue<
        CrossSectionKnot,
        std::vector<CrossSectionKnot>,
        CompareCrossSectionKnots> next_knots;

    std::vector<Accumulator> process_slopes(processes.size(), 0.0L);
    Accumulator total_sigma = 0.0L;
    Accumulator total_slope = 0.0L;

    for (int ip = 0; ip < static_cast<int>(processes.size()); ++ip)
    {
        auto const& energies = processes[ip].getEnergyGrid();
        auto const& sigmas = processes[ip].getCrossSectionGrid();

        total_sigma += static_cast<Accumulator>(sigmas[0]);

        int next_knot = 0;
        if (energies[0] == 0.0)
        {
            process_slopes[ip] =
                (static_cast<Accumulator>(sigmas[1]) - sigmas[0]) /
                (static_cast<Accumulator>(energies[1]) - energies[0]);
            total_slope += process_slopes[ip];
            next_knot = 1;
        }
        next_knots.push({energies[next_knot], ip, next_knot});
    }

    auto const collision_speed = [this] (Accumulator const energy)
    {
        if (energy <= 0.0L) { return 0.0L; }

        auto const mass = static_cast<Accumulator>(m_mass1);
        auto const q_e = static_cast<Accumulator>(PhysConst::q_e_v<double>);
        if (m_use_relativistic_electron_kinematics)
        {
            auto const c = static_cast<Accumulator>(PhysConst::c_v<double>);
            auto const c2 = static_cast<Accumulator>(PhysConst::c2_v<double>);
            auto const tau = energy*q_e/(mass*c2);
            return c*std::sqrt(tau*(tau + 2.0L))/(tau + 1.0L);
        }
        return std::sqrt(2.0L*q_e*energy/mass);
    };

    Accumulator left_energy = 0.0L;
    Accumulator left_sigma = total_sigma;
    Accumulator max_sigma_v = 0.0L;

    while (!next_knots.empty())
    {
        auto const right_energy =
            static_cast<Accumulator>(next_knots.top().energy);
        auto right_sigma = left_sigma + total_slope*(right_energy - left_energy);
        right_sigma = std::max(0.0L, right_sigma);

        max_sigma_v = std::max(
            max_sigma_v,
            std::max(left_sigma, right_sigma)*collision_speed(right_energy));

        left_energy = right_energy;
        left_sigma = right_sigma;

        while (!next_knots.empty() &&
               static_cast<Accumulator>(next_knots.top().energy) == right_energy)
        {
            auto const event = next_knots.top();
            next_knots.pop();

            auto const& energies = processes[event.process_index].getEnergyGrid();
            auto const& sigmas = processes[event.process_index].getCrossSectionGrid();

            total_slope -= process_slopes[event.process_index];
            if (event.knot_index + 1 < static_cast<int>(energies.size()))
            {
                auto const j = event.knot_index;
                process_slopes[event.process_index] =
                    (static_cast<Accumulator>(sigmas[j+1]) - sigmas[j]) /
                    (static_cast<Accumulator>(energies[j+1]) - energies[j]);
                next_knots.push({energies[j+1], event.process_index, j+1});
            }
            else
            {
                process_slopes[event.process_index] = 0.0L;
            }
            total_slope += process_slopes[event.process_index];
        }
    }

    // Electron speed is bounded by c, so this also covers the constant
    // high-energy extrapolation used by ScatteringProcess::getCrossSection.
    if (m_use_relativistic_electron_kinematics)
    {
        auto const c = static_cast<Accumulator>(PhysConst::c_v<double>);
        max_sigma_v = std::max(max_sigma_v, left_sigma*c);
    }

    auto nu_max = static_cast<Accumulator>(m_max_background_density)*max_sigma_v;
    if (nu_max <= 0.0L) { return 0; }

    // Round the analytically conservative bound upward after the long-double
    // accumulation. This also absorbs small interpolation roundoff differences.
    nu_max *= 1.0L + 16.0L*static_cast<Accumulator>(
        std::numeric_limits<amrex::ParticleReal>::epsilon());
    auto result = static_cast<amrex::ParticleReal>(nu_max);
    return std::nextafter(
        result, std::numeric_limits<amrex::ParticleReal>::infinity());
}

void
BackgroundMCCCollision::doCollisions (
    amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc)
{
    ABLASTR_PROFILE("BackgroundMCCCollision::doCollisions()");
    using namespace amrex::literals;

    auto& species1 = mypc->GetParticleContainerFromName(m_species_names[0]);
    auto& species2 = (m_species_names.size() == 2)
        ? mypc->GetParticleContainerFromName(m_species_names[1])
        : species1;

    bool const first_call = !init_flag;
    if (!init_flag)
    {
        m_mass1 = species1.getMass();
        m_use_relativistic_electron_kinematics =
            species1.AmIA<PhysicalSpecies::electron>();

        if (m_background_mass < 0.0_prt)
        {
            if (ionization_flag) {
                m_background_mass = species2.getMass() + PhysConst::m_e;
            } else {
                m_background_mass = species1.getMass();
            }
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_background_mass > 0.0_prt,
            "The background neutral mass must be greater than 0."
        );

        if (!m_user_nu_max) {
            m_nu_max = get_nu_max(m_processes);
        }
        init_flag = true;
    }

    auto const coll_n = m_nu_max*dt;
    m_total_collision_prob = static_cast<amrex::ParticleReal>(
        -std::expm1(-static_cast<double>(coll_n)));

    if (coll_n > 0.1_prt && !m_warned_large_dt)
    {
        ablastr::warn_manager::WMRecordWarning(
            "BackgroundMCC Collisions",
            "nu_max*dt = " + std::to_string(coll_n) +
            " is greater than 0.1. Use collision subcycling for converged "
            "one-event-per-substep results."
        );
        m_warned_large_dt = true;
    }

    if (first_call)
    {
        amrex::Print() << Utils::TextMsg::Info(
            "Setting up Monte-Carlo collisions for " + m_species_names[0] + " with:\n"
            + "     nu_max: " + std::to_string(m_nu_max)
            + (m_user_nu_max ? " (user supplied)" : " (automatic)")
            + "\n     total collision probability: "
            + std::to_string(m_total_collision_prob)
        );
    }

    if (m_processes.empty() || m_total_collision_prob <= 0.0_prt) { return; }

    auto const finest_level = species1.finestLevel();

    if (!ionization_flag)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            auto* cost = WarpX::getCosts(lev);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti)
            {
                if (cost && WarpX::load_balance_costs_update_algo ==
                            LoadBalanceCostsUpdateAlgo::Timers)
                {
                    amrex::Gpu::synchronize();
                }
                auto wt = static_cast<amrex::Real>(amrex::second());

                doBackgroundCollisionsWithinTile(pti, cur_time, nullptr);

                if (cost && WarpX::load_balance_costs_update_algo ==
                            LoadBalanceCostsUpdateAlgo::Timers)
                {
                    amrex::Gpu::synchronize();
                    wt = static_cast<amrex::Real>(amrex::second()) - wt;
                    amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
                }
            }
        }
        return;
    }

    const SmartCopyFactory copy_factory_elec(species1, species1);
    const SmartCopyFactory copy_factory_ion(species1, species2);
    const auto CopyElec = copy_factory_elec.getSmartCopy();
    const auto CopyIon = copy_factory_ion.getSmartCopy();

    const amrex::ParticleReal sqrt_kb_m =
        std::sqrt(PhysConst::kb/m_background_mass);
    const auto Transform = ImpactIonizationTransformFunc(
        m_processes[m_ionization_process_index].getEnergyPenalty(),
        m_mass1,
        sqrt_kb_m,
        m_background_temperature_func,
        cur_time);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto* cost = WarpX::getCosts(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(species1, lev); pti.isValid(); ++pti)
        {
            if (cost && WarpX::load_balance_costs_update_algo ==
                        LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            auto wt = static_cast<amrex::Real>(amrex::second());

            auto& elec_tile = species1.ParticlesAt(lev, pti);
            auto& ion_tile = species2.ParticlesAt(lev, pti);
            const amrex::Long np_elec = elec_tile.numParticles();
            const amrex::Long np_ion = ion_tile.numParticles();

            amrex::Gpu::DeviceVector<amrex::Long> ionization_mask(np_elec);
            doBackgroundCollisionsWithinTile(
                pti, cur_time, ionization_mask.dataPtr());

            const auto num_added = filterCopyTransformParticles<1>(
                species1,
                species2,
                elec_tile,
                ion_tile,
                elec_tile,
                ionization_mask.dataPtr(),
                np_elec,
                np_ion,
                CopyElec,
                CopyIon,
                Transform);

            setNewParticleIDs(elec_tile, np_elec, num_added);
            setNewParticleIDs(ion_tile, np_ion, num_added);

            if (cost && WarpX::load_balance_costs_update_algo ==
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
BackgroundMCCCollision::doBackgroundCollisionsWithinTile (
    WarpXParIter& pti, amrex::Real t, amrex::Long* ionization_mask)
{
    using namespace amrex::literals;
    using std::sqrt;

    const long np = pti.numParticles();
    auto n_a_func = m_background_density_func;
    auto T_a_func = m_background_temperature_func;

    auto* processes = m_processes_exe.data();
    auto const process_count = static_cast<int>(m_processes_exe.size());
    auto const ionization_process_index = m_ionization_process_index;
    auto const total_collision_prob = m_total_collision_prob;
    auto const nu_max = m_nu_max;
    auto const use_relativistic_electron_kinematics =
        m_use_relativistic_electron_kinematics;

    auto const m = m_mass1;
    auto const M = m_background_mass;
    auto GetPosition = GetParticlePosition<PIdx>(pti);

    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    amrex::ParallelForRNG(
        np,
        [=] AMREX_GPU_HOST_DEVICE (
            long ip, amrex::RandomEngine const& engine)
        {
            if (ionization_mask != nullptr) { ionization_mask[ip] = 0; }
            if (amrex::Random(engine) > total_collision_prob) { return; }

            amrex::ParticleReal x, y, z;
            GetPosition(ip, x, y, z);

            const amrex::ParticleReal n_a = n_a_func(x, y, z, t);
            const amrex::ParticleReal T_a = T_a_func(x, y, z, t);

            const auto vel_std = sqrt(PhysConst::kb*T_a/M);
            const amrex::ParticleReal ua_x =
                vel_std*amrex::RandomNormal(0_prt, 1.0_prt, engine);
            const amrex::ParticleReal ua_y =
                vel_std*amrex::RandomNormal(0_prt, 1.0_prt, engine);
            const amrex::ParticleReal ua_z =
                vel_std*amrex::RandomNormal(0_prt, 1.0_prt, engine);

            amrex::ParticleReal E_coll;
            amrex::ParticleReal v_coll;
            if (use_relativistic_electron_kinematics)
            {
                BackgroundMCCUtils::getElectronNeutralCollisionParameters(
                    ux[ip], uy[ip], uz[ip], ua_x, ua_y, ua_z, m, E_coll, v_coll);
            }
            else
            {
                const amrex::ParticleReal vx = ux[ip] - ua_x;
                const amrex::ParticleReal vy = uy[ip] - ua_y;
                const amrex::ParticleReal vz = uz[ip] - ua_z;
                const amrex::ParticleReal v_coll2 = vx*vx + vy*vy + vz*vz;
                double gamma;
                double E_coll_double;
                v_coll = sqrt(v_coll2);
                ParticleUtils::getCollisionEnergy(
                    v_coll2, m, M, gamma, E_coll_double);
                E_coll = static_cast<amrex::ParticleReal>(E_coll_double);
            }

            const amrex::ParticleReal process_draw = amrex::Random(engine);
            amrex::ParticleReal cumulative_probability = 0.0_prt;

            for (int i = 0; i < process_count; ++i)
            {
                auto const& process = processes[i];
                const auto sigma = process.getCrossSection(E_coll);
                cumulative_probability += n_a*sigma*v_coll/nu_max;
                if (process_draw > cumulative_probability) { continue; }

                if (i == ionization_process_index)
                {
                    ionization_mask[ip] = 1;
                    break;
                }

                amrex::ParticleReal u1x_out, u1y_out, u1z_out;
                amrex::ParticleReal u2x_out, u2y_out, u2z_out;
                TwoProductComputeProductMomenta(
                    ux[ip], uy[ip], uz[ip], m,
                    ua_x, ua_y, ua_z, M,
                    u1x_out, u1y_out, u1z_out, m,
                    u2x_out, u2y_out, u2z_out, M,
                    -process.m_energy_penalty*PhysConst::q_e,
                    process.m_scattering_angle_model,
                    engine);

                ux[ip] = u1x_out;
                uy[ip] = u1y_out;
                uz[ip] = u1z_out;
                break;
            }
        });
}
