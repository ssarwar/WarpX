/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCCollision.H"

#include "BackgroundMCCElasticKinematics.H"
#include "BackgroundMCCProducts.H"
#include "BackgroundMCCUtils.H"
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
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
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
            std::isfinite(static_cast<double>(background_density)) &&
                background_density > 0,
            "The background density must be finite and greater than 0."
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
            std::isfinite(static_cast<double>(background_temperature)) &&
                background_temperature >= 0,
            "The background temperature must be finite and non-negative."
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

    auto const has_max_background_density = utils::parser::queryWithParser(
        pp_collision_name, "max_background_density", m_max_background_density);
    if (has_max_background_density)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(m_max_background_density)) &&
                m_max_background_density > 0.0_prt,
            "The maximum background density must be finite and greater than 0."
        );
    }
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
    amrex::Vector<int> process_product_group;
    process_product_group.reserve(processes.size());
    amrex::Vector<IonizationEnergySharingModel> ionization_energy_models;
    amrex::Vector<BackgroundMCCIonizationTarget> ionization_targets;
    amrex::Vector<std::string> differential_cross_sections;
    ionization_energy_models.reserve(processes.size());
    ionization_targets.reserve(processes.size());
    differential_cross_sections.reserve(processes.size());
    amrex::ParticleReal n2_maximum_energy = 0.0_prt;
    amrex::ParticleReal o2_maximum_energy = 0.0_prt;

    for (auto& process : processes)
    {
        auto const process_type = process.type();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            process_type != ScatteringProcessType::INVALID,
            "Cannot add an unknown scattering process type."
        );

        auto energy_sharing_model = IonizationEnergySharingModel::Equal;
        auto ionization_target = BackgroundMCCIonizationTarget::None;
        if (process_type == ScatteringProcessType::IONIZATION)
        {
            pp_collision_name.query_enum_case_insensitive(
                process.name() + "_energy_sharing_model", energy_sharing_model);

            if (energy_sharing_model == IonizationEnergySharingModel::RBEQ)
            {
                std::string target_name;
                pp_collision_name.get(
                    process.name() + "_rbeq_target", target_name);
                ionization_target =
                    BackgroundMCCIonizationModel::parseTarget(target_name);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    ionization_target != BackgroundMCCIonizationTarget::None,
                    "RBEQ ionization target must be N2 or O2."
                );

                auto const outer_binding =
                    BackgroundMCCIonizationModel::outerBindingEnergy(
                        ionization_target);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::abs(process.getEnergyPenalty() - outer_binding) <= 0.05_prt,
                    "RBEQ ionization energy must match the target's outer-shell "
                    "binding energy to within 0.05 eV."
                );

                if (ionization_target == BackgroundMCCIonizationTarget::N2)
                {
                    n2_maximum_energy = std::max(
                        n2_maximum_energy, process.getMaxEnergyInput());
                }
                else
                {
                    o2_maximum_energy = std::max(
                        o2_maximum_energy, process.getMaxEnergyInput());
                }
            }
        }

        std::string differential_cross_section;
        auto const has_differential_cross_section = pp_collision_name.query(
            process.name() + "_differential_cross_section", differential_cross_section);
        if (process.scatteringAngleModel() == ScatteringAngleModel::IAA)
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                process_type == ScatteringProcessType::ELASTIC ||
                    process_type == ScatteringProcessType::EXCITATION ||
                    process_type == ScatteringProcessType::IONIZATION,
                "The IAA scattering-angle model applies only to Background MCC "
                "elastic, excitation and ionization processes."
            );
            if (process_type == ScatteringProcessType::ELASTIC ||
                process_type == ScatteringProcessType::EXCITATION)
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    has_differential_cross_section,
                    "IAA elastic or excitation scattering requires a "
                    "<process>_differential_cross_section file."
                );
                m_has_iaa_differential_processes = true;
            }
            else
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    !has_differential_cross_section,
                    "A differential-cross-section file is not valid for "
                    "IAA ionization scattering."
                );
            }
        }
        else
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                !has_differential_cross_section,
                "<process>_differential_cross_section requires "
                "<process>_scattering_angle_model = IAA."
            );
        }

        if (process_type == ScatteringProcessType::ATTACHMENT)
        {
            auto const units_key = process.name() + "_cross_section_units";
            std::string cross_section_units;
            auto const has_cross_section_units =
                pp_collision_name.query(units_key, cross_section_units);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                has_cross_section_units,
                "Every attachment process must specify <process>_cross_section_units "
                "as either m2 or m5."
            );
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                cross_section_units == "m2" || cross_section_units == "m5",
                "Attachment cross_section_units must be either m2 or m5."
            );

            amrex::ParticleReal third_body_density = 0.0_prt;
            auto const third_body_density_key =
                process.name() + "_third_body_density";
            auto const has_third_body_density = utils::parser::queryWithParser(
                pp_collision_name,
                third_body_density_key.c_str(),
                third_body_density);

            if (cross_section_units == "m5")
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    has_third_body_density,
                    "Attachment cross sections in m5 require a positive "
                    "<process>_third_body_density in m^-3."
                );
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(static_cast<double>(third_body_density)) &&
                    third_body_density > 0.0_prt,
                    "Attachment third_body_density must be finite and greater than 0."
                );
                process.setCrossSectionMultiplier(
                    static_cast<double>(third_body_density));
            }
            else
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    !has_third_body_density,
                    "Attachment third_body_density is only valid when "
                    "cross_section_units = m5."
                );
            }
        }

        int product_group = -1;
        if (process_type == ScatteringProcessType::IONIZATION ||
            process_type == ScatteringProcessType::ATTACHMENT)
        {
            m_has_ionization_processes = m_has_ionization_processes ||
                process_type == ScatteringProcessType::IONIZATION;
            m_has_attachment_processes = m_has_attachment_processes ||
                process_type == ScatteringProcessType::ATTACHMENT;

            std::string product_species;
            pp_collision_name.get(process.name() + "_species", product_species);

            for (int i = 0; i < static_cast<int>(m_product_groups.size()); ++i)
            {
                if (m_product_groups[i].type == process_type &&
                    m_product_groups[i].species_name == product_species)
                {
                    product_group = i;
                    break;
                }
            }
            if (product_group < 0)
            {
                product_group = static_cast<int>(m_product_groups.size());
                m_product_groups.push_back({process_type, product_species});
            }

            if (std::find(
                    m_species_names.begin(), m_species_names.end(), product_species) ==
                m_species_names.end())
            {
                m_species_names.push_back(product_species);
            }
        }

        process_product_group.push_back(product_group);
        ionization_energy_models.push_back(energy_sharing_model);
        ionization_targets.push_back(ionization_target);
        differential_cross_sections.push_back(
            std::move(differential_cross_section));
        m_processes.push_back(std::move(process));
    }

    m_process_selector =
        std::make_unique<BackgroundMCCProcessSelector>(m_processes);

    amrex::Gpu::HostVector<BackgroundMCCElasticScatteringModel::Executor>
        host_differential_scattering_processes;
    host_differential_scattering_processes.reserve(m_processes.size());
    m_differential_scattering_models.reserve(m_processes.size());
    std::unordered_map<std::string, std::size_t> differential_model_indices;
    differential_model_indices.reserve(m_processes.size());
    for (int i = 0; i < static_cast<int>(m_processes.size()); ++i)
    {
        BackgroundMCCElasticScatteringModel::Executor executor;
        if (!differential_cross_sections[i].empty())
        {
            auto const& file_name = differential_cross_sections[i];
            auto model = differential_model_indices.find(file_name);
            if (model == differential_model_indices.end())
            {
                auto const model_index = m_differential_scattering_models.size();
                m_differential_scattering_models.push_back(
                    std::make_unique<BackgroundMCCElasticScatteringModel>(file_name));
                model = differential_model_indices.emplace(file_name, model_index).first;
            }
            executor = m_differential_scattering_models[model->second]->executor();
        }
        host_differential_scattering_processes.push_back(executor);
    }

    if (n2_maximum_energy > 0.0_prt)
    {
        m_n2_ionization_model = std::make_unique<BackgroundMCCIonizationModel>(
            BackgroundMCCIonizationTarget::N2, n2_maximum_energy);
    }
    if (o2_maximum_energy > 0.0_prt)
    {
        m_o2_ionization_model = std::make_unique<BackgroundMCCIonizationModel>(
            BackgroundMCCIonizationTarget::O2, o2_maximum_energy);
    }

    amrex::Gpu::HostVector<BackgroundMCCIonizationModel::Executor>
        host_ionization_processes;
    host_ionization_processes.reserve(m_processes.size());
    for (int i = 0; i < static_cast<int>(m_processes.size()); ++i)
    {
        BackgroundMCCIonizationModel::Executor executor;
        executor.m_model = ionization_energy_models[i];
        if (ionization_targets[i] == BackgroundMCCIonizationTarget::N2)
        {
            executor = m_n2_ionization_model->executor();
        }
        else if (ionization_targets[i] == BackgroundMCCIonizationTarget::O2)
        {
            executor = m_o2_ionization_model->executor();
        }
        host_ionization_processes.push_back(executor);
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

    m_ionization_processes_exe.resize(host_ionization_processes.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice,
        host_ionization_processes.begin(),
        host_ionization_processes.end(),
        m_ionization_processes_exe.begin());

    m_differential_scattering_processes_exe.resize(
        host_differential_scattering_processes.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice,
        host_differential_scattering_processes.begin(),
        host_differential_scattering_processes.end(),
        m_differential_scattering_processes_exe.begin());

    m_process_product_group.resize(process_product_group.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice,
        process_product_group.begin(),
        process_product_group.end(),
        m_process_product_group.begin());
    amrex::Gpu::streamSynchronize();
#else
    for (auto const& process : m_processes) {
        m_processes_exe.push_back(process.executor());
    }
    for (auto const& ionization_process : host_ionization_processes) {
        m_ionization_processes_exe.push_back(ionization_process);
    }
    for (auto const& differential_scattering_process :
         host_differential_scattering_processes) {
        m_differential_scattering_processes_exe.push_back(
            differential_scattering_process);
    }
    for (auto const product_group : process_product_group) {
        m_process_product_group.push_back(product_group);
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

    bool const first_call = !init_flag;
    if (!init_flag)
    {
        m_mass1 = species1.getMass();
        m_use_relativistic_electron_kinematics =
            species1.AmIA<PhysicalSpecies::electron>();

        if (!m_product_groups.empty())
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_use_relativistic_electron_kinematics,
                "Background MCC ionization and attachment require an electron "
                "incident species."
            );

            amrex::ParticleReal inferred_background_mass = -1.0_prt;
            auto const charge_tolerance = 100.0_prt*
                std::numeric_limits<amrex::ParticleReal>::epsilon()*PhysConst::q_e;
            auto const mass_tolerance = 100.0_prt*
                std::numeric_limits<amrex::ParticleReal>::epsilon();

            for (auto const& product_group : m_product_groups)
            {
                auto& product = mypc->GetParticleContainerFromName(
                    product_group.species_name);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    product_group.species_name != m_species_names[0],
                    "Background MCC product species must differ from the incident "
                    "electron species."
                );

                if (product_group.type == ScatteringProcessType::IONIZATION)
                {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        std::abs(product.getCharge() - PhysConst::q_e) <=
                            charge_tolerance,
                        "Background MCC ionization product species must have charge +q_e."
                    );

                    auto const candidate_mass = product.getMass() + PhysConst::m_e;
                    if (inferred_background_mass < 0.0_prt)
                    {
                        inferred_background_mass = candidate_mass;
                    }
                    else
                    {
                        auto const mass_scale = std::max(
                            std::abs(inferred_background_mass),
                            std::abs(candidate_mass));
                        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                            std::abs(candidate_mass - inferred_background_mass) <=
                                mass_tolerance*mass_scale,
                            "Ionization product masses imply different neutral masses. "
                            "Specify background_mass explicitly for dissociative or "
                            "multi-target channels."
                        );
                    }
                }
                else
                {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        std::abs(product.getCharge() + PhysConst::q_e) <=
                            charge_tolerance,
                        "Background MCC attachment product species must have charge -q_e."
                    );
                }
            }

            if (m_background_mass < 0.0_prt)
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    inferred_background_mass > 0.0_prt,
                    "background_mass is required for attachment-only Background MCC "
                    "blocks because dissociative product mass does not determine the "
                    "target-neutral mass."
                );
                m_background_mass = inferred_background_mass;
            }
        }
        else if (m_background_mass < 0.0_prt)
        {
            m_background_mass = species1.getMass();
        }

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(m_background_mass)) &&
                m_background_mass > 0.0_prt,
            "The background neutral mass must be finite and greater than 0."
        );
        if (m_has_iaa_differential_processes)
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_use_relativistic_electron_kinematics,
                "IAA differential scattering requires an electron incident species."
            );
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_background_mass > m_mass1,
                "IAA differential scattering requires a neutral target heavier than "
                "the incident electron."
            );
        }

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
            + "\n     product groups: " + std::to_string(m_product_groups.size())
        );
    }

    if (m_processes.empty() || m_total_collision_prob <= 0.0_prt) { return; }

    auto const finest_level = species1.finestLevel();

    if (m_product_groups.empty())
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

                doBackgroundCollisionsWithinTile(
                    pti, cur_time, nullptr, nullptr, nullptr, nullptr);

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
    const auto CopyElec = copy_factory_elec.getSmartCopy();

    std::vector<WarpXParticleContainer*> product_species;
    std::vector<std::unique_ptr<SmartCopyFactory>> product_copy_factories;
    std::vector<SmartCopy> product_copies;
    product_species.reserve(m_product_groups.size());
    product_copy_factories.reserve(m_product_groups.size());
    product_copies.reserve(m_product_groups.size());

    for (auto const& product_group : m_product_groups)
    {
        auto& product = mypc->GetParticleContainerFromName(product_group.species_name);
        product_species.push_back(&product);
        product_copy_factories.push_back(
            std::make_unique<SmartCopyFactory>(species1, product));
        product_copies.push_back(product_copy_factories.back()->getSmartCopy());
    }

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

            auto& source_tile = species1.ParticlesAt(lev, pti);
            amrex::Long const np_source = source_tile.numParticles();
            if (np_source > 0)
            {
                amrex::Long const metadata_size = m_has_ionization_processes
                    ? 2*np_source : np_source;
                amrex::Gpu::DeviceVector<int> selected_process(metadata_size);
                amrex::Gpu::DeviceVector<amrex::ParticleReal> neutral_vx(np_source);
                amrex::Gpu::DeviceVector<amrex::ParticleReal> neutral_vy(np_source);
                amrex::Gpu::DeviceVector<amrex::ParticleReal> neutral_vz(np_source);

                if (metadata_size > np_source)
                {
                    int* const selected = selected_process.dataPtr();
                    amrex::ParallelFor(
                        np_source,
                        [=] AMREX_GPU_HOST_DEVICE (amrex::Long i) noexcept
                        {
                            selected[np_source + i] = -1;
                        });
                }

                doBackgroundCollisionsWithinTile(
                    pti,
                    cur_time,
                    selected_process.dataPtr(),
                    neutral_vx.dataPtr(),
                    neutral_vy.dataPtr(),
                    neutral_vz.dataPtr());

                ABLASTR_PROFILE_VAR(
                    "BackgroundMCCCollision::createProducts()", prof_create_products);
                for (int group_index = 0;
                     group_index < static_cast<int>(m_product_groups.size());
                     ++group_index)
                {
                    auto& product = *product_species[group_index];
                    auto& product_tile = product.ParticlesAt(lev, pti);
                    amrex::Long const old_product_count =
                        product_tile.numParticles();

                    const auto Filter = BackgroundMCCProductFilterFunc(
                        group_index,
                        selected_process.dataPtr(),
                        m_process_product_group.dataPtr());

                    amrex::Long num_added = 0;
                    if (m_product_groups[group_index].type ==
                        ScatteringProcessType::ATTACHMENT)
                    {
                        const auto Transform = AttachmentTransformFunc(
                            neutral_vx.dataPtr(),
                            neutral_vy.dataPtr(),
                            neutral_vz.dataPtr());
                        num_added = filterCopyTransformParticles<1>(
                            product,
                            product_tile,
                            source_tile,
                            old_product_count,
                            Filter,
                            product_copies[group_index],
                            Transform);
                        setNewParticleIDs(
                            product_tile, old_product_count, num_added);
                    }
                    else
                    {
                        amrex::Long const old_electron_count =
                            source_tile.numParticles();
                        const auto Transform = ImpactIonizationTransformFunc(
                            m_processes_exe.data(),
                            m_ionization_processes_exe.data(),
                            selected_process.dataPtr(),
                            neutral_vx.dataPtr(),
                            neutral_vy.dataPtr(),
                            neutral_vz.dataPtr(),
                            m_mass1,
                            product.getMass());
                        num_added = filterCopyTransformParticles<1>(
                            species1,
                            product,
                            source_tile,
                            product_tile,
                            source_tile,
                            old_electron_count,
                            old_product_count,
                            Filter,
                            CopyElec,
                            product_copies[group_index],
                            Transform);
                        setNewParticleIDs(
                            source_tile, old_electron_count, num_added);
                        setNewParticleIDs(
                            product_tile, old_product_count, num_added);

                        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                            source_tile.numParticles() <= metadata_size,
                            "Background MCC created more than one electron per "
                            "original incident particle."
                        );
                    }
                }
                ABLASTR_PROFILE_VAR_STOP(prof_create_products);
            }

            if (cost && WarpX::load_balance_costs_update_algo ==
                        LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = static_cast<amrex::Real>(amrex::second()) - wt;
                amrex::HostDevice::Atomic::Add(&(*cost)[pti.index()], wt);
            }
        }
    }

    if (m_has_attachment_processes)
    {
        ABLASTR_PROFILE("BackgroundMCCCollision::compactAttachedElectrons()");
        species1.deleteInvalidParticles();
    }
}

void
BackgroundMCCCollision::doBackgroundCollisionsWithinTile (
    WarpXParIter& pti,
    amrex::Real t,
    int* selected_process,
    amrex::ParticleReal* neutral_vx,
    amrex::ParticleReal* neutral_vy,
    amrex::ParticleReal* neutral_vz)
{
    ABLASTR_PROFILE("BackgroundMCCCollision::selectAndScatter()");
    using namespace amrex::literals;
    using std::sqrt;

    const long np = pti.numParticles();
    auto n_a_func = m_background_density_func;
    auto T_a_func = m_background_temperature_func;

    auto* processes = m_processes_exe.data();
    auto const process_count = static_cast<int>(m_processes_exe.size());
    auto const process_selector = m_process_selector->executor();
    auto const* differential_scattering_processes =
        m_differential_scattering_processes_exe.data();
    auto const* process_product_group = m_process_product_group.data();
    auto const total_collision_prob = m_total_collision_prob;
    auto const nu_max = m_nu_max;
    auto const user_nu_max = m_user_nu_max;
    auto const max_background_density = m_max_background_density;
    auto const use_relativistic_electron_kinematics =
        m_use_relativistic_electron_kinematics;

    auto const m = m_mass1;
    auto const M = m_background_mass;
    auto GetPosition = GetParticlePosition<PIdx>(pti);

    auto& attribs = pti.GetAttribs();
    amrex::ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    amrex::ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

    auto const tolerance = 64.0_prt*
        std::numeric_limits<amrex::ParticleReal>::epsilon();

    amrex::ParallelForRNG(
        np,
        [=] AMREX_GPU_HOST_DEVICE (
            long ip, amrex::RandomEngine const& engine)
        {
            if (selected_process != nullptr) { selected_process[ip] = -1; }
            if (amrex::Random(engine) > total_collision_prob) { return; }

            amrex::ParticleReal x, y, z;
            GetPosition(ip, x, y, z);

            const amrex::ParticleReal n_a = n_a_func(x, y, z, t);
            const amrex::ParticleReal T_a = T_a_func(x, y, z, t);

            bool const valid_density = n_a >= 0.0_prt &&
                n_a <= std::numeric_limits<amrex::ParticleReal>::max() &&
                (user_nu_max ||
                 n_a <= max_background_density*(1.0_prt + tolerance));
            bool const valid_temperature = T_a >= 0.0_prt &&
                T_a <= std::numeric_limits<amrex::ParticleReal>::max();
            AMREX_IF_ON_DEVICE((
                AMREX_DEVICE_ASSERT(valid_density);
                AMREX_DEVICE_ASSERT(valid_temperature);
            ))
            AMREX_IF_ON_HOST((
                if (!valid_density) {
                    amrex::Abort(
                        "Background MCC density is negative or exceeds "
                        "max_background_density.");
                }
                if (!valid_temperature) {
                    amrex::Abort("Background MCC temperature is negative.");
                }
            ))
            if (n_a == 0.0_prt) { return; }

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
            if (v_coll <= 0.0_prt || nu_max <= 0.0_prt) { return; }

            const amrex::ParticleReal process_draw = amrex::Random(engine);
            amrex::ParticleReal cumulative_probability = 0.0_prt;
            int chosen_process = -1;

            if (process_selector.enabled())
            {
                auto const probability_per_cross_section = n_a * v_coll / nu_max;
                auto const cross_section_draw =
                    process_draw / probability_per_cross_section;
                amrex::ParticleReal total_cross_section;
                process_selector.select(
                    E_coll, cross_section_draw, chosen_process, total_cross_section);
                cumulative_probability =
                    probability_per_cross_section * total_cross_section;
            }
            else
            {
                for (int i = 0; i < process_count; ++i)
                {
                    auto const& process = processes[i];
                    const auto sigma = process.getCrossSection(E_coll);
                    cumulative_probability += n_a*sigma*v_coll/nu_max;
                    if (chosen_process < 0 &&
                        process_draw <= cumulative_probability)
                    {
                        chosen_process = i;
                        if (!user_nu_max) { break; }
                    }
                }
            }

            if (user_nu_max)
            {
                bool const valid_majorant =
                    cumulative_probability <= 1.0_prt + tolerance;
                AMREX_IF_ON_DEVICE((
                    AMREX_DEVICE_ASSERT(valid_majorant);
                ))
                AMREX_IF_ON_HOST((
                    if (!valid_majorant) {
                        amrex::Abort(
                            "User-specified Background MCC nu_max is smaller "
                            "than the local total collision frequency.");
                    }
                ))
            }

            if (chosen_process < 0) { return; }

            if (process_product_group[chosen_process] >= 0)
            {
                selected_process[ip] = chosen_process;
                neutral_vx[ip] = ua_x;
                neutral_vy[ip] = ua_y;
                neutral_vz[ip] = ua_z;
                return;
            }

            auto const& process = processes[chosen_process];
            amrex::ParticleReal u1x_out, u1y_out, u1z_out;
            amrex::ParticleReal u2x_out, u2y_out, u2z_out;
            if (process.m_scattering_angle_model == ScatteringAngleModel::IAA)
            {
                auto const cosine =
                    differential_scattering_processes[chosen_process].sampleCosine(
                        E_coll, amrex::Random(engine));
                if (process.m_type == ScatteringProcessType::ELASTIC)
                {
                    BackgroundMCCElasticKinematics::compute(
                        ux[ip], uy[ip], uz[ip], ua_x, ua_y, ua_z, m, M, cosine, engine,
                        u1x_out, u1y_out, u1z_out, u2x_out, u2y_out, u2z_out);
                }
                else
                {
                    BackgroundMCCElasticKinematics::computeExcitation(
                        ux[ip], uy[ip], uz[ip], ua_x, ua_y, ua_z, m, M,
                        process.m_energy_penalty, cosine, engine,
                        u1x_out, u1y_out, u1z_out, u2x_out, u2y_out, u2z_out);
                }
            }
            else
            {
                TwoProductComputeProductMomenta(
                    ux[ip], uy[ip], uz[ip], m,
                    ua_x, ua_y, ua_z, M,
                    u1x_out, u1y_out, u1z_out, m,
                    u2x_out, u2y_out, u2z_out, M,
                    -process.m_energy_penalty*PhysConst::q_e,
                    process.m_scattering_angle_model,
                    engine);
            }

            ux[ip] = u1x_out;
            uy[ip] = u1y_out;
            uz[ip] = u1z_out;
        });
}
