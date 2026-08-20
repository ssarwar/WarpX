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

BackgroundMCCCollision::BackgroundMCCCollision (std::string const& collision_name)
    : CollisionBase(collision_name)
{
    using namespace amrex::literals;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_species_names.size() == 1,
        "Background MCC must have exactly one incident species.");

    const amrex::ParmParse pp_collision_name(collision_name);

    amrex::ParticleReal background_density = 0.0_prt;
    if (utils::parser::queryWithParser(
            pp_collision_name, "background_density", background_density))
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            background_density > 0.0_prt,
            "The background density must be greater than zero.");
        m_background_density_parser = utils::parser::makeParser(
            std::to_string(background_density), {"x", "y", "z", "t"});
    }
    else
    {
        std::string background_density_str;
        utils::parser::Store_parserString(
            pp_collision_name, "background_density(x,y,z,t)",
            background_density_str);
        m_background_density_parser = utils::parser::makeParser(
            background_density_str, {"x", "y", "z", "t"});
    }

    amrex::ParticleReal background_temperature = -1.0_prt;
    if (utils::parser::queryWithParser(
            pp_collision_name, "background_temperature", background_temperature))
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            background_temperature >= 0.0_prt,
            "The background temperature must be non-negative.");
        m_background_temperature_parser = utils::parser::makeParser(
            std::to_string(background_temperature), {"x", "y", "z", "t"});
    }
    else
    {
        std::string background_temperature_str;
        utils::parser::Store_parserString(
            pp_collision_name, "background_temperature(x,y,z,t)",
            background_temperature_str);
        m_background_temperature_parser = utils::parser::makeParser(
            background_temperature_str, {"x", "y", "z", "t"});
    }

    m_background_density_func = m_background_density_parser.compile<4>();
    m_background_temperature_func = m_background_temperature_parser.compile<4>();

    utils::parser::queryWithParser(
        pp_collision_name, "max_background_density", m_max_background_density);
    if (m_max_background_density == 0.0_prt && background_density > 0.0_prt) {
        m_max_background_density = background_density;
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_max_background_density > 0.0_prt,
        "The maximum background density must be greater than zero.");

    utils::parser::queryWithParser(
        pp_collision_name, "background_mass", m_background_mass);
    utils::parser::queryWithParser(
        pp_collision_name, "nu_max_safety_factor", m_nu_max_safety_factor);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_nu_max_safety_factor >= 1.0_prt,
        "The background-MCC nu_max_safety_factor must be at least one.");

    pp_collision_name.queryarr("scattering_processes", m_process_names);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_process_names.empty(),
        "Background MCC requires at least one scattering process.");

    m_processes = BinaryCollisionUtils::parse_scattering_processes(collision_name);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_processes.size() == m_process_names.size(),
        "Internal mismatch while parsing background-MCC processes.");

    m_product_species_names.resize(m_processes.size());
    m_rate_multipliers_h.resize(m_processes.size(), 1.0_prt);

    int ionization_process_count = 0;
    for (int i = 0; i < static_cast<int>(m_processes.size()); ++i) {
        const auto process_type = m_processes[i].type();
        const auto& process_name = m_process_names[i];

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            process_type != ScatteringProcessType::INVALID,
            "Cannot add an unknown background-MCC process type: " + process_name);

        if (process_type == ScatteringProcessType::IONIZATION ||
            process_type == ScatteringProcessType::ATTACHMENT)
        {
            const std::string species_keyword = process_name + "_species";
            pp_collision_name.get(
                species_keyword.c_str(), m_product_species_names[i]);
        }

        if (process_type == ScatteringProcessType::IONIZATION) {
            ++ionization_process_count;
        }

        bool zero_below_min = false;
        bool zero_above_max = false;
        const std::string zero_below_keyword =
            process_name + "_zero_below_min";
        const std::string zero_above_keyword =
            process_name + "_zero_above_max";
        pp_collision_name.query(
            zero_below_keyword.c_str(), zero_below_min);
        pp_collision_name.query(
            zero_above_keyword.c_str(), zero_above_max);
        m_processes[i].setCrossSectionExtrapolation(
            zero_below_min, zero_above_max);

        amrex::ParticleReal third_body_density = 0.0_prt;
        const std::string third_body_density_keyword =
            process_name + "_third_body_density";
        if (utils::parser::queryWithParser(
                pp_collision_name,
                third_body_density_keyword.c_str(),
                third_body_density))
        {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                process_type == ScatteringProcessType::ATTACHMENT,
                "third_body_density is currently supported only for attachment channels.");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                third_body_density > 0.0_prt,
                "The third-body density must be greater than zero.");
            m_rate_multipliers_h[i] = third_body_density;
        }
    }

    // A single ionization channel preserves the existing one-secondary-electron
    // SmartCopy implementation. Different neutral components should use separate
    // MCC objects, e.g. one object for N2 and one for O2.
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        ionization_process_count <= 1,
        "Background MCC currently supports at most one ionization channel per object.");

#ifdef AMREX_USE_GPU
    amrex::Gpu::HostVector<ScatteringProcess::Executor> h_processes_exe;
    h_processes_exe.reserve(m_processes.size());
    for (auto const& process : m_processes) {
        h_processes_exe.push_back(process.executor());
    }

    m_processes_exe.resize(h_processes_exe.size());
    m_rate_multipliers_d.resize(m_rate_multipliers_h.size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice,
        h_processes_exe.begin(), h_processes_exe.end(),
        m_processes_exe.begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice,
        m_rate_multipliers_h.begin(), m_rate_multipliers_h.end(),
        m_rate_multipliers_d.begin());
    amrex::Gpu::streamSynchronize();
#else
    for (auto const& process : m_processes) {
        m_processes_exe.push_back(process.executor());
    }
    for (auto const multiplier : m_rate_multipliers_h) {
        m_rate_multipliers_d.push_back(multiplier);
    }
#endif
}
