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
BackgroundMCCCollision::initialize (MultiParticleContainer* mypc)
{
    using namespace amrex::literals;

    auto& source_species =
        mypc->GetParticleContainerFromName(m_species_names[0]);
    m_projectile_mass = source_species.getMass();

    bool has_attachment = false;
    if (m_background_mass <= 0.0) {
        for (int i = 0; i < static_cast<int>(m_processes.size()); ++i) {
            if (m_processes[i].type() == ScatteringProcessType::IONIZATION) {
                auto& product_species = mypc->GetParticleContainerFromName(
                    m_product_species_names[i]);
                m_background_mass = product_species.getMass();
                break;
            }
            if (m_processes[i].type() == ScatteringProcessType::ATTACHMENT) {
                has_attachment = true;
            }
        }
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !(has_attachment && m_background_mass <= 0.0),
        "background_mass must be specified for an MCC object containing attachment. "
        "For dissociative attachment, the negative-ion product mass need not equal the "
        "incident neutral mass.");

    if (m_background_mass <= 0.0) {
        m_background_mass = m_projectile_mass;
    }

    const amrex::ParticleReal source_charge = source_species.getCharge();
    const amrex::ParticleReal charge_tolerance =
        1.0e-8_prt*PhysConst::q_e;
    for (int i = 0; i < static_cast<int>(m_processes.size()); ++i) {
        const auto process_type = m_processes[i].type();
        if (process_type != ScatteringProcessType::IONIZATION &&
            process_type != ScatteringProcessType::ATTACHMENT)
        {
            continue;
        }

        auto& product_species = mypc->GetParticleContainerFromName(
            m_product_species_names[i]);
        const amrex::ParticleReal product_charge = product_species.getCharge();

        if (process_type == ScatteringProcessType::IONIZATION) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::abs(product_charge + source_charge) <= charge_tolerance,
                "For impact ionization, the new charged product must balance the "
                "charge of the newly created copy of the incident species.");
        }
        else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::abs(product_charge - source_charge) <= charge_tolerance,
                "For attachment, the negative-ion product must have the same charge "
                "as the incident electron that is removed.");
        }
    }

    m_nu_max = get_nu_max();

    amrex::Print() << Utils::TextMsg::Info(
        "Setting up unified background MCC for " + m_species_names[0]
        + " with " + std::to_string(m_processes.size())
        + " competing channels and nu_max = " + std::to_string(m_nu_max)
        + " s^-1.");

    m_initialized = true;
}

amrex::ParticleReal
BackgroundMCCCollision::get_nu_max () const
{
    using namespace amrex::literals;

    std::vector<amrex::ParticleReal> energy_samples;
    for (auto const& process : m_processes) {
        const int grid_size = process.getEnergyGridSize();
        for (int i = 0; i < grid_size; ++i) {
            energy_samples.push_back(process.getEnergyInput(i));
            if (i + 1 < grid_size) {
                const amrex::ParticleReal energy_lo = process.getEnergyInput(i);
                const amrex::ParticleReal energy_hi = process.getEnergyInput(i + 1);
                const amrex::ParticleReal interval = energy_hi - energy_lo;
                energy_samples.push_back(energy_lo + 0.25_prt*interval);
                energy_samples.push_back(energy_lo + 0.50_prt*interval);
                energy_samples.push_back(energy_lo + 0.75_prt*interval);
            }
        }
    }

    // Include a nearly ultrarelativistic point so the majorant remains valid
    // when a process uses constant high-energy extrapolation. Processes with
    // zero_above_max enabled contribute zero at this point.
    const amrex::ParticleReal rest_energy_eV = static_cast<amrex::ParticleReal>(
        m_projectile_mass*PhysConst::c2/PhysConst::q_e);
    energy_samples.push_back(100.0_prt*rest_energy_eV);

    std::sort(energy_samples.begin(), energy_samples.end());
    energy_samples.erase(
        std::unique(energy_samples.begin(), energy_samples.end()),
        energy_samples.end());

    amrex::ParticleReal nu_max = 0.0_prt;
    for (const auto energy : energy_samples) {
        amrex::ParticleReal effective_cross_section = 0.0_prt;
        for (int i = 0; i < static_cast<int>(m_processes.size()); ++i) {
            effective_cross_section +=
                m_rate_multipliers_h[i]*m_processes[i].getCrossSection(energy);
        }

        const amrex::ParticleReal speed =
            relativeSpeedFromCollisionEnergy(
                energy, m_projectile_mass, m_background_mass);
        const amrex::ParticleReal frequency =
            m_max_background_density*effective_cross_section*speed;
        nu_max = std::max(nu_max, frequency);
    }

    return m_nu_max_safety_factor*nu_max;
}
