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
BackgroundMCCCollision::createReactionProductsWithinTile (
    int lev,
    WarpXParIter& pti,
    long source_particle_count,
    amrex::Gpu::DeviceVector<int> const& event_process,
    WarpXParticleContainer& source_species,
    MultiParticleContainer* mypc,
    amrex::Real t) const
{
    using namespace amrex::literals;

    if (source_particle_count == 0) {
        return;
    }

    const int* event_process_ptr = event_process.dataPtr();
    auto& source_tile = source_species.ParticlesAt(lev, pti);

    // Attachment channels are processed first. They do not resize the source
    // electron tile, so all masks retain the original source-particle extent.
    for (int process_index = 0;
         process_index < static_cast<int>(m_processes.size());
         ++process_index)
    {
        if (m_processes[process_index].type() !=
            ScatteringProcessType::ATTACHMENT)
        {
            continue;
        }

        amrex::Gpu::DeviceVector<amrex::Long> mask(source_particle_count);
        amrex::Long* mask_ptr = mask.dataPtr();
        amrex::ParallelFor(
            source_particle_count,
            [=] AMREX_GPU_DEVICE (long ip) noexcept
            {
                mask_ptr[ip] =
                    (event_process_ptr[ip] == process_index) ? 1L : 0L;
            });

        auto& product_species = mypc->GetParticleContainerFromName(
            m_product_species_names[process_index]);
        auto& product_tile = product_species.ParticlesAt(lev, pti);
        const amrex::Long old_product_count = product_tile.numParticles();

        const SmartCopyFactory copy_factory(source_species, product_species);
        const SmartCopy copy_product = copy_factory.getSmartCopy();
        const amrex::ParticleReal sqrt_kb_m_product =
            std::sqrt(PhysConst::kb/product_species.getMass());
        const AttachmentTransformFunc transform(
            sqrt_kb_m_product, m_background_temperature_func, t);

        const amrex::Long num_added = filterCopyTransformParticles<1>(
            product_species,
            product_tile,
            source_tile,
            mask.dataPtr(),
            old_product_count,
            copy_product,
            transform);
        setNewParticleIDs(product_tile, old_product_count, num_added);
    }

    // The existing ionization transform creates one new copy of the incident
    // species and one positive ion. It is deliberately processed after all
    // attachments because it resizes the source tile.
    for (int process_index = 0;
         process_index < static_cast<int>(m_processes.size());
         ++process_index)
    {
        if (m_processes[process_index].type() !=
            ScatteringProcessType::IONIZATION)
        {
            continue;
        }

        amrex::Gpu::DeviceVector<amrex::Long> mask(source_particle_count);
        amrex::Long* mask_ptr = mask.dataPtr();
        amrex::ParallelFor(
            source_particle_count,
            [=] AMREX_GPU_DEVICE (long ip) noexcept
            {
                mask_ptr[ip] =
                    (event_process_ptr[ip] == process_index) ? 1L : 0L;
            });

        auto& product_species = mypc->GetParticleContainerFromName(
            m_product_species_names[process_index]);
        auto& product_tile = product_species.ParticlesAt(lev, pti);
        const amrex::Long old_source_count = source_tile.numParticles();
        const amrex::Long old_product_count = product_tile.numParticles();

        const SmartCopyFactory electron_copy_factory(
            source_species, source_species);
        const SmartCopyFactory product_copy_factory(
            source_species, product_species);
        const SmartCopy copy_electron =
            electron_copy_factory.getSmartCopy();
        const SmartCopy copy_product =
            product_copy_factory.getSmartCopy();

        const amrex::ParticleReal sqrt_kb_m =
            std::sqrt(PhysConst::kb/m_background_mass);
        const ImpactIonizationTransformFunc transform(
            m_processes[process_index].getEnergyPenalty(),
            m_projectile_mass,
            sqrt_kb_m,
            m_background_temperature_func,
            t);

        const amrex::Long num_added = filterCopyTransformParticles<1>(
            source_species,
            product_species,
            source_tile,
            product_tile,
            source_tile,
            mask.dataPtr(),
            old_source_count,
            old_product_count,
            copy_electron,
            copy_product,
            transform);

        setNewParticleIDs(source_tile, old_source_count, num_added);
        setNewParticleIDs(product_tile, old_product_count, num_added);
    }
}
