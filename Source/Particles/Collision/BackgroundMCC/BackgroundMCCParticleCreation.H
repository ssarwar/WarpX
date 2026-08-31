/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_PARTICLE_CREATION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_PARTICLE_CREATION_H_

#include "BackgroundMCCProducts.H"

#include "Particles/ParticleCreation/DefaultInitialization.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/TextMsg.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Particle.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

/** Data retained for one selected Background MCC product-changing event. */
struct BackgroundMCCProductEvent
{
    amrex::ParticleReal m_neutral_vx;
    amrex::ParticleReal m_neutral_vy;
    amrex::ParticleReal m_neutral_vz;
    int m_source_index;
    int m_process;
    int m_group_offset;
};

/** Create every Background MCC product selected in one source particle tile.
 *
 * The collision-selection kernel assigns each product event an offset within
 * its destination group. This class consequently resizes every destination
 * once and creates all secondary electrons and ions in one device kernel,
 * independent of the number of configured product channels or groups.
 */
class BackgroundMCCParticleCreation
{
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = ParticleTileType::ParticleTileDataType;

    struct ProductTileMetadata
    {
        ParticleTileDataType m_data;
        amrex::Long m_first_id;
        int m_old_count;
        int m_electron_offset;
    };

public:
    BackgroundMCCParticleCreation (
        amrex::Vector<ScatteringProcessType> const& product_group_types,
        std::vector<WarpXParticleContainer*> const& product_species,
        SmartCopy const copy_electron,
        std::vector<SmartCopy> const& copy_products,
        ScatteringProcess::Executor const* processes,
        BackgroundMCCIonizationModel::Executor const* ionization_processes,
        int const* process_product_group,
        double const electron_mass)
        : m_product_group_types_h(product_group_types.size()),
          m_product_masses_h(product_group_types.size()),
          m_copy_electron(copy_electron),
          m_processes(processes),
          m_ionization_processes(ionization_processes),
          m_process_product_group(process_product_group),
          m_electron_mass(electron_mass)
    {
        auto const group_count = static_cast<std::size_t>(product_group_types.size());
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(product_species.size() == group_count,
            "Background MCC product species and group counts must match.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(copy_products.size() == group_count,
            "Background MCC product copy and group counts must match.");

        for (std::size_t group = 0; group < group_count; ++group)
        {
            m_product_group_types_h[group] = static_cast<int>(product_group_types[group]);
            m_product_masses_h[group] = static_cast<double>(product_species[group]->getMass());
        }

#ifdef AMREX_USE_GPU
        m_product_group_types_d.resize(group_count);
        m_product_masses_d.resize(group_count);
        m_copy_products_d.resize(group_count);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            m_product_group_types_h.begin(),
            m_product_group_types_h.end(),
            m_product_group_types_d.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            m_product_masses_h.begin(),
            m_product_masses_h.end(),
            m_product_masses_d.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            copy_products.begin(),
            copy_products.end(),
            m_copy_products_d.begin());
        amrex::Gpu::streamSynchronize();
#else
        m_copy_products_h.resize(group_count);
        std::copy(copy_products.begin(), copy_products.end(), m_copy_products_h.begin());
#endif
    }

    /** Create products and return the number of attached source electrons.
     *
     * @param source_species source electron particle container
     * @param source_tile source electron particle tile
     * @param product_species product particle containers, ordered by group
     * @param product_tiles product particle tiles, ordered by group
     * @param product_events compact selected product-changing events
     * @param product_counts number of selected events in each product group
     */
    [[nodiscard]] amrex::Long
    operator() (
        WarpXParticleContainer& source_species,
        ParticleTileType& source_tile,
        std::vector<WarpXParticleContainer*> const& product_species,
        amrex::Vector<ParticleTileType*> const& product_tiles,
        BackgroundMCCProductEvent const* product_events,
        amrex::Vector<int> const& product_counts) const
    {
        auto const group_count = static_cast<std::size_t>(product_counts.size());
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(product_species.size() == group_count,
            "Background MCC product species and group counts must match.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<std::size_t>(product_tiles.size()) == group_count,
            "Background MCC product tile and group counts must match.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_product_group_types_h.size() == group_count,
            "Background MCC product metadata and group counts must match.");

        int total_events = 0;
        int total_ionization_events = 0;
        amrex::Long attachment_events = 0;
        amrex::Gpu::HostVector<int> electron_group_offsets_h(group_count, -1);
        for (std::size_t group = 0; group < group_count; ++group)
        {
            auto const count = product_counts[group];
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                count >= 0, "Background MCC product counts must be non-negative.");
            total_events += count;
            auto const type = static_cast<ScatteringProcessType>(m_product_group_types_h[group]);
            if (type == ScatteringProcessType::IONIZATION)
            {
                electron_group_offsets_h[group] = total_ionization_events;
                total_ionization_events += count;
            }
            else
            {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(type == ScatteringProcessType::ATTACHMENT,
                    "Background MCC product groups must be ionization or "
                    "attachment.");
                attachment_events += count;
            }
        }
        if (total_events == 0)
        {
            return 0;
        }

        auto const old_source_count = source_tile.numParticles();
        auto const new_source_count = old_source_count + total_ionization_events;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(new_source_count <= std::numeric_limits<int>::max(),
            "Background MCC particle tiles must contain fewer than INT_MAX "
            "particles.");
        source_tile.resize(new_source_count);
        amrex::Long const first_electron_id = total_ionization_events > 0
            ? reserveParticleIDs(total_ionization_events)
            : -1;

        amrex::Gpu::HostVector<ProductTileMetadata> product_tile_metadata_h(group_count);
        for (std::size_t group = 0; group < group_count; ++group)
        {
            auto& product_tile = *product_tiles[group];
            auto const old_count = product_tile.numParticles();
            auto const count = product_counts[group];
            auto const new_count = old_count + count;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(new_count <= std::numeric_limits<int>::max(),
                "Background MCC particle tiles must contain fewer than INT_MAX "
                "particles.");
            product_tile.resize(new_count);
            product_tile_metadata_h[group] = {
                product_tile.getParticleTileData(),
                count > 0 ? reserveParticleIDs(count) : -1,
                static_cast<int>(old_count),
                electron_group_offsets_h[group]};
        }

        int const* product_group_types = m_product_group_types_h.data();
        double const* product_masses = m_product_masses_h.data();
        SmartCopy const* copy_products = m_copy_products_h.data();
        ProductTileMetadata const* product_tile_metadata =
            product_tile_metadata_h.data();

#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<ProductTileMetadata> product_tile_metadata_d(group_count);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            product_tile_metadata_h.begin(),
            product_tile_metadata_h.end(),
            product_tile_metadata_d.begin());
        product_group_types = m_product_group_types_d.data();
        product_masses = m_product_masses_d.data();
        copy_products = m_copy_products_d.data();
        product_tile_metadata = product_tile_metadata_d.data();
#endif

        auto source_data = source_tile.getParticleTileData();
        auto const copy_electron = m_copy_electron;
        auto const* processes = m_processes;
        auto const* ionization_processes = m_ionization_processes;
        auto const* process_product_group = m_process_product_group;
        auto const electron_mass = m_electron_mass;
        auto const old_source_index = static_cast<int>(old_source_count);
        int const cpu = amrex::ParallelDescriptor::MyProc();

        auto const ionization_transform =
            ImpactIonizationTransformFunc(processes, ionization_processes, electron_mass);
        AttachmentTransformFunc const attachment_transform;

        amrex::ParallelForRNG(total_events,
            [=] AMREX_GPU_DEVICE(int const event_index, amrex::RandomEngine const& engine) noexcept
            {
                auto const& event = product_events[event_index];
                int const source_index = event.m_source_index;
                int const process = event.m_process;
                int const group = process_product_group[process];
                int const group_offset = event.m_group_offset;
                auto const& group_metadata = product_tile_metadata[group];
                int const product_index = group_metadata.m_old_count + group_offset;
                auto source = source_data;
                auto product = group_metadata.m_data;
                copy_products[group](product, source, source_index, product_index, engine);

                auto const type = static_cast<ScatteringProcessType>(product_group_types[group]);
                if (type == ScatteringProcessType::ATTACHMENT)
                {
                    attachment_transform(product,
                        source,
                        source_index,
                        product_index,
                        event.m_neutral_vx,
                        event.m_neutral_vy,
                        event.m_neutral_vz,
                        engine);
                }
                else
                {
                    int const electron_index =
                        old_source_index + group_metadata.m_electron_offset + group_offset;
                    copy_electron(source, source, source_index, electron_index, engine);
                    ionization_transform(source,
                        product,
                        source,
                        source_index,
                        electron_index,
                        product_index,
                        process,
                        event.m_neutral_vx,
                        event.m_neutral_vy,
                        event.m_neutral_vz,
                        product_masses[group],
                        engine);
                    source.m_idcpu[electron_index] = amrex::SetParticleIDandCPU(
                        first_electron_id + group_metadata.m_electron_offset +
                            group_offset,
                        cpu);
                }
                product.m_idcpu[product_index] = amrex::SetParticleIDandCPU(
                    group_metadata.m_first_id + group_offset, cpu);
            });

        if (total_ionization_events > 0)
        {
            ParticleCreation::DefaultInitializeRuntimeAttributes(source_tile,
                source_species,
                static_cast<int>(old_source_count),
                static_cast<int>(new_source_count));
        }
        for (std::size_t group = 0; group < group_count; ++group)
        {
            auto const count = product_counts[group];
            if (count == 0)
            {
                continue;
            }
            auto const old_count = product_tile_metadata_h[group].m_old_count;
            auto const new_count = old_count + count;
            ParticleCreation::DefaultInitializeRuntimeAttributes(
                *product_tiles[group], *product_species[group], old_count, new_count);
        }

        // Particle-tile data and group offsets are local temporaries. Complete
        // every dependent kernel before their device storage leaves scope.
        amrex::Gpu::synchronize();
        return attachment_events;
    }

private:
    /** Reserve a contiguous range of particle IDs and return its first ID. */
    [[nodiscard]] static amrex::Long
    reserveParticleIDs (amrex::Long const count)
    {
        amrex::Long first_id;
#ifdef AMREX_USE_OMP
#pragma omp critical (ionization_nextid)
#endif
        {
            first_id = ParticleTileType::ParticleType::NextID();
            ParticleTileType::ParticleType::NextID(first_id + count);
        }
        return first_id;
    }

    amrex::Gpu::HostVector<int> m_product_group_types_h;
    amrex::Gpu::HostVector<double> m_product_masses_h;
    amrex::Gpu::HostVector<SmartCopy> m_copy_products_h;
#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<int> m_product_group_types_d;
    amrex::Gpu::DeviceVector<double> m_product_masses_d;
    amrex::Gpu::DeviceVector<SmartCopy> m_copy_products_d;
#endif
    SmartCopy m_copy_electron;
    ScatteringProcess::Executor const* m_processes;
    BackgroundMCCIonizationModel::Executor const* m_ionization_processes;
    int const* m_process_product_group;
    double m_electron_mass;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_PARTICLE_CREATION_H_
