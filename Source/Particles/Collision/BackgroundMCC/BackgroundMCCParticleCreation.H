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
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/TextMsg.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

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
     * @param selected_process selected process index for each original electron
     * @param product_offsets destination-group offset for each product event
     * @param neutral_vx sampled neutral x velocity for each product event
     * @param neutral_vy sampled neutral y velocity for each product event
     * @param neutral_vz sampled neutral z velocity for each product event
     * @param product_counts number of selected events in each product group
     */
    [[nodiscard]] amrex::Long
    operator() (
        WarpXParticleContainer& source_species,
        ParticleTileType& source_tile,
        std::vector<WarpXParticleContainer*> const& product_species,
        amrex::Vector<ParticleTileType*> const& product_tiles,
        int const* selected_process,
        int const* product_offsets,
        amrex::ParticleReal const* neutral_vx,
        amrex::ParticleReal const* neutral_vy,
        amrex::ParticleReal const* neutral_vz,
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

        amrex::Gpu::HostVector<int> old_product_counts_h(group_count);
        amrex::Gpu::HostVector<ParticleTileDataType> product_data_h(group_count);
        for (std::size_t group = 0; group < group_count; ++group)
        {
            auto& product_tile = *product_tiles[group];
            auto const old_count = product_tile.numParticles();
            auto const new_count = old_count + product_counts[group];
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(new_count <= std::numeric_limits<int>::max(),
                "Background MCC particle tiles must contain fewer than INT_MAX "
                "particles.");
            old_product_counts_h[group] = static_cast<int>(old_count);
            product_tile.resize(new_count);
            product_data_h[group] = product_tile.getParticleTileData();
        }

        int const* product_group_types = m_product_group_types_h.data();
        double const* product_masses = m_product_masses_h.data();
        SmartCopy const* copy_products = m_copy_products_h.data();
        int const* old_product_counts = old_product_counts_h.data();
        int const* electron_group_offsets = electron_group_offsets_h.data();
        ParticleTileDataType* product_data = product_data_h.data();

#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<int> old_product_counts_d(group_count);
        amrex::Gpu::DeviceVector<int> electron_group_offsets_d(group_count);
        amrex::Gpu::DeviceVector<ParticleTileDataType> product_data_d(group_count);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            old_product_counts_h.begin(),
            old_product_counts_h.end(),
            old_product_counts_d.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            electron_group_offsets_h.begin(),
            electron_group_offsets_h.end(),
            electron_group_offsets_d.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
            product_data_h.begin(),
            product_data_h.end(),
            product_data_d.begin());
        product_group_types = m_product_group_types_d.data();
        product_masses = m_product_masses_d.data();
        copy_products = m_copy_products_d.data();
        old_product_counts = old_product_counts_d.data();
        electron_group_offsets = electron_group_offsets_d.data();
        product_data = product_data_d.data();
#endif

        auto source_data = source_tile.getParticleTileData();
        auto const copy_electron = m_copy_electron;
        auto const* processes = m_processes;
        auto const* ionization_processes = m_ionization_processes;
        auto const* process_product_group = m_process_product_group;
        auto const electron_mass = m_electron_mass;
        auto const old_source_index = static_cast<int>(old_source_count);

        auto const ionization_transform = ImpactIonizationTransformFunc(processes,
            ionization_processes,
            selected_process,
            neutral_vx,
            neutral_vy,
            neutral_vz,
            electron_mass);
        auto const attachment_transform =
            AttachmentTransformFunc(neutral_vx, neutral_vy, neutral_vz);

        amrex::ParallelForRNG(old_source_count,
            [=] AMREX_GPU_DEVICE(int const source_index, amrex::RandomEngine const& engine) noexcept
            {
                int const process = selected_process[source_index];
                if (process < 0)
                {
                    return;
                }

                int const group = process_product_group[process];
                int const group_offset = product_offsets[source_index];
                int const product_index = old_product_counts[group] + group_offset;
                auto source = source_data;
                auto product = product_data[group];
                copy_products[group](product, source, source_index, product_index, engine);

                auto const type = static_cast<ScatteringProcessType>(product_group_types[group]);
                if (type == ScatteringProcessType::ATTACHMENT)
                {
                    attachment_transform(product, source, source_index, product_index, engine);
                }
                else
                {
                    int const electron_index =
                        old_source_index + electron_group_offsets[group] + group_offset;
                    copy_electron(source, source, source_index, electron_index, engine);
                    ionization_transform(source,
                        product,
                        source,
                        source_index,
                        electron_index,
                        product_index,
                        product_masses[group],
                        engine);
                }
            });

        if (total_ionization_events > 0)
        {
            ParticleCreation::DefaultInitializeRuntimeAttributes(source_tile,
                source_species,
                static_cast<int>(old_source_count),
                static_cast<int>(new_source_count));
            setNewParticleIDs(source_tile, old_source_count, total_ionization_events);
        }
        for (std::size_t group = 0; group < group_count; ++group)
        {
            auto const count = product_counts[group];
            if (count == 0)
            {
                continue;
            }
            auto const old_count = old_product_counts_h[group];
            auto const new_count = old_count + count;
            ParticleCreation::DefaultInitializeRuntimeAttributes(
                *product_tiles[group], *product_species[group], old_count, new_count);
            setNewParticleIDs(*product_tiles[group], old_count, count);
        }

        // Particle-tile data and group offsets are local temporaries. Complete
        // every dependent kernel before their device storage leaves scope.
        amrex::Gpu::synchronize();
        return attachment_events;
    }

private:
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
