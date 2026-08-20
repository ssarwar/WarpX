/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_

#include "Particles/Collision/CollisionBase.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>

class BackgroundMCCCollision final
    : public CollisionBase
{
public:
    explicit BackgroundMCCCollision (std::string const& collision_name);

    ~BackgroundMCCCollision () override = default;

    BackgroundMCCCollision (BackgroundMCCCollision const&) = delete;
    BackgroundMCCCollision& operator= (BackgroundMCCCollision const&) = delete;
    BackgroundMCCCollision (BackgroundMCCCollision&&) = delete;
    BackgroundMCCCollision& operator= (BackgroundMCCCollision&&) = delete;

    /** Return a conservative null-collision majorant for all channels. */
    [[nodiscard]] amrex::ParticleReal get_nu_max () const;

    /** Perform all background-MCC channels using one competing-process draw. */
    void doCollisions (
        amrex::Real cur_time, amrex::Real dt,
        MultiParticleContainer* mypc) override;

private:
    void initialize (MultiParticleContainer* mypc);

    void selectCollisionEventsWithinTile (
        WarpXParIter& pti,
        amrex::Real t,
        amrex::ParticleReal total_collision_probability,
        amrex::Gpu::DeviceVector<int>& event_process) const;

    void createReactionProductsWithinTile (
        int lev,
        WarpXParIter& pti,
        long source_particle_count,
        amrex::Gpu::DeviceVector<int> const& event_process,
        WarpXParticleContainer& source_species,
        MultiParticleContainer* mypc,
        amrex::Real t) const;

    amrex::Vector<ScatteringProcess> m_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor> m_processes_exe;

    // One process name, optional product species, and rate multiplier per channel.
    // A multiplier of n_M turns an m^5 three-body attachment table into an
    // effective m^2 cross section.
    amrex::Vector<std::string> m_process_names;
    amrex::Vector<std::string> m_product_species_names;
    amrex::Vector<amrex::ParticleReal> m_rate_multipliers_h;
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_rate_multipliers_d;

    bool m_initialized = false;
    bool m_warned_large_step = false;

    double m_projectile_mass = -1.0;
    double m_background_mass = -1.0;

    amrex::ParticleReal m_max_background_density = 0.0;
    amrex::ParticleReal m_nu_max = 0.0;
    amrex::ParticleReal m_nu_max_safety_factor = 1.05;

    amrex::Parser m_background_density_parser;
    amrex::Parser m_background_temperature_parser;
    amrex::ParserExecutor<4> m_background_density_func;
    amrex::ParserExecutor<4> m_background_temperature_func;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_
