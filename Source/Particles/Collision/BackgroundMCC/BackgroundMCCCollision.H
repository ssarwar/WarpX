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

    /** Return a null-collision majorant for all supplied processes. */
    [[nodiscard]] amrex::ParticleReal get_nu_max (
        amrex::Vector<ScatteringProcess> const& processes) const;

    /** Perform the collisions.
     *
     * @param cur_time Current time
     * @param dt Collision time step
     * @param mypc Container of species involved
     */
    void doCollisions (
        amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc) override;

    /** Select and apply at most one Background MCC process per particle.
     *
     * Product-producing ionization events are recorded in ionization_mask and
     * transformed in a later particle-creation pass. Pure scattering events are
     * applied in place.
     *
     * @param pti particle iterator
     * @param t current time
     * @param ionization_mask optional per-particle mask for ionization events
     */
    void doBackgroundCollisionsWithinTile (
        WarpXParIter& pti, amrex::Real t, amrex::Long* ionization_mask);

private:
    amrex::Vector<ScatteringProcess> m_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor> m_processes_exe;

    bool init_flag = false;
    bool ionization_flag = false;
    bool m_use_relativistic_electron_kinematics = false;
    bool m_user_nu_max = false;
    bool m_warned_large_dt = false;

    int m_ionization_process_index = -1;

    amrex::ParticleReal m_mass1 = 0;
    amrex::ParticleReal m_max_background_density = 0;
    amrex::ParticleReal m_background_mass = -1;
    amrex::ParticleReal m_total_collision_prob = 0;
    amrex::ParticleReal m_nu_max = 0;

    amrex::Parser m_background_density_parser;
    amrex::Parser m_background_temperature_parser;

    amrex::ParserExecutor<4> m_background_density_func;
    amrex::ParserExecutor<4> m_background_temperature_func;
};

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_
