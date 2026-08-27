/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCCCOLLISION_H_

#include "BackgroundMCCElasticScattering.H"
#include "BackgroundMCCIonization.H"
#include "Particles/Collision/CollisionBase.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Particles/MultiParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>
#include <vector>

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
     * Product-changing events are recorded for a later grouped particle-creation
     * pass. Pure scattering events are applied in place.
     *
     * @param pti particle iterator
     * @param t current time
     * @param selected_process optional selected process index per particle
     * @param neutral_vx optional sampled neutral x velocity per particle
     * @param neutral_vy optional sampled neutral y velocity per particle
     * @param neutral_vz optional sampled neutral z velocity per particle
     */
    void doBackgroundCollisionsWithinTile (
        WarpXParIter& pti,
        amrex::Real t,
        int* selected_process,
        amrex::ParticleReal* neutral_vx,
        amrex::ParticleReal* neutral_vy,
        amrex::ParticleReal* neutral_vz);

private:
    struct ProductGroup
    {
        ScatteringProcessType type = ScatteringProcessType::INVALID;
        std::string species_name;
    };

    amrex::Vector<ScatteringProcess> m_processes;
    amrex::Gpu::DeviceVector<ScatteringProcess::Executor> m_processes_exe;
    amrex::Gpu::DeviceVector<BackgroundMCCElasticScatteringModel::Executor>
        m_elastic_scattering_processes_exe;
    amrex::Gpu::DeviceVector<BackgroundMCCIonizationModel::Executor>
        m_ionization_processes_exe;
    std::vector<std::unique_ptr<BackgroundMCCElasticScatteringModel>>
        m_elastic_scattering_models;
    std::unique_ptr<BackgroundMCCIonizationModel> m_n2_ionization_model;
    std::unique_ptr<BackgroundMCCIonizationModel> m_o2_ionization_model;
    amrex::Vector<ProductGroup> m_product_groups;
    amrex::Gpu::DeviceVector<int> m_process_product_group;

    bool init_flag = false;
    bool m_has_iaa_elastic_processes = false;
    bool m_has_ionization_processes = false;
    bool m_has_attachment_processes = false;
    bool m_use_relativistic_electron_kinematics = false;
    bool m_user_nu_max = false;
    bool m_warned_large_dt = false;

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
