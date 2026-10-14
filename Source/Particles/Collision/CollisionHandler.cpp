/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "CollisionHandler.H"

#include "Particles/Collision/BackgroundMCC/BackgroundMCCCollision.H"
#include "Particles/Collision/BackgroundStopping/BackgroundStopping.H"
#include "Particles/Collision/BinaryCollision/BinaryCollision.H"
#include "Particles/Collision/BinaryCollision/Bremsstrahlung/BremsstrahlungFunc.H"
#include "Particles/Collision/BinaryCollision/Bremsstrahlung/PhotonCreationFunc.H"
#include "Particles/Collision/BinaryCollision/Coulomb/PairWiseCoulombCollisionFunc.H"
#include "Particles/Collision/BinaryCollision/DSMC/DSMCFunc.H"
#include "Particles/Collision/BinaryCollision/DSMC/SplitAndScatterFunc.H"
#include "Particles/Collision/BinaryCollision/NuclearFusion/NuclearFusionFunc.H"
#include "Particles/Collision/BinaryCollision/LinearBreitWheeler/LinearBreitWheelerCollisionFunc.H"
#include "Particles/Collision/BinaryCollision/LinearCompton/LinearComptonCollisionFunc.H"
#include "Particles/Collision/BinaryCollision/ParticleCreationFunc.H"
#ifdef WARPX_QED
#include "Particles/Collision/BinaryCollision/VirtualPhotonCreation.H"
#endif
#include "Particles/Collision/HybridResistiveDrag/HybridResistiveDrag.H"
#include "Particles/Collision/InverseBremsstrahlung/InverseBremsstrahlung.H"
#include "Particles/Collision/ProtonImpactIonization/ProtonImpactIonization.H"
#include "Particles/Collision/PulsedDecay/PulsedDecay.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

std::string
CollisionHandler::CheckpointConfiguration () const
{
    std::string configuration;
    for (auto const& collision : allcollisions) { configuration += collision->CheckpointConfiguration(); }
    return configuration;
}

void
CollisionHandler::WriteCheckpoint (std::string const& directory) const
{
    auto const configuration = CheckpointConfiguration();
    if (configuration.empty() || !amrex::ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream output(directory+"/PrescribedSources");
    output << "WarpX prescribed sources 1\n" << configuration;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(output.good(), "Cannot write prescribed-source checkpoint metadata.");
}

void
CollisionHandler::ValidateRestart (std::string const& directory) const
{
    auto const configuration = CheckpointConfiguration();
    auto const path = directory+"/PrescribedSources";
    if (configuration.empty() && !std::filesystem::exists(path)) { return; }
    amrex::Vector<char> contents;
    amrex::ParallelDescriptor::ReadAndBcastFile(path, contents);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::string(contents.data()) ==
        "WarpX prescribed sources 1\n"+configuration,
        "Prescribed collision sources or their immutable physics/sampling configuration changed.");
    auto const& warpx = WarpX::GetInstance();
    for (auto const& collision : allcollisions) {
        for (auto const& field : collision->CheckpointFields()) {
            auto const field_path = directory+"/Level_0/"+warpx.m_fields.mf_name(field, 0);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::VisMF::Exist(field_path),
                "Missing required prescribed-source checkpoint state: "+field_path);
        }
    }
}

void
CollisionHandler::ValidateRestartState () const
{
    for (auto const& collision : allcollisions) { collision->ValidateRestartState(); }
}

CollisionHandler::CollisionHandler(MultiParticleContainer const * const mypc)
{

    // Read in collision input
    const amrex::ParmParse pp_collisions("collisions");
    pp_collisions.queryarr("collision_names", collision_names);

    // Create instances based on the collision type
    auto const ncollisions = collision_names.size();
    collision_types.resize(ncollisions);
    allcollisions.resize(ncollisions);
    m_schedule.reserve(ncollisions);
    for (int i = 0; i < static_cast<int>(ncollisions); ++i) {
        const amrex::ParmParse pp_collision_name(collision_names[i]);

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::n_rz_azimuthal_modes==1,
        "RZ mode `warpx.n_rz_azimuthal_modes` must be 1 when using the binary collision module.");

        // For legacy, pairwisecoulomb is the default
        std::string type = "pairwisecoulomb";

        pp_collision_name.query("type", type);
        collision_types[i] = type;

        if (type == "pairwisecoulomb") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<PairWiseCoulombCollisionFunc>>(
                    collision_names[i], mypc
                );
            m_use_global_debye_length |= allcollisions[i]->use_global_debye_length();
        }
        else if (type == "background_mcc") {
            allcollisions[i] = std::make_unique<BackgroundMCCCollision>(collision_names[i]);
        }
        else if (type == "pulsed_decay") {
            allcollisions[i] = std::make_unique<PulsedDecay>(collision_names[i], mypc);
        }
        else if (type == "proton_impact_ionization") {
            allcollisions[i] =
                std::make_unique<ProtonImpactIonizationCollision>(collision_names[i], mypc);
        }
        else if (type == "background_stopping") {
            allcollisions[i] = std::make_unique<BackgroundStopping>(collision_names[i]);
        }
        else if (type == "hybrid_resistive_drag") {
            allcollisions[i] = std::make_unique<HybridResistiveDrag>(collision_names[i]);
        }
        else if (type == "dsmc") {
            allcollisions[i] =
                std::make_unique<BinaryCollision<DSMCFunc, SplitAndScatterFunc>>(
                    collision_names[i], mypc
                );
        }
        else if (type == "nuclearfusion") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<NuclearFusionFunc, ParticleCreationFunc>>(
                    collision_names[i], mypc
                );
        }
        else if (type == "bremsstrahlung") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<BremsstrahlungFunc, PhotonCreationFunc>>(
                    collision_names[i], mypc
                );
        }
        else if (type == "inverse_bremsstrahlung") {
            allcollisions[i] = std::make_unique<InverseBremsstrahlung>(collision_names[i], mypc);
            m_use_global_debye_length = true;
        }
        else if (type == "linear_breit_wheeler") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<LinearBreitWheelerCollisionFunc, ParticleCreationFunc>>(
                    collision_names[i], mypc
               );
        }
        else if (type == "linear_compton") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<LinearComptonCollisionFunc, ParticleCreationFunc>>(
                    collision_names[i], mypc
               );
        }
        else{
            WARPX_ABORT_WITH_MESSAGE("Unknown collision type.");
        }

    }

}

/* \brief Allocate any data needed for the collision */
void CollisionHandler::AllocData ()
{
    for (auto& collision : allcollisions) {
        collision->AllocData();
    }
}

/** Perform all collisions
 *
 * @param step Current iteration
 * @param cur_time Current time
 * @param dt Time step
 * @param mypc MultiParticleContainer calling this method
 *
 */
void CollisionHandler::doCollisions ( int step, amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc)
{

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    /* In RZ and RCYLINDER geometry, macroparticles can collide with other macroparticles
     * in the same *cylindrical* cell, or in RSPHERE the same *spherical* shell.
     * Because of this, the colliding macroparticles would not nessecarily be spatially
     * near each other. This would violate the underlying assumptions that particles within the
     * same cylindrical or spherical cell represent a cylindrically- or spherically-symmetric
     * momentum distribution function and are spatially local. Therefore, we temporarily rotate
     * the momentum of the macroparticles to the curvilinear frame, equivalent to the x-axis.
     * (This is only valid if we use only the m=0 azimuthal mode in the simulation;
     * there is a corresponding assert statement at initialization.) */
    mypc->TransformMomentumToCurvilinear(/*forward*/true);
#endif

#ifdef WARPX_QED
    // For QED incoherent processes (e.g. Bethe-Heitler, Landau-Lifschitz), the process is mediated by virtual photons.
    // The virtual photons are newly generated here and participate in the collisions.
    // Here, the virtual photons are regenerated from scratch, i.e. they are overwritten by new ones at each time step.
    if(mypc->nSpecies() > 0) {
        collision::binarycollision::virtualphotons::GenerateVirtualPhotons(mypc);
    }
#endif

    if (m_use_global_debye_length) {
        // This will calculate the temperature, Vbar, and particle number that are needed by
        // the various collision algorithms
        mypc->GenerateGlobalDebyeLength();
    }

    bool const after_push = WarpX::GetInstance().evolve_scheme != EvolveScheme::Explicit;
    auto const start_time = after_push ? cur_time-dt : cur_time;
    m_schedule.clear();
    for (int i = 0; i < static_cast<int>(allcollisions.size()); ++i) {
        auto const& collision = allcollisions[i];
        // Skip collisions before their start step
        const int start_step = collision->get_start_step();
        if (step < start_step) { continue; }

        const int ndt = collision->get_ndt();
        if (collision->get_collision_stepping_mode() == CollisionSteppingMode::Subcycle) {
            m_schedule.push_back({i, 0, ndt});
        } else if ((step - start_step) % ndt == 0) {
            m_schedule.push_back({i, 0, 1});
        }
    }

    // Interleave coupled operators at their substep endpoints. Advancing every
    // substep of one gas before the next retains a full-PIC-step splitting error.
    // Integer fractions preserve coincident endpoints and input order, even for
    // different subcycle counts. The heap stores one entry per collision object
    // and reuses its allocation across PIC steps.
    auto later = [after_push] (CollisionStep const& a, CollisionStep const& b) {
        auto const left = (std::int64_t(a.m_substep) + int(after_push)) * b.m_subcycles;
        auto const right = (std::int64_t(b.m_substep) + int(after_push)) * a.m_subcycles;
        return left != right ? left > right : a.m_collision > b.m_collision;
    };
    std::make_heap(m_schedule.begin(), m_schedule.end(), later);
    while (!m_schedule.empty()) {
        std::pop_heap(m_schedule.begin(), m_schedule.end(), later);
        auto next = m_schedule.back();
        m_schedule.pop_back();
        auto& collision = allcollisions[next.m_collision];
        const int ndt = collision->get_ndt();
        if (collision->get_collision_stepping_mode() == CollisionSteppingMode::Subcycle) {
            const amrex::Real dt_sub = dt / ndt;
            // Explicit collisions use left endpoints; implicit collisions use
            // right endpoints. Both advance the same physical interval.
            int const offset = after_push ? next.m_substep + 1 - ndt : next.m_substep;
            collision->doCollisionsInInterval(cur_time + offset*dt_sub,
                start_time + next.m_substep*dt_sub, dt_sub, mypc);
        } else {
            collision->doCollisionsInInterval(cur_time, start_time, dt*ndt, mypc);
        }
        if (++next.m_substep < next.m_subcycles) {
            m_schedule.push_back(next);
            std::push_heap(m_schedule.begin(), m_schedule.end(), later);
        }
    }

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    // Undo the rotation above
    mypc->TransformMomentumToCurvilinear(/*forward*/false);
#endif

}
