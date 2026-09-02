/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_H_
#define WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_H_

#include "PJGModel.H"
#include "Particles/Collision/CollisionBase.H"

#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

/** Rigid-beam proton and bare-ion impact ionization of a neutral background.
 *
 * The incident projectile is not modified. Each sampled event creates one
 * kinetic electron and one singly charged molecular ion.
 */
class ProtonImpactIonizationCollision final : public CollisionBase
{
public:
    ProtonImpactIonizationCollision (std::string const& collision_name,
                                     MultiParticleContainer const* mypc);

    ~ProtonImpactIonizationCollision () override = default;

    ProtonImpactIonizationCollision (ProtonImpactIonizationCollision const&) = delete;
    ProtonImpactIonizationCollision& operator= (ProtonImpactIonizationCollision const&) = delete;
    ProtonImpactIonizationCollision (ProtonImpactIonizationCollision&&) = delete;
    ProtonImpactIonizationCollision& operator= (ProtonImpactIonizationCollision&&) = delete;

    void AllocData () override;

    void doCollisions (amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc) override;

private:
    amrex::Vector<std::string> m_product_species;
    ProtonImpactIonization::PJGTarget m_target = ProtonImpactIonization::PJGTarget::N2;
    std::unique_ptr<ProtonImpactIonization::PJGModel> m_pjg_model;

    amrex::Parser m_background_density_parser;
    amrex::Parser m_background_temperature_parser;
    amrex::ParserExecutor<4> m_background_density_func;
    amrex::ParserExecutor<4> m_background_temperature_func;

    std::string m_remainder_field_name;
    amrex::ParticleReal m_fixed_product_weight = 0.0;
    amrex::ParticleReal m_projectile_charge_squared = 0.0;
    int m_max_products_per_cell = 64;
};

#endif // WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_H_
