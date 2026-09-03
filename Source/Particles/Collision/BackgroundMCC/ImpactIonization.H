/* Copyright 2021 Modern Electron
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_IMPACT_IONIZATION_H_
#define WARPX_PARTICLES_COLLISION_IMPACT_IONIZATION_H_

#include "Particles/Algorithms/KineticEnergy.H"
#include "Utils/ParticleUtils.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Random.H>
#include <AMReX_REAL.H>

/** Transform functor for the legacy Background MCC ionization outcome.
 *
 * Process selection is performed by BackgroundMCCCollision. This functor
 * retains the existing equal-energy-sharing outcome until configurable
 * ionization energy-sharing and angular models are added separately.
 */
class ImpactIonizationTransformFunc
{
public:
    /**
     * @param[in] energy_cost energy cost of ionization in eV
     * @param[in] mass1 mass of the colliding species
     * @param[in] sqrt_kb_m value of sqrt(kB/M), where M is neutral mass
     * @param[in] T_a_func background temperature in K as a function of space and time
     * @param[in] t current simulation time
     */
    ImpactIonizationTransformFunc (
        amrex::ParticleReal energy_cost,
        double mass1,
        amrex::ParticleReal sqrt_kb_m,
        amrex::ParserExecutor<4> const& T_a_func,
        amrex::Real t)
        : m_energy_cost(energy_cost),
          m_mass1(mass1),
          m_sqrt_kb_m(sqrt_kb_m),
          m_T_a_func(T_a_func),
          m_t(t)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (
        DstData& dst1,
        DstData& dst2,
        SrcData& src,
        int const i_src,
        int const i_dst1,
        int const i_dst2,
        amrex::RandomEngine const& engine) const noexcept
    {
        using namespace amrex;
        using std::sqrt;

        const auto& p = src.getSuperParticle(i_src);
        ParticleReal x, y, z;
        get_particle_position(p, x, y, z);

        const ParticleReal ion_vel_std =
            m_sqrt_kb_m*sqrt(m_T_a_func(x, y, z, m_t));

        auto& ux = src.m_rdata[PIdx::ux][i_src];
        auto& uy = src.m_rdata[PIdx::uy][i_src];
        auto& uz = src.m_rdata[PIdx::uz][i_src];

        auto& e_ux = dst1.m_rdata[PIdx::ux][i_dst1];
        auto& e_uy = dst1.m_rdata[PIdx::uy][i_dst1];
        auto& e_uz = dst1.m_rdata[PIdx::uz][i_dst1];
        auto& i_ux = dst2.m_rdata[PIdx::ux][i_dst2];
        auto& i_uy = dst2.m_rdata[PIdx::uy][i_dst2];
        auto& i_uz = dst2.m_rdata[PIdx::uz][i_dst2];

        constexpr auto eV = PhysConst::q_e;
        const double E_coll =
            Algorithms::KineticEnergy<double>(ux, uy, uz, m_mass1)/eV;

        const auto E_out = static_cast<ParticleReal>(
            (E_coll - m_energy_cost)/2.0_prt*eV);

        constexpr auto c2 = PhysConst::c2;
        const auto mc2 = m_mass1*c2;
        const ParticleReal up =
            sqrt(E_out*(E_out + 2.0_prt*mc2)/c2)/m_mass1;

        ParticleUtils::RandomizeVelocity(ux, uy, uz, up, engine);
        ParticleUtils::RandomizeVelocity(e_ux, e_uy, e_uz, up, engine);

        i_ux = ion_vel_std*RandomNormal(0_prt, 1.0_prt, engine);
        i_uy = ion_vel_std*RandomNormal(0_prt, 1.0_prt, engine);
        i_uz = ion_vel_std*RandomNormal(0_prt, 1.0_prt, engine);
    }

private:
    amrex::ParticleReal m_energy_cost;
    double m_mass1;
    amrex::ParticleReal m_sqrt_kb_m;
    amrex::ParserExecutor<4> m_T_a_func;
    amrex::Real m_t;
};

#endif // WARPX_PARTICLES_COLLISION_IMPACT_IONIZATION_H_
