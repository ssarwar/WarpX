/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PJG_TEST_UTILS_H_
#define WARPX_PJG_TEST_UTILS_H_

#include "Source/Particles/Collision/ProtonImpactIonization/PJGModel.H"

#include <AMReX_Gpu.H>
#include <AMReX_Vector.H>

using ProtonImpactIonization::PJGModel;
template <typename Real> struct CopiedTable
{
    amrex::Gpu::DeviceVector<Real> m_total;
    amrex::Gpu::DeviceVector<Real> m_energy;
    amrex::Gpu::DeviceVector<Real> m_binding;
    PJGModel::ExecutorT<Real> m_executor;

    explicit CopiedTable (PJGModel::Executor const& source)
    {
        auto copy = [] (auto const* pointer, int count, auto& destination) {
            amrex::Vector<amrex::ParticleReal> host(count);
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, pointer, pointer + count, host.begin());
            amrex::Vector<Real> converted(host.begin(), host.end());
            destination.resize(count);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, converted.begin(), converted.end(),
                             destination.begin());
        };
        copy(source.m_cross_section, PJGModel::table_energy_points, m_total);
        auto const count = PJGModel::table_energy_points * PJGModel::table_quantile_points;
        copy(source.m_log_secondary_energy, count, m_energy);
        copy(source.m_binding_energy, count, m_binding);
        m_executor = {m_total.data(),
                      m_energy.data(),
                      m_binding.data(),
                      Real(source.m_log_projectile_energy_min),
                      Real(source.m_inv_log_projectile_energy_step),
                      Real(source.m_projectile_energy_min),
                      Real(source.m_projectile_energy_max),
                      Real(source.m_projectile_rest_energy),
                      Real(source.m_neutral_rest_energy),
                      Real(source.m_minimum_binding_energy)};
    }
};

#endif // WARPX_PJG_TEST_UTILS_H_
