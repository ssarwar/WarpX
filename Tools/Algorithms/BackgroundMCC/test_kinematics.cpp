/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCElasticKinematics.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCIonizationKinematics.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
    struct Case
    {
        double m_energy, m_loss, m_target_mass, m_cosine, m_share;
        int m_kind;
        ScatteringAngleModel m_angle;
    };

    // Independent on-shell kinetic energy, avoiding subtraction of rest energies.
    long double
    kineticEnergy (amrex::ParticleReal const* u, long double mass)
    {
        long double u2 = 0;
        for (int d = 0; d < 3; ++d) {
            u2 += static_cast<long double>(u[d]) * u[d];
        }
        auto const c2 = static_cast<long double>(PhysConst::c2_v<double>);
        return mass * u2 / ((std::sqrt(1 + u2 / c2) + 1) * PhysConst::q_e_v<double>);
    }
} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr double me = PhysConst::m_e_v<double>;
        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double qe = PhysConst::q_e_v<double>;
        amrex::Vector<Case> cases;
        for (auto const mass_number : {28.0134, 31.9988}) {
            auto const mass = mass_number * PhysConst::m_u_v<double>;
            auto const rest = mass * c2 / qe;
            for (auto const loss : {0.0, 8.0, 15.0}) {
                auto const threshold = loss * (1 + me / mass) + loss * loss / (2 * rest);
                for (auto const excess : {1.e-10, 1.e-8, 1.e-5, 1.0, 1.e3, 1.e9}) {
                    for (auto const cosine : {-1.0, -.5, 0.0, .5, 1.0}) {
                        cases.push_back({threshold + excess, loss, mass, cosine, 0,
                                         loss == 0 ? 0 : 1, ScatteringAngleModel::IAA});
                    }
                }
            }
            for (auto const excess : {.001, .1, 1., 100., 1.e6, 1.e9, 1.e11}) {
                for (auto const share : {0., 1.e-8, .01, .25, .5}) {
                    for (auto const angle :
                         {ScatteringAngleModel::IAA, ScatteringAngleModel::Forward,
                          ScatteringAngleModel::Backward, ScatteringAngleModel::Isotropic}) {
                        cases.push_back({15.59 + excess, 15.59, mass - me, 0, share, 2, angle});
                    }
                }
            }
            // Above a shell binding, but below its recoil-shifted threshold,
            // no three-product state exists. The source must remain unchanged.
            for (auto const energy : {15.0, 15.59, 15.5902}) {
                cases.push_back(
                    {energy, 15.59, mass - me, 0, .5, 3, ScatteringAngleModel::Backward});
            }
        }
        amrex::Gpu::DeviceVector<Case> device_cases(cases.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(),
                         device_cases.begin());
        amrex::Gpu::DeviceVector<amrex::ParticleReal> device_results(12 * cases.size());
        auto const* inputs = device_cases.data();
        auto* results = device_results.data();
        amrex::ParallelForRNG(
            static_cast<int>(cases.size()),
            [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const& engine) noexcept {
                auto const state = inputs[i];
                auto* r = results + 12 * i;
                r[0] = 0;
                r[1] = 0;
                r[2] = static_cast<amrex::ParticleReal>(
                    BackgroundMCCKinematics::properSpeedFromKineticEnergy(state.m_energy, me));
                for (int d = 3; d < 12; ++d) {
                    r[d] = 0;
                }
                if (state.m_kind == 0) {
                    BackgroundMCCElasticKinematics::compute(
                        r[0], r[1], r[2], 0, 0, 0, me, state.m_target_mass, state.m_cosine, engine,
                        r[3], r[4], r[5], r[9], r[10], r[11]);
                } else if (state.m_kind == 1) {
                    BackgroundMCCElasticKinematics::computeExcitation(
                        r[0], r[1], r[2], 0, 0, 0, me, state.m_target_mass,
                        static_cast<amrex::ParticleReal>(state.m_loss), state.m_cosine, engine,
                        r[3], r[4], r[5], r[9], r[10], r[11]);
                } else {
                    auto const accepted = BackgroundMCCIonizationKinematics::compute(
                        r[0], r[1], r[2], 0, 0, 0, me, state.m_target_mass,
                        static_cast<amrex::ParticleReal>(state.m_loss),
                        static_cast<amrex::ParticleReal>(state.m_share), state.m_angle, engine,
                        r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11]);
                    if (state.m_kind == 3) {
                        r[6] = accepted ? 1 : 0;
                    } else if (!accepted) {
                        r[3] = std::numeric_limits<amrex::ParticleReal>::quiet_NaN();
                    }
                }
            });
        amrex::Vector<amrex::ParticleReal> host(device_results.size());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_results.begin(), device_results.end(),
                         host.begin());
        std::array<double, 3> max_errors{};
        int failures = 0;
        auto const tolerance = 512 * std::numeric_limits<amrex::ParticleReal>::epsilon();
        for (std::size_t i = 0; i < cases.size(); ++i) {
            auto const& state = cases[i];
            auto const* r = host.data() + 12 * i;
            if (state.m_kind == 3) {
                for (int d = 3; d < 12; ++d) {
                    if (r[d] != 0) {
                        throw std::runtime_error("Inadmissible ionization changed outputs");
                    }
                }
                continue;
            }
            auto const mass =
                state.m_target_mass + (state.m_kind == 1 ? state.m_loss * qe / c2 : 0);
            auto const initial = kineticEnergy(r, me);
            auto const final = kineticEnergy(r + 3, me) + kineticEnergy(r + 6, me) +
                               kineticEnergy(r + 9, mass) + state.m_loss;
            auto const error = static_cast<double>(std::abs(final - initial) / initial);
            max_errors[state.m_kind] = std::max(max_errors[state.m_kind], error);
            bool valid = std::isfinite(error) && error < tolerance;
            for (int d = 0; d < 3; ++d) {
                auto const residual = me * (r[d] - r[3 + d] - r[6 + d]) - mass * r[9 + d];
                valid = valid && std::abs(residual) < tolerance * me * r[2];
            }
            if (!valid) {
                if (failures < 12 || state.m_kind == 2) {
                    std::cout << "Failed kind=" << state.m_kind << " E=" << state.m_energy
                              << " loss=" << state.m_loss << " mu=" << state.m_cosine
                              << " share=" << state.m_share << " angle=" << int(state.m_angle)
                              << " relative energy error=" << error << '\n';
                }
                ++failures;
            }
        }
        std::cout << "Elastic/excitation/ionization maximum relative energy errors: "
                  << max_errors[0] << ' ' << max_errors[1] << ' ' << max_errors[2] << '\n';
        if (failures != 0) {
            throw std::runtime_error(std::to_string(failures) + " MCC conservation checks failed");
        }
    }
    amrex::Finalize();
}
