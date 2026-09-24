/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/ProtonImpactIonization/ProtonImpactIonizationKinematics.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
    using BackgroundMCCKinematics::Vector3;
    using ProtonImpactIonization::Products;

    double
    kinetic (Vector3 u, double rest)
    {
        double const u2 = (u.x * u.x + u.y * u.y + u.z * u.z) / PhysConst::c2_v<double>;
        return rest * u2 / (std::sqrt(1 + u2) + 1);
    }

    template <typename Real>
    void
    checkRecoil ()
    {
        constexpr int count = 257;
        constexpr double c = PhysConst::c_v<double>;
        constexpr double conversion = PhysConst::c2_v<double> / PhysConst::q_e_v<double>;
        constexpr double mass = PhysConst::m_p_v<double>;
        constexpr double projectile = mass * conversion;
        constexpr double electron = 510998.95069;
        double worst_energy = 0, worst_momentum = 0;
        amrex::Gpu::DeviceVector<Products> device(count);
        amrex::Vector<Products> host(count);
        auto* output = device.data();
        for (double a : {28.0134, 31.9988}) {
            double const neutral_mass = a * 1.66053906660e-27;
            double const neutral = neutral_mass * conversion;
            for (double energy : {5000., 5.e4, 8.e8, 1.e10}) {
                double const speed = c * std::sqrt(energy * (energy + 2 * projectile)) / projectile;
                Vector3 const incoming{.36 * speed, -.48 * speed, .8 * speed};
                for (double temperature_case : {0., 1.}) {
                    Vector3 const velocity{300 * temperature_case, -150 * temperature_case,
                                           600 * temperature_case};
                    auto const frame =
                        ProtonImpactIonization::neutralFrame(incoming, velocity, mass);
                    for (double binding : {0., 15.59, 543.5}) {
                        double const maximum =
                            ProtonImpactIonization::molecularMaximumSecondaryEnergy(
                                frame.m_energy, projectile, neutral, binding);
                        double const free = ProtonImpactIonization::relativisticMaximumTransfer(
                            frame.m_energy, projectile);
                        for (double secondary :
                             {0., 1.e-8, .1 * free, .9 * free, 1.1 * free, .5 * maximum,
                              .999999 * maximum, double(Real(maximum))}) {
                            amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int j) noexcept {
                                auto products = ProtonImpactIonization::compute(
                                    frame, secondary, binding, mass, neutral_mass, (j + .5) / count,
                                    2 * double(MathConst::pi) * j / count,
                                    std::numeric_limits<Real>::epsilon());
                                // Exercise the stored particle precision, not only double helpers.
                                products.m_electron = {Real(products.m_electron.x),
                                                       Real(products.m_electron.y),
                                                       Real(products.m_electron.z)};
                                products.m_projectile = {Real(products.m_projectile.x),
                                                         Real(products.m_projectile.y),
                                                         Real(products.m_projectile.z)};
                                products.m_ion = {Real(products.m_ion.x), Real(products.m_ion.y),
                                                  Real(products.m_ion.z)};
                                output[j] = products;
                            });
                            amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(),
                                             host.begin());
                            double const beta2 =
                                (velocity.x * velocity.x + velocity.y * velocity.y +
                                 velocity.z * velocity.z) /
                                (c * c);
                            double const gamma = 1 / std::sqrt(1 - beta2);
                            double const initial_energy =
                                energy + neutral * gamma * gamma * beta2 / (gamma + 1);
                            Vector3 const initial_pc{
                                projectile * incoming.x / c + neutral * gamma * velocity.x / c,
                                projectile * incoming.y / c + neutral * gamma * velocity.y / c,
                                projectile * incoming.z / c + neutral * gamma * velocity.z / c};
                            for (auto const& result : host) {
                                if (!result.m_valid) {
                                    std::cerr << "Invalid state E=" << energy << " T=" << secondary
                                              << " B=" << binding << '\n';
                                    throw std::runtime_error(
                                        "No recoil state for allowed electron");
                                }
                                double const ion = result.m_ion_rest_energy;
                                double const final_energy =
                                    kinetic(result.m_electron, electron) +
                                    kinetic(result.m_projectile, projectile) +
                                    kinetic(result.m_ion, ion) + binding;
                                double const energy_error =
                                    std::abs(final_energy - initial_energy) / initial_energy;
                                Vector3 const difference{
                                    (electron * result.m_electron.x +
                                     projectile * result.m_projectile.x + ion * result.m_ion.x) /
                                            c -
                                        initial_pc.x,
                                    (electron * result.m_electron.y +
                                     projectile * result.m_projectile.y + ion * result.m_ion.y) /
                                            c -
                                        initial_pc.y,
                                    (electron * result.m_electron.z +
                                     projectile * result.m_projectile.z + ion * result.m_ion.z) /
                                            c -
                                        initial_pc.z};
                                double const momentum_error = std::sqrt(
                                    (difference.x * difference.x + difference.y * difference.y +
                                     difference.z * difference.z) /
                                    (initial_pc.x * initial_pc.x + initial_pc.y * initial_pc.y +
                                     initial_pc.z * initial_pc.z));
                                worst_energy = std::max(worst_energy, energy_error);
                                worst_momentum = std::max(worst_momentum, momentum_error);
                                double const tolerance = 32 * std::numeric_limits<Real>::epsilon();
                                if (energy_error > tolerance || momentum_error > tolerance) {
                                    std::cerr << "E=" << energy << " T=" << secondary
                                              << " B=" << binding
                                              << " relative energy=" << energy_error
                                              << " momentum=" << momentum_error << '\n';
                                    throw std::runtime_error(
                                        "Recoil violates four-momentum conservation");
                                }
                                if (binding == 0 && secondary <= free && temperature_case == 0 &&
                                    kinetic(result.m_ion, ion) > 1.e-14 * energy) {
                                    throw std::runtime_error(
                                        "Free binary collision recoils the spectator ion");
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "Recoil maximum relative energy/momentum errors (" << sizeof(Real)
                  << " bytes): " << worst_energy << " / " << worst_momentum << '\n';
    }
} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    checkRecoil<float>();
    checkRecoil<double>();
    amrex::Finalize();
}
