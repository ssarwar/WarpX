/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCKinematics.H"
#include "Particles/Collision/ProtonImpactIonization/ProtonImpactIonizationKinematics.H"
#include "Particles/Collision/ProtonImpactIonization/RelativisticIonization.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
    void
    require (bool condition, char const* message)
    {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    void
    accumulate (double value, double& sum, double& correction)
    {
        auto const adjusted = value - correction;
        auto const next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }

    template <typename Real>
    void
    checkAngles ()
    {
        constexpr int count = 4096;
        constexpr double m = 510998.95069;
        constexpr double mass = 938272088.16;
        constexpr auto tol = 32 * std::numeric_limits<Real>::epsilon();
        amrex::Gpu::DeviceVector<Real> output(4 * count);
        auto* data = output.data();
        amrex::Vector<Real> host(output.size());
        for (auto const energy : {5.e3, 5.e4, 8.e8, 1.e10}) {
            auto const maximum =
                ProtonImpactIonization::relativisticMaximumTransfer(Real(energy), Real(mass));
            for (auto const fraction : {0., 1.e-6, .1, .9, 1., 2.}) {
                auto const secondary = Real(fraction * maximum);
                for (auto const binding : {Real(0), Real(12.07), Real(15.59), Real(543.5)}) {
                    amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int j) noexcept {
                        auto const u = (Real(j) + Real(.5)) / Real(count);
                        auto const mu =
                            ProtonImpactIonization::polarCosine(secondary, binding, maximum, u);
                        data[4 * j] = mu;
                        using namespace BackgroundMCCKinematics;
                        // Alternate axes exercise both transverse-basis branches.
                        Vector3 const axis =
                            j % 2 == 0 ? Vector3{0., 0., -1.} : Vector3{.36, -.48, .8};
                        Vector3 e1, e2;
                        transverseDirections(axis, e1, e2);
                        auto const direction = directionFromPolarAngle(
                            axis, e1, e2, mu,
                            2.0 * MathConst::pi * double((j * 2371) % count) / count);
                        data[4 * j + 1] = Real(dot(direction, direction));
                        data[4 * j + 2] = Real(dot(direction, axis));
                        data[4 * j + 3] =
                            ProtonImpactIonization::freeElectronCosine(secondary, maximum);
                    });
                    amrex::Gpu::copy(amrex::Gpu::deviceToHost, output.begin(), output.end(),
                                     host.begin());
                    auto const t = std::min(double(secondary), double(maximum));
                    auto const mu_free =
                        std::sqrt(t * (maximum + 2 * m) / (double(maximum) * (t + 2 * m)));
                    auto const denom = double(secondary) + binding;
                    auto const a = denom > 0 ? mu_free * (secondary + .5 * binding) / denom : 0;
                    auto const b = denom > 0 ? binding / denom : 1;
                    auto const lo = std::max(-1., a - b);
                    auto const hi = std::min(1., a + b);
                    double mean = 0, second = 0, mean_correction = 0, second_correction = 0;
                    int forward = 0;
                    for (int j = 0; j < count; ++j) {
                        auto const mu = host[4 * j];
                        auto const expected = lo + (hi - lo) * (j + .5) / count;
                        require(std::isfinite(mu) && mu >= -1 && mu <= 1, "Invalid polar cosine");
                        require(std::abs(mu - expected) < tol, "Angular CDF/quantile mismatch");
                        require(std::abs(host[4 * j + 1] - 1) < tol, "Non-unit emitted direction");
                        require(std::abs(host[4 * j + 2] - mu) < tol, "Axis-dependent polar angle");
                        require(std::abs(host[4 * j + 3] - mu_free) < tol, "Wrong binary angle");
                        forward += mu == 1;
                        accumulate(double(mu) / count, mean, mean_correction);
                        accumulate(double(mu) * mu / count, second, second_correction);
                    }
                    require(std::abs(mean - (lo + hi) / 2) < tol, "Wrong angular mean");
                    require(std::abs(second - (lo * lo + lo * hi + hi * hi) / 3) <
                                tol + 1. / (count * count),
                            "Wrong angular second moment");
                    if (hi - lo > .01) {
                        require(forward == 0, "Artificial forward point mass");
                    }
                    // Free scattering independently leaves the projectile on shell:
                    // p_projectile * p_e * mu = (E_total + m_e) * T.
                    if (secondary > 0 && secondary <= maximum) {
                        auto const p = std::sqrt(energy * (energy + 2 * mass));
                        auto const pe = std::sqrt(double(secondary) * (secondary + 2 * m));
                        require(std::abs(p * pe * host[3] / ((energy + mass + m) * secondary) - 1) <
                                    tol,
                                "Binary cone violates projectile mass shell");
                    }
                }
            }
        }
    }
} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    checkAngles<float>();
    checkAngles<double>();
    std::cout << "PASS: isotropic, binary, bound-tail and rotated angular "
                 "distributions\n";
    amrex::Finalize();
}
