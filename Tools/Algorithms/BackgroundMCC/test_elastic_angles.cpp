/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCElasticScattering.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCUtils.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <cmath>
#include <iostream>
#include <limits>

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr int count = 65536;
        constexpr double width = 36.0 / count;
        amrex::Gpu::DeviceVector<double> device(count);
        amrex::Vector<double> host(count);
        auto* output = device.data();
        amrex::ParallelForRNG(count,
            [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const& engine) noexcept {
                double const draw = BackgroundMCCUtils::uniformDouble(engine);
                double const grid = draw * 16777216.0;
                output[i] = draw >= .5 && draw < 1 && grid != std::floor(grid) ? 1 : 0;
            });
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), host.begin());
        double resolved = 0;
        for (double value : host) { resolved += value; }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(resolved > .4 * count,
            "IAA angular uniforms must resolve probabilities below the binary32 grid.");

        for (double radius : {0.6052, 0.5677}) {
            BackgroundMCCElasticScatteringModel::Executor angular;
            angular.m_screening_radius = radius;
            for (amrex::ParticleReal energy : {amrex::ParticleReal(1e4),
                    amrex::ParticleReal(2.5e6), amrex::ParticleReal(1e9)}) {
                // Logarithmic strata in 1-u resolve the rare large-angle tail.
                // Each stratum retains its actual uniform-CDF probability mass.
                amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int i) noexcept {
                    double const u = -std::expm1(-(i + .5) * width);
                    double const weight = std::exp(-i * width) * -std::expm1(-width);
                    output[i] = weight * (1 - angular.sampleCosine(energy, u));
                });
                amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                    device.begin(), device.end(), host.begin());
                long double moment = 0;
                for (double value : host) { moment += value; }

                long double const tau = static_cast<long double>(energy) / 510998.95069L;
                long double const alpha = 7.2973525693e-3L;
                long double const eta = alpha * alpha /
                    (4 * radius * radius * tau * (tau + 2));
                // Integral of d*f(d), f(d)=2*eta*(1+eta)/(d+2*eta)^2,
                // for d=1-cos(theta) in [0,2]. This is sigma_MT/sigma_total.
                long double const exact = 2 * eta * ((1 + eta) * std::log1p(1 / eta) - 1);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::abs(moment - exact) < 2e-5L * exact +
                        4 * std::numeric_limits<double>::epsilon(),
                    "Screened-Rutherford momentum transfer lost its rare angular tail.");

                // These quantiles straddle the analytically known backward tail.
                // Both round to one in binary32 at one GeV.
                AMREX_ALWAYS_ASSERT(angular.sampleCosine(energy, double(1 - 2 * eta)) > .3);
                AMREX_ALWAYS_ASSERT(angular.sampleCosine(energy, double(1 - eta / 2)) < -.3);
            }
        }
        std::cout << "PASS: IAA angular precision and exact momentum-transfer moments\n";
    }
    amrex::Finalize();
}
