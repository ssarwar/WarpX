/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCIonizationKinematics.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
template <typename Real>
void
checkAngles ()
{
    constexpr int count = 1024;
    constexpr Real rest = Real(510998.95069);
    constexpr double tolerance = 16 * std::numeric_limits<Real>::epsilon();
    amrex::Gpu::DeviceVector<Real> device(count + 2);
    amrex::Vector<Real> host(count + 2);
    auto* output = device.data();
    for (Real nominal_available : {Real(1.e-3), Real(50), Real(1.e3), Real(2.5e6), Real(1.e9)}) {
        for (Real binding : {Real(0), Real(15.58), Real(543.8)}) {
            Real const energy = nominal_available + binding;
            Real const available = energy - binding;
            for (Real fraction : {Real(0), Real(1.e-6), Real(.5), Real(1)}) {
                Real const secondary = fraction * available;
                amrex::ParallelFor(count + 2, [=] AMREX_GPU_DEVICE(int i) noexcept {
                    Real const u = i < count ? (Real(i) + Real(.5)) / Real(count) : Real(i - count);
                    output[i] = BackgroundMCCIonizationKinematics::secondaryCosine(
                        energy, binding, secondary, rest, u);
                });
                amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(),
                                 host.begin());
                // Independently impose the outgoing primary mass shell:
                // p_in * p_secondary * mu_binary = (available + 2*m) * secondary.
                long double const a = available, t = secondary, m = rest;
                long double const pin = std::sqrt(a * (a + 2 * m));
                long double const pout = std::sqrt(t * (t + 2 * m));
                double const binary = t > 0 ? double((a + 2 * m) * t / (pin * pout)) : 0;
                double const weight =
                    secondary + binding > 0 ? double(secondary) / (double(secondary) + binding) : 0;
                double const center = weight * binary, width = 1 - weight;
                double mean = 0, second = 0, mean_correction = 0, second_correction = 0;
                for (int i = 0; i < count + 2; ++i) {
                    double const value = host[i];
                    if (!std::isfinite(value) || value < -1 || value > 1) {
                        throw std::runtime_error("Secondary cosine requires clipping");
                    }
                    double const u = i < count ? (i + .5) / count : i - count;
                    if (std::abs(value - (center + width * (1 - 2 * u))) > tolerance) {
                        throw std::runtime_error("Secondary angular quantile mismatch");
                    }
                    if (i < count) {
                        auto accumulate = [] (double v, double& sum, double& correction) {
                            double const adjusted = v - correction;
                            double const next = sum + adjusted;
                            correction = (next - sum) - adjusted;
                            sum = next;
                        };
                        accumulate(value / count, mean, mean_correction);
                        accumulate(value * value / count, second, second_correction);
                    }
                }
                if (std::abs(mean - center) > tolerance ||
                    std::abs(second - (center * center + width * width / 3)) >
                        tolerance + 1.0 / (count * count)) {
                    throw std::runtime_error("Incorrect secondary angular moments");
                }
                if (binding == 0 && available < .01 &&
                    std::abs(mean - std::sqrt(double(fraction))) > 1.e-8 + tolerance) {
                    throw std::runtime_error("Incorrect nonrelativistic binary limit");
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
    std::cout
        << "PASS: secondary angular mass shell, isotropy, nonrelativistic limit and moments\n";
    amrex::Finalize();
}
