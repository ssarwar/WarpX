/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Fluids/RigidBeam.H"
#include "Particles/ShapeFactors.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
void
require (bool condition, char const* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

// Second antiderivative of the truncated unit Gaussian, anchored at -cutoff.
// This closed form is independent of the production space-time quadrature.
long double
gaussianPrimitive (long double x, long double cutoff) {
    if (x <= -cutoff) {
        return 0.0L;
    }
    auto const upper = std::min(x, cutoff);
    auto const lower_exp = std::exp(-cutoff * cutoff / 2);
    auto const integral =
        std::sqrt(std::acos(-1.0L) / 2) * (std::erf(upper / std::sqrt(2.0L)) +
                                           std::erf(cutoff / std::sqrt(2.0L)));
    return x * integral + std::exp(-upper * upper / 2) - lower_exp;
}

void
compareSourceIntegrals () {
    // Overlapping pulses, unequal amplitudes and a nonzero reference plane.
    std::vector<double> const times{-0.7, 0.2, 0.35};
    std::vector<double> const amplitudes{0.5, 0.0, 1.5};
    amrex::Gpu::DeviceVector<double> times_d(times.size()),
        amplitudes_d(amplitudes.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, times.begin(), times.end(),
                     times_d.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, amplitudes.begin(),
                     amplitudes.end(), amplitudes_d.begin());
    for (double cutoff : {2.0, 8.0, std::numeric_limits<double>::infinity()}) {
        for (double velocity : {-1.31, 1.31}) {
            RigidBeam::Executor profile;
            profile.m_velocity = velocity;
            profile.m_sigma_r = 0.71;
            profile.m_sigma_z = 0.43;
            profile.m_peak_density = 2.3;
            profile.m_z_reference = 0.17;
            profile.m_cutoff_r = profile.m_cutoff_z = cutoff;
            profile.m_times = times_d.data();
            profile.m_amplitudes = amplitudes_d.data();
            profile.m_count = static_cast<int>(times.size());
            amrex::Gpu::DeviceVector<double> values(48);
            auto* results = values.data();
            amrex::ParallelFor(24, [=] AMREX_GPU_DEVICE(int i) noexcept {
                double const zlo = -0.9 + (i % 4) * 0.5;
                double const width = i % 2 == 0 ? 0.07 : 1.4;
                double const start = -2.0 + (i / 4) * 0.55;
                double const duration = i % 3 == 0 ? 5.3 : 0.17;
                results[2 * i] = profile.integratedLongitudinal(
                    zlo, zlo + width, start, duration);
                results[2 * i + 1] =
                    profile.radialIntegral(i * 0.13, (i + 1) * 0.13);
            });
            std::vector<double> host(values.size());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, values.begin(),
                             values.end(), host.begin());
            for (int i = 0; i < 24; ++i) {
                long double const zlo = -0.9 + (i % 4) * 0.5;
                long double const width = i % 2 == 0 ? 0.07 : 1.4;
                long double const start = -2.0 + (i / 4) * 0.55;
                long double const duration = i % 3 == 0 ? 5.3 : 0.17;
                long double const sigma = profile.m_sigma_z;
                long double expected = 0;
                for (int pulse = 0; pulse < profile.m_count; ++pulse) {
                    auto const a = (zlo - profile.m_z_reference -
                                    velocity * (start - times[pulse])) /
                                   sigma;
                    auto const b = a - velocity * duration / sigma;
                    expected += amplitudes[pulse] * sigma * sigma / velocity *
                                (gaussianPrimitive(a + width / sigma, cutoff) -
                                 gaussianPrimitive(a, cutoff) -
                                 gaussianPrimitive(b + width / sigma, cutoff) +
                                 gaussianPrimitive(b, cutoff));
                }
                expected *= profile.m_peak_density;
                double const scale =
                    profile.m_peak_density * double(width * duration);
                require(std::abs(host[2 * i] - expected) < 5e-13 * scale,
                        "Beam space-time source differs from the Gaussian "
                        "antiderivative");
                auto const lo = i * 0.13;
                auto const hi =
                    std::min((i + 1) * 0.13, cutoff * profile.m_sigma_r);
                auto const sigma2 = profile.m_sigma_r * profile.m_sigma_r;
                double const radial =
                    hi > lo ? 2 * std::acos(-1.0) * sigma2 *
                                  (std::exp(-lo * lo / (2 * sigma2)) -
                                   std::exp(-hi * hi / (2 * sigma2)))
                            : 0;
                require(std::abs(host[2 * i + 1] - radial) < 2e-15 * sigma2,
                        "Beam annular source differs from the cylindrical "
                        "Gaussian integral");
            }
        }
    }
}

template <int Order>
void
compareParticles (double dx, double cutoff) {
    using warpx::fluid::projectedIntegral;
    constexpr double sigma = 0.71;
    double const center = 0.17 * dx;
    int const cells =
        static_cast<int>(std::ceil((cutoff * sigma + 2.5 * dx) / dx));
    std::vector<double> reference(2 * cells + 1, 0.0);
    // Independent composite-midpoint particle quadrature. Weights sample
    // the Gaussian within each sampling bin; deposition uses WarpX's actual
    // kinetic shape routine, not the analytic fluid projection.
    constexpr int samples = 1000000;
    double const spacing = 2 * cutoff * sigma / samples;
    Compute_shape_factor<Order> const compute;
    for (int p = 0; p < samples; ++p) {
        double const z = -cutoff * sigma + (p + 0.5) * spacing;
        double const weight = spacing * std::exp(-z * z / (2 * sigma * sigma));
        amrex::Real values[Order + 1];
        int const left = compute(values, (z + center) / dx + cells);
        for (int i = 0; i <= Order; ++i) {
            reference[left + i] += weight * values[i] / dx;
        }
    }
    double error = 0.0;
    double total = 0.0;
    for (int node = -cells; node <= cells; ++node) {
        auto const lo = (node - 0.5) * dx - center;
        auto const hi = (node + 0.5) * dx - center;
        auto const density =
            projectedIntegral(Order - 1, lo, hi, dx, sigma, cutoff) / dx;
        error = std::max(error, std::abs(density - reference[node + cells]));
        auto const instantaneous = warpx::fluid::projectedDensity(
            Order, node * dx - center, dx, sigma, cutoff);
        error =
            std::max(error, std::abs(instantaneous - reference[node + cells]));
        total += density * dx;
    }
    amrex::Print() << "order=" << Order << " dx=" << dx << " cutoff=" << cutoff
                   << " particle error=" << error << '\n';
    require(
        error < 2e-8,
        "Fluid projection differs from resolved kinetic particle deposition");
    double const expected = std::sqrt(2 * std::acos(-1.0)) * sigma *
                            std::erf(cutoff / std::sqrt(2.0));
    require(std::abs(total / expected - 1) < 3e-14,
            "Projected pulse has incorrect normalization");

    // Test the device path at both velocity signs and near the support edge.
    amrex::Gpu::DeviceVector<double> residual(200);
    auto* result = residual.data();
    amrex::ParallelFor(200, [=] AMREX_GPU_DEVICE(int i) noexcept {
        double const z = (i - 100) * dx / 8;
        double const dt = 0.037;
        double const v = i % 2 == 0 ? 1.31 : -1.31;
        double const old_rho =
            projectedIntegral(Order - 1, z - dx / 2, z + dx / 2, dx, sigma,
                              cutoff) /
            dx;
        double const new_rho =
            projectedIntegral(Order - 1, z - dx / 2 - v * dt,
                              z + dx / 2 - v * dt, dx, sigma, cutoff) /
            dx;
        double const a = std::min(z - dx / 2, z - dx / 2 - v * dt);
        double const b = std::max(z - dx / 2, z - dx / 2 - v * dt);
        double const sign = v > 0 ? 1.0 : -1.0;
        double const left_j =
            sign * projectedIntegral(Order - 1, a, b, dx, sigma, cutoff) / dt;
        double const right_j =
            sign *
            projectedIntegral(Order - 1, a + dx, b + dx, dx, sigma, cutoff) /
            dt;
        result[i] = (new_rho - old_rho) + dt / dx * (right_j - left_j);
    });
    std::vector<double> host(200);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, residual.begin(), residual.end(),
                     host.begin());
    for (double value : host) {
        require(std::abs(value) < 2e-14,
                "Prescribed current violates discrete continuity");
    }
}
} // namespace

int
main (int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    {
        compareSourceIntegrals();
        for (double cutoff : {2.0, 8.0}) {
            for (double dx : {0.1, 0.7, 1.4}) {
                compareParticles<1>(dx, cutoff);
                compareParticles<2>(dx, cutoff);
                compareParticles<3>(dx, cutoff);
                compareParticles<4>(dx, cutoff);
            }
        }
        amrex::Print() << "PASS: exact source integrals, normalization, "
                          "kinetic shapes 1--4 "
                       << "and discrete continuity\n";
    }
    amrex::Finalize();
}
