/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Fluids/RigidBeamProfile.H"
#include "Particles/ShapeFactors.H"

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
    void require (bool condition, char const* message)
    {
        if (!condition) { throw std::runtime_error(message); }
    }

    template <int Order>
    void compareParticles (double dx, double cutoff)
    {
        using warpx::fluid::projectedIntegral;
        constexpr double sigma = 0.71;
        double const center = 0.17*dx;
        int const cells = static_cast<int>(std::ceil((cutoff*sigma+2.5*dx)/dx));
        std::vector<double> reference(2*cells+1, 0.0);
        // Independent composite-midpoint particle quadrature. Weights sample
        // the Gaussian within each sampling bin; deposition uses WarpX's actual
        // kinetic shape routine, not the analytic fluid projection.
        constexpr int samples = 1000000;
        double const spacing = 2*cutoff*sigma/samples;
        Compute_shape_factor<Order> const compute;
        for (int p = 0; p < samples; ++p) {
            double const z = -cutoff*sigma + (p+0.5)*spacing;
            double const weight = spacing*std::exp(-z*z/(2*sigma*sigma));
            amrex::Real values[Order+1];
            int const left = compute(values, (z+center)/dx+cells);
            for (int i = 0; i <= Order; ++i) { reference[left+i] += weight*values[i]/dx; }
        }
        double error = 0.0;
        double total = 0.0;
        for (int node = -cells; node <= cells; ++node) {
            auto const lo = (node-0.5)*dx-center;
            auto const hi = (node+0.5)*dx-center;
            auto const density = projectedIntegral(Order-1, lo, hi, dx, sigma, cutoff)/dx;
            error = std::max(error, std::abs(density-reference[node+cells]));
            auto const instantaneous = warpx::fluid::projectedDensity(
                Order, node*dx-center, dx, sigma, cutoff);
            error = std::max(error, std::abs(instantaneous-reference[node+cells]));
            total += density*dx;
        }
        amrex::Print() << "order=" << Order << " dx=" << dx << " cutoff=" << cutoff
                       << " particle error=" << error << '\n';
        require(error < 2e-8, "Fluid projection differs from resolved kinetic particle deposition");
        double const expected = std::sqrt(2*std::acos(-1.0))*sigma*std::erf(cutoff/std::sqrt(2.0));
        require(std::abs(total/expected-1) < 3e-14, "Projected pulse has incorrect normalization");

        // Test the device path at both velocity signs and near the support edge.
        amrex::Gpu::DeviceVector<double> residual(200);
        auto* result = residual.data();
        amrex::ParallelFor(200, [=] AMREX_GPU_DEVICE(int i) noexcept {
            double const z = (i-100)*dx/8;
            double const dt = 0.037;
            double const v = i%2 == 0 ? 1.31 : -1.31;
            double const old_rho = projectedIntegral(Order-1, z-dx/2, z+dx/2,
                                                       dx, sigma, cutoff)/dx;
            double const new_rho = projectedIntegral(Order-1, z-dx/2-v*dt, z+dx/2-v*dt,
                                                       dx, sigma, cutoff)/dx;
            double const a = std::min(z-dx/2, z-dx/2-v*dt);
            double const b = std::max(z-dx/2, z-dx/2-v*dt);
            double const sign = v > 0 ? 1.0 : -1.0;
            double const left_j = sign*projectedIntegral(Order-1, a, b, dx, sigma, cutoff)/dt;
            double const right_j = sign*projectedIntegral(Order-1, a+dx, b+dx,
                                                          dx, sigma, cutoff)/dt;
            result[i] = (new_rho-old_rho) + dt/dx*(right_j-left_j);
        });
        std::vector<double> host(200);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, residual.begin(), residual.end(), host.begin());
        for (double value : host) {
            require(std::abs(value) < 2e-14, "Prescribed current violates discrete continuity");
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        for (double cutoff : {2.0, 8.0}) {
            for (double dx : {0.1, 0.7, 1.4}) {
                compareParticles<1>(dx, cutoff);
                compareParticles<2>(dx, cutoff);
                compareParticles<3>(dx, cutoff);
                compareParticles<4>(dx, cutoff);
            }
        }
        amrex::Print() << "PASS: normalization, kinetic shapes 1--4 and discrete continuity\n";
    }
    amrex::Finalize();
}
