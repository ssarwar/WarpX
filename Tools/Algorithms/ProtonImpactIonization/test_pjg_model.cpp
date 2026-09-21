/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PJGReferenceData.H"
#include "PJGTestUtils.H"
#include "Source/Particles/Collision/ProtonImpactIonization/PJGModel.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

using ProtonImpactIonization::PJGModel;
using ProtonImpactIonization::PJGTarget;

namespace
{
    constexpr double proton_mass = 938272088.16;

    void
    require (bool condition, char const* message)
    {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    template <typename Real>
    void
    checkTables (PJGModel const& model, PJGTarget target)
    {
        CopiedTable<Real> table(model.executor());
        auto const exec = table.m_executor;
        // Integrate the actual accelerator executor with positive Gauss weights
        // in the clustered quantile coordinate, resolving rare hard electrons.
        constexpr int cells = 4096;
        constexpr int nodes = 4;
        constexpr int count = cells * nodes;
        constexpr std::array<double, nodes> gx{-.8611363115940526, -.3399810435848563,
                                               .3399810435848563, .8611363115940526};
        constexpr std::array<double, nodes> gw{.3478548451374538, .6521451548625461,
                                               .6521451548625461, .3478548451374538};
        amrex::Vector<double> probabilities(count), weights(count);
        for (int i = 0; i < count; ++i) {
            auto const x = (i / nodes + .5 * (gx[i % nodes] + 1)) / cells;
            auto const a = std::pow(x, 4), b = std::pow(1 - x, 4);
            probabilities[i] = a / (a + b);
            weights[i] =
                4 * std::pow(x * (1 - x), 3) / ((a + b) * (a + b)) * gw[i % nodes] / (2 * cells);
        }
        amrex::Gpu::DeviceVector<double> device_q(count);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, probabilities.begin(), probabilities.end(),
                         device_q.begin());
        amrex::Gpu::DeviceVector<Real> sampled(6 * count);
        using SamplingState = typename PJGModel::ExecutorT<Real>::SamplingState;
        amrex::Gpu::DeviceVector<SamplingState> prepared(1);
        auto* prepared_device = prepared.data();
        auto* output = sampled.data();
        auto const* q = device_q.data();
        amrex::Vector<Real> host(6 * count);
        double max_total = 0, max_mean = 0, max_second = 0, max_binding = 0;
        double max_host_device_difference = 0;
        for (int row = 0; row <= 48; ++row) {
            auto const e = static_cast<Real>(5e3 * std::pow(2e6, row / 48.0));
            auto const reference = PJGModel::integratedMoments(target, e, proton_mass);
            auto const state = exec.prepareSampling(e);
            amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
                prepared_device[0] = exec.prepareSampling(e);
            });
            amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int i) noexcept {
                exec.sample(state, q[i], output[i], output[count + i]);
                exec.sample(prepared_device[0], q[i], output[2 * count + i], output[3 * count + i]);
                exec.sample(e, q[i], output[4 * count + i], output[5 * count + i]);
            });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, sampled.begin(), sampled.end(),
                             host.begin());
            for (int i = 0; i < count; ++i) {
                // Compare caching on the same backend. Host and device log/exp
                // implementations need not yield identical floating-point bits.
                require(host[2 * count + i] == host[4 * count + i] &&
                            host[3 * count + i] == host[5 * count + i],
                        "Device-prepared sampling changed the secondary distribution");
                max_host_device_difference = std::max(max_host_device_difference,
                    std::abs(double(host[i])-host[2 * count + i]) /
                    std::max(1.0, double(host[2 * count + i])));
            }
            // Host-cached (rigid beam), device-cached and per-event (particles)
            // paths all satisfy the same independent physical moment bounds.
            for (int variant = 0; variant < 3; ++variant) {
                int const offset = 2 * variant * count;
                double mean = 0, second = 0, binding = 0;
                for (int i = 0; i < count; ++i) {
                    double const energy = host[offset + i];
                    require(std::isfinite(energy) && energy >= 0, "Invalid sampled energy");
                    require(i == 0 || energy >= host[offset + i - 1], "Nonmonotone inverse CDF");
                    mean += weights[i] * energy;
                    second += weights[i] * energy * energy;
                    binding += weights[i] * host[offset + count + i];
                }
                max_mean = std::max(max_mean,
                    std::abs(mean * reference.m_total / reference.m_kinetic - 1));
                max_second = std::max(max_second,
                    std::abs(second * reference.m_total / reference.m_kinetic_second - 1));
                max_binding = std::max(max_binding,
                    std::abs(binding * reference.m_total / reference.m_binding - 1));
            }
            // Check the total with the same device execution path, not a host dereference.
            amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
                output[0] = exec.crossSection(e);
                exec.sample(e, 1.0, output[1], output[2]);
                output[3] = exec.crossSection(Real(1));
                exec.sample(e, -0.1, output[4], output[5]);
                output[6] = exec.prepareSampling(e).m_secondary_endpoint;
                output[7] = exec.m_minimum_binding_energy;
            });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, sampled.begin(), sampled.begin() + 8,
                             host.begin());
            require(host[3] == 0 && host[4] == 0 && host[5] == 0, "Invalid input not rejected");
            auto const endpoint = ProtonImpactIonization::molecularMaximumSecondaryEnergy(
                e, exec.m_projectile_rest_energy, exec.m_neutral_rest_energy,
                exec.m_minimum_binding_energy);
            if (!(std::abs(host[1] / endpoint - 1) < 3e-6)) {
                SamplingState device_state;
                amrex::Gpu::copy(amrex::Gpu::deviceToHost, prepared.begin(), prepared.end(),
                                 &device_state);
                amrex::Print().SetPrecision(17)
                    << PJGModel::targetName(target) << " endpoint mismatch at E=" << e
                    << ": sample=" << host[1] << " host=" << endpoint
                    << " device=" << device_state.m_secondary_endpoint << '\n';
            }
            require(std::abs(host[1] / endpoint - 1) < 3e-6, "Molecular endpoint not preserved");
            require(host[1] == host[6] && host[2] == host[7],
                    "Unit quantile did not return the exact molecular endpoint");
            max_total = std::max(max_total, std::abs(host[0] / reference.m_total - 1));
        }
        amrex::Print() << PJGModel::targetName(target)
                       << (std::is_same_v<Real, float> ? " float" : " double")
                       << " table max relative errors: total=" << max_total << " mean=" << max_mean
                       << " second=" << max_second << " binding=" << max_binding
                       << " host/device sample difference=" << max_host_device_difference << '\n';
        require(max_total < 1e-3, "Total-table interpolation error");
        require(max_mean < 1e-3, "Mean-energy interpolation error");
        require(max_second < 2e-3, "Second-moment interpolation error");
        require(max_binding < 1e-3, "Binding-table interpolation error");
    }

    void
    checkMonoenergetic (PJGModel const& full, PJGTarget target)
    {
        auto const reference = full.executor();
        for (auto energy : {5e3, 1.5e5, 8e8, 1e10}) {
            PJGModel const mono(target, proton_mass, 5e3, 1e10, energy);
            auto const sample = mono.executor();
            require(sample.m_rows == 2, "Monoenergetic table constructed unused rows");
            constexpr int n = 4097;
            amrex::Gpu::DeviceVector<amrex::ParticleReal> output(4*n+2);
            auto* values = output.data();
            amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE(int i) noexcept {
                // Include both endpoints and clustered quantiles in the rare tail.
                double const x = static_cast<double>(i)/(n-1);
                double const a = x*x*x*x, b = (1-x)*(1-x)*(1-x)*(1-x);
                reference.sample(energy, a/(a+b), values[4*i], values[4*i+1]);
                sample.sample(energy, a/(a+b), values[4*i+2], values[4*i+3]);
                if (i == 0) {
                    values[4*n] = reference.crossSection(energy);
                    values[4*n+1] = sample.crossSection(energy);
                }
            });
            amrex::Vector<amrex::ParticleReal> host(output.size());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, output.begin(), output.end(), host.begin());
            for (int i = 0; i < n; ++i) {
                require(host[4*i] == host[4*i+2] && host[4*i+1] == host[4*i+3],
                        "Monoenergetic tables changed the PJG samples");
            }
            require(host[4*n] == host[4*n+1] &&
                        host[4*n] == mono.monoenergeticCrossSection(),
                    "Monoenergetic tables changed the total cross section");
        }
    }

    void
    checkReference ()
    {
        for (auto const invalid :
             {0.0, -1.0, 1000.0, 2.0e10, std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::quiet_NaN()}) {
            for (auto const target : {PJGTarget::N2, PJGTarget::O2}) {
                require(PJGModel::integratedCrossSection(target, invalid, proton_mass) == 0.0,
                        "Host total did not reject an unsupported incident energy");
                require(PJGModel::differentialCrossSection(target, invalid, 10.0, proton_mass) ==
                            0.0,
                        "Host SDCS did not reject an unsupported incident energy");
            }
        }
        for (auto const& row : reference_rows) {
            auto const target = row.m_target == 0 ? PJGTarget::N2 : PJGTarget::O2;
            auto const e = row.m_moments[0];
            auto const m = PJGModel::integratedMoments(target, e, proton_mass);
            std::array<double, 5> const values{m.m_total, m.m_kinetic, m.m_kinetic_second,
                                               m.m_binding, m.m_above_free};
            for (int i = 0; i < 5; ++i) {
                require(std::abs(values[i] / row.m_moments[i + 1] - 1) < 2e-7,
                        "C++/independent Python moment mismatch");
            }
            for (int i = 0; i < 6; ++i) {
                auto const value =
                    PJGModel::differentialCrossSection(target, e, row.m_secondary[i], proton_mass);
                require(std::abs(value / row.m_sdcs[i] - 1) < 2e-12,
                        "C++/independent Python SDCS mismatch");
            }
        }
    }
} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        checkReference();
        for (auto const target : {PJGTarget::N2, PJGTarget::O2}) {
            auto const start = std::chrono::steady_clock::now();
            PJGModel const model(target, proton_mass, 5e3, 1e10);
            amrex::Print()
                << PJGModel::targetName(target) << " table initialization seconds="
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                << '\n';
            checkTables<float>(model, target);
            checkTables<double>(model, target);
            checkMonoenergetic(model, target);
        }
    }
    amrex::Finalize();
}
