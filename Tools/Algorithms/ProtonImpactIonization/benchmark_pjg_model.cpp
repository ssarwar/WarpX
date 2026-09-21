/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PJGTestUtils.H"
#include "Source/Particles/Collision/ProtonImpactIonization/IonizationSampling.H"

#include <AMReX.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <type_traits>

namespace
{
    template <typename Real>
    void
    benchmark (PJGModel const& model)
    {
        CopiedTable<Real> table(model.executor());
        auto const exec = table.m_executor;
        constexpr int count = 1 << 20;
        amrex::Gpu::DeviceVector<Real> result(count);
        auto* output = result.data();
        std::array<double, 7> elapsed{};
        for (int repeat = -1; repeat < 7; ++repeat) {
            amrex::Gpu::streamSynchronize();
            auto const start = std::chrono::steady_clock::now();
            amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int i) noexcept {
                // Span the full table with interleaved parent energies and
                // independent low-discrepancy electron probabilities.
                auto const parent = (static_cast<std::uint32_t>(i) * 2654435761u) >> 24;
                auto const energy = Real(5e3) * std::exp(Real(14.508657738524219) *
                                                         (Real(parent) + Real(.5)) / Real(256));
                auto const probability = ProtonImpactIonization::shiftedRadicalInverse(
                    static_cast<std::uint32_t>(i), .3717626816462);
                Real secondary, binding;
                exec.sample(energy, probability, secondary, binding);
                output[i] = secondary + binding + Real(1e22) * exec.crossSection(energy);
            });
            amrex::Gpu::streamSynchronize();
            auto const duration =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (repeat >= 0) {
                elapsed[repeat] = duration;
            }
        }
        amrex::Vector<Real> host(count);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, result.begin(), result.end(), host.begin());
        double checksum = 0;
        for (auto const value : host) {
            checksum += value;
        }
        std::sort(elapsed.begin(), elapsed.end());
        amrex::Print() << (std::is_same_v<Real, float> ? "float" : "double")
                       << " median ns/(total+sample)=" << 1e9 * elapsed[3] / count
                       << " range=" << 1e9 * elapsed.front() / count << ':'
                       << 1e9 * elapsed.back() / count << " checksum=" << checksum << '\n';
    }

    template <typename Real, bool prepared>
    void
    benchmarkSharedParent (PJGModel const& model)
    {
        CopiedTable<Real> table(model.executor());
        auto const exec = table.m_executor;
        constexpr int groups = 1 << 14;
        constexpr int parents = 8;
        constexpr int products = 64;
        amrex::Gpu::DeviceVector<Real> parent_data(2 * groups * parents);
        auto* energies = parent_data.data();
        auto* scores = energies + groups * parents;
        amrex::ParallelFor(groups * parents, [=] AMREX_GPU_DEVICE(int i) noexcept {
            energies[i] = Real(8e8) * (Real(1) + Real(i % 7) / Real(100));
            scores[i] = Real(1) + Real(i % 5) / Real(10);
        });
        amrex::Gpu::DeviceVector<Real> result(groups * products);
        auto* output = result.data();
        std::array<double, 7> elapsed{};
        for (int repeat = -1; repeat < 7; ++repeat) {
            amrex::Gpu::streamSynchronize();
            auto const start = std::chrono::steady_clock::now();
            amrex::ParallelFor(groups, [=] AMREX_GPU_DEVICE(int group) noexcept {
                Real total_score = 0;
                for (int parent = 0; parent < parents; ++parent) {
                    total_score += scores[group * parents + parent];
                }
                auto const spacing = total_score / Real(products);
                Real cumulative_score = 0;
                Real energy = 0;
                int next_parent = 0;
                int selected_parent = -1;
                int previous_parent = -1;
                typename PJGModel::ExecutorT<Real>::SamplingState state;
                for (int product = 0; product < products; ++product) {
                    // Match the source's ordered, weighted parent selection.
                    // Its parent energy can change inside the product loop.
                    auto const score_target = (Real(product) + Real(.5)) * spacing;
                    while (next_parent < parents && cumulative_score <= score_target) {
                        selected_parent = next_parent++;
                        auto const parent_index = group * parents + selected_parent;
                        energy = energies[parent_index];
                        cumulative_score += scores[parent_index];
                    }
                    int const index = group * products + product;
                    auto const probability = ProtonImpactIonization::shiftedRadicalInverse(
                        static_cast<std::uint32_t>(index), .3717626816462);
                    Real secondary, binding;
                    if constexpr (prepared) {
                        if (selected_parent != previous_parent) {
                            state = exec.prepareSampling(energy);
                            previous_parent = selected_parent;
                        }
                        exec.sample(state, probability, secondary, binding);
                    } else {
                        exec.sample(energy, probability, secondary, binding);
                    }
                    output[index] = secondary + binding;
                }
            });
            amrex::Gpu::streamSynchronize();
            if (repeat >= 0) {
                elapsed[repeat] =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            }
        }
        amrex::Vector<Real> host(result.size());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, result.begin(), result.end(), host.begin());
        double checksum = 0.0;
        for (auto const value : host) { checksum += value; }
        std::sort(elapsed.begin(), elapsed.end());
        amrex::Print() << (std::is_same_v<Real, float> ? "float" : "double")
                       << (prepared ? " prepared" : " per-event")
                       << " 64 samples/8 parents median ns/sample="
                       << 1e9 * elapsed[3] / (groups * products)
                       << " checksum=" << checksum << '\n';
    }
} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        for (auto const target :
             {ProtonImpactIonization::PJGTarget::N2, ProtonImpactIonization::PJGTarget::O2}) {
            PJGModel const model(target, 938272088.16, 5e3, 1e10);
            amrex::Print() << PJGModel::targetName(target) << '\n';
            benchmark<float>(model);
            benchmark<double>(model);
            benchmarkSharedParent<float, false>(model);
            benchmarkSharedParent<float, true>(model);
            benchmarkSharedParent<double, false>(model);
            benchmarkSharedParent<double, true>(model);
        }
    }
    amrex::Finalize();
}
