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
        }
    }
    amrex::Finalize();
}
