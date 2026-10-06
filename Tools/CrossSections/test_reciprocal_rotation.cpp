/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCReciprocalRotation.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCUtils.H"
#include "Particles/Collision/ScatteringProcess.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace
{
struct Result
{
    amrex::GpuArray<double, 18> m_values;
    bool m_valid;
};
} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::ParmParse pp;
        std::string file, cross_section, output;
        pp.get("file", file);
        pp.get("cross_section", cross_section);
        pp.get("output", output);
        int samples = 32768;
        int timing_repetitions = 0;
        bool lookup_check = false;
        bool cumulative = false;
        double temperature = 300;
        pp.query("samples", samples);
        pp.query("timing_repetitions", timing_repetitions);
        pp.query("lookup_check", lookup_check);
        pp.query("cumulative", cumulative);
        pp.query("temperature", temperature);
        auto const free_before = amrex::Gpu::Device::freeMemAvailable();
        auto const start = std::chrono::steady_clock::now();
        auto model = BackgroundMCCReciprocalRotation::get(file, temperature, cumulative);
        amrex::Gpu::synchronize();
        double const load_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        auto const free_after = amrex::Gpu::Device::freeMemAvailable();
        amrex::Print() << "TABLE_LOAD_SECONDS " << load_seconds << " DEVICE_ALLOCATION_BYTES "
                       << free_before - free_after << '\n';
        auto shared = BackgroundMCCReciprocalRotation::get(file, temperature, cumulative);
        AMREX_ALWAYS_ASSERT(model == shared);
        ScatteringProcess source("elastic", cross_section, 0, ScatteringAngleModel::IAA);
        model->checkInclusiveRate(source);
        auto const executor = model->executor();
        auto const& host = model->hostExecutor();
        AMREX_ALWAYS_ASSERT(host.inRange(host.m_maximum_energy));
        auto const endpoint = static_cast<amrex::ParticleReal>(host.m_maximum_energy);
        AMREX_ALWAYS_ASSERT(host.inRange(std::nextafter(
            endpoint, std::numeric_limits<amrex::ParticleReal>::infinity())));
        AMREX_ALWAYS_ASSERT(!host.inRange(host.m_maximum_energy *
            (1 + 16 * std::numeric_limits<amrex::ParticleReal>::epsilon())));
        AMREX_ALWAYS_ASSERT(!host.inRange(1.001 * host.m_maximum_energy));
        AMREX_ALWAYS_ASSERT(!host.inRange(-1));
        AMREX_ALWAYS_ASSERT(!host.inRange(std::numeric_limits<double>::quiet_NaN()));
        std::vector<double> energy{0,     1e-9,  0.0001, 0.001, 0.0015, 0.003,
                                   0.01,  0.025, 0.1,    1,     1.25,   2.22,
                                   2.47,  10,    20,     100,   205,    211.04623958760578,
                                   500,   999,   1000,   6000,  8000,   9000,
                                   10000, 1e6,   2.5e6,  1e9};
        amrex::Vector<Result> result(samples);
        amrex::Gpu::DeviceVector<Result> device(samples);
        auto* values = device.data();
        // Binary32 uniforms in [1/2,1) all lie on the 2^-24 grid. Verify that
        // the actual device RNG used for rare acceptance and outcome tails
        // resolves that interval more finely, including in all-single builds.
        amrex::ParallelForRNG(samples, [=] AMREX_GPU_DEVICE (
            int i, amrex::RandomEngine const& rng) noexcept {
            double const draw = BackgroundMCCUtils::uniformDouble(rng);
            double const coordinate = draw * 16777216.0;
            values[i].m_valid = draw >= 0 && draw < 1;
            values[i].m_values[0] = draw >= .5 && coordinate != std::floor(coordinate) ? 1 : 0;
        });
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), result.begin());
        int resolved = 0;
        for (auto const& value : result) {
            AMREX_ALWAYS_ASSERT(value.m_valid);
            resolved += int(value.m_values[0]);
        }
        if (samples >= 4096) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(resolved > .4 * samples,
                "Reciprocal sampling uniforms lost resolution in the upper half interval.");
        }
        std::ofstream stream(output);
        stream << std::setprecision(17);
        for (double e : energy) {
            if (!host.inRange(e)) {
                continue;
            }
            auto draw_samples = [&] {
                amrex::ParallelForRNG(samples, [=] AMREX_GPU_DEVICE(
                                                   int i, amrex::RandomEngine const& rng) {
                    auto const state = executor.interpolate(e);
                    auto const draw = executor.sample(state, BackgroundMCCUtils::uniformDouble(rng),
                                                      BackgroundMCCUtils::uniformDouble(rng));
                    double const loss = draw.m_outcome.m_loss;
                    double const base[6]{loss == 0 ? 1.0 : 0.0,  loss > 0 ? 1.0 : 0.0,
                                         loss < 0 ? 1.0 : 0.0,   amrex::max(loss, 0.0),
                                         amrex::max(-loss, 0.0), loss * loss};
                    for (int k = 0; k < 6; ++k) {
                        values[i].m_values[3 * k] = base[k];
                        values[i].m_values[3 * k + 1] = base[k] * draw.m_deflection;
                        values[i].m_values[3 * k + 2] =
                            base[k] * draw.m_deflection * draw.m_deflection;
                    }
                    values[i].m_valid =
                        draw.m_valid && draw.m_deflection >= 0 && draw.m_deflection <= 2;
                });
            };
            draw_samples();
            if (timing_repetitions > 0) {
                amrex::Gpu::synchronize();
                auto const begin = std::chrono::steady_clock::now();
                for (int repeat = 0; repeat < timing_repetitions; ++repeat) {
                    draw_samples();
                }
                amrex::Gpu::synchronize();
                double const seconds =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
                amrex::Print() << "SAMPLER_SECONDS " << seconds << " ENERGY_EV " << e << " SAMPLES "
                               << double(samples) * timing_repetitions << '\n';
            }
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(),
                             result.begin());
            std::array<double, 18> sum{}, square{};
            for (auto const& r : result) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(r.m_valid, "Invalid reciprocal sample.");
                for (int k = 0; k < 18; ++k) {
                    sum[k] += r.m_values[k];
                    square[k] += r.m_values[k] * r.m_values[k];
                }
            }
            stream << e << ' ' << host.interpolate(e).m_rate;
            for (int k = 0; k < 18; ++k) {
                double const mean = sum[k] / samples;
                double const variance = std::max(0.0, square[k] / samples - mean * mean);
                stream << ' ' << mean << ' ' << std::sqrt(variance / samples);
            }
            stream << '\n';
            if (lookup_check && (executor.m_angular_lookup || executor.m_conditional_lookup)) {
                auto reference = executor;
                reference.m_angular_lookup = nullptr;
                reference.m_conditional_lookup = nullptr;
                amrex::ParallelForRNG(samples, [=] AMREX_GPU_DEVICE (
                    int i, amrex::RandomEngine const& rng) noexcept {
                    auto const state = executor.interpolate(e);
                    double const angle = BackgroundMCCUtils::uniformDouble(rng);
                    double const outcome = BackgroundMCCUtils::uniformDouble(rng);
                    auto const fast = executor.sample(state, angle, outcome);
                    auto const slow = reference.sample(state, angle, outcome);
                    values[i].m_valid = fast.m_valid == slow.m_valid &&
                        fast.m_deflection == slow.m_deflection &&
                        fast.m_outcome.m_loss == slow.m_outcome.m_loss &&
                        fast.m_outcome.m_initial_energy == slow.m_outcome.m_initial_energy;
                });
                amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), result.begin());
                for (auto const& value : result) {
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(value.m_valid,
                        "Quantile lookup changed a sampled event.");
                }
            }
        }
        if (lookup_check) {
            amrex::Print() << "PASS: indexed and full searches agree for identical uniforms.\n";
        }
        amrex::Print() << "Reciprocal tables: " << model->tableBytes()
                       << " bytes; valid samples and shared storage verified.\n";
    }
    amrex::Finalize();
}
