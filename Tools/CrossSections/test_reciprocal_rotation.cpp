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

int main (int argc, char *argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::ParmParse pp;
        std::string file, cross_section, output;
        pp.get("file", file);
        pp.get("cross_section", cross_section);
        pp.get("output", output);
        int samples = 32768;
        bool cumulative = false;
        double temperature = 300;
        pp.query("samples", samples);
        pp.query("cumulative", cumulative);
        pp.query("temperature", temperature);
        auto model = BackgroundMCCReciprocalRotation::get(file, temperature, cumulative);
        auto shared = BackgroundMCCReciprocalRotation::get(file, temperature, cumulative);
        AMREX_ALWAYS_ASSERT(model == shared);
        ScatteringProcess source("elastic", cross_section, 0, ScatteringAngleModel::IAA);
        model->checkInclusiveRate(source);
        auto const executor = model->executor();
        auto const &host = model->hostExecutor();
        AMREX_ALWAYS_ASSERT(host.inRange(host.m_maximum_energy));
        AMREX_ALWAYS_ASSERT(!host.inRange(1.001 * host.m_maximum_energy));
        AMREX_ALWAYS_ASSERT(!host.inRange(-1));
        std::vector<double> energy{0,     1e-9,  0.0001, 0.001, 0.0015, 0.003,
                                   0.01,  0.025, 0.1,    1,     1.25,   2.22,
                                   2.47,  10,    20,     100,   205,    211.04623958760578,
                                   500,   999,   1000,   6000,  8000,   9000,
                                   10000, 1e6,   2.5e6,  1e9};
        amrex::Vector<Result> result(samples);
        amrex::Gpu::DeviceVector<Result> device(samples);
        auto *values = device.data();
        std::ofstream stream(output);
        stream << std::setprecision(17);
        for (double e : energy) {
            if (!host.inRange(e)) {
                continue;
            }
            amrex::ParallelForRNG(samples, [=] AMREX_GPU_DEVICE(int i,
                                                                amrex::RandomEngine const &rng) {
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
                    values[i].m_values[3 * k + 2] = base[k] * draw.m_deflection * draw.m_deflection;
                }
                values[i].m_valid =
                    draw.m_valid && draw.m_deflection >= 0 && draw.m_deflection <= 2;
            });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(),
                             result.begin());
            std::array<double, 18> sum{}, square{};
            for (auto const &r : result) {
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
        }
        amrex::Print() << "Reciprocal tables: " << model->tableBytes()
                       << " bytes; valid samples and shared storage verified.\n";
    }
    amrex::Finalize();
}
