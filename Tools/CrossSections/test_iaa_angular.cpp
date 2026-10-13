/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCElasticKinematics.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCElasticScattering.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCKinematics.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCProcessSelector.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCUtils.H"
#include "Particles/Collision/ScatteringProcess.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::ParmParse pp;
        std::string file, reference;
        pp.get("file", file);
        pp.get("reference", reference);
        int count = 65536, repeats = 3, relaxation = 0;
        pp.query("samples", count);
        pp.query("repeats", repeats);
        pp.query("relaxation", relaxation);
        AMREX_ALWAYS_ASSERT(count > 0 && repeats > 0);
        double const start = amrex::second();
        auto model = BackgroundMCCElasticScatteringModel::get(file);
        auto const angular = model->executor();
        amrex::Gpu::synchronize();
        double const initialization = amrex::second() - start;
        amrex::Gpu::DeviceVector<double> device(5 * count);
        amrex::Vector<double> host(5 * count);
        auto* output = device.data();
        std::string native_grid;
        if (pp.query("native_grid", native_grid)) {
            ScatteringProcess attachment("attachment", native_grid, 0,
                                         ScatteringAngleModel::Isotropic);
            attachment.setCrossSectionMultiplier(2e25);
            auto const process = attachment.executor();
            std::string large, small;
            pp.get("large_rate", large);
            pp.get("small_rate", small);
            amrex::Vector<ScatteringProcess> processes;
            processes.emplace_back("excitation_large", large, 0, ScatteringAngleModel::Isotropic);
            processes.emplace_back("excitation_small", small, 0, ScatteringAngleModel::Isotropic);
            BackgroundMCCProcessSelector selector(processes);
            auto const select = selector.executor();
            double const draw = double(processes[0].getCrossSectionGrid()[0]) +
                                0.5 * double(processes[1].getCrossSectionGrid()[0]);
            amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
                output[0] = process.getCrossSection(amrex::ParticleReal(10));
                output[1] = select.select(select.interpolate(amrex::ParticleReal(10)), draw);
            });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), host.begin());
            AMREX_ALWAYS_ASSERT(std::abs(host[0] / 4e-25 - 1) <
                                (sizeof(amrex::ParticleReal) == 4 ? 5e-7 : 1e-12));
            AMREX_ALWAYS_ASSERT(host[1] == 1);
        }
        std::ifstream input(reference);
        AMREX_ALWAYS_ASSERT(input.good());
        std::string line;
        double seconds = 0.0;
        int rows = 0;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            std::istringstream stream(line);
            double energy;
            std::array<double, 5> expected{};
            AMREX_ALWAYS_ASSERT(stream >> energy);
            for (auto& value : expected) {
                AMREX_ALWAYS_ASSERT(stream >> value);
            }
            std::vector<double> expected_cdf;
            double cumulative;
            while (stream >> cumulative) {
                expected_cdf.push_back(cumulative);
            }
            for (int repeat = -1; repeat < repeats; ++repeat) {
                double const begin = amrex::second();
                amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int i) noexcept {
                    double const u = (i + 0.5) / count;
                    double const mu =
                        1 - angular.sampleDeflection(static_cast<amrex::ParticleReal>(energy), u);
                    output[5 * i] = mu;
                    output[5 * i + 1] = (3 * mu * mu - 1) / 2;
                    output[5 * i + 2] = (5 * mu * mu * mu - 3 * mu) / 2;
                    output[5 * i + 3] = (35 * mu * mu * mu * mu - 30 * mu * mu + 3) / 8;
                    output[5 * i + 4] = mu < 0 ? 1.0 : 0.0;
                });
                amrex::Gpu::synchronize();
                if (repeat >= 0) {
                    seconds += amrex::second() - begin;
                }
            }
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), host.begin());
            std::array<long double, 5> measured{};
            for (int i = 0; i < count; ++i) {
                for (int j = 0; j < 5; ++j) {
                    measured[j] += host[5 * i + j] / count;
                }
            }
            for (int j = 0; j < 5; ++j) {
                if (std::abs(measured[j] - expected[j]) >= 1e-3) {
                    std::cerr << "Angular mismatch E=" << energy << " moment=" << j
                              << " measured=" << measured[j] << " expected=" << expected[j] << '\n';
                    amrex::Abort("IAA angular sampler failed independent quadrature.");
                }
            }
            if (!expected_cdf.empty()) {
                AMREX_ALWAYS_ASSERT(expected_cdf.size() >= 2);
                std::vector<double> cosines(count);
                for (int i = 0; i < count; ++i) {
                    cosines[i] = host[5 * i];
                }
                std::sort(cosines.begin(), cosines.end());
                for (std::size_t j = 0; j < expected_cdf.size(); ++j) {
                    double const cosine =
                        std::cos(double(j) * MathConst::pi / double(expected_cdf.size() - 1));
                    double const observed =
                        double(cosines.end() -
                               std::lower_bound(cosines.begin(), cosines.end(), cosine)) /
                        count;
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(observed - expected_cdf[j]) < 1e-3,
                                                     "IAA angular CDF differs from the independent "
                                                     "reference.");
                }
            }
            ++rows;
        }
        AMREX_ALWAYS_ASSERT(rows > 0);
        if (relaxation) {
            constexpr int steps = 24;
            constexpr double probability = 0.125;
            amrex::ParallelForRNG(
                count, [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const& engine) noexcept {
                    using namespace BackgroundMCCKinematics;
                    Vector3 direction{0, 0, 1};
                    for (int step = 0; step < steps; ++step) {
                        if (BackgroundMCCUtils::uniformDouble(engine) < probability) {
                            double const d = angular.sampleDeflection(
                                amrex::ParticleReal(2), BackgroundMCCUtils::uniformDouble(engine));
                            Vector3 a, b;
                            transverseDirections(direction, a, b);
                            direction = directionFromPolarAngle(
                                direction, a, b, 1 - d,
                                2 * MathConst::pi * BackgroundMCCUtils::uniformDouble(engine), d);
                        }
                    }
                    output[5 * i] = direction.z;
                    output[5 * i + 1] = (3 * direction.z * direction.z - 1) / 2;
                });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), host.begin());
            double const a1 = std::pow(1 - probability, steps);
            double const a2 = std::pow(1 - probability * (1 - 1.0 / 7), steps);
            double const a4 = std::pow(1 - probability * (1 - 8.0 / 63), steps);
            double const variance1 = 1.0 / 3 + 2 * a2 / 3 - a1 * a1;
            double const variance2 = 1.0 / 5 + 2 * a2 / 7 + 18 * a4 / 35 - a2 * a2;
            double m1 = 0, m2 = 0;
            for (int i = 0; i < count; ++i) {
                m1 += host[5 * i] / count;
                m2 += host[5 * i + 1] / count;
            }
            AMREX_ALWAYS_ASSERT(std::abs(m1 - a1) < 6 * std::sqrt(variance1 / count) + 1e-3);
            AMREX_ALWAYS_ASSERT(std::abs(m2 - a2) < 6 * std::sqrt(variance2 / count) + 1e-3);
        }
        // Preserve a deflection far below binary32 cosine resolution through
        // the actual excitation kinematics, including at relativistic energies.
        amrex::ParallelForRNG(
            1, [=] AMREX_GPU_DEVICE(int, amrex::RandomEngine const& engine) noexcept {
                constexpr double d = 1e-20;
                double const incoming =
                    BackgroundMCCKinematics::properSpeedFromKineticEnergy(1e9, PhysConst::m_e);
                amrex::ParticleReal x, y, z, rx, ry, rz;
                BackgroundMCCElasticKinematics::computeExcitation(
                    0, 0, static_cast<amrex::ParticleReal>(incoming), 0, 0, 0, PhysConst::m_e,
                    28.0134 * 1.66053906660e-27, 7, 1 - d, engine, x, y, z, rx, ry, rz, d);
                output[0] = (double(x) * x + double(y) * y) / (double(z) * z);
            });
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device.begin(), device.end(), host.begin());
        AMREX_ALWAYS_ASSERT(std::abs(host[0] / 2e-20 - 1) < 1e-5);
        std::cout << std::setprecision(12)
                  << "IAA_ANGULAR initialization_seconds=" << initialization
                  << " table_bytes=" << model->tableBytes() << " sample_seconds=" << seconds
                  << " samples=" << rows * double(count) * repeats << " PASS\n";
    }
    amrex::Finalize();
}
