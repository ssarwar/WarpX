/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCElasticKinematics.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCThermalRotation.H"
#include "Utils/WarpXConst.H"

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
#include <memory>
#include <string>

struct Sample
{
    double loss, cosine, energy_residual, momentum_residual;
};

/** Portable sampler, conservation and throughput driver; independent references
 * are Python generated. */
int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        amrex::ParmParse pp;
        std::string reference_file, output_file;
        pp.get("reference", reference_file);
        pp.get("output", output_file);
        int count = 131072;
        pp.query("samples", count);
        bool cumulative = false;
        pp.query("cumulative", cumulative);
        AMREX_ALWAYS_ASSERT(count > 0);
        std::ifstream reference(reference_file);
        std::ofstream output(output_file);
        AMREX_ALWAYS_ASSERT(reference && output);
        output << "# T E K startup_s bytes sampling_s loss loss2 cos cos_loss "
                  "p0 p_loss p_gain cos2\n"
               << std::setprecision(17);
        std::string file;
        double temperature, energy, expected_rate;
        std::array<double, 8> expected{};
        std::array<double, 8> expected_square{};
        std::shared_ptr<BackgroundMCCThermalRotation> model;
        amrex::Gpu::DeviceVector<Sample> samples(count);
        amrex::Gpu::HostVector<Sample> host(count);
        int cases = 0;
        while (reference >> std::quoted(file) >> temperature >> energy >> expected_rate) {
            ++cases;
            for (auto& value : expected) {
                reference >> value;
            }
            for (auto& value : expected_square) {
                reference >> value;
            }
            AMREX_ALWAYS_ASSERT(reference);
            auto const start = std::chrono::steady_clock::now();
            model =
                BackgroundMCCThermalRotation::get(file, "analytic_test", temperature, cumulative);
            auto const initialized = std::chrono::steady_clock::now();
            auto const executor = model->executor();
            auto const state_host =
                model->hostExecutor().interpolate(static_cast<amrex::ParticleReal>(energy));
            AMREX_ALWAYS_ASSERT(std::abs(state_host.m_rate / expected_rate - 1) < 0.002);
            auto* data = samples.data();
            double const mass = model->neutralMass();
            auto const sample_start = std::chrono::steady_clock::now();
            amrex::ParallelForRNG(count, [=] AMREX_GPU_DEVICE(int i,
                                                              amrex::RandomEngine const& engine) {
                double const cosine = 1 - 2 * amrex::Random(engine);
                auto const state = executor.interpolate(static_cast<amrex::ParticleReal>(energy));
                auto const outcome =
                    executor.sample(state, cosine, amrex::Random(engine), amrex::Random(engine));
                data[i] = {outcome.m_loss, cosine, 0, 0};
            });
            amrex::Gpu::streamSynchronize();
            auto const sampled = std::chrono::steady_clock::now();
            // Separate timing keeps recoil work out of the lookup
            // microbenchmark.
            amrex::ParallelForRNG(
                count, [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const& engine) {
                    constexpr double c = PhysConst::c_v<double>;
                    constexpr double c2 = PhysConst::c2_v<double>;
                    constexpr double qe = PhysConst::q_e_v<double>;
                    constexpr double me = PhysConst::m_e_v<double>;
                    constexpr double rest = me * c2 / qe;
                    // Exercise the signed solver through beam energies independently
                    // of the intentionally limited analytic rotational source bundle.
                    int const high_case = (i / 16) % 3;
                    double const kinetic_energy = i % 16 == 0 ? (high_case == 0   ? 1.0e4
                                                                 : high_case == 1 ? 2.5e6
                                                                                  : 1.0e9)
                                                              : energy;
                    auto const u = static_cast<amrex::ParticleReal>(
                        c * std::sqrt(kinetic_energy * (kinetic_energy + 2 * rest)) / rest);
                    amrex::ParticleReal ex, ey, ez, nx, ny, nz;
                    double const vx = i % 3 == 0 ? 0 : 300;
                    double const vy = i % 3 == 0 ? 0 : -150;
                    double const vz = i % 3 == 0 ? 0 : 600;
                    double const v2 = vx * vx + vy * vy + vz * vz;
                    double const neutral_gamma = 1 / std::sqrt(1 - v2 / c2);
                    auto const change = data[i].loss;
                    bool const ok = BackgroundMCCElasticKinematics::computeInternalEnergyChangeLab(
                        0, 0, u, vx, vy, vz, me, mass, change, data[i].cosine, engine, ex, ey, ez,
                        nx, ny, nz);
                    if (!ok) {
                        double const incident_total = rest * std::sqrt(1 + double(u) * u / c2);
                        double const relative_energy =
                            neutral_gamma * (incident_total - vz * rest * u / c2) - rest;
                        double const target_rest = mass * c2 / qe;
                        double const threshold =
                            change * (1 + rest / target_rest) + change * change / (2 * target_rest);
                        data[i].energy_residual = relative_energy > threshold + 1.0e-8 ? 1 : 0;
                        return;
                    }
                    double const mf = mass + change * qe / c2;
                    double const ue2 = double(ex) * ex + double(ey) * ey + double(ez) * ez;
                    double const un2 = double(nx) * nx + double(ny) * ny + double(nz) * nz;
                    double const ke = me * ue2 / (qe * (1 + std::sqrt(1 + ue2 / c2)));
                    double const kn = mf * un2 / (qe * (1 + std::sqrt(1 + un2 / c2)));
                    double const incoming_ke =
                        me * double(u) * u / (qe * (1 + std::sqrt(1 + double(u) * u / c2)));
                    // Scale by the total available kinetic/internal energy,
                    // including energy released to a cold electron in a
                    // superelastic event.
                    double const neutral_ke =
                        mass * neutral_gamma * neutral_gamma * v2 / (qe * (neutral_gamma + 1));
                    data[i].energy_residual =
                        (ke + kn + change - incoming_ke - neutral_ke) /
                        amrex::max(kinetic_energy + neutral_ke + std::abs(change), .001);
                    double const px = me * ex + mf * nx - mass * neutral_gamma * vx;
                    double const py = me * ey + mf * ny - mass * neutral_gamma * vy;
                    double const pz = me * (ez - u) + mf * nz - mass * neutral_gamma * vz;
                    data[i].momentum_residual = std::sqrt(px * px + py * py + pz * pz) /
                                                (me * amrex::max(double(u), std::sqrt(ue2)) +
                                                 mass * neutral_gamma * std::sqrt(v2));
                });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, samples.begin(), samples.end(),
                             host.begin());
            std::array<double, 8> sum{}, square{};
            double max_energy_error = 0, max_momentum_error = 0;
            for (auto const& sample : host) {
                double const loss = sample.loss, mu = sample.cosine;
                std::array<double, 8> values{loss,
                                             loss * loss,
                                             mu,
                                             mu * loss,
                                             loss == 0 ? 1.0 : 0.0,
                                             loss > 0 ? 1.0 : 0.0,
                                             loss < 0 ? 1.0 : 0.0,
                                             mu * mu};
                for (int j = 0; j < 8; ++j) {
                    sum[j] += values[j];
                    square[j] += values[j] * values[j];
                }
                double const tolerance = sizeof(amrex::ParticleReal) == 4 ? 2.0e-6 : 1.0e-10;
                max_energy_error = std::max(max_energy_error, std::abs(sample.energy_residual));
                max_momentum_error = std::max(max_momentum_error, sample.momentum_residual);
                if (!std::isfinite(sample.energy_residual) ||
                    !std::isfinite(sample.momentum_residual) ||
                    std::abs(sample.energy_residual) > tolerance ||
                    sample.momentum_residual > tolerance) {
                    status = 1;
                }
            }
            if (status != 0) {
                amrex::Print() << "Conservation residuals " << temperature << ' ' << energy << ' '
                               << max_energy_error << ' ' << max_momentum_error << '\n';
            }
            output << temperature << ' ' << energy << ' ' << state_host.m_rate << ' '
                   << std::chrono::duration<double>(initialized - start).count() << ' '
                   << model->tableBytes() << ' '
                   << std::chrono::duration<double>(sampled - sample_start).count();
            for (int j = 0; j < 8; ++j) {
                double const mean = sum[j] / count;
                double const variance = std::max(square[j] / count - mean * mean,
                                                 expected_square[j] - expected[j] * expected[j]);
                double const error = std::sqrt(std::max(0.0, variance) / count);
                double const tolerance =
                    7 * error + 0.002 * std::max(std::abs(expected[j]), 1.0e-10);
                if (std::abs(mean - expected[j]) > tolerance) {
                    amrex::Print() << "moment failure " << temperature << ' ' << energy << ' ' << j
                                   << ' ' << mean << ' ' << expected[j] << ' ' << tolerance << '\n';
                    status = 1;
                }
                output << ' ' << mean;
            }
            output << '\n';
        }
        AMREX_ALWAYS_ASSERT(cases > 0 && reference.eof());
        amrex::Print() << "Thermal rotation sampler and signed recoil: "
                       << (status == 0 ? "PASS" : "FAIL") << '\n';
    }
    amrex::Finalize();
    return status;
}
