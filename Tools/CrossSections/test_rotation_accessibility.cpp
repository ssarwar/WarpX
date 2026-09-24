/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCElasticKinematics.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCThermalRotation.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
    struct Case
    {
        double m_energy, m_loss, m_cosine;
    };
    struct Result
    {
        double m_energy_error, m_angle_error;
        bool m_valid;
    };

    // Independent minimum final kinetic energy at the specified laboratory angle.
    // Numerical minimization of the on-shell energy budget does not use the
    // implementation's analytic discriminant or maximum internal energy formula.
    double
    minimumCost (Case value, double m, double mass)
    {
        auto cost = [=] (double t) {
            double const p = std::sqrt(value.m_energy * (value.m_energy + 2 * m));
            double const q = std::sqrt(t * (t + 2 * m));
            double const recoil2 = p * p + q * q - 2 * p * q * value.m_cosine;
            double const final_mass = mass + value.m_loss;
            return t + recoil2 / (std::sqrt(final_mass * final_mass + recoil2) + final_mass);
        };
        double lo = 0, hi = std::max(1., value.m_energy + std::abs(value.m_loss));
        constexpr double ratio = .6180339887498948482;
        double x = hi - ratio * (hi - lo), y = lo + ratio * (hi - lo);
        for (int i = 0; i < 160; ++i) {
            if (cost(x) < cost(y)) {
                hi = y;
                y = x;
                x = hi - ratio * (hi - lo);
            } else {
                lo = x;
                x = y;
                y = lo + ratio * (hi - lo);
            }
        }
        return std::min(cost(0), std::min(cost(x), cost(y)));
    }
    void
    checkConditionalSampler ()
    {
        using Model = BackgroundMCCThermalRotation;
        constexpr int draws = 8192;
        constexpr double c2 = PhysConst::c2_v<double>, qe = PhysConst::q_e_v<double>;
        constexpr double m = PhysConst::m_e_v<double> * c2 / qe;
        constexpr double mass = 28.0134 * 1.66053906660e-27 * c2 / qe;
        amrex::ParticleReal const loss = amrex::ParticleReal(.0014863225708172047);
        amrex::Vector<Model::Alias> ah{{.3, 1, 0}, {1, 1, 1},  {.9, 1, 2},
                                       {.9, 2, 0}, {.9, 2, 1}, {1, 2, 2}};
        amrex::Vector<double> ch{.1, .7, 1, .3, .6, 1};
        amrex::Vector<Model::Cell> sh{{0, 3}, {3, 3}};
        amrex::Vector<Model::Outcome> oh{{-loss, loss}, {0, 0}, {loss, 0}};
        amrex::Vector<amrex::ParticleReal> rh{1, 2};
        amrex::Gpu::DeviceVector<Model::Alias> aliases(ah.size());
        amrex::Gpu::DeviceVector<double> cdf(ch.size());
        amrex::Gpu::DeviceVector<Model::Cell> cells(sh.size());
        amrex::Gpu::DeviceVector<Model::Outcome> outcomes(oh.size());
        amrex::Gpu::DeviceVector<amrex::ParticleReal> rates(rh.size());
        auto upload = [] (auto const& host, auto& device) {
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
        };
        upload(ah, aliases);
        upload(ch, cdf);
        upload(sh, cells);
        upload(oh, outcomes);
        upload(rh, rates);
        Model::Executor executor;
        executor.m_loss_alias = aliases.data();
        executor.m_cdf = cdf.data();
        executor.m_cells = cells.data();
        executor.m_outcomes = outcomes.data();
        executor.m_rates = rates.data();
        executor.m_neutral_rest_energy = mass;
        double const threshold = double(loss) * (1 + m / mass) + double(loss) * loss / (2 * mass);
        amrex::Vector<Case> cases;
        for (double relative : {-1e-8, 2e-11, 1e-10, 1e-8}) {
            for (double cosine : {-1., .3, .8, 1.}) {
                cases.push_back({threshold * (1 + relative), loss, cosine});
            }
        }
        amrex::Gpu::DeviceVector<Case> inputs(cases.size());
        upload(cases, inputs);
        amrex::Gpu::DeviceVector<int> samples(draws * cases.size());
        auto const* input = inputs.data();
        auto* output = samples.data();
        amrex::Vector<int> host(samples.size());
        for (bool cumulative : {false, true}) {
            executor.m_cumulative = cumulative;
            amrex::ParallelFor(
                static_cast<int>(samples.size()), [=] AMREX_GPU_DEVICE(int i) noexcept {
                    int const j = i % draws;
                    auto const value = input[i / draws];
                    Model::Executor::Interpolation const state{0, .25, 1.25, value.m_energy};
                    unsigned reversed = 0;
                    for (int bit = 0; bit < 13; ++bit) {
                        reversed = (reversed << 1) | ((unsigned(j) >> bit) & 1u);
                    }
                    // Use Hammersley quadrature for the two independent uniforms.
                    auto const result = executor.sample(state, value.m_cosine, (j + .5) / draws,
                                                        (reversed + .5) / draws);
                    output[i] = result.m_loss < 0 ? 0 : result.m_loss == 0 ? 1 : 2;
                });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, samples.begin(), samples.end(),
                             host.begin());
            for (std::size_t k = 0; k < cases.size(); ++k) {
                bool const allowed =
                    cases[k].m_energy - cases[k].m_loss > minimumCost(cases[k], m, mass);
                double expected[3]{.75 * .1 + .5 * .3, .75 * .6 + .5 * .3,
                                   allowed ? .75 * .3 + .5 * .4 : 0};
                double const sum = expected[0] + expected[1] + expected[2];
                int count[3]{};
                for (int j = 0; j < draws; ++j) {
                    ++count[host[k * draws + j]];
                }
                for (int j = 0; j < 3; ++j) {
                    if (std::abs(double(count[j]) / draws - expected[j] / sum) > 8. / draws) {
                        std::cerr << "Case " << k << " cumulative=" << cumulative << " label=" << j
                                  << " allowed=" << allowed
                                  << " actual=" << double(count[j]) / draws
                                  << " expected=" << expected[j] / sum << '\n';
                        throw std::runtime_error("Incorrect accessibility-conditioned row mixture");
                    }
                }
                if (!allowed && count[2] != 0) {
                    throw std::runtime_error("Sampled an inaccessible rotational excitation");
                }
            }
        }
        std::cout << "PASS: alias/prefix mixtures conditioned on independent accessibility\n";
    }

} // namespace

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr double c = PhysConst::c_v<double>;
        constexpr double qe = PhysConst::q_e_v<double>;
        constexpr double me = PhysConst::m_e_v<double>;
        constexpr double mass = 28.0134 * 1.66053906660e-27;
        constexpr double m = me * c * c / qe, target = mass * c * c / qe;
        amrex::Vector<Case> cases;
        for (double loss : {.0014863225708172047, .02, .1}) {
            double const threshold = loss * (1 + m / target) + loss * loss / (2 * target);
            for (double shift : {-1e-5, -1e-8, -1e-12, 0., 1e-12, 1e-10, 1e-8, 1e-5}) {
                for (double cosine : {-1., -.3, 0., .3, .8, 1.}) {
                    // Test what the actual particle precision represents.
                    amrex::ParticleReal const u = static_cast<amrex::ParticleReal>(
                        c * std::sqrt(threshold * (1 + shift) * (threshold * (1 + shift) + 2 * m)) /
                        m);
                    double const actual = m * (double(u) * u / (c * c)) /
                                          (std::sqrt(1 + double(u) * u / (c * c)) + 1);
                    cases.push_back({actual, loss, cosine});
                }
            }
        }
        for (double energy : {0., 1e-12, .01, 1e9}) {
            for (double cosine : {-1., 0., 1.}) {
                cases.push_back({energy, -.01, cosine});
            }
        }
        amrex::Gpu::DeviceVector<Case> inputs(cases.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), inputs.begin());
        amrex::Gpu::DeviceVector<Result> outputs(cases.size());
        auto const* in = inputs.data();
        auto* out = outputs.data();
        amrex::ParallelForRNG(
            static_cast<int>(cases.size()),
            [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const& engine) noexcept {
                auto const state = in[i];
                auto const u = static_cast<amrex::ParticleReal>(
                    c * std::sqrt(state.m_energy * (state.m_energy + 2 * m)) / m);
                amrex::ParticleReal ex, ey, ez, ix, iy, iz;
                bool const valid = BackgroundMCCElasticKinematics::computeInternalEnergyChangeLab(
                    0, 0, u, 0, 0, 0, me, mass, state.m_loss, state.m_cosine, engine, ex, ey, ez,
                    ix, iy, iz);
                out[i] = {0, 0, valid};
                if (valid) {
                    double const e2 = double(ex) * ex + double(ey) * ey + double(ez) * ez;
                    double const i2 = double(ix) * ix + double(iy) * iy + double(iz) * iz;
                    double const ke = m * e2 / (c * c * (1 + std::sqrt(1 + e2 / (c * c))));
                    double const ki =
                        (target + state.m_loss) * i2 / (c * c * (1 + std::sqrt(1 + i2 / (c * c))));
                    out[i].m_energy_error =
                        std::abs(ke + ki + state.m_loss - state.m_energy) /
                        amrex::max(state.m_energy + std::abs(state.m_loss), .001);
                    out[i].m_angle_error =
                        e2 > 0 ? std::abs(double(ez) / std::sqrt(e2) - state.m_cosine) : 0;
                }
            });
        amrex::Vector<Result> result(cases.size());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, outputs.begin(), outputs.end(), result.begin());
        double const tolerance = sizeof(amrex::ParticleReal) == 4 ? 2e-6 : 1e-10;
        int restricted = 0;
        for (std::size_t i = 0; i < cases.size(); ++i) {
            auto const value = cases[i];
            double const gap = value.m_energy - value.m_loss - minimumCost(value, m, target);
            double const boundary = 1e-14 * std::max(std::abs(value.m_loss), value.m_energy);
            if (std::abs(gap) > boundary && result[i].m_valid != (gap > 0)) {
                std::cerr << "E=" << value.m_energy << " Q=" << value.m_loss
                          << " mu=" << value.m_cosine << " independent gap=" << gap
                          << " valid=" << result[i].m_valid << '\n';
                throw std::runtime_error("Incorrect angular accessibility");
            }
            if (result[i].m_energy_error > tolerance || result[i].m_angle_error > tolerance) {
                throw std::runtime_error("Laboratory-angle signed recoil changed angle or energy");
            }
            restricted += !result[i].m_valid;
        }
        if (restricted == 0) {
            throw std::runtime_error("Accessibility test did not exercise forbidden states");
        }
        std::cout
            << "PASS: numerical energy-budget accessibility, signed recoil and preserved angle\n";
    }
    checkConditionalSampler();
    amrex::Finalize();
}
