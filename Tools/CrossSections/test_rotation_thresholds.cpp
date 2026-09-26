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
        double m_mass = 28.0134 * 1.66053906660e-27;
    };
    struct Result
    {
        double m_energy_error, m_angle_error, m_outgoing, m_momentum_error;
        bool m_valid;
    };

    // Independent scalar energy-budget solve, with no implementation quadratic.
    // Above the continuation band, the larger physical root lies in [0,E-loss].
    double
    outgoingReference (Case value, double m, double mass)
    {
        double const p = std::sqrt(value.m_energy * (value.m_energy + 2 * m));
        double const final_mass = mass + value.m_loss;
        double lo = 0, hi = value.m_energy - value.m_loss;
        for (int i = 0; i < 100; ++i) {
            double const t = (lo + hi) / 2;
            double const q = std::sqrt(t * (t + 2 * m));
            double const recoil2 = p * p + q * q - 2 * p * q * value.m_cosine;
            double const recoil =
                recoil2 / (std::sqrt(final_mass * final_mass + recoil2) + final_mass);
            if (t + recoil < value.m_energy - value.m_loss) {
                lo = t;
            } else {
                hi = t;
            }
        }
        return (lo + hi) / 2;
    }
    void
    checkIndependentSampler ()
    {
        using Model = BackgroundMCCThermalRotation;
        constexpr int draws = 8192;
        double const loss = .0014863225708172047;
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
        auto upload = [] (auto const &host, auto &device) {
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

        amrex::Vector<Case> cases;
        for (double relative : {-1e-8, 2e-11, 1e-10, 1e-8}) {
            for (double cosine : {-1., .3, .8, 1.}) {
                cases.push_back({loss * (1 + relative), loss, cosine});
            }
        }
        amrex::Gpu::DeviceVector<Case> inputs(cases.size());
        upload(cases, inputs);
        amrex::Gpu::DeviceVector<int> samples(draws * cases.size());
        auto const *input = inputs.data();
        auto *output = samples.data();
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
                    auto const result =
                        executor.sample(state, (j + .5) / draws, (reversed + .5) / draws);
                    output[i] = result.m_loss < 0 ? 0 : result.m_loss == 0 ? 1 : 2;
                });
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, samples.begin(), samples.end(),
                             host.begin());
            for (std::size_t k = 0; k < cases.size(); ++k) {
                bool const allowed = cases[k].m_energy >= cases[k].m_loss;
                double expected[3]{.75 * .1 + .5 * .3, .75 * .6 + .5 * .3,
                                   allowed ? .75 * .3 + .5 * .4 : 0};
                if (!allowed) {
                    expected[1] += .75 * .3 + .5 * .4;
                }
                double const sum = 1.25;
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
                        throw std::runtime_error("Incorrect energy-only threshold or row mixture");
                    }
                }
                if (!allowed && count[2] != 0) {
                    throw std::runtime_error("Sampled a subthreshold rotational excitation");
                }
            }
        }
        std::cout
            << "PASS: independent aliases, strict thresholds and unchanged gain probabilities\n";
    }

} // namespace

int
main (int argc, char *argv[])
{
    amrex::Initialize(argc, argv);
    {
        constexpr double c = PhysConst::c_v<double>;
        constexpr double qe = PhysConst::q_e_v<double>;
        constexpr double me = PhysConst::m_e_v<double>;
        constexpr double m = me * c * c / qe;
        amrex::Vector<Case> cases;
        for (double mass : {28.0134 * 1.66053906660e-27, 31.9988 * 1.66053906660e-27}) {
            for (double loss : {.0014863225708172047, .0017828927734694198, .02, .1}) {
                for (double shift :
                     {-1e-5, -1e-8, -1e-12, 0., 1e-12, 1e-8, 1e-5, 2e-5, 3e-5, 1e-4, 1e-2}) {
                    for (double cosine : {-1., -.3, 0., .3, .8, 1.}) {
                        // Compare the actual energy represented by particle precision.
                        amrex::ParticleReal const u = static_cast<amrex::ParticleReal>(
                            c * std::sqrt(loss * (1 + shift) * (loss * (1 + shift) + 2 * m)) / m);
                        double const actual = m * (double(u) * u / (c * c)) /
                                              (std::sqrt(1 + double(u) * u / (c * c)) + 1);
                        cases.push_back({actual, loss, cosine, mass});
                    }
                }
            }
            for (double energy : {0., 1e-12, .01, 1e4, 2.5e6, 1e9}) {
                for (double loss : {-.01, 0.}) {
                    for (double cosine : {-1., 0., 1.}) {
                        cases.push_back({energy, loss, cosine, mass});
                    }
                }
            }
        }
        amrex::Gpu::DeviceVector<Case> inputs(cases.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, cases.begin(), cases.end(), inputs.begin());
        amrex::Gpu::DeviceVector<Result> outputs(cases.size());
        auto const *in = inputs.data();
        auto *out = outputs.data();
        amrex::ParallelForRNG(
            static_cast<int>(cases.size()),
            [=] AMREX_GPU_DEVICE(int i, amrex::RandomEngine const &engine) noexcept {
                auto const state = in[i];
                auto const u = static_cast<amrex::ParticleReal>(
                    c * std::sqrt(state.m_energy * (state.m_energy + 2 * m)) / m);
                amrex::ParticleReal ex, ey, ez, ix, iy, iz;
                double const target = state.m_mass * c * c / qe;
                bool const valid = BackgroundMCCElasticKinematics::computeRotation(
                    0, 0, u, 0, 0, 0, me, state.m_mass, state.m_loss, state.m_cosine, engine, ex,
                    ey, ez, ix, iy, iz);
                out[i] = {0, 0, 0, 0, valid};
                if (valid) {
                    double const e2 = double(ex) * ex + double(ey) * ey + double(ez) * ez;
                    double const i2 = double(ix) * ix + double(iy) * iy + double(iz) * iz;
                    double const ke = m * e2 / (c * c * (1 + std::sqrt(1 + e2 / (c * c))));
                    double const ki =
                        (target + state.m_loss) * i2 / (c * c * (1 + std::sqrt(1 + i2 / (c * c))));
                    out[i].m_outgoing = ke;
                    double const factor = (target + state.m_loss) / m;
                    double const px = ex + factor * ix, py = ey + factor * iy;
                    double const pz = ez + factor * iz - u;
                    out[i].m_momentum_error = std::sqrt(px * px + py * py + pz * pz) /
                                              amrex::max(double(u) + std::sqrt(e2), 1.0);
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
        int restricted = 0, continued = 0;
        for (std::size_t i = 0; i < cases.size(); ++i) {
            auto const value = cases[i];
            double const target = value.m_mass * c * c / qe;
            double const available = value.m_energy - value.m_loss;
            // Avoid classifying the double-precision velocity reconstruction's
            // final ulp as a threshold violation; both neighboring sides are tested.
            if (std::abs(available) > 8 * std::numeric_limits<double>::epsilon() *
                                          std::max(value.m_energy, std::abs(value.m_loss)) &&
                result[i].m_valid != (available >= 0)) {
                throw std::runtime_error("Incorrect nominal rotational threshold");
            }
            restricted += !result[i].m_valid;
            if (!result[i].m_valid) {
                continue;
            }
            bool const continuation =
                value.m_loss > 0 && available <= 2 * m / target * value.m_energy;
            double const expected =
                continuation ? std::max(available, 0.) : outgoingReference(value, m, target);
            double const scale = std::max(value.m_energy + std::abs(value.m_loss), .001);
            double const energy_tolerance = continuation ? 2 * m / target + tolerance : tolerance;
            if (!std::isfinite(result[i].m_outgoing) || result[i].m_outgoing < 0 ||
                std::abs(result[i].m_outgoing - expected) / scale > tolerance ||
                result[i].m_energy_error > energy_tolerance ||
                result[i].m_momentum_error > tolerance || result[i].m_angle_error > tolerance) {
                std::cerr << "E=" << value.m_energy << " loss=" << value.m_loss
                          << " cosine=" << value.m_cosine << " continued=" << continuation
                          << " outgoing=" << result[i].m_outgoing << " expected=" << expected
                          << " energy_error=" << result[i].m_energy_error
                          << " angle_error=" << result[i].m_angle_error << '\n';
                throw std::runtime_error("Independent-angle recoil reference mismatch");
            }
            continued += continuation;
        }
        if (restricted == 0 || continued == 0) {
            throw std::runtime_error("Threshold test missed forbidden or continuation cases");
        }
        std::cout
            << "PASS: nominal thresholds, independent angles, recoil and bounded continuation\n";
    }
    checkIndependentSampler();
    amrex::Finalize();
}
