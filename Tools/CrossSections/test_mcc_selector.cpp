/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCProcessSelector.H"
#include "Particles/Collision/ScatteringProcess.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        using Real = amrex::ParticleReal;
        int points = 4096, repeats = 0;
        amrex::ParmParse pp;
        pp.query("points", points);
        pp.query("repeats", repeats);
        AMREX_ALWAYS_ASSERT(points >= 64 && repeats >= 0);
        amrex::Vector<ScatteringProcess> processes;
        for (int channel = 0; channel < 3; ++channel) {
            auto const file = "selector-grid-" + std::to_string(channel) + ".txt";
            std::ofstream table(file);
            table << std::setprecision(std::numeric_limits<Real>::max_digits10);
            table << "0 " << (channel == 0 ? 1.e-20 : 0.) << '\n';
            for (int i = 0; i < points; ++i) {
                // Offset grids, abrupt zero intervals, a very rare channel,
                // and a long range of exponents exercise the union table.
                Real const energy = Real(std::exp2(-20. + (50. * i + .125 * channel) / points));
                double const sigma = channel == 0 ? 1.e-20 / (1 + energy)
                    : channel == 1 ? (i % 31 == 0 ? 1.e-21 : 0.) : 1.e-30;
                table << energy << ' ' << sigma << '\n';
            }
            table.close();
            processes.emplace_back("elastic_" + std::to_string(channel), file, 0,
                                    ScatteringAngleModel::Isotropic);
        }
        BackgroundMCCProcessSelector selector(processes);
        auto const fast = selector.executor();
        AMREX_ALWAYS_ASSERT(fast.m_energy_lookup != nullptr);
        AMREX_ALWAYS_ASSERT(selector.tableBytes() <= selector.maximum_table_bytes);
        auto reference = fast;
        reference.m_energy_lookup = nullptr;

        auto const& grid = selector.energyGrid();
        std::vector<Real> queries;
        auto add_neighbors = [&queries] (Real energy) {
            if (std::isfinite(energy) && energy >= 0) {
                queries.push_back(energy);
                queries.push_back(std::nextafter(energy, Real(0)));
                queries.push_back(std::nextafter(energy, std::numeric_limits<Real>::infinity()));
            }
        };
        for (auto energy : grid) { add_neighbors(energy); }
        // Include every search-bin boundary and both adjacent floating-point values.
        for (int exponent = -24; exponent <= 32; ++exponent) {
            for (int bin = 0; bin < fast.energy_lookup_bins; ++bin) {
                add_neighbors(std::ldexp(Real(1) + Real(bin) / fast.energy_lookup_bins, exponent));
            }
        }
        add_neighbors(std::numeric_limits<Real>::min());
        add_neighbors(std::numeric_limits<Real>::denorm_min());
        add_neighbors(std::numeric_limits<Real>::max());
        std::mt19937 generator(37829);
        for (int i = 0; i < 65536; ++i) {
            queries.push_back(Real(std::exp2(-24. + 60. * generator() / generator.max())));
        }
        // The standard-library reference does not use the production bisection.
        std::vector<int> expected(queries.size());
        for (std::size_t i = 0; i < queries.size(); ++i) {
            expected[i] = std::clamp(int(std::upper_bound(grid.begin(), grid.end(), queries[i]) -
                                         grid.begin()) - 1, 0, int(grid.size()) - 2);
        }
        amrex::Gpu::DeviceVector<Real> device_queries(queries.size());
        amrex::Gpu::DeviceVector<int> device_expected(expected.size()), errors(queries.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, queries.begin(), queries.end(),
                         device_queries.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, expected.begin(), expected.end(),
                         device_expected.begin());
        auto const* energy = device_queries.data();
        auto const* interval = device_expected.data();
        auto* error = errors.data();
        int const count = static_cast<int>(queries.size());
        amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (int i) noexcept {
            auto const actual = fast.interpolate(energy[i]);
            auto const original = reference.interpolate(energy[i]);
            bool matches = actual.m_energy_index == interval[i] &&
                actual.m_energy_index == original.m_energy_index &&
                actual.m_energy_fraction == original.m_energy_fraction &&
                actual.m_total == original.m_total;
            double const fractions[]{0., 1.e-10, .25, .5, .9999999999, 1.};
            for (double fraction : fractions) {
                double const draw = fraction * original.m_total;
                matches = matches && fast.select(actual, draw) == reference.select(original, draw);
            }
            error[i] = matches ? 0 : 1;
        });
        std::vector<int> host_errors(queries.size());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                         errors.begin(), errors.end(), host_errors.begin());
        for (std::size_t i = 0; i < queries.size(); ++i) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(host_errors[i] == 0,
                "Bounded MCC search changed an interval, interpolation or selected channel.");
        }
        BackgroundMCCProcessSelector disabled(processes, 0);
        AMREX_ALWAYS_ASSERT(!disabled.enabled());
        std::cout << "PASS: " << queries.size()
                  << " MCC queries identical with bounded and full searches; "
                  << "table bytes=" << selector.tableBytes() << '\n';

        // Optional timings compare the same query distribution, outputs and
        // production executors. They are not CTest acceptance criteria.
        amrex::Gpu::DeviceVector<double> values(queries.size());
        auto* value = values.data();
        for (int trial = 0; trial < 6 && repeats > 0; ++trial) {
            bool const bounded = trial % 2 != 0;
            auto const executor = bounded ? fast : reference;
            amrex::Gpu::synchronize();
            double const start = amrex::second();
            for (int repeat = 0; repeat < repeats; ++repeat) {
                amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE (int i) noexcept {
                    auto const state = executor.interpolate(energy[i]);
                    value[i] = state.m_total + executor.select(state, .37 * state.m_total);
                });
            }
            amrex::Gpu::synchronize();
            std::cout << "selector_benchmark bounded=" << bounded << " trial=" << trial
                      << " seconds=" << amrex::second() - start << '\n';
        }
    }
    amrex::Finalize();
}
