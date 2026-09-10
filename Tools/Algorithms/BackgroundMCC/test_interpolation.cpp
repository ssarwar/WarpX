/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Particles/Collision/BackgroundMCC/BackgroundMCCIonization.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCProcessSelector.H"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        using Real = amrex::ParticleReal;
        // A decreasing table must retain a positive terminal value, including
        // constant extrapolation and knots where lo + f*(hi-lo) cancels.
        amrex::Gpu::DeviceVector<Real> table(4), output(12);
        amrex::Vector<Real> host{Real(1), Real(2), Real(1.e-20), Real(1.e-37)};
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, host.begin(), host.end(), table.begin());
        ScatteringProcess::Executor process{};
        process.m_energies_data = table.data();
        process.m_sigmas_data = table.data() + 2;
        process.m_grid_size = 2;
        process.m_energy_lo = 1;
        process.m_energy_hi = 2;
        process.m_sigma_lo = host[2];
        process.m_sigma_hi = host[3];
        BackgroundMCCProcessSelector::Executor selector;
        selector.m_energies = table.data();
        selector.m_cumulative_cross_sections = table.data() + 2;
        selector.m_energy_count = 2;
        selector.m_process_count = 1;
        selector.m_energy_lo = 1;
        selector.m_energy_hi = 2;
        auto* result = output.data();
        amrex::ParallelFor(4, [=] AMREX_GPU_DEVICE(int i) noexcept {
            auto const energy = Real(i + 1);
            int selected;
            Real total;
            selector.select(energy, Real(0), selected, total);
            result[3 * i] = process.getCrossSection(energy);
            result[3 * i + 1] = total;
            result[3 * i + 2] = Real(selected);
        });
        amrex::Vector<Real> results(12);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, output.begin(), output.end(), results.begin());
        for (int i = 0; i < 4; ++i) {
            auto const expected = i == 0 ? host[2] : host[3];
            if (results[3 * i] != expected || results[3 * i + 1] != expected ||
                results[3 * i + 2] != 0) {
                throw std::runtime_error("MCC interpolation lost a positive endpoint");
            }
        }
        // A manufactured inverse table linear in the transformed coordinate
        // isolates probability precision from the empirical RBEQ spectrum.
        // Run this test with AMReX_PARTICLES_PRECISION=SINGLE as well as
        // DOUBLE.
        BackgroundMCCIonizationModel::Executor energy_model;
        energy_model.m_model = IonizationEnergySharingModel::RBEQ;
        energy_model.m_shell_count = 1;
        energy_model.m_energy_grid_size = 2;
        energy_model.m_shell_energy_grid_size = 2;
        energy_model.m_quantile_grid_size = 3;
        energy_model.m_energy_min = 1;
        energy_model.m_inverse_log_energy_step = Real(1 / std::log(2.));
        energy_model.m_binding_energies[0] = 1;
        energy_model.m_uniform_threshold_coordinates[0] = -1;
        constexpr int shells = BackgroundMCCIonizationModel::max_shell_count;
        amrex::Vector<Real> weights(2 * shells, Real(1)), inverse(6 * shells);
        for (int i = 0; i < 6 * shells; ++i) {
            inverse[i] = Real(i % 3) / 2;
        }
        amrex::Gpu::DeviceVector<Real> device_weights(weights.size()),
            device_inverse(inverse.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, weights.begin(), weights.end(),
                         device_weights.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, inverse.begin(), inverse.end(),
                         device_inverse.begin());
        energy_model.m_shell_cross_sections = device_weights.data();
        energy_model.m_inverse_cdf = device_inverse.data();
        amrex::ParallelFor(4, [=] AMREX_GPU_DEVICE(int i) noexcept {
            auto const q = 1.0 - std::pow(10.0, -8.0 - i);
            energy_model.sampleEnergy(2, 1, 0, q, result[2 * i], result[2 * i + 1]);
        });
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, output.begin(), output.end(), results.begin());
        for (int i = 0; i < 4; ++i) {
            auto const q = 1.0 - std::pow(10.0, -8.0 - i);
            auto const expected = .5 / (1 + std::pow((1 - q) / q, .25));
            if (std::abs(results[2 * i + 1] - expected) >
                    8 * std::numeric_limits<Real>::epsilon() ||
                !(results[2 * i + 1] < Real(.5))) {
                throw std::runtime_error("RBEQ inverse-CDF probability lost its hard tail");
            }
        }
        std::cout << "PASS: positive endpoints and constant extrapolation in "
                     "both MCC selectors\n";
    }
    amrex::Finalize();
}
