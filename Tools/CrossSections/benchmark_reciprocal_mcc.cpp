/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Initialization/WarpXInit.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCCollision.H"
#include "Particles/Collision/BinaryCollision/BinaryCollisionUtils.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"

#include <ablastr/utils/timer/Timer.H>

#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cmath>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

int main (int argc, char* argv[])
{
    warpx::initialization::initialize_external_libraries(argc, argv);
    {
        auto& simulation = WarpX::GetInstance();
        simulation.InitData();
        auto& particles = simulation.GetPartContainer();
        amrex::ParmParse pp("benchmark");
        int steps = 64, warmup = 4;
        bool broad = false;
        bool linear_selector = false;
        bool print_majorant = false;
        double direct_query = 0;
        pp.query("steps", steps);
        pp.query("warmup", warmup);
        pp.query("broad", broad);
        pp.query("linear_selector", linear_selector);
        pp.query("print_majorant", print_majorant);
        bool const query_directly = pp.query("direct_query", direct_query);
        if (broad) {
            auto& electrons = particles.GetParticleContainerFromName("electrons");
            for (WarpXParIter tile(electrons, 0); tile.isValid(); ++tile) {
                auto& data = tile.GetAttribs();
                auto* ux = data[PIdx::ux].dataPtr();
                auto* uy = data[PIdx::uy].dataPtr();
                auto* uz = data[PIdx::uz].dataPtr();
                amrex::ParallelForRNG(tile.numParticles(), [=] AMREX_GPU_DEVICE(
                                                               int i,
                                                               amrex::RandomEngine const& rng) {
                    constexpr double rest = 510998.95069;
                    double const energy =
                        std::exp(std::log(0.002) + amrex::Random(rng) * std::log(3e6 / 0.002));
                    ux[i] = 0;
                    uy[i] = 0;
                    uz[i] = static_cast<amrex::ParticleReal>(
                        PhysConst::c_v<double> * std::sqrt(energy * (energy + 2 * rest)) / rest);
                });
            }
        }
        amrex::Vector<std::string> names;
        amrex::ParmParse("collisions").getarr("collision_names", names);
        std::vector<std::unique_ptr<BackgroundMCCCollision>> operators;
        for (auto const& name : names) {
            if (query_directly) {
                auto processes = BinaryCollisionUtils::parse_scattering_processes(name);
                amrex::Print() << "DIRECT_CROSS_SECTION "
                               << processes.front().getCrossSection(
                                      static_cast<amrex::ParticleReal>(direct_query))
                               << '\n';
            }
            operators.push_back(std::make_unique<BackgroundMCCCollision>(
                name, linear_selector ? 0 : BackgroundMCCProcessSelector::maximum_table_bytes));
        }
        auto const dt = simulation.getdt(0);
        auto advance = [&] (int step) {
            for (auto& collision : operators) {
                collision->doCollisions(step * dt, dt, &particles);
            }
        };
        for (int step = 0; step < warmup; ++step) {
            advance(step);
        }
        if (print_majorant) {
            for (std::size_t i = 0; i < names.size(); ++i) {
                auto processes = BinaryCollisionUtils::parse_scattering_processes(names[i]);
                amrex::Print() << std::setprecision(17) << "MCC_MAJORANT " << names[i] << ' '
                               << operators[i]->get_nu_max(processes) << '\n';
            }
        }
        amrex::Gpu::synchronize();
        ablastr::utils::timer::Timer timer;
        timer.record_start_time();
        for (int step = 0; step < steps; ++step) {
            advance(step + warmup);
        }
        amrex::Gpu::synchronize();
        timer.record_stop_time();
        amrex::Print() << "MCC_OPERATOR_SECONDS " << timer.get_global_duration() << " STEPS "
                       << steps << '\n';
        std::size_t device_bytes = 0;
#ifdef AMREX_USE_GPU
        device_bytes =
            amrex::Gpu::Device::totalGlobalMem() - amrex::Gpu::Device::freeMemAvailable();
#endif
        amrex::Print() << "DEVICE_USED_BYTES " << device_bytes << '\n';
        operators.clear();
        WarpX::Finalize();
    }
    warpx::initialization::finalize_external_libraries();
}
