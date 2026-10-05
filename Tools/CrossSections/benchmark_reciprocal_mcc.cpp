/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Initialization/WarpXInit.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCCollision.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"

#include <ablastr/utils/timer/Timer.H>

#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

int main (int argc, char *argv[])
{
    warpx::initialization::initialize_external_libraries(argc, argv);
    {
        auto &simulation = WarpX::GetInstance();
        simulation.InitData();
        auto &particles = simulation.GetPartContainer();
        amrex::ParmParse pp("benchmark");
        int steps = 64, warmup = 4;
        bool broad = false;
        pp.query("steps", steps);
        pp.query("warmup", warmup);
        pp.query("broad", broad);
        if (broad) {
            auto &electrons = particles.GetParticleContainerFromName("electrons");
            for (WarpXParIter tile(electrons, 0); tile.isValid(); ++tile) {
                auto &data = tile.GetAttribs();
                auto *ux = data[PIdx::ux].dataPtr();
                auto *uy = data[PIdx::uy].dataPtr();
                auto *uz = data[PIdx::uz].dataPtr();
                amrex::ParallelForRNG(tile.numParticles(), [=] AMREX_GPU_DEVICE(
                                                               int i,
                                                               amrex::RandomEngine const &rng) {
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
        for (auto const &name : names) {
            operators.push_back(std::make_unique<BackgroundMCCCollision>(name));
        }
        auto const dt = simulation.getdt(0);
        auto advance = [&] (int step) {
            for (auto &collision : operators) {
                collision->doCollisions(step * dt, dt, &particles);
            }
        };
        for (int step = 0; step < warmup; ++step) {
            advance(step);
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
        operators.clear();
        WarpX::Finalize();
    }
    warpx::initialization::finalize_external_libraries();
}
