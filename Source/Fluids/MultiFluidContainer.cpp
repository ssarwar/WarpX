/* Copyright 2023 Grant Johnson, Remi Lehe
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

using namespace amrex;

MultiFluidContainer::MultiFluidContainer ()
{
    const ParmParse pp_fluids("fluids");
    pp_fluids.queryarr("species_names", species_names);

    const int nspecies = static_cast<int>(species_names.size());

    allcontainers.resize(nspecies);
    for (int i = 0; i < nspecies; ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::count(species_names.begin(), species_names.end(), species_names[i]) == 1,
            "Fluid species names must be unique.");
        allcontainers[i] = std::make_unique<WarpXFluidContainer>(i, species_names[i]);
    }
}

WarpXFluidContainer*
MultiFluidContainer::FindSpecies (std::string const& name) const
{
    for (auto const& fluid : allcontainers) {
        if (fluid->getName() == name) { return fluid.get(); }
    }
    return nullptr;
}

std::string
MultiFluidContainer::CheckpointConfiguration () const
{
    std::ostringstream config;
    config << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    for (auto const& fluid : allcontainers) {
        if (fluid->isPrescribed()) {
            config << fluid->getName() << ' ' << static_cast<int>(fluid->getModel()) << ' '
                   << fluid->getMass() << ' ' << fluid->getCharge() << ' ' << WarpX::nox << '\n';
        }
    }
    return config.str();
}

void
MultiFluidContainer::WriteCheckpoint (std::string const& directory) const
{
    auto const config = CheckpointConfiguration();
    if (config.empty() || !amrex::ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream output(directory + "/FluidModels");
    output << "WarpX prescribed fluids 1\n" << config;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(output.good(), "Cannot write fluid checkpoint metadata.");
}

void
MultiFluidContainer::ValidateRestart (
    std::string const& directory, ablastr::fields::MultiFabRegister const& fields) const
{
    auto const config = CheckpointConfiguration();
    if (config.empty()) { return; }
    amrex::Vector<char> contents;
    amrex::ParallelDescriptor::ReadAndBcastFile(directory + "/FluidModels", contents);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::string(contents.data()) == "WarpX prescribed fluids 1\n" + config,
        "The checkpoint's prescribed-fluid species, models, masses, charges or shapes changed.");
    for (auto const& fluid : allcontainers) {
        if (!fluid->isPrescribed()) { continue; }
        auto const path = directory + "/Level_0/" + fields.mf_name(fluid->name_mf_N, 0);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::VisMF::Exist(path),
            "Missing required prescribed-fluid checkpoint density: " + path);
    }
}

void
MultiFluidContainer::CommitDensityIncrements (ablastr::fields::MultiFabRegister& fields)
{
    for (auto const& fluid : allcontainers) {
        if (fluid->getModel() == FluidModel::Immobile) {
            fluid->CommitDensityIncrement(fields, 0);
        }
    }
}

void
MultiFluidContainer::AllocateLevelMFs (ablastr::fields::MultiFabRegister& m_fields, const BoxArray& ba, const DistributionMapping& dm, int lev)
{
    for (auto& fl : allcontainers) {
        fl->AllocateLevelMFs(m_fields, ba, dm, lev);
    }
}

void
MultiFluidContainer::InitData (
    ablastr::fields::MultiFabRegister& m_fields, amrex::Box init_box, amrex::Real cur_time, int lev,
    const amrex::Geometry& geom_lev, const amrex::Real gamma_boost, const amrex::Real beta_boost)
{
    for (auto& fl : allcontainers) {
        fl->InitData(m_fields, init_box, cur_time, lev, geom_lev, gamma_boost, beta_boost);
    }
}


void
MultiFluidContainer::DepositCharge (ablastr::fields::MultiFabRegister& m_fields, amrex::MultiFab &rho, int lev)
{
    for (auto& fl : allcontainers) {
        fl->DepositCharge(m_fields,rho,lev);
    }
}

void
MultiFluidContainer::DepositCurrent (ablastr::fields::MultiFabRegister& m_fields,
    amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, int lev)
{
    for (auto& fl : allcontainers) {
        fl->DepositCurrent(m_fields,jx,jy,jz,lev);
    }
}

void
MultiFluidContainer::Evolve (ablastr::fields::MultiFabRegister& fields,
                            int lev,
                            std::string const& current_fp_string,
                            amrex::Real cur_time,
                            bool skip_deposition)
{
    for (auto& fl : allcontainers) {
        fl->Evolve(fields, lev, current_fp_string, cur_time, skip_deposition);
    }
}
