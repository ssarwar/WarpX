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
#include <filesystem>
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
            if (fluid->getRigidBeam()) { config << fluid->getRigidBeam()->configuration() << '\n'; }
        }
    }
    if (!config.str().empty()) {
        bool correction = true;
        amrex::ParmParse("boundary").query("verboncoeur_axis_correction", correction);
        auto const& warpx = WarpX::GetInstance();
        config << "discretization " << static_cast<int>(WarpX::electromagnetic_solver_id) << ' '
            << static_cast<int>(warpx.evolve_scheme) << ' ' << correction << '\n';
    }
    return config.str();
}

bool
MultiFluidContainer::InitializeSelfFields () const
{
    return std::any_of(allcontainers.begin(), allcontainers.end(),
                       [](auto const& fluid) { return fluid->InitializeSelfFields(); });
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
    if (config.empty() && !std::filesystem::exists(directory+"/FluidModels")) { return; }
    amrex::Vector<char> contents;
    amrex::ParallelDescriptor::ReadAndBcastFile(directory + "/FluidModels", contents);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::string(contents.data()) == "WarpX prescribed fluids 1\n" + config,
        "The checkpoint's prescribed-fluid species, models, masses, charges or shapes or discretization changed.");
    for (auto const& fluid : allcontainers) {
        if (!fluid->isPrescribed()) { continue; }
        auto const path = directory + "/Level_0/" + fields.mf_name(fluid->name_mf_N, 0);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::VisMF::Exist(path),
            "Missing required prescribed-fluid checkpoint density: " + path);
    }
}

void
MultiFluidContainer::ValidateRestartState (ablastr::fields::MultiFabRegister& fields) const
{
    for (auto const& fluid : allcontainers) {
        if (!fluid->isPrescribed()) { continue; }
        auto const& density = *fields.get(fluid->name_mf_N, 0);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(density.min(0, density.nGrow()) >= 0.0 &&
            !density.contains_nan(0, 1, density.nGrow()) &&
            !density.contains_inf(0, 1, density.nGrow()),
            "Invalid prescribed-fluid checkpoint density: "+fluid->getName());
        if (auto* beam = fluid->getRigidBeam()) {
            auto const& warpx = WarpX::GetInstance();
            beam->UpdateCurrentDiagnostic(*fields.get("fluid_current_"+fluid->getName(),
                ablastr::fields::Direction{2}, 0), warpx.Geom(0), warpx.gett_new(0));
        }
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
MultiFluidContainer::UpdatePrescribedDensities (
    ablastr::fields::MultiFabRegister& fields, amrex::Real time)
{
    for (auto& fluid : allcontainers) {
        if (auto* beam = fluid->getRigidBeam()) {
            beam->UpdateDensity(*fields.get(fluid->name_mf_N, 0), WarpX::GetInstance().Geom(0), time);
            beam->UpdateCurrentDiagnostic(*fields.get("fluid_current_"+fluid->getName(),
                ablastr::fields::Direction{2}, 0), WarpX::GetInstance().Geom(0), time);
        }
    }
}

void
MultiFluidContainer::DepositPrescribedSources (
    ablastr::fields::MultiFabRegister& fields, amrex::Real start, amrex::Real dt,
    bool deposit_charge)
{
    using ablastr::fields::Direction;
    using warpx::fields::FieldType;
    auto const& geom = WarpX::GetInstance().Geom(0);
    for (auto& fluid : allcontainers) {
        if (!fluid->isPrescribed()) { continue; }
        auto* beam = fluid->getRigidBeam();
        if (!fluid->do_not_deposit && deposit_charge && fields.has(FieldType::rho_fp, 0)) {
            auto& rho = *fields.get(FieldType::rho_fp, 0);
            for (int comp = 0; comp < rho.nComp(); ++comp) {
                if (beam) {
                    beam->UpdateDensity(*fields.get(fluid->name_mf_N, 0), geom, start+comp*dt);
                }
                fluid->DepositCharge(fields, rho, 0, comp);
            }
        }
        if (beam) {
            // Mid-step diagnostics see the same physical time as the implicit fields.
            beam->UpdateDensity(*fields.get(fluid->name_mf_N, 0), geom, start+0.5*dt);
            if (!fluid->do_not_deposit) {
                beam->DepositCurrent(*fields.get(FieldType::current_fp, Direction{2}, 0),
                                     geom, start, dt);
            }
        }
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
