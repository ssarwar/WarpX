/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PrescribedSourceBudget.H"

#include "Fluids/MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "WarpX.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <array>
#include <fstream>

PrescribedSourceBudget::PrescribedSourceBudget (std::string const& name) : ReducedDiags(name)
{
    auto& warpx = WarpX::GetInstance();
    amrex::Vector<std::string> collisions;
    amrex::ParmParse("collisions").queryarr("collision_names", collisions);
    for (auto const& collision : collisions) {
        amrex::ParmParse pp(collision);
        std::string type;
        pp.query("type", type);
        if (type != "proton_impact_ionization" || !warpx.DoFluidSpecies()) {
            continue;
        }
        amrex::Vector<std::string> species;
        pp.getarr("species", species);
        auto const* fluid = warpx.GetFluidContainer().FindSpecies(species.at(0));
        if (fluid && fluid->getRigidBeam()) {
            m_sources.push_back(collision);
        }
    }
    m_data.resize(5 * m_sources.size(), 0.0);
    if (amrex::ParallelDescriptor::IOProcessor() && m_write_header) {
        std::ofstream output(m_path + m_rd_name + "." + m_extension);
        output << "#[0]step()" << m_sep << "[1]time(s)";
        int column = 2;
        std::array<std::string, 5> const names{"pending()", "emitted()", "electron_energy(J)",
                                               "binding_energy(J)", "discarded_ion_energy(J)"};
        for (auto const& source : m_sources) {
            for (auto const& quantity : names) {
                output << m_sep << '[' << column++ << ']' << source << '_' << quantity;
            }
        }
        output << '\n';
    }
}

void
PrescribedSourceBudget::ComputeDiags (int step)
{
    if (!m_intervals.contains(step + 1)) {
        return;
    }
    auto const& fields = WarpX::GetInstance().m_fields;
    for (std::size_t source = 0; source < m_sources.size(); ++source) {
        m_data[5 * source] = fields.get(m_sources[source] + "_product_weight_remainder", 0)->sum(0);
        auto const& budget = *fields.get(m_sources[source] + "_source_budget", 0);
        for (int comp = 0; comp < 4; ++comp) {
            // These cell-centered states contain extensive numbers and energies,
            // already integrated over each cylindrical cell, not densities.
            m_data[5 * source + comp + 1] = budget.sum(comp);
        }
    }
}
