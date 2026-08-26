/* Copyright 2021-2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Modern Electron, Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ScatteringProcess.H"

#include "Utils/TextMsg.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

ScatteringProcess::ScatteringProcess (
                        const std::string& scattering_process,
                        const std::string& cross_section_file,
                        const amrex::ParticleReal energy,
                        const ScatteringAngleModel scattering_angle_model )
{
    // Retain the unscaled values in double precision. Three-body attachment
    // tables in m^5 can be much smaller than the single-precision range.
    readCrossSectionFileRaw(cross_section_file, m_energies, m_sigmas_unscaled);

    init(scattering_process, energy, scattering_angle_model);
}

template <typename InputVector>
ScatteringProcess::ScatteringProcess (
                        const std::string& scattering_process,
                        const InputVector&& energies,
                        const InputVector&& sigmas,
                        const amrex::ParticleReal energy,
                        const ScatteringAngleModel scattering_angle_model )
{
    m_energies.reserve(energies.size());
    m_sigmas_unscaled.reserve(sigmas.size());
    for (auto const value : energies) {
        m_energies.push_back(static_cast<amrex::ParticleReal>(value));
    }
    for (auto const value : sigmas) {
        m_sigmas_unscaled.push_back(static_cast<double>(value));
    }

    init(scattering_process, energy, scattering_angle_model);
}

void
ScatteringProcess::init (const std::string& scattering_process, const amrex::ParticleReal energy,
                         const ScatteringAngleModel scattering_angle_model)
{
    using namespace amrex::literals;

    m_name = scattering_process;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_energies.size() == static_cast<amrex::Long>(m_sigmas_unscaled.size()),
        "Cross-section energy and value arrays must have the same length."
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_energies.size() >= 2u,
        "Cross-section tables must contain at least two points."
    );
    sanityCheckEnergyGrid(m_energies);
    for (auto const sigma : m_sigmas_unscaled) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(sigma) && sigma >= 0.0,
            "Cross-section values must be finite and non-negative."
        );
    }

    m_sigmas_h.resize(m_sigmas_unscaled.size());
    m_exe_h.m_energies_data = m_energies.data();
    m_exe_h.m_sigmas_data = m_sigmas_h.data();

    // save energy grid parameters for easy use
    const int grid_size = static_cast<int>(m_energies.size());
    m_exe_h.m_grid_size = grid_size;
    m_exe_h.m_energy_lo = m_energies[0];
    m_exe_h.m_energy_hi = m_energies[grid_size-1];
    // The energy grid does not need to be evenly spaced; `m_dE` is only used as a
    // representative energy step. Use the smallest spacing so that finely resolved
    // regions of a non-uniform grid are not skipped over.
    m_exe_h.m_dE = m_energies[grid_size-1] - m_energies[0];
    for (int i = 1; i < grid_size; i++) {
        m_exe_h.m_dE = std::min(m_exe_h.m_dE, m_energies[i] - m_energies[i-1]);
    }
    m_exe_h.m_energy_penalty = energy;
    m_exe_h.m_type = parseProcessType(scattering_process);
    m_exe_h.m_scattering_angle_model = scattering_angle_model;
    m_exe_h.m_produces_products = (
        m_exe_h.m_type == ScatteringProcessType::IONIZATION ||
        m_exe_h.m_type == ScatteringProcessType::ATTACHMENT ||
        m_exe_h.m_type == ScatteringProcessType::TWOPRODUCT_REACTION ||
        m_exe_h.m_type == ScatteringProcessType::CHARGE_EXCHANGE);

    setCrossSectionMultiplier(1.0);

    // check that the cross-section is 0 at the energy cost if the energy
    // cost is > 0 - this is to prevent the possibility of negative left
    // over energy after a collision event
    if (m_exe_h.m_energy_penalty > 0) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            (getCrossSection(m_exe_h.m_energy_penalty) == 0),
            "Cross-section > 0 at energy cost for collision."
        );
    }
}

void
ScatteringProcess::setCrossSectionMultiplier (double const multiplier)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(multiplier) && multiplier > 0.0,
        "Cross-section multiplier must be finite and greater than 0."
    );

    auto const max_value = static_cast<long double>(
        std::numeric_limits<amrex::ParticleReal>::max());
    for (std::size_t i = 0; i < m_sigmas_unscaled.size(); ++i)
    {
        auto const scaled =
            static_cast<long double>(m_sigmas_unscaled[i]) * multiplier;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(scaled) && scaled <= max_value,
            "Scaled cross-section is not representable as ParticleReal."
        );
        m_sigmas_h[i] = static_cast<amrex::ParticleReal>(scaled);
    }

    m_exe_h.m_sigmas_data = m_sigmas_h.data();
    m_exe_h.m_sigma_lo = m_sigmas_h.front();
    m_exe_h.m_sigma_hi = m_sigmas_h.back();

#ifdef AMREX_USE_GPU
    m_exe_d = m_exe_h;
    m_energies_d.resize(m_energies.size());
    m_sigmas_d.resize(m_sigmas_h.size());
    m_exe_d.m_energies_data = m_energies_d.data();
    m_exe_d.m_sigmas_data = m_sigmas_d.data();
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_energies.begin(), m_energies.end(),
                          m_energies_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, m_sigmas_h.begin(), m_sigmas_h.end(),
                          m_sigmas_d.begin());
    amrex::Gpu::streamSynchronize();
#endif
}

ScatteringProcessType
ScatteringProcess::parseProcessType(const std::string& scattering_process)
{
    auto const starts_with = [&scattering_process] (std::string const& prefix)
    {
        return scattering_process.rfind(prefix, 0) == 0;
    };

    if (starts_with("elastic")) {
        // Prefix matching allows several uniquely named channels in one collision.
        return ScatteringProcessType::ELASTIC;
    } else if (scattering_process == "charge_exchange") {
        return ScatteringProcessType::CHARGE_EXCHANGE;
    } else if (scattering_process == "two_product_reaction") {
        return ScatteringProcessType::TWOPRODUCT_REACTION;
    } else if (starts_with("ionization")) {
        return ScatteringProcessType::IONIZATION;
    } else if (starts_with("attachment")) {
        return ScatteringProcessType::ATTACHMENT;
    } else if (starts_with("excitation")) {
        return ScatteringProcessType::EXCITATION;
    } else {
        return ScatteringProcessType::INVALID;
    }
}

void
ScatteringProcess::readCrossSectionFileRaw (
    const std::string& cross_section_file,
    amrex::Vector<amrex::ParticleReal>& energies,
    std::vector<double>& sigmas)
{
    std::ifstream infile(cross_section_file);
    if (!infile.is_open()) {
        WARPX_ABORT_WITH_MESSAGE("Failed to open cross-section data file");
    }

    double energy;
    double sigma;
    while (infile >> energy >> sigma) {
        energies.push_back(static_cast<amrex::ParticleReal>(energy));
        sigmas.push_back(sigma);
    }
    if (infile.bad()) {
        WARPX_ABORT_WITH_MESSAGE("Failed to read cross-section data from file.");
    }
}

void
ScatteringProcess::readCrossSectionFile (
    const std::string& cross_section_file,
    amrex::Vector<amrex::ParticleReal>& energies,
    amrex::Gpu::HostVector<amrex::ParticleReal>& sigmas)
{
    std::vector<double> raw_sigmas;
    readCrossSectionFileRaw(cross_section_file, energies, raw_sigmas);
    sigmas.reserve(raw_sigmas.size());
    for (auto const sigma : raw_sigmas) {
        sigmas.push_back(static_cast<amrex::ParticleReal>(sigma));
    }
}

void
ScatteringProcess::sanityCheckEnergyGrid (
                                   const amrex::Vector<amrex::ParticleReal>& energies
                                   )
{
    // The energy grid does not need to be evenly spaced, but it must be finite,
    // non-negative and sorted in strictly increasing order for interpolation.
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(static_cast<double>(energies[0])) && energies[0] >= 0.0,
        "Cross-section energies must be finite and non-negative."
    );
    for (unsigned i = 1; i < energies.size(); i++) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(energies[i])) && energies[i] >= 0.0,
            "Cross-section energies must be finite and non-negative."
        );
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                                         (energies[i] > energies[i-1]),
                                         "Cross-section energy grid must be sorted in "
                                         "strictly increasing order."
                                         );
    }
}
