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
#include <sstream>

ScatteringProcess::ScatteringProcess (
                        const std::string& scattering_process,
                        const std::string& cross_section_file,
                        const amrex::ParticleReal energy,
                        const ScatteringAngleModel scattering_angle_model )
{
    // Retain the unscaled values in double precision. Three-body attachment
    // tables in m^5 can be much smaller than the single-precision range.
    readCrossSectionFileRaw(cross_section_file, m_energies, m_sigmas_unscaled, &m_metadata);

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
        std::isfinite(static_cast<double>(energy)) && energy >= 0.0_prt,
        "The scattering-process energy cost must be finite and non-negative."
    );

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
    bool const has_range = m_metadata.count("outside_energy_range") != 0;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        has_range || (m_metadata.count("energy_min_eV") == 0 &&
                      m_metadata.count("energy_max_eV") == 0),
        "Cross-section energy limits require outside_energy_range = error.");
    if (has_range) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_metadata.at("outside_energy_range") == "error" &&
                m_metadata.count("energy_min_eV") == 1 &&
                m_metadata.count("energy_max_eV") == 1,
            "Bounded cross sections require energy_min_eV, energy_max_eV and "
            "outside_energy_range = error.");
        auto parse_limit = [this] (std::string const& name) {
            std::istringstream input(m_metadata.at(name));
            double value = -1;
            std::string trailing;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                input >> value && !(input >> trailing) &&
                    std::isfinite(value) && value >= 0,
                "Invalid cross-section limit: " + name);
            return value;
        };
        m_exe_h.m_valid_energy_lo = parse_limit("energy_min_eV");
        m_exe_h.m_valid_energy_hi = parse_limit("energy_max_eV");
        m_exe_h.m_bounded_energy = true;
        constexpr double range_tolerance =
            8 * std::numeric_limits<amrex::ParticleReal>::epsilon();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_exe_h.m_valid_energy_lo < m_exe_h.m_valid_energy_hi &&
                m_exe_h.m_valid_energy_hi <= double(m_energies.back()) * (1 + range_tolerance) &&
                (m_exe_h.m_valid_energy_lo >= double(m_energies.front()) * (1 - range_tolerance) ||
                 m_sigmas_unscaled.front() == 0),
            "Cross-section table does not cover its declared energy interval.");
    }
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
        m_exe_h.m_type == ScatteringProcessType::TWO_PRODUCT_REACTION ||
        m_exe_h.m_type == ScatteringProcessType::CHARGE_EXCHANGE);

    setCrossSectionMultiplier(1.0);

    // A zero at the threshold alone does not exclude an unphysical island
    // below it. Check the unscaled data too, before an m^5 table can underflow
    // in particle precision. Linear interpolation then preserves the cutoff.
    if (m_exe_h.m_energy_penalty > 0) {
        for (std::size_t i = 0; i < m_sigmas_unscaled.size(); ++i) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_energies[i] > energy || m_sigmas_unscaled[i] == 0.0,
                "Cross-section must be zero at and below the energy cost for collision."
            );
        }
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

void
ScatteringProcess::useZeroRateGrid (amrex::Gpu::HostVector<amrex::ParticleReal> const& energies)
{
    m_energies.assign(energies.begin(), energies.end());
    m_sigmas_unscaled.assign(energies.size(), 0.0);
    init(m_name, m_exe_h.m_energy_penalty, m_exe_h.m_scattering_angle_model);
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
        return ScatteringProcessType::TWO_PRODUCT_REACTION;
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
ScatteringProcess::readCrossSectionFileRaw (const std::string& cross_section_file,
                                            amrex::Vector<amrex::ParticleReal>& energies,
                                            std::vector<double>& sigmas,
                                            std::map<std::string, std::string>* metadata)
{
    std::ifstream infile(cross_section_file);
    if (!infile.is_open()) {
        WARPX_ABORT_WITH_MESSAGE("Failed to open cross-section data file");
    }

    std::string line;
    int line_number = 0;
    bool resample_native = false;
    std::vector<double> source_energies;
    while (std::getline(infile, line))
    {
        ++line_number;
        std::istringstream row(line);
        row >> std::ws;
        if (row.eof()) {
            continue;
        }
        if (row.peek() == '#') {
            row.get();
            std::string key, equal, value;
            if (row >> key &&
                (key == "rbeq_model" || key == "rbeq_normalization" || key == "rbeq_target" ||
                 key == "outside_energy_range" || key == "energy_min_eV" ||
                 key == "energy_max_eV" || key == "native_energy_grid")) {
                std::string trailing;
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    row >> equal >> value && equal == "=" && !(row >> trailing),
                    "Invalid cross-section metadata: " + key);
                if (key == "native_energy_grid") {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(value == "resample",
                                                     "native_energy_grid must be resample.");
                    resample_native = true;
                }
                if (metadata) {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(metadata->count(key) == 0,
                                                     "Duplicate cross-section metadata: " + key);
                    metadata->emplace(key, value);
                }
            }
            continue;
        }

        double energy;
        double sigma;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<bool>(row >> energy >> sigma),
            "Each cross-section table row must contain an energy and a cross-section. "
            "Invalid row " + std::to_string(line_number) + " in " + cross_section_file + "."
        );

        std::string trailing;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !(row >> trailing) || trailing.rfind('#', 0) == 0,
            "Cross-section tables must contain exactly two columns. Invalid row " +
            std::to_string(line_number) + " in " + cross_section_file + "."
        );

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(energy) && energy >= 0.0 &&
                (source_energies.empty() || energy > source_energies.back()),
            "Source cross-section energies must be finite, non-negative and "
            "strictly increasing.");
        source_energies.push_back(energy);
        energies.push_back(static_cast<amrex::ParticleReal>(energy));
        sigmas.push_back(sigma);
    }
    if (infile.bad()) {
        WARPX_ABORT_WITH_MESSAGE("Failed to read cross-section data from file.");
    }
    if (resample_native && sizeof(amrex::ParticleReal) < sizeof(double) && !energies.empty()) {
        // Resolve the source interpolant at representable particle energies.
        // Keeping the first/last sigma in each rounded bin would bias onsets.
        // Bracket source zeros with native neighbors. Otherwise rounding an
        // onset upward removes the zero prefix, or rounding a terminal zero
        // downward creates a spurious positive tail across a long zero interval.
        std::vector<amrex::ParticleReal> native_zeros;
        for (std::size_t i = 0; i < source_energies.size(); ++i) {
            if (sigmas[i] == 0.0) {
                auto const zero = static_cast<amrex::ParticleReal>(source_energies[i]);
                native_zeros.push_back(zero);
                if (zero > 0) {
                    energies.push_back(std::nextafter(
                        zero, -std::numeric_limits<amrex::ParticleReal>::infinity()));
                    energies.push_back(
                        std::nextafter(zero, std::numeric_limits<amrex::ParticleReal>::infinity()));
                }
            }
        }
        std::sort(energies.begin(), energies.end());
        energies.erase(std::unique(energies.begin(), energies.end()), energies.end());
        std::vector<double> resampled;
        resampled.reserve(energies.size());
        std::size_t interval = 0;
        for (auto const native_energy : energies) {
            double const energy = native_energy;
            if (std::binary_search(native_zeros.begin(), native_zeros.end(), native_energy)) {
                resampled.push_back(0.0);
            } else if (energy <= source_energies.front()) {
                resampled.push_back(sigmas.front());
            } else if (energy >= source_energies.back()) {
                resampled.push_back(sigmas.back());
            } else {
                while (source_energies[interval + 1] < energy) {
                    ++interval;
                }
                double const fraction = (energy - source_energies[interval]) /
                                        (source_energies[interval + 1] - source_energies[interval]);
                resampled.push_back((1 - fraction) * sigmas[interval] +
                                    fraction * sigmas[interval + 1]);
            }
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::any_of(resampled.begin(), resampled.end(),
                        [] (double value) { return value > 0; }) ||
                std::none_of(sigmas.begin(), sigmas.end(), [] (double value) { return value > 0; }),
            "A cross-section channel is unresolved in particle precision; use double precision.");
        sigmas = std::move(resampled);
    }
}

void
ScatteringProcess::readCrossSectionFile (const std::string& cross_section_file,
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
