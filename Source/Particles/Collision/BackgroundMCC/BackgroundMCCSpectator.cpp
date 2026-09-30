/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCSpectator.H"

#include "Utils/TextMsg.H"

#include <AMReX_Gpu.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <sstream>
#include <vector>

namespace
{
    constexpr std::size_t maximum_bytes = 512u * 1024u * 1024u;

    template <typename T>
    void readArray (std::ifstream& input, std::size_t count, amrex::Gpu::HostVector<T>& values)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(count <= maximum_bytes / sizeof(T),
                                         "Spectator table exceeds the 512 MiB limit.");
        values.resize(count);
        input.read(reinterpret_cast<char*>(values.data()),
                   static_cast<std::streamsize>(count * sizeof(T)));
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(input.good(), "Truncated spectator data.");
    }
} // namespace

struct BackgroundMCCSpectator::Data
{
    explicit Data (std::string const& file);
    amrex::Gpu::HostVector<amrex::ParticleReal> m_energies_h;
    amrex::Gpu::HostVector<int> m_angle_offsets_h;
    amrex::Gpu::HostVector<float> m_angles_h, m_basis_h, m_weights_h;
    Executor m_executor_h;
#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amrex::ParticleReal> m_energies_d;
    amrex::Gpu::DeviceVector<int> m_angle_offsets_d;
    amrex::Gpu::DeviceVector<float> m_angles_d, m_basis_d, m_weights_d;
    Executor m_executor_d;
#endif
    std::size_t m_bytes = 0;
};

std::shared_ptr<BackgroundMCCSpectator>
BackgroundMCCSpectator::get (std::string const& file, double temperature, bool cumulative)
{
    static std::mutex mutex;
    static std::map<std::string, std::weak_ptr<BackgroundMCCSpectator>> cache;
    std::ostringstream key;
    key << std::filesystem::canonical(file).string() << ':' << std::hexfloat << temperature << ':'
        << cumulative;
    std::lock_guard<std::mutex> lock(mutex);
    auto& entry = cache[key.str()];
    if (auto existing = entry.lock()) {
        return existing;
    }
    auto result = std::make_shared<BackgroundMCCSpectator>(file, temperature, cumulative);
    entry = result;
    return result;
}

BackgroundMCCSpectator::Data::Data (std::string const& file)
{
    std::ifstream input(file, std::ios::binary);
    std::string magic, target, model;
    int maximum_j = 0, energy_count = 0, angle_count = 0, reserved = 0;
    input >> magic >> target >> model >> maximum_j >> reserved >> energy_count >> angle_count;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        input && magic == "WARPX_THERMAL_ROTATION_V5" && target == "N2" &&
            model == "iaa_spectator" && reserved == 0 && maximum_j >= 8 && maximum_j <= 4096 &&
            energy_count >= 2 && energy_count <= 100000 && angle_count >= 2 * energy_count,
        "Invalid IAA spectator bundle header.");
    std::uint32_t const endian = 1;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(*reinterpret_cast<unsigned char const*>(&endian) == 1 &&
                                         sizeof(double) == 8 && sizeof(float) == 4 &&
                                         sizeof(int) == 4,
                                     "Spectator bundles require little-endian IEEE hosts.");
    auto const weight_count = static_cast<std::size_t>(energy_count) * (maximum_j + 1) * 28;
    auto const byte_count = std::size_t(energy_count) * 8 + std::size_t(energy_count + 1) * 4 +
                            std::size_t(angle_count) * 20 + weight_count * 4;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(byte_count <= maximum_bytes,
                                     "Spectator data exceed the 512 MiB limit.");
    char newline = 0;
    input.get(newline);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(newline == '\n', "Invalid spectator payload boundary.");
    amrex::Gpu::HostVector<double> energies;
    readArray(input, energy_count, energies);
    readArray(input, energy_count + 1, m_angle_offsets_h);
    readArray(input, angle_count, m_angles_h);
    readArray(input, 4u * angle_count, m_basis_h);
    readArray(input, weight_count, m_weights_h);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(input.peek() == std::ifstream::traits_type::eof(),
                                     "Trailing spectator data.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(energies.front() == 0 && m_angle_offsets_h.front() == 0 &&
                                         m_angle_offsets_h.back() == angle_count,
                                     "Invalid spectator grid endpoints.");
    amrex::ParticleReal previous = -1;
    constexpr double b = 0.0002477204284695341;
    for (int e = 0; e < energy_count; ++e) {
        auto const energy = static_cast<amrex::ParticleReal>(energies[e]);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(energies[e]) && std::isfinite(energy) &&
                                             energy > previous,
                                         "Spectator energies must remain finite and increasing.");
        m_energies_h.push_back(energy);
        previous = energy;
        int const first = m_angle_offsets_h[e];
        int const last = m_angle_offsets_h[e + 1] - 1;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(first >= 0 && last > first && last < angle_count &&
                                             m_angles_h[first] == 0 && m_angles_h[last] == 1,
                                         "Invalid spectator angular row.");
        for (int a = first + 1; a <= last; ++a) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(m_angles_h[a]) &&
                                                 m_angles_h[a] > m_angles_h[a - 1],
                                             "Spectator angles must remain finite and increasing.");
        }
        for (int j = 0; j <= maximum_j; ++j) {
            for (int choice = 0; choice < 7; ++choice) {
                int const final = j + 2 * (choice - 3);
                double const loss = b * (final * (final + 1) - j * (j + 1));
                for (int rank = 0; rank < 4; ++rank) {
                    float const value =
                        m_weights_h[((e * (maximum_j + 1) + j) * 7 + choice) * 4 + rank];
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        std::isfinite(value) && value >= 0 && value <= 1 &&
                            ((final >= 0 && loss <= energies[e]) || value == 0),
                        "Invalid, subthreshold, or negative-state spectator "
                        "weight.");
                }
            }
        }
    }
    for (auto value : m_basis_h) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(value) && value >= 0,
                                         "Invalid spectator angular basis.");
    }

    m_bytes = m_energies_h.size() * sizeof(amrex::ParticleReal) +
              m_angle_offsets_h.size() * sizeof(int) +
              (m_angles_h.size() + m_basis_h.size() + m_weights_h.size()) * sizeof(float);
    m_executor_h = {m_energies_h.data(),
                    m_angle_offsets_h.data(),
                    m_angles_h.data(),
                    m_basis_h.data(),
                    m_weights_h.data(),
                    nullptr,
                    nullptr,
                    energy_count,
                    maximum_j + 1,
                    0,
                    energies.back(),
                    false};
#ifdef AMREX_USE_GPU
    auto upload = [] (auto const& host, auto& device) {
        device.resize(host.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
    };
    upload(m_energies_h, m_energies_d);
    upload(m_angle_offsets_h, m_angle_offsets_d);
    upload(m_angles_h, m_angles_d);
    upload(m_basis_h, m_basis_d);
    upload(m_weights_h, m_weights_d);
    m_executor_d = {m_energies_d.data(),
                    m_angle_offsets_d.data(),
                    m_angles_d.data(),
                    m_basis_d.data(),
                    m_weights_d.data(),
                    nullptr,
                    nullptr,
                    energy_count,
                    maximum_j + 1,
                    0,
                    energies.back(),
                    false};
    amrex::Gpu::streamSynchronize();
#endif
}

BackgroundMCCSpectator::BackgroundMCCSpectator (std::string const& file, double temperature,
                                                bool cumulative)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(temperature) && temperature >= 0,
                                     "Rotational temperature must be finite and nonnegative.");
    // Differential weights are temperature independent. Share their host and
    // device storage across all populations and sampling reference modes.
    static std::mutex mutex;
    static std::map<std::string, std::weak_ptr<Data const>> cache;
    {
        std::lock_guard<std::mutex> lock(mutex);
        auto& entry = cache[std::filesystem::canonical(file).string()];
        m_data = entry.lock();
        if (!m_data) {
            m_data = std::make_shared<Data>(file);
            entry = m_data;
        }
    }
    int const maximum_j = m_data->m_executor_h.m_state_count - 1;
    constexpr double b = 0.0002477204284695341;

    // This is only a distribution of virtual collision partners. No rotational
    // state is attached to a neutral particle or updated after an event.
    std::vector<double> population;
    std::vector<int> levels;
    double partition = 0;
    constexpr double kb = 8.617333262145e-5;
    for (int j = 0; j <= maximum_j; ++j) {
        double const value = temperature == 0 ? (j == 0 ? 1.0 : 0.0)
                                              : (2 * j + 1) * (j % 2 == 0 ? 6 : 3) *
                                                    std::exp(-b * j * (j + 1) / (kb * temperature));
        if (value > 0) {
            population.push_back(value);
            levels.push_back(j);
            partition += value;
        }
    }
    if (temperature > 0) {
        double const a = b / (kb * temperature);
        double const first = maximum_j + 1;
        double const tail = 6 * (2 * first + 1 + 1 / a) * std::exp(-a * first * (first + 1));
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(a * (2 * first + 1) * (2 * first + 1) >= 2 &&
                                             tail / partition < 1.0e-10,
                                         "Spectator data omit too much Boltzmann population.");
    }
    int const count = static_cast<int>(population.size());
    std::vector<int> small, large;
    double prefix = 0;
    for (int i = 0; i < count; ++i) {
        prefix += population[i] / partition;
        if (cumulative) {
            m_cdf_h.push_back(std::min(prefix, 1.0));
        }
        population[i] *= count / partition;
        m_alias_h.push_back({1, i, levels[i]});
        (population[i] < 1 ? small : large).push_back(i);
    }
    if (cumulative) {
        m_cdf_h.back() = 1;
    }
    while (!small.empty() && !large.empty()) {
        int const lo = small.back();
        small.pop_back();
        int const hi = large.back();
        large.pop_back();
        m_alias_h[lo].m_probability = static_cast<amrex::ParticleReal>(population[lo]);
        m_alias_h[lo].m_alternate = hi;
        population[hi] -= 1 - population[lo];
        (population[hi] < 1 ? small : large).push_back(hi);
    }
    m_table_bytes =
        m_data->m_bytes + m_alias_h.size() * sizeof(Alias) + m_cdf_h.size() * sizeof(double);
    m_executor_h = m_data->m_executor_h;
    m_executor_h.m_initial_alias = m_alias_h.data();
    m_executor_h.m_initial_cdf = m_cdf_h.data();
    m_executor_h.m_population_count = count;
    m_executor_h.m_cumulative = cumulative;
#ifdef AMREX_USE_GPU
    auto upload = [] (auto const& host, auto& device) {
        device.resize(host.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
    };
    upload(m_alias_h, m_alias_d);
    upload(m_cdf_h, m_cdf_d);
    m_executor_d = m_data->m_executor_d;
    m_executor_d.m_initial_alias = m_alias_d.data();
    m_executor_d.m_initial_cdf = m_cdf_d.data();
    m_executor_d.m_population_count = count;
    m_executor_d.m_cumulative = cumulative;
    amrex::Gpu::streamSynchronize();
#endif
}
