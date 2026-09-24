/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCThermalRotation.H"

#include "Particles/Collision/ScatteringProcess.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Gpu.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>

namespace
{
constexpr std::size_t maximum_bytes = 512u * 1024u * 1024u;

template <typename T>
std::vector<T>
readArray (std::ifstream& input, std::size_t count)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(count <= maximum_bytes / sizeof(T),
                                     "Thermal-rotation reference array exceeds 512 MiB.");
    std::vector<T> result(count);
    input.read(reinterpret_cast<char*>(result.data()),
               static_cast<std::streamsize>(count * sizeof(T)));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<bool>(input),
                                     "Truncated thermal-rotation bundle.");
    return result;
}

std::vector<double>
populations (bool nitrogen, double b, int maximum_j, double temperature)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(temperature) && temperature >= 0,
                                     "Rotational temperature must be finite and nonnegative.");
    int const ground = nitrogen ? 0 : 1;
    std::vector<double> result(maximum_j + 1, 0);
    if (temperature == 0) {
        result[ground] = 1;
        return result;
    }
    constexpr double kb = 8.617333262145e-5;
    double const a = b / (kb * temperature);
    double partition = 0;
    for (int j = 0; j <= maximum_j; ++j) {
        int const spin = nitrogen ? (j % 2 == 0 ? 6 : 3) : j % 2;
        result[j] = (2 * j + 1) * spin * std::exp(-a * (j * (j + 1) - ground * (ground + 1)));
        partition += result[j];
    }
    double const first = maximum_j + 1;
    double const tail =
        6 * (2 * first + 1 + 1 / a) * std::exp(-a * (first * (first + 1) - ground * (ground + 1)));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a * (2 * first + 1) * (2 * first + 1) >= 2 && tail / partition < 1.0e-10,
        "Thermal-rotation bundle omits too much Boltzmann population.");
    for (auto& value : result) {
        value /= partition;
    }
    return result;
}

void
appendAlias (std::vector<double> const& weights, std::vector<int> const& outcomes,
             amrex::Gpu::HostVector<BackgroundMCCThermalRotation::Alias>& table)
{
    using Alias = BackgroundMCCThermalRotation::Alias;
    double const sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        sum > 0 && weights.size() == outcomes.size(),
        "Cannot construct an empty thermal-rotation alias distribution.");
    auto const count = static_cast<int>(weights.size());
    std::vector<double> scaled(weights);
    std::vector<int> small, large;
    std::vector<Alias> row(count);
    for (int i = 0; i < count; ++i) {
        scaled[i] *= count / sum;
        row[i].m_outcome = outcomes[i];
        row[i].m_alternate = i;
        (scaled[i] < 1 ? small : large).push_back(i);
    }
    while (!small.empty() && !large.empty()) {
        int const low = small.back();
        small.pop_back();
        int const high = large.back();
        large.pop_back();
        row[low].m_probability = static_cast<amrex::ParticleReal>(scaled[low]);
        row[low].m_alternate = high;
        scaled[high] -= 1 - scaled[low];
        (scaled[high] < 1 ? small : large).push_back(high);
    }
    table.insert(table.end(), row.data(), row.data() + row.size());
}
} // namespace

std::shared_ptr<BackgroundMCCThermalRotation>
BackgroundMCCThermalRotation::get (std::string const& file, std::string const& model,
                                   double temperature, bool cumulative)
{
    static std::mutex mutex;
    static std::map<std::string, std::weak_ptr<BackgroundMCCThermalRotation>> cache;
    std::ostringstream key;
    key << std::filesystem::canonical(file).string() << ':' << model << ':' << std::hexfloat
        << temperature << ":" << cumulative;
    std::lock_guard<std::mutex> lock(mutex);
    auto& entry = cache[key.str()];
    if (auto existing = entry.lock()) {
        return existing;
    }
    auto result =
        std::make_shared<BackgroundMCCThermalRotation>(file, model, temperature, cumulative);
    entry = result;
    return result;
}

BackgroundMCCThermalRotation::BackgroundMCCThermalRotation (std::string const& file,
                                                            std::string const& model,
                                                            double temperature, bool cumulative)
{
    std::ifstream input(file, std::ios::binary);
    std::string magic, target, source_model;
    int maximum_j = 0, energy_count = 0, transition_count = 0;
    double reference_temperature = 0;
    input >> magic >> target >> source_model >> maximum_j >> reference_temperature >>
        energy_count >> transition_count;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        input && magic == "WARPX_THERMAL_ROTATION_V2" && (target == "N2" || target == "O2") &&
            source_model == model && (model == "elastic_dcs" || model == "analytic_test"),
        "Invalid thermal-rotation bundle header or model mismatch.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(maximum_j >= 8 && maximum_j <= 4096 && energy_count >= 2 &&
                                         energy_count <= 100000 && transition_count >= 1 &&
                                         transition_count <= 100000,
                                     "Invalid thermal-rotation bundle dimensions.");
    std::uint32_t const endian = 1;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        *reinterpret_cast<unsigned char const*>(&endian) == 1 && sizeof(double) == 8,
        "Thermal-rotation bundles require little-endian binary64 hosts.");
    char newline = 0;
    input.get(newline);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(newline == '\n', "Invalid thermal-rotation payload boundary.");

    auto const energies = readArray<double>(input, energy_count);
    auto const transitions = readArray<std::int32_t>(input, 2u * transition_count);
    auto const component_count = static_cast<std::size_t>(transition_count) + 1;
    auto const rates =
        readArray<double>(input, static_cast<std::size_t>(energy_count) * component_count);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(input.peek() == std::ifstream::traits_type::eof(),
                                     "Trailing data in thermal-rotation bundle.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(energies.front() == 0,
                                     "Thermal-rotation grid must start at zero energy.");
    auto copy_grid = [] (auto const& source, auto& destination) {
        double previous = -std::numeric_limits<double>::infinity();
        for (auto value : source) {
            auto const narrowed = static_cast<amrex::ParticleReal>(value);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(value) && std::isfinite(static_cast<double>(narrowed)) &&
                    narrowed > previous,
                "Thermal-rotation grids must remain finite and increasing in "
                "particle precision.");
            destination.push_back(narrowed);
            previous = narrowed;
        }
    };
    copy_grid(energies, m_energies_h);
    bool const nitrogen = target == "N2";
    double const b = nitrogen ? 0.0002477204284695341 : 0.00017828927734694198;
    m_neutral_mass = (nitrogen ? 28.0134 : 31.9988) * 1.66053906660e-27;
    auto const population = populations(nitrogen, b, maximum_j, temperature);
    auto const reference_population = populations(nitrogen, b, maximum_j, reference_temperature);
    std::vector<double> factors(component_count, 1), reference_factors(component_count, 1);
    std::vector<double> thresholds(component_count, 0);
    m_outcomes_h.push_back({});
    for (int i = 0; i < transition_count; ++i) {
        int const initial = transitions[2 * i], final = transitions[2 * i + 1];
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(initial >= 0 && final >= 0 && initial <= maximum_j &&
                                             final <= maximum_j && initial != final &&
                                             (initial - final) % 2 == 0 &&
                                             (nitrogen || (initial % 2 == 1 && final % 2 == 1)),
                                         "Invalid homonuclear rotational transition.");
        double const initial_energy = b * initial * (initial + 1);
        double const loss = b * (final * (final + 1) - initial * (initial + 1));
        double const mass =
            m_neutral_mass * PhysConst::c2_v<double> / PhysConst::q_e_v<double> + initial_energy;
        thresholds[i + 1] = loss * (1 + 510998.95069 / mass) + loss * loss / (2 * mass);
        factors[i + 1] = population[initial];
        reference_factors[i + 1] = reference_population[initial];
        m_outcomes_h.push_back({static_cast<amrex::ParticleReal>(loss),
                                static_cast<amrex::ParticleReal>(initial_energy)});
        if (loss > 0 && thresholds[i + 1] < energies.back()) {
            auto const knot = std::find_if(energies.begin(), energies.end(), [&] (double e) {
                return static_cast<float>(e) == static_cast<float>(thresholds[i + 1]);
            });
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(knot != energies.end(),
                                             "Thermal-rotation bundle is missing a "
                                             "physical excitation threshold knot.");
        }
    }
    std::vector<int> order(component_count);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&] (int a, int b_index) {
        auto const av = m_outcomes_h[a], bv = m_outcomes_h[b_index];
        return av.m_loss == bv.m_loss ? av.m_initial_energy > bv.m_initial_energy
                                      : av.m_loss < bv.m_loss;
    });
    double const neutral_rest = m_neutral_mass * PhysConst::c2_v<double> / PhysConst::q_e_v<double>;
    double const largest_loss = m_outcomes_h[order.back()].m_loss;
    m_maximum_loss = largest_loss;
    double const electron_rest =
        PhysConst::m_e_v<double> * PhysConst::c2_v<double> / PhysConst::q_e_v<double>;
    // Backward scattering sets the most restrictive laboratory threshold.
    // This factor bounds it for every retained loss and initial rotational mass.
    double const all_angle_factor =
        (neutral_rest + largest_loss / 2) / (neutral_rest - electron_rest) +
        8 * std::numeric_limits<double>::epsilon();
    double const mass_span = b * maximum_j * (maximum_j + 1);
    // Accessibility changes at E ~ loss. Bound its variation with the initial
    // rotational mass below every distinct canonical level gap.
    double const ordering_bound = 16 * mass_span * largest_loss * (510998.95069 + largest_loss) /
                                  (neutral_rest * neutral_rest);
    double previous_loss = 0;
    for (int index : order) {
        double const loss = m_outcomes_h[index].m_loss;
        if (loss > previous_loss) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                loss - previous_loss > ordering_bound,
                "Rotational levels do not admit ordered kinematic accessibility.");
            previous_loss = loss;
        }
    }
    m_rates_h.resize(energy_count);
    m_reference_rates.assign(energy_count, 0);
    m_cells_h.reserve(energy_count);
    for (int e = 0; e < energy_count; ++e) {
        std::vector<double> weights;
        std::vector<int> outcomes;
        for (int index : order) {
            auto const rate = rates[static_cast<std::size_t>(e) * component_count + index];
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(rate) && rate >= 0 &&
                    (thresholds[index] <= 0 ||
                     energies[e] >
                         thresholds[index] * (1 - 32 * std::numeric_limits<double>::epsilon()) ||
                     rate == 0),
                "Negative, nonfinite or sub-threshold thermal-rotation reference rate.");
            m_reference_rates[e] += reference_factors[index] * rate;
            auto const weighted = factors[index] * rate;
            if (weighted > 0) {
                weights.push_back(weighted);
                outcomes.push_back(index);
            }
        }
        double const total = std::accumulate(weights.begin(), weights.end(), 0.0);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(total) && total <= std::numeric_limits<amrex::ParticleReal>::max(),
            "Invalid aggregate rotational rate.");
        m_rates_h[e] = static_cast<amrex::ParticleReal>(total);
        m_maximum_rate = std::max(m_maximum_rate, total);
        if (weights.empty()) {
            weights.push_back(1);
            outcomes.push_back(0);
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_outcomes_h[outcomes.front()].m_loss <= 0,
            "Every thermal-rotation row needs an unchanged or de-excitation outcome.");
        m_cells_h.push_back(
            {static_cast<int>(m_loss_alias_h.size()), static_cast<int>(weights.size())});
        appendAlias(weights, outcomes, m_loss_alias_h);
        double prefix = 0;
        double const normalization = total > 0 ? total : 1;
        for (double weight : weights) {
            prefix += weight / normalization;
            m_cdf_h.push_back(std::min(prefix, 1.0));
        }
        m_cdf_h.back() = 1;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_loss_alias_h.size() * (sizeof(Alias) + sizeof(double)) <=
                                             maximum_bytes,
                                         "Thermal-rotation sampling tables exceed 512 MiB.");
    }
    m_table_bytes = (m_energies_h.size() + m_rates_h.size()) * sizeof(amrex::ParticleReal) +
                    m_cdf_h.size() * sizeof(double) + m_loss_alias_h.size() * sizeof(Alias) +
                    m_cells_h.size() * sizeof(Cell) + m_outcomes_h.size() * sizeof(Outcome);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_table_bytes <= maximum_bytes && m_maximum_rate > 0,
                                     "Empty or oversized thermal-rotation sampling tables.");
    m_executor_h = {m_energies_h.data(),   m_rates_h.data(),    m_cdf_h.data(),
                    m_loss_alias_h.data(), m_cells_h.data(),    m_outcomes_h.data(),
                    energy_count,          m_energies_h.back(), neutral_rest,
                    electron_rest,         all_angle_factor,    cumulative};
#ifdef AMREX_USE_GPU
    auto upload = [] (auto const& host, auto& device) {
        device.resize(host.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
    };
    upload(m_energies_h, m_energies_d);
    upload(m_rates_h, m_rates_d);
    upload(m_cdf_h, m_cdf_d);
    upload(m_loss_alias_h, m_loss_alias_d);
    upload(m_cells_h, m_cells_d);
    upload(m_outcomes_h, m_outcomes_d);
    m_executor_d = {m_energies_d.data(),   m_rates_d.data(),    m_cdf_d.data(),
                    m_loss_alias_d.data(), m_cells_d.data(),    m_outcomes_d.data(),
                    energy_count,          m_energies_h.back(), neutral_rest,
                    electron_rest,         all_angle_factor,    cumulative};
    amrex::Gpu::streamSynchronize();
#endif
}

void
BackgroundMCCThermalRotation::checkInclusiveRate (ScatteringProcess const& process) const
{
    constexpr double rest = 510998.95069;
    for (std::size_t i = 1; i < m_energies_h.size(); ++i) {
        double const e = m_energies_h[i];
        double const velocity = PhysConst::c_v<double> * std::sqrt(e * (e + 2 * rest)) / (e + rest);
        double const expected = velocity * process.getCrossSection(m_energies_h[i]);
        double const actual = m_reference_rates[i];
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs(actual - expected) <= 0.002 * std::max(actual, expected),
            "Thermal-rotation bundle does not reconstruct the supplied "
            "inclusive elastic rate.");
    }
}
