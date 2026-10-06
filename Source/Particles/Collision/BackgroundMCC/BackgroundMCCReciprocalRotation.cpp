/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BackgroundMCCReciprocalRotation.H"

#include "BackgroundMCCReciprocalText.H"

#include "Particles/Collision/ScatteringProcess.H"
#include "Utils/TextMsg.H"

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <utility>

namespace {
constexpr std::size_t maximum_bytes = 1024u * 1024u * 1024u;
constexpr std::size_t part_bytes = 32u * 1024u * 1024u;
struct Array {
    std::string m_type;
    std::size_t m_count = 0, m_offset = 0;
    std::string m_file;
};

template <typename T>
std::vector<T>
readArray (std::filesystem::path const& directory,
           std::vector<std::string> const& parts,
           std::map<std::string, Array>& arrays, std::string const& name,
           std::string const& type, std::size_t total) {
    auto const it = arrays.find(name);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(it != arrays.end(),
                                     "Missing V6 array: " + name);
    auto const entry = it->second;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        entry.m_type == type && entry.m_offset % 8 == 0 &&
            entry.m_count <= maximum_bytes / sizeof(T) &&
            entry.m_offset <= total &&
            entry.m_count * sizeof(T) <= total - entry.m_offset,
        "Invalid V6 array: " + name);
    if (!entry.m_file.empty()) {
        if constexpr (!std::is_same_v<T,
                                      BackgroundMCCReciprocalRotation::Alias>) {
            auto values = BackgroundMCCReciprocalText::readArray<T>(
                (directory / entry.m_file).string(), entry.m_count);
            arrays.erase(it);
            return values;
        }
    }
    std::vector<T> values(entry.m_count);
    auto* output = reinterpret_cast<char*>(values.data());
    auto offset = entry.m_offset;
    auto remaining = entry.m_count * sizeof(T);
    while (remaining > 0) {
        auto const part = offset / part_bytes;
        auto const position = offset % part_bytes;
        auto const count = std::min(remaining, part_bytes - position);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(part < parts.size(),
                                         "Truncated V6 part list.");
        std::ifstream input(directory / parts[part], std::ios::binary);
        input.seekg(static_cast<std::streamoff>(position));
        input.read(output, static_cast<std::streamsize>(count));
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(input.good(),
                                         "Truncated V6 array: " + name);
        output += count;
        offset += count;
        remaining -= count;
    }
    arrays.erase(it);
    return values;
}

template <typename T>
T const*
upload (std::vector<T> const& input, amrex::Gpu::DeviceVector<T>& output) {
    output.resize(input.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, input.begin(), input.end(),
                     output.begin());
    return output.data();
}

template <typename T>
void
checkIncreasing (std::vector<T> const& values, std::string const& name) {
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(values.size() >= 2,
                                     "V6 grid needs two nodes: " + name);
    double previous = -1;
    for (auto value : values) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(value) &&
                                             value > previous,
                                         "Invalid V6 grid: " + name);
        previous = value;
    }
}
} // namespace

struct BackgroundMCCReciprocalRotation::Data {
    amrex::Gpu::DeviceVector<double> m_energies, m_rates, m_angular_u,
        m_deflection, m_conditional_u, m_high_edges, m_cdf;
    amrex::Gpu::DeviceVector<std::uint32_t> m_coordinates, m_angular_offsets,
        m_conditional_offsets, m_conditional_cells, m_cell_offsets,
        m_high_cells, m_angular_lookup, m_conditional_lookup;
    amrex::Gpu::DeviceVector<float> m_changing;
    amrex::Gpu::DeviceVector<Alias> m_aliases;
    amrex::Gpu::DeviceVector<std::uint16_t> m_outcome_ids;
    amrex::Gpu::DeviceVector<Outcome> m_outcomes;
};

BackgroundMCCReciprocalRotation::~BackgroundMCCReciprocalRotation () = default;

std::shared_ptr<BackgroundMCCReciprocalRotation>
BackgroundMCCReciprocalRotation::get (std::string const& file,
                                      double temperature, bool cumulative) {
    static std::mutex mutex;
    static std::map<std::string, std::weak_ptr<BackgroundMCCReciprocalRotation>>
        cache;
    std::ostringstream key;
    key << std::filesystem::canonical(file).string() << ':' << std::hexfloat
        << temperature << ':' << cumulative;
    std::lock_guard<std::mutex> lock(mutex);
    auto& slot = cache[key.str()];
    if (auto existing = slot.lock()) {
        return existing;
    }
    auto result = std::make_shared<BackgroundMCCReciprocalRotation>(
        std::filesystem::canonical(file).string(), temperature, cumulative);
    std::size_t bytes = result->tableBytes();
    for (auto const& item : cache) {
        if (auto existing = item.second.lock()) {
            bytes += existing->tableBytes();
        }
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(bytes <= maximum_bytes,
                                     "Loaded reciprocal rotational tables "
                                     "exceed the 1 GiB per-process budget.");
    slot = result;
    return result;
}

BackgroundMCCReciprocalRotation::BackgroundMCCReciprocalRotation (
    std::string const& file, double temperature, bool cumulative)
    : m_data(std::make_unique<Data>()) {
    std::ifstream index(file);
    std::string magic, target, model, sampling;
    double table_temperature;
    auto& e = m_executor;
    index >> magic >> target >> model >> sampling >> table_temperature >>
        m_neutral_mass >> e.m_maximum_energy >> e.m_high_energy >>
        e.m_rutherford_energy >> e.m_separation >> e.m_screening_radius;
    bool const text_format = magic == "WARPX_THERMAL_ROTATION_V7";
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        index && (text_format || magic == "WARPX_THERMAL_ROTATION_V6") &&
            (target == "N2" || target == "O2") &&
            model == "reciprocal_hybrid" &&
            sampling == (text_format ? "probabilities"
                                     : (cumulative ? "cumulative" : "alias")),
        "Invalid V6 model/representation. Prepare the requested sampling "
        "tables offline.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(temperature) && temperature >= 0 &&
            std::abs(table_temperature - temperature) <=
                1e-10 * std::max(1.0, temperature),
        "V6 rotational_temperature does not match the fixed-temperature "
        "bundle.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(e.m_maximum_energy) && e.m_maximum_energy > 0 &&
            e.m_maximum_energy <= 1e9 && std::isfinite(e.m_rutherford_energy) &&
            e.m_rutherford_energy > e.m_high_energy && e.m_high_energy > 0 &&
            e.m_separation > 0 && e.m_screening_radius > 0 &&
            std::isfinite(m_neutral_mass) && m_neutral_mass > 0,
        "Invalid V6 physical metadata.");
    double const expected_mass =
        (target == "N2" ? 28.0134 : 31.9988) * 1.66053906660e-27;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(m_neutral_mass / expected_mass - 1) < 1e-12 &&
            std::isfinite(e.m_separation) &&
            std::isfinite(e.m_screening_radius),
        "V6 molecular constants do not match the target.");
    std::uint32_t endian = 1;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (text_format || *reinterpret_cast<unsigned char*>(&endian) == 1) && sizeof(double) == 8,
        "V6 files require a little-endian binary64 host.");
    std::size_t array_count = 0, part_count = 0, total = 0;
    std::map<std::string, Array> arrays;
    auto const directory = std::filesystem::path(file).parent_path();
    std::vector<std::string> parts;
    if (text_format) {
        index >> array_count;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(index && array_count > 0 &&
                                             array_count < 32,
                                         "Invalid V7 array count.");
        for (std::size_t i = 0; i <= array_count; ++i) {
            std::string name;
            Array value;
            index >> name >> value.m_type >> value.m_count >> value.m_file;
            if (i == array_count) {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    name == "probabilities" && value.m_type == "alias",
                    "Missing V7 probability table.");
                name = "aliases";
            }
            std::size_t const width =
                value.m_type == "f64" || value.m_type == "alias" ? 8
                : value.m_type == "f32" || value.m_type == "u32" ? 4
                : value.m_type == "u16"                          ? 2
                                                                 : 0;
            value.m_offset = (total + 7) / 8 * 8;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                index && width > 0 && value.m_offset <= maximum_bytes &&
                    value.m_count <= (maximum_bytes - value.m_offset) / width &&
                    std::filesystem::path(value.m_file).filename() ==
                        value.m_file &&
                    value.m_file != "." && value.m_file != ".." &&
                    arrays.emplace(name, value).second,
                "Invalid V7 array directory.");
            total = value.m_offset + width * value.m_count;
        }
        std::string trailing;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!(index >> trailing),
                                         "Unexpected V7 index fields.");
    } else {
        index >> array_count >> part_count >> total;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            index && array_count <= 32 && part_count > 0 &&
                part_count <= maximum_bytes / part_bytes &&
                total <= maximum_bytes,
            "Invalid V6 bundle dimensions.");
        std::size_t previous_end = 0;
        for (std::size_t i = 0; i < array_count; ++i) {
            std::string name;
            Array value;
            index >> name >> value.m_type >> value.m_count >> value.m_offset;
            std::size_t const width =
                value.m_type == "f64" || value.m_type == "alias" ? 8
                : value.m_type == "f32" || value.m_type == "u32" ? 4
                : value.m_type == "u16"                          ? 2
                                                                 : 0;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                index && width > 0 && value.m_offset >= previous_end &&
                    value.m_offset % 8 == 0 &&
                    value.m_count <= maximum_bytes / width &&
                    value.m_offset <= total &&
                    width * value.m_count <= total - value.m_offset &&
                    arrays.emplace(name, value).second,
                "Invalid V6 array directory.");
            previous_end = value.m_offset + width * value.m_count;
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(previous_end == total,
                                         "Unexpected V6 trailing bytes.");
        std::size_t part_total = 0;
        for (std::size_t i = 0; i < part_count; ++i) {
            std::string name;
            std::size_t bytes;
            index >> name >> bytes;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                index && std::filesystem::path(name).filename() == name &&
                    bytes > 0 && bytes <= part_bytes &&
                    (i + 1 == part_count || bytes == part_bytes) &&
                    std::filesystem::file_size(directory / name) == bytes,
                "Missing or invalid V6 binary part.");
            parts.push_back(name);
            part_total += bytes;
        }
        std::string trailing;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            part_total == total && !(index >> trailing),
            "V6 index has inconsistent length or trailing fields.");
    }
    m_table_bytes = total;
    auto read_double = [&] (std::string const& name) {
        return readArray<double>(directory, parts, arrays, name, "f64", total);
    };
    auto read_uint = [&] (std::string const& name) {
        return readArray<std::uint32_t>(directory, parts, arrays, name, "u32",
                                        total);
    };
    m_energies = read_double("energies");
    m_rates = read_double("rates");
    m_coordinates = read_uint("coordinates");
    checkIncreasing(m_energies, "energies");
    auto const n = m_energies.size();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        n <= 100000 && m_energies.front() == 0 &&
            m_energies.back() == e.m_maximum_energy && m_rates.size() == n &&
            m_coordinates.size() == n,
        "V6 rate grids do not match the supported energy interval.");
    for (std::size_t i = 0; i < n; ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(m_rates[i]) && m_rates[i] >= 0 &&
                m_coordinates[i] <= 1,
            "Invalid V6 rate or interpolation coordinate.");
    }
    auto ao = read_uint("angular_offsets");
    auto au = read_double("angular_u");
    auto deflection = read_double("deflection");
    auto changing =
        readArray<float>(directory, parts, arrays, "changing", "f32", total);
    auto co = read_uint("conditional_offsets");
    auto cu = read_double("conditional_u");
    auto cc = read_uint("conditional_cells");
    auto cells = read_uint("cell_offsets");
    auto outcomes_raw = read_double("outcomes");
    auto high_edges = read_double("high_edges");
    auto high_cells = read_uint("high_cells");
    auto offsets_valid = [] (auto const& offsets, std::size_t rows,
                             std::size_t values) {
        return offsets.size() == rows + 1 && offsets.front() == 0 &&
               offsets.back() == values &&
               std::is_sorted(offsets.begin(), offsets.end());
    };
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        offsets_valid(ao, n, au.size()) && au.size() == deflection.size() &&
            au.size() == changing.size() && offsets_valid(co, n, cu.size()) &&
            cu.size() == cc.size() && cells.size() >= 2 && cells.front() == 0 &&
            std::is_sorted(cells.begin(), cells.end()) &&
            outcomes_raw.size() >= 2 && outcomes_raw.size() % 2 == 0 &&
            outcomes_raw.size() / 2 <= 65535,
        "Invalid V6 sampling dimensions.");
    for (std::size_t row = 0; row < n; ++row) {
        for (int kind = 0; kind < 2; ++kind) {
            auto const& offsets = kind == 0 ? ao : co;
            auto const& grid = kind == 0 ? au : cu;
            auto const first = offsets[row], last = offsets[row + 1];
            if (last == first) {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    kind != 0 || m_energies[row] > e.m_rutherford_energy,
                    "V6 angular row is missing.");
                continue;
            }
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                last - first >= 2 && grid[first] == 0 && grid[last - 1] == 1,
                "V6 quantile row is missing its endpoints.");
            for (auto i = first; i < last; ++i) {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(grid[i]) &&
                        (i == first || grid[i] > grid[i - 1]),
                    "Invalid V6 quantile row.");
            }
        }
        for (auto i = ao[row]; i < ao[row + 1]; ++i) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                std::isfinite(deflection[i]) && deflection[i] >= 0 &&
                    deflection[i] <= 2 &&
                    (i == ao[row] || deflection[i] >= deflection[i - 1]) &&
                    std::isfinite(changing[i]) && changing[i] >= 0 &&
                    changing[i] <= 1 &&
                    (changing[i] == 0 || co[row + 1] > co[row]),
                "Invalid V6 angular probabilities.");
        }
        if (ao[row + 1] > ao[row]) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                deflection[ao[row]] == 0 && deflection[ao[row + 1] - 1] == 2,
                "V6 angular support is incomplete.");
        }
    }
    std::vector<std::uint32_t> angular_lookup, conditional_lookup;
    auto read_lookup = [&] (std::string const& name, auto const& offsets,
                            auto const& grid) {
        std::vector<std::uint32_t> lookup;
        constexpr int bins = Executor::lookup_bins;
        if (text_format) {
            lookup.resize(n * (bins + 1), 0);
            for (std::size_t row = 0; row < n; ++row) {
                auto const first = offsets[row], end = offsets[row + 1];
                if (first == end) {
                    continue;
                }
                for (int bin = 0; bin <= bins; ++bin) {
                    double const quantile = double(bin) / bins;
                    auto const upper = std::upper_bound(
                        grid.begin() + first, grid.begin() + end, quantile);
                    lookup[row * (bins + 1) + bin] =
                        static_cast<std::uint32_t>(std::min<std::size_t>(
                            upper - grid.begin() - 1, end - 2));
                }
            }
            m_table_bytes += lookup.size() * sizeof(std::uint32_t);
        } else if (arrays.count(name) != 0) {
            lookup = read_uint(name);
        } else {
            return lookup;
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            lookup.size() == n * (bins + 1),
            "Invalid V6 quantile lookup dimensions.");
        for (std::size_t row = 0; row < n; ++row) {
            auto const first = offsets[row], end = offsets[row + 1];
            for (int bin = 0; bin <= bins; ++bin) {
                auto const i = lookup[row * (bins + 1) + bin];
                if (first == end) {
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        i == 0, "Invalid empty V6 lookup row.");
                    continue;
                }
                double const quantile = double(bin) / bins;
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    i >= first && i < end - 1 && grid[i] <= quantile &&
                        (bin == bins ? i == end - 2 : grid[i + 1] > quantile),
                    "V6 lookup does not bracket its quantile bin.");
            }
        }
        return lookup;
    };
    angular_lookup = read_lookup("angular_lookup", ao, au);
    conditional_lookup = read_lookup("conditional_lookup", co, cu);
    std::vector<Outcome> outcomes(outcomes_raw.size() / 2);
    double const rotational_constant =
        target == "N2" ? 0.0002477204284695341 : 0.00017828927734694198;
    int const ground = target == "N2" ? 0 : 1;
    double const dissociation = target == "N2" ? 9.759 : 5.116;
    auto level = [&] (double energy) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(energy) && energy >= -1e-12 && energy < dissociation,
            "V6 rotational level is outside the bound manifold.");
        double const value =
            energy / rotational_constant + ground * (ground + 1);
        int const j =
            static_cast<int>(std::lround((std::sqrt(1 + 4 * value) - 1) / 2));
        double const canonical =
            rotational_constant * (j * (j + 1) - ground * (ground + 1));
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            j >= ground && (target == "N2" || j % 2 == 1) &&
                std::abs(energy - canonical) <= 1e-12 && energy < dissociation,
            "V6 outcome is not a canonical bound rotational level.");
        return j;
    };
    for (std::size_t i = 0; i < outcomes.size(); ++i) {
        outcomes[i] = {outcomes_raw[2 * i], outcomes_raw[2 * i + 1]};
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(outcomes[i].m_loss) &&
                std::isfinite(outcomes[i].m_initial_energy) &&
                outcomes[i].m_initial_energy >= 0 &&
                outcomes[i].m_initial_energy + outcomes[i].m_loss >= -1e-12,
            "Invalid V6 rotational energy change.");
        int const initial = level(outcomes[i].m_initial_energy);
        int const final =
            level(outcomes[i].m_initial_energy + outcomes[i].m_loss);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            (initial - final) % 2 == 0,
            "V6 outcome changes homonuclear rotational parity.");
    }
    checkIncreasing(high_edges, "momentum transfer");
    constexpr double electron_rest_energy = 510998.95069;
    constexpr double fine_structure = 7.2973525693e-3;
    double const maximum_q =
        e.m_separation *
        std::sqrt(e.m_maximum_energy *
                  (e.m_maximum_energy + 2 * electron_rest_energy)) /
        (fine_structure * electron_rest_energy);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        high_edges.back() >= maximum_q,
        "V6 momentum-transfer tables do not cover the supported energy range.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        high_edges.front() == 0 && high_cells.size() + 1 == high_edges.size(),
        "Invalid V6 momentum-transfer cells.");
    for (auto cell : cc) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(cell + 1 < cells.size(),
                                         "Invalid V6 conditional cell.");
    }
    for (auto cell : high_cells) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(cell + 1 < cells.size(),
                                         "Invalid V6 high-energy cell.");
    }
    auto& d = *m_data;
    std::vector<double> largest_loss(cells.size() - 1, 0);
    std::vector<Alias> prepared;
    if (text_format) {
        auto const item = arrays.find("aliases");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            item != arrays.end() && item->second.m_count == cells.back(),
            "V7 probability length does not match its cells.");
        prepared = BackgroundMCCReciprocalText::readProbabilities(
            (directory / item->second.m_file).string(), cells, outcomes.size());
        arrays.erase(item);
    }
    if (cumulative) {
        std::vector<double> cdf;
        std::vector<std::uint16_t> ids;
        if (text_format) {
            cdf.resize(prepared.size(), 0);
            ids.resize(prepared.size());
            for (std::size_t cell = 0; cell + 1 < cells.size(); ++cell) {
                auto const first = cells[cell], end = cells[cell + 1];
                double const count = end - first;
                for (auto i = first; i < end; ++i) {
                    auto const entry = prepared[i];
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        entry.m_alternate < end - first,
                        "Invalid prepared alias index.");
                    cdf[i] += entry.m_probability / count;
                    cdf[first + entry.m_alternate] +=
                        (1.0 - entry.m_probability) / count;
                    ids[i] = entry.m_outcome;
                }
                long double sum = 0;
                for (auto i = first; i < end; ++i) {
                    sum += cdf[i];
                }
                long double prefix = 0;
                for (auto i = first; i < end; ++i) {
                    prefix += cdf[i];
                    cdf[i] = static_cast<double>(prefix / sum);
                }
                cdf[end - 1] = 1;
            }
            m_table_bytes += ids.size() * sizeof(std::uint16_t);
            std::vector<Alias>().swap(prepared);
        } else {
            cdf = read_double("cdf");
            ids = readArray<std::uint16_t>(directory, parts, arrays,
                                           "outcome_ids", "u16", total);
        }
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(cdf.size() == ids.size() &&
                                             cells.back() == cdf.size(),
                                         "Invalid V6 cumulative table.");
        for (std::size_t c = 0; c + 1 < cells.size(); ++c) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                cells[c + 1] > cells[c] && cells[c + 1] - cells[c] <= 65535 &&
                    cdf[cells[c + 1] - 1] == 1,
                "Invalid V6 cumulative cell.");
            double previous = 0;
            for (auto i = cells[c]; i < cells[c + 1]; ++i) {
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(cdf[i]) && cdf[i] >= previous &&
                        cdf[i] <= 1 && ids[i] < outcomes.size(),
                    "Invalid V6 cumulative entry.");
                previous = cdf[i];
                largest_loss[c] =
                    std::max(largest_loss[c], outcomes[ids[i]].m_loss);
            }
        }
        e.m_cdf = upload(cdf, d.m_cdf);
        e.m_outcome_ids = upload(ids, d.m_outcome_ids);
    } else {
        auto aliases = text_format
                           ? std::move(prepared)
                           : readArray<Alias>(directory, parts, arrays,
                                              "aliases", "alias", total);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(cells.back() == aliases.size(),
                                         "Invalid V6 alias length.");
        for (std::size_t c = 0; c + 1 < cells.size(); ++c) {
            auto const count = cells[c + 1] - cells[c];
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(count > 0 && count <= 65535,
                                             "Invalid V6 alias cell size.");
            for (auto i = cells[c]; i < cells[c + 1]; ++i) {
                auto const entry = aliases[i];
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(entry.m_probability) &&
                        entry.m_probability >= 0 && entry.m_probability <= 1 &&
                        entry.m_alternate < count &&
                        entry.m_outcome < outcomes.size(),
                    "Invalid V6 alias entry.");
                largest_loss[c] =
                    std::max(largest_loss[c], outcomes[entry.m_outcome].m_loss);
            }
        }
        e.m_aliases = upload(aliases, d.m_aliases);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(arrays.empty(), "Unknown V6 arrays.");
    for (std::size_t row = 0; row < n && m_energies[row] <= e.m_high_energy;
         ++row) {
        double const minimum_energy = row == 0 ? 0 : m_energies[row - 1];
        for (auto j = co[row]; j < co[row + 1]; ++j) {
            // An upper row may be selected anywhere in its preceding interval.
            // Its support must therefore be open at that interval's lower end.
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                largest_loss[cc[j]] <= minimum_energy,
                "V6 row mixtures would permit a subthreshold rotational "
                "excitation.");
        }
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_table_bytes <= maximum_bytes,
        "Prepared rotational tables exceed the 1 GiB budget.");
    e.m_energies = upload(m_energies, d.m_energies);
    e.m_rates = upload(m_rates, d.m_rates);
    e.m_coordinates = upload(m_coordinates, d.m_coordinates);
    e.m_angular_offsets = upload(ao, d.m_angular_offsets);
    e.m_angular_u = upload(au, d.m_angular_u);
    e.m_deflection = upload(deflection, d.m_deflection);
    e.m_changing = upload(changing, d.m_changing);
    e.m_conditional_offsets = upload(co, d.m_conditional_offsets);
    e.m_conditional_u = upload(cu, d.m_conditional_u);
    e.m_conditional_cells = upload(cc, d.m_conditional_cells);
    if (!angular_lookup.empty()) {
        e.m_angular_lookup = upload(angular_lookup, d.m_angular_lookup);
    }
    if (!conditional_lookup.empty()) {
        e.m_conditional_lookup =
            upload(conditional_lookup, d.m_conditional_lookup);
    }
    e.m_cell_offsets = upload(cells, d.m_cell_offsets);
    e.m_outcomes = upload(outcomes, d.m_outcomes);
    e.m_high_edges = upload(high_edges, d.m_high_edges);
    e.m_high_cells = upload(high_cells, d.m_high_cells);
    e.m_energy_count = static_cast<int>(n);
    e.m_high_count = static_cast<int>(high_cells.size());
    e.m_cumulative = cumulative;
    m_host_executor = {};
    m_host_executor.m_energies = m_energies.data();
    m_host_executor.m_rates = m_rates.data();
    m_host_executor.m_coordinates = m_coordinates.data();
    m_host_executor.m_energy_count = e.m_energy_count;
    m_host_executor.m_maximum_energy = e.m_maximum_energy;
}

void
BackgroundMCCReciprocalRotation::checkInclusiveRate (
    ScatteringProcess const& process) const {
    constexpr double rest = 510998.95069;
    constexpr double c = 299792458.0;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(process.getMinEnergyInput() == 0 &&
                                         process.getMaxEnergyInput() >=
                                             m_energies.back(),
                                     "The elastic source table does not cover "
                                     "the reciprocal bundle's energy range.");
    for (std::size_t i = 0; i < m_energies.size(); ++i) {
        double const energy = m_energies[i];
        if (energy < 0.001) {
            continue;
        }
        double const expected =
            c * std::sqrt(energy * (energy + 2 * rest)) / (energy + rest) *
            process.getCrossSection(static_cast<amrex::ParticleReal>(energy));
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            expected > 0 && std::abs(m_rates[i] / expected - 1) < 0.002,
            "Reciprocal bundle does not match the inclusive elastic source "
            "rate.");
    }
}
