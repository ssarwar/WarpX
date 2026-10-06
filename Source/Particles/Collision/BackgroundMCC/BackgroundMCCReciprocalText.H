/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_RECIPROCAL_TEXT_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_RECIPROCAL_TEXT_H_

#include "BackgroundMCCReciprocalRotation.H"
#include "Utils/TextMsg.H"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

/** Host-only decoding of readable probability tables. No source-model fitting
 * or per-particle work occurs here. See
 * theory/multiphysics/rotational_scattering.
 */
namespace BackgroundMCCReciprocalText {
class Input {
  public:
    explicit Input (std::string const& file) : m_stream(file), m_file(file) {}

    template <typename T>
    T
    next () {
        skipComments();
        T value{};
        m_stream >> value;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            bool(m_stream), "Invalid rotational text data: " + m_file);
        return value;
    }

    void
    finish () {
        skipComments();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_stream.peek() == std::char_traits<char>::eof(),
            "Trailing rotational text data: " + m_file);
    }

  private:
    void
    skipComments () {
        m_stream >> std::ws;
        while (m_stream.peek() == '#') {
            m_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            m_stream >> std::ws;
        }
    }
    std::ifstream m_stream;
    std::string m_file;
};

template <typename T>
std::vector<T>
readArray (std::string const& file, std::size_t size) {
    Input input(file);
    std::vector<T> result(size);
    std::size_t position = 0;
    while (position < size) {
        auto const command = input.next<std::string>();
        auto const count = input.next<std::size_t>();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            count > 0 && count <= size - position,
            "Invalid rotational text array length.");
        if (command == "copy") {
            auto const first = input.next<std::size_t>();
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                first <= position && count <= position - first,
                "Rotational text copy must reference earlier values.");
            std::copy_n(result.begin() + first, count,
                        result.begin() + position);
        } else {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                command == "values" || command == "repeat",
                "Unknown rotational text array command.");
            T repeated{};
            for (std::size_t j = 0; j < count; ++j) {
                if (command == "values" || j == 0) {
                    double const value = input.next<double>();
                    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                        std::isfinite(value) &&
                            value >= static_cast<double>(
                                         std::numeric_limits<T>::lowest()) &&
                            value <= static_cast<double>(
                                         std::numeric_limits<T>::max()) &&
                            (!std::numeric_limits<T>::is_integer ||
                             value == std::floor(value)),
                        "Invalid rotational text array value.");
                    repeated = static_cast<T>(value);
                }
                result[position + j] = repeated;
            }
        }
        position += count;
    }
    input.finish();
    return result;
}

/** Expand a numerical matrix factorization, preserving the exact support.
 * Signed factors are not physical cross sections. Negative reconstruction
 * roundoff is projected to zero, as checked by the offline compression audit.
 * An invalid physical residual is never accepted by that source calculation.
 */
inline std::vector<BackgroundMCCReciprocalRotation::Alias>
readProbabilities (std::string const& file,
                   std::vector<std::uint32_t> const& offsets,
                   std::size_t outcome_count) {
    using Alias = BackgroundMCCReciprocalRotation::Alias;
    Input input(file);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        input.next<std::string>() == "WARPX_PROBABILITY_FACTORS_V1",
        "Invalid rotational probability encoding.");
    auto const cells = input.next<std::size_t>();
    auto const entries = input.next<std::size_t>();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        cells + 1 == offsets.size() && entries == offsets.back() &&
            entries <= (1024u * 1024u * 1024u) / sizeof(Alias),
        "Rotational probability dimensions do not match.");
    std::vector<Alias> result(entries);
    std::size_t position = 0;
    while (position < cells) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            input.next<std::string>() == "block",
            "Missing rotational probability block.");
        auto const first = input.next<std::size_t>();
        auto const count = input.next<std::size_t>();
        auto const rows = input.next<std::size_t>();
        auto const rank = input.next<std::size_t>();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            first == position && count > 0 && count <= 256 &&
                count <= cells - first && rows > 0 && rows <= outcome_count &&
                rank > 0 && rank <= std::min(rows, count),
            "Invalid rotational factor dimensions.");
        std::vector<std::uint16_t> ids(rows);
        std::vector<unsigned char> support(rows * count, 0);
        std::vector<double> left(rows * rank);
        for (std::size_t i = 0; i < rows; ++i) {
            auto const id = input.next<std::size_t>();
            auto const intervals = input.next<std::size_t>();
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                id < outcome_count && intervals > 0 && intervals <= count &&
                    (i == 0 || id > ids[i - 1]),
                "Invalid rotational support identifier.");
            ids[i] = static_cast<std::uint16_t>(id);
            std::size_t end = 0;
            for (std::size_t k = 0; k < intervals; ++k) {
                auto const lo = input.next<std::size_t>();
                auto const hi = input.next<std::size_t>();
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    lo >= end && hi > lo && hi <= count,
                    "Invalid rotational support interval.");
                std::fill(support.begin() + i * count + lo,
                          support.begin() + i * count + hi, 1);
                end = hi;
            }
            for (std::size_t k = 0; k < rank; ++k) {
                left[i * rank + k] = input.next<double>();
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(left[i * rank + k]),
                    "Nonfinite rotational factor.");
            }
        }
        std::vector<double> right(rank), probability;
        std::vector<std::uint16_t> outcomes;
        std::vector<std::size_t> small, large;
        probability.reserve(rows);
        outcomes.reserve(rows);
        small.reserve(rows);
        large.reserve(rows);
        for (std::size_t j = 0; j < count; ++j) {
            for (auto& value : right) {
                value = input.next<double>();
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    std::isfinite(value), "Nonfinite rotational factor.");
            }
            probability.clear();
            outcomes.clear();
            double negative = 0;
            for (std::size_t i = 0; i < rows; ++i) {
                if (support[i * count + j] == 0) {
                    continue;
                }
                double value = 0;
                for (std::size_t k = 0; k < rank; ++k) {
                    value += left[i * rank + k] * right[k];
                }
                negative += std::max(-value, 0.0);
                probability.push_back(std::max(value, 0.0));
                outcomes.push_back(ids[i]);
            }
            double const sum =
                std::accumulate(probability.begin(), probability.end(), 0.0);
            auto const start = offsets[first + j];
            auto const n = probability.size();
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                n > 0 && n <= 65535 && n == offsets[first + j + 1] - start &&
                    std::isfinite(sum) && std::abs(sum - 1) < 2e-6 &&
                    negative < 2e-6,
                "Rotational factors do not reconstruct a probability "
                "distribution.");
            small.clear();
            large.clear();
            for (std::size_t k = 0; k < n; ++k) {
                probability[k] *= double(n) / sum;
                result[start + k] = {1, static_cast<std::uint16_t>(k),
                                     outcomes[k]};
                (probability[k] < 1 ? small : large).push_back(k);
            }
            // Vose's alias construction is linear in cell size and runs
            // once on the host. Every device event still costs one lookup.
            while (!small.empty() && !large.empty()) {
                auto const lo = small.back(), hi = large.back();
                small.pop_back();
                large.pop_back();
                result[start + lo].m_probability =
                    static_cast<float>(probability[lo]);
                result[start + lo].m_alternate = static_cast<std::uint16_t>(hi);
                probability[hi] -= 1 - probability[lo];
                (probability[hi] < 1 ? small : large).push_back(hi);
            }
        }
        position += count;
    }
    input.finish();
    return result;
}
} // namespace BackgroundMCCReciprocalText
#endif
