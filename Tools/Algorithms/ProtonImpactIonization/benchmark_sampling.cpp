/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Source/Particles/Collision/ProtonImpactIonization/IonizationSampling.H"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <utility>

namespace
{
    template <typename Real, bool Shifted>
    std::pair<double, double>
    measure ()
    {
        constexpr std::uint32_t iterations = 1u << 22;
        double checksum = 0;
        auto const start = std::chrono::steady_clock::now();
        for (std::uint32_t iteration = 0; iteration < iterations; ++iteration) {
            auto const product = iteration & 4095u;
            auto const shift = (static_cast<Real>(iteration >> 12) + Real(.5)) / Real(1024);
            Real quantile;
            if constexpr (Shifted) {
                quantile = ProtonImpactIonization::shiftedRadicalInverse(product, shift);
            } else {
                // Timing baseline only: the ordered sequence biases mixed parents.
                quantile = (static_cast<Real>(product) + shift) / Real(4096);
            }
            checksum += static_cast<double>(quantile * quantile);
        }
        auto const stop = std::chrono::steady_clock::now();
        auto const nanoseconds = std::chrono::duration<double, std::nano>(stop - start).count();
        return {nanoseconds / iterations, checksum};
    }

    template <typename Real>
    void
    benchmark (char const* precision)
    {
        std::array<double, 7> ordered{};
        std::array<double, 7> shifted{};
        double checksum = 0;
        for (std::size_t repeat = 0; repeat < ordered.size(); ++repeat) {
            auto const first = measure<Real, false>();
            auto const second = measure<Real, true>();
            ordered[repeat] = first.first;
            shifted[repeat] = second.first;
            checksum += first.second + second.second;
        }
        std::sort(ordered.begin(), ordered.end());
        std::sort(shifted.begin(), shifted.end());
        std::cout << precision << ": median ns/quantile, ordered=" << ordered[3]
                  << ", shifted=" << shifted[3] << ", checksum=" << checksum << '\n';
    }
} // namespace

int
main ()
{
    // This isolates host quantile generation, not the full source or GPU cost.
    benchmark<float>("float");
    benchmark<double>("double");
}
