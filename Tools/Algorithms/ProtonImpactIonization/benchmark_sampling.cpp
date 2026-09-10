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
#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>

namespace
{
    template <typename Real, bool FixedPoint>
    std::pair<double, double>
    measurePhase ()
    {
        constexpr std::uint32_t iterations = 1u << 22;
        double second_moment = 0;
        auto const start = std::chrono::steady_clock::now();
        for (std::uint32_t j = 0; j < iterations; ++j) {
            auto const index = (1u << 28) + j;
            Real value;
            if constexpr (FixedPoint) {
                value =
                    ProtonImpactIonization::shiftedKronecker<Real>(index, 0x12345678u, 0x6a09e667u);
            } else {
                // Historical float arithmetic loses every fractional bit at
                // these indices. Report the moment as well as the timing.
                auto const phase = Real(.071111111) + Real(.4142135623730950488) *
                                                          (static_cast<Real>(index) + Real(.5));
                value = phase - std::floor(phase);
            }
            second_moment += double(value) * value;
        }
        auto const elapsed =
            std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start)
                .count();
        return {elapsed / iterations, second_moment / iterations};
    }

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
        std::array<double, 7> floating_phase{}, fixed_phase{};
        double old_moment = 0, new_moment = 0;
        double checksum = 0;
        for (std::size_t repeat = 0; repeat < ordered.size(); ++repeat) {
            auto const first = measure<Real, false>();
            auto const second = measure<Real, true>();
            ordered[repeat] = first.first;
            shifted[repeat] = second.first;
            checksum += first.second + second.second;
            auto const old_phase = measurePhase<Real, false>();
            auto const new_phase = measurePhase<Real, true>();
            floating_phase[repeat] = old_phase.first;
            fixed_phase[repeat] = new_phase.first;
            old_moment = old_phase.second;
            new_moment = new_phase.second;
        }
        std::sort(ordered.begin(), ordered.end());
        std::sort(shifted.begin(), shifted.end());
        std::sort(floating_phase.begin(), floating_phase.end());
        std::sort(fixed_phase.begin(), fixed_phase.end());
        std::cout << precision << ": median ns/quantile, ordered=" << ordered[3]
                  << ", shifted=" << shifted[3] << ", checksum=" << checksum << '\n';
        std::cout << precision << ": median ns/phase, floating=" << floating_phase[3]
                  << ", fixed=" << fixed_phase[3] << "; second moments=" << old_moment << ", "
                  << new_moment << " (uniform reference: 1/3)\n";
    }
} // namespace

int
main ()
{
    // This isolates host quantile generation, not the full source or GPU cost.
    benchmark<float>("float");
    benchmark<double>("double");
}
