/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Source/Particles/Collision/ProtonImpactIonization/IonizationSampling.H"

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace
{
    void
    require (bool const condition, char const* message)
    {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    template <typename Real>
    void
    checkMixture ()
    {
        constexpr std::uint32_t count = 4096;
        for (auto const shift : {Real(.00001), Real(.17), Real(.54321), Real(.99999)}) {
            for (bool const reversed : {false, true}) {
                double first_moment = 0;
                double second_moment = 0;
                std::array<int, count / 2> strata{};
                for (std::uint32_t product = 0; product < count; ++product) {
                    // An analytic two-parent test: equal source weights and
                    // conditional uniform spectra on [0,1] and [0,10].
                    // Their exact mixture moments are 2.75 and 101/6.
                    auto const low_energy_parent = (product < count / 2) != reversed;
                    auto const scale = low_energy_parent ? 1.0 : 10.0;
                    auto const quantile = ProtonImpactIonization::shiftedRadicalInverse(
                        product, shift);
                    require(quantile >= Real(0) && quantile < Real(1), "Invalid quantile");
                    auto const energy = scale * static_cast<double>(quantile);
                    first_moment += energy / count;
                    second_moment += energy * energy / count;
                    if (product < count / 2) {
                        auto const bin = static_cast<std::uint32_t>(quantile * Real(count / 2));
                        ++strata[bin];
                    }
                }
                require(std::abs(first_moment - 2.75) < 0.0014, "Biased mixture mean");
                require(std::abs(second_moment - 101.0 / 6.0) < 0.025, "Biased mixture moment");
                for (auto const population : strata) {
                    require(population == 1, "Lost stratification within a parent block");
                }
            }
        }

        // Independently integrate the random shift for one fixed parent/index.
        // The conditional quantile must be uniform, not confined to a stratum
        // determined by that parent's position in the source-weight scan.
        for (auto const index : {0u, 1u, 12345u, 0xffffffffu}) {
            double mean = 0;
            for (std::uint32_t sample = 0; sample < count; ++sample) {
                auto const shift = (static_cast<Real>(sample) + Real(.5)) / Real(count);
                mean += static_cast<double>(
                    ProtonImpactIonization::shiftedRadicalInverse(index, shift)) / count;
            }
            require(std::abs(mean - 0.5) < 0.00013, "Conditional quantile is not uniform");
        }
    }
} // namespace

int
main ()
{
    checkMixture<float>();
    checkMixture<double>();
    std::cout << "PASS: independently shifted energy sampling, float and double\n";
}
