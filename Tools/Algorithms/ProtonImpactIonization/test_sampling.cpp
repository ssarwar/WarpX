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
    checkKronecker ()
    {
        constexpr std::uint32_t count = 65536;
        for (auto const increment :
             {0x6a09e667u, 0xbb67ae85u, 0x3c6ef373u, 0xa54ff53bu, 0x510e527fu, 0x7311c281u}) {
            for (auto const start : {0u, 1u << 20, 1u << 28, 0xffff0000u}) {
                double mean = 0, second = 0;
                std::array<int, 64> bins{};
                for (std::uint32_t j = 0; j < count; ++j) {
                    auto const q = ProtonImpactIonization::shiftedKronecker<Real>(
                        start + j, 0x87654321u, increment);
                    require(q > Real(0) && q < Real(1), "Closed Kronecker endpoint");
                    mean += double(q) / count;
                    second += double(q) * q / count;
                    ++bins[static_cast<int>(q * Real(bins.size()))];
                }
                require(std::abs(mean - .5) < 1.e-4, "Biased large-index angular/thermal sequence");
                require(std::abs(second - 1. / 3.) < 1.e-4, "Biased sequence second moment");
                for (int const population : bins) {
                    require(std::abs(population - 1024) <= 12,
                            "Large-index sequence lost coverage");
                }
            }
        }
        for (auto const shift : {0u, 1u, 0xfffffffeu, 0xffffffffu}) {
            auto const q = ProtonImpactIonization::shiftedKronecker<Real>(0u, shift, 1u);
            require(q > Real(0) && q < Real(1), "Extreme phase is not in open unit interval");
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
    checkKronecker<float>();
    checkKronecker<double>();
    std::cout << "PASS: independently shifted energy sampling, float and double\n";
}
