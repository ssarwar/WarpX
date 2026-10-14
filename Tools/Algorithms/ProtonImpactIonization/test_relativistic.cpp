/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Source/Particles/Collision/ProtonImpactIonization/RelativisticIonization.H"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>

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
    checkPrecision ()
    {
        constexpr long double electron_mass = 510998.95069L;
        constexpr long double proton_mass = 938272088.16L;
        constexpr auto tolerance = std::is_same_v<Real, float> ? 7.e-7L : 2.e-14L;
        constexpr auto factor_tolerance = std::is_same_v<Real, float> ? 8.e-7L : 5.e-11L;
        constexpr std::array<long double, 4> mass_numbers{1.L, 4.L, 16.L, 238.L};

        for (auto const mass_number : mass_numbers) {
            auto const mass = static_cast<Real>(mass_number * proton_mass);
            auto const m = static_cast<long double>(mass);
            for (int index = 0; index <= 360; ++index) {
                auto const energy = static_cast<Real>(std::pow(10.L, -3.L + index / 20.L));
                auto const e = static_cast<long double>(energy);
                // Independent invariant result, without constructing gamma.
                auto const reference = 2.L * electron_mass * e * (e + 2.L * m) /
                    ((m + electron_mass) * (m + electron_mass) + 2.L * electron_mass * e);
                auto const maximum =
                    ProtonImpactIonization::relativisticMaximumTransfer(energy, mass);
                require(std::isfinite(maximum) && maximum > Real(0), "Invalid Tmax");
                require(std::abs(maximum / reference - 1.L) < tolerance, "Tmax mismatch");
                require(maximum <= energy, "Transfer exceeds incident kinetic energy");

                auto const beta_squared =
                    ProtonImpactIonization::relativisticBetaSquared(energy, mass);
                auto const reference_beta = (e / (e + m)) * ((e + 2.L * m) / (e + m));
                require(std::abs(beta_squared / reference_beta - 1.L) < tolerance,
                        "beta squared mismatch");
                require(beta_squared > Real(0) && beta_squared <= Real(1), "Invalid beta squared");

                for (auto const fraction : {Real(0), Real(.1), Real(.5), Real(.9), Real(1)}) {
                    auto const secondary = fraction * maximum;
                    auto const factor = ProtonImpactIonization::bhabhaSpinHalfFactor(
                        energy, mass, secondary);
                    // Tree-level spin-averaged trace for distinguishable Dirac particles:
                    // [a^2+b^2+2*t*(M^2+m_e^2)]/(8*m_e^2*E_total^2).
                    // a=s-M^2-m_e^2; b=u-M^2-m_e^2; t=-2*m_e*T.
                    auto const me = static_cast<long double>(static_cast<Real>(electron_mass));
                    auto const transfer = static_cast<long double>(secondary);
                    auto const total = e + m;
                    auto const a = 2.L * me * total;
                    auto const b = -2.L * me * (total - transfer);
                    auto const t = -2.L * me * transfer;
                    auto const trace = (a * a + b * b + 2.L * t * (m * m + me * me)) /
                                       (8.L * me * me * total * total);
                    require(std::isfinite(factor) && factor > Real(0),
                            "Non-positive Bhabha factor");
                    // A float-rounded Tmax shifts the reference endpoint very slightly.
                    auto const endpoint_error = std::abs(maximum / reference - 1.L);
                    require(std::abs(factor - trace) < factor_tolerance + endpoint_error,
                            "Bhabha trace mismatch");
                }
            }
        }

        auto const massive = static_cast<Real>(238.L * proton_mass);
        auto const slow = ProtonImpactIonization::relativisticMaximumTransfer(Real(1000), massive);
        require(slow > Real(.009) && slow < Real(.010), "Slow heavy-ion transfer lost to rounding");
        require(ProtonImpactIonization::relativisticMaximumTransfer(Real(0), massive) == Real(0),
                "Nonzero transfer at rest");
    }

    template <typename Real>
    void
    checkMolecularEndpoint ()
    {
        constexpr auto tolerance = std::is_same_v<Real, float> ? 3.e-6L : 3.e-13L;
        constexpr long double proton_mass = 938272088.16L;
        constexpr long double atomic_mass = 931494103.72L;
        auto const me = static_cast<long double>(Real(510998.95069));

        for (auto const mass_number : {1.L, 4.L, 238.L}) {
            auto const mass = static_cast<Real>(mass_number * proton_mass);
            for (auto const target_number : {28.0134L, 31.9988L}) {
                auto const neutral = static_cast<Real>(target_number * atomic_mass);
                for (auto const binding : {Real(12.07), Real(15.59), Real(543.5)}) {
                    auto const m = static_cast<long double>(mass);
                    auto const a = static_cast<long double>(neutral);
                    auto const i = static_cast<long double>(binding);
                    auto const threshold = i * (1.L + m / a) + i * i / (2.L * a);
                    require(ProtonImpactIonization::molecularMaximumSecondaryEnergy(
                                Real(.99L * threshold), mass, neutral, binding) == Real(0),
                            "Nonzero molecular endpoint below reaction threshold");
                    for (int index = 0; index < 100; ++index) {
                        auto const energy = static_cast<Real>(
                            threshold * (1.001L + std::pow(10.L, index / 10.L)));
                        auto const e = static_cast<long double>(energy);
                        auto const maximum =
                            ProtonImpactIonization::molecularMaximumSecondaryEnergy(
                                energy, mass, neutral, binding);
                        // Independently solve the invariant recoil-mass quadratic.
                        // K*T - p*sqrt(T*(T+2*m_e)) = C for forward emission.
                        auto const k = e + m + a;
                        auto const p2 = e * (e + 2.L * m);
                        auto const s = (m + a) * (m + a) + 2.L * a * e;
                        auto const c = e * (a - me) - i * (m + a - me) - .5L * i * i;
                        auto const reference = (k * c + p2 * me +
                            std::sqrt(p2 * (c * c + 2.L * k * c * me + p2 * me * me))) / s;
                        require(std::isfinite(maximum) && maximum > Real(0),
                                "Invalid molecular endpoint");
                        require(std::abs(maximum / reference - 1.L) < tolerance,
                                "Molecular endpoint disagrees with invariant quadratic");
                        require(maximum <= (e - i) * (1.L + tolerance),
                                "Molecular endpoint violates energy conservation");
                    }
                }
            }
        }

        auto const proton = static_cast<Real>(proton_mass);
        auto const nitrogen = static_cast<Real>(28.0134L * atomic_mass);
        auto const maximum = ProtonImpactIonization::molecularMaximumSecondaryEnergy(
            Real(5000), proton, nitrogen, Real(15.59));
        auto const free_maximum =
            ProtonImpactIonization::relativisticMaximumTransfer(Real(5000), proton);
        require(maximum > Real(4800) && maximum < Real(5000), "Incorrect N2 molecular endpoint");
        require(free_maximum > Real(10) && free_maximum < Real(11), "Incorrect binary endpoint");

        // Removing the residual ion and binding recovers the free-target limit.
        for (auto const energy : {Real(1000), Real(8.e8), Real(1.e12)}) {
            auto const molecular = ProtonImpactIonization::molecularMaximumSecondaryEnergy(
                energy, proton, Real(510998.95069), Real(0));
            auto const free = ProtonImpactIonization::relativisticMaximumTransfer(energy, proton);
            require(std::abs(static_cast<long double>(molecular) / free - 1.L) < tolerance,
                    "Molecular endpoint does not recover the free-target limit");
        }
    }
} // namespace

int
main ()
{
    // These host checks do not initialize MPI or require a GPU. Backend builds
    // also compile the host/device functions; device execution is a separate test.
    checkPrecision<float>();
    checkPrecision<double>();
    checkMolecularEndpoint<float>();
    checkMolecularEndpoint<double>();
    std::cout << "PASS: 2888 free-collision and 3600 molecular states, float and double\n";
}
