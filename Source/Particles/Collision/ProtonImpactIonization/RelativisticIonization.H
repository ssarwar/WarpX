/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_RELATIVISTIC_IONIZATION_H_
#define WARPX_PARTICLES_COLLISION_RELATIVISTIC_IONIZATION_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

namespace ProtonImpactIonization
{
    /** Squared projectile speed in units of c, without subtracting nearby numbers.
     *
     * Both arguments are energies in the same units. The kinetic energy must
     * be nonnegative and the rest energy positive. The template permits
     * independent tests of both particle precisions in a single build.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    relativisticBetaSquared (Real const kinetic_energy, Real const rest_energy) noexcept
    {
        auto const u = kinetic_energy / rest_energy;
        auto const inverse_gamma = Real(1) / (Real(1) + u);
        if (u > Real(1)) {
            // No cancellation here; this also prevents a rounded beta^2 > 1.
            return Real(1) - inverse_gamma * inverse_gamma;
        }
        return (u * inverse_gamma) * ((u + Real(2)) * inverse_gamma);
    }

    /** Maximum kinetic energy transferred to a free stationary electron, in eV.
     *
     * The projectile kinetic energy is nonnegative and its rest energy is
     * positive, both in eV. This is not the maximum energy of an electron
     * emitted from a molecule with a recoiling residual ion.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    relativisticMaximumTransfer (Real const kinetic_energy, Real const rest_energy) noexcept
    {
        constexpr auto electron_rest_energy = Real(510998.95069);
        auto const u = kinetic_energy / rest_energy;
        auto const mass_ratio = electron_rest_energy / rest_energy;
        auto const denominator = Real(1) + Real(2) * (Real(1) + u) * mass_ratio +
                                 mass_ratio * mass_ratio;
        // beta^2 gamma^2 = u(u+2) also resolves slow, massive projectiles in float.
        return (Real(2) * electron_rest_energy * u) * ((u + Real(2)) / denominator);
    }

    /** Maximum secondary energy for ionization of a stationary molecule, in eV.
     *
     * All arguments are energies in eV. The positive projectile and neutral
     * rest energies and ionization energy specify the three-body final state:
     * projectile + electron + residual ion. The residual ion rest energy is
     * neutral_rest_energy - m_e*c^2 + ionization_energy. This bound assumes
     * that residual state and does not include additional excitation energy.
     * It is a kinematic bound, not a prediction of the bound-electron tail.
     * Below the reaction threshold the function returns zero.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    molecularMaximumSecondaryEnergy (Real const kinetic_energy, Real const rest_energy,
                                     Real const neutral_rest_energy,
                                     Real const ionization_energy) noexcept
    {
        using std::sqrt;

        constexpr auto electron_rest_energy = Real(510998.95069);
        auto const initial_rest_energy = rest_energy + neutral_rest_energy;
        auto const threshold = ionization_energy *
            (Real(1) + rest_energy / neutral_rest_energy) +
            ionization_energy * (ionization_energy / (Real(2) * neutral_rest_energy));
        if (kinetic_energy < threshold) {
            return Real(0);
        }

        auto const invariant_mass = sqrt(initial_rest_energy * initial_rest_energy +
                                         Real(2) * neutral_rest_energy * kinetic_energy);
        // Rationalize sqrt(s) - (M+A+I), resolving energies close to threshold.
        auto const excess = Real(2) * neutral_rest_energy * (kinetic_energy - threshold) /
            (invariant_mass + initial_rest_energy + ionization_energy);
        // The other two final particles co-move at the maximum electron energy.
        auto const cm_kinetic = excess *
            (Real(1) - (electron_rest_energy + Real(0.5) * excess) / invariant_mass);
        auto const cm_momentum = sqrt(cm_kinetic *
                                      (cm_kinetic + Real(2) * electron_rest_energy));
        auto const momentum_squared = kinetic_energy * (kinetic_energy + Real(2) * rest_energy);
        auto const gamma_cm_minus_one = momentum_squared /
            (invariant_mass * (kinetic_energy + initial_rest_energy + invariant_mass));
        return cm_kinetic + gamma_cm_minus_one * (cm_kinetic + electron_rest_energy) +
               sqrt(momentum_squared) / invariant_mass * cm_momentum;
    }

    /** Dimensionless Bhabha factor for a pointlike spin-1/2 projectile.
     *
     * Multiply by 2*pi*r_e^2*m_e*c^2/(beta^2*T^2) for the per-electron
     * SDCS. Energies are in eV, with E > 0, M*c^2 > 0, and 0 <= T <= Tmax.
     * The factored evaluation remains positive at the relativistic endpoint.
     * This is a tree-level free-electron result, not a molecular SDCS.
     */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    bhabhaSpinHalfFactor (Real const kinetic_energy, Real const rest_energy,
                          Real const secondary_energy) noexcept
    {
        auto const maximum = relativisticMaximumTransfer(kinetic_energy, rest_energy);
        auto const fraction = secondary_energy / maximum;
        auto const inverse_gamma = rest_energy / (kinetic_energy + rest_energy);
        auto const transfer_fraction = secondary_energy / (kinetic_energy + rest_energy);
        return (Real(1) - fraction) + fraction * inverse_gamma * inverse_gamma +
               Real(0.5) * transfer_fraction * transfer_fraction;
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_RELATIVISTIC_IONIZATION_H_
