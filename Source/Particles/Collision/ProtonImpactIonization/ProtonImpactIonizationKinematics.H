/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_

#include "Particles/Collision/BackgroundMCC/BackgroundMCCKinematics.H"
#include "RelativisticIonization.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>
#include <limits>

namespace ProtonImpactIonization
{
    /** Free-electron binary-encounter polar cosine for a heavy projectile. */
    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    freeElectronCosine (Real const secondary_energy, Real const maximum_transfer) noexcept
    {
        using std::sqrt;

        if (secondary_energy <= Real(0) || maximum_transfer <= Real(0)) {
            return Real(0);
        }
        constexpr auto electron_rest_energy = Real(510998.95069);
        // The energy model has a molecular tail. Its free reference direction
        // is forward above Tmax; do not evaluate free kinematics out of domain.
        auto const free_energy = amrex::min(secondary_energy, maximum_transfer);
        return sqrt(amrex::min(
            Real(1), free_energy * (maximum_transfer + Real(2) * electron_rest_energy) /
                         (maximum_transfer * (free_energy + Real(2) * electron_rest_energy))));
    }

    struct AngularInterval
    {
        double m_lower = -1;
        double m_upper = 1;
        double m_minimum = -1;
        bool m_valid = false;
    };

    /** Heavy-projectile extension of the bound/free interpolation in Eq. (2.66).
     *
     * All masses are rest energies and all energies are in eV. Binding is paid
     * before constructing the free-binary reference. Above its endpoint the
     * reference is forward; the molecular SDCS tail remains intact. Condition
     * the angular interval on the exact projectile-plus-ion invariant mass.
     * This phenomenological closure does not constitute a molecular DDCS.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE AngularInterval
    angularInterval (double secondary, double binding, double incident, double projectile,
                     double neutral) noexcept
    {
        constexpr double electron = 510998.95069;
        double const threshold =
            binding * (1 + projectile / neutral) + binding * binding / (2 * neutral);
        if (incident < threshold || secondary < 0 || incident <= binding) {
            return {};
        }
        double const maximum = relativisticMaximumTransfer(incident - binding, projectile);
        double const denominator = secondary + binding;
        double const weight = denominator > 0 ? secondary / denominator : 0;
        double const center = weight * freeElectronCosine(secondary, maximum);
        double const width = 1 - weight;
        double const p = std::sqrt(incident * (incident + 2 * projectile));
        double const pe = std::sqrt(secondary * (secondary + 2 * electron));
        double const total = incident + projectile + neutral;
        double const excess = incident * (neutral - electron) -
                              binding * (projectile + neutral - electron) - binding * binding / 2;
        double const minimum = pe > 0 ? (total * secondary - excess) / (p * pe) : -1;
        double const lower = amrex::max(center - width, minimum);
        double const upper = center + width;
        return {lower, upper, minimum, lower <= upper && lower >= -1 && upper <= 1};
    }

    template <typename Real>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
    polarCosine (Real secondary, Real binding, Real incident, Real projectile, Real neutral,
                 Real uniform) noexcept
    {
        auto const interval = angularInterval(secondary, binding, incident, projectile, neutral);
        return interval.m_valid ? Real((1 - double(uniform)) * interval.m_lower +
                                       double(uniform) * interval.m_upper)
                                : std::numeric_limits<Real>::quiet_NaN();
    }

    using BackgroundMCCKinematics::Vector3;

    struct NeutralFrame
    {
        Vector3 m_velocity{};
        Vector3 m_direction{0, 0, 1};
        Vector3 m_transverse1{}, m_transverse2{};
        double m_gamma = 1;
        double m_boost = 0;
        double m_energy = 0;
    };

    /** Transform the actual incoming projectile to the sampled neutral frame. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE NeutralFrame
    neutralFrame (Vector3 incident, Vector3 velocity, double projectile_mass) noexcept
    {
        using namespace BackgroundMCCKinematics;
        NeutralFrame frame;
        frame.m_velocity = velocity;
        frame.m_gamma = 1 / std::sqrt(1 - dot(velocity, velocity) / PhysConst::c2_v<double>);
        frame.m_boost =
            frame.m_gamma * frame.m_gamma / ((frame.m_gamma + 1) * PhysConst::c2_v<double>);
        auto const incoming = boostToRest(incident, velocity, frame.m_gamma, frame.m_boost);
        auto const norm = std::sqrt(dot(incoming, incoming));
        if (norm > 0) {
            frame.m_direction = scale(incoming, 1 / norm);
        }
        transverseDirections(frame.m_direction, frame.m_transverse1, frame.m_transverse2);
        frame.m_energy = kineticEnergyFromProperVelocity(incoming, projectile_mass);
        return frame;
    }

    struct Products
    {
        Vector3 m_electron{}, m_projectile{}, m_ion{};
        double m_ion_rest_energy = 0;
        bool m_valid = false;
    };

    /** Construct an on-shell electron, projectile and ion before reservoir discard.
     *
     * The sampled electron energy is preserved apart from an input-precision endpoint guard. The
     * remaining two-body system is partitioned at minimum ion kinetic energy in the neutral frame.
     * This selects one recoil closure without another random draw or iterative solve. The binding
     * defect is included in the ion rest energy. A rigid projectile may discard its outgoing
     * momentum after this construction.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Products
    compute (
        NeutralFrame const& frame, double secondary, double binding, double projectile_mass,
        double neutral_mass, double uniform, double azimuth,
        double relative_precision = std::numeric_limits<amrex::ParticleReal>::epsilon()) noexcept
    {
        using namespace BackgroundMCCKinematics;
        constexpr double c = PhysConst::c_v<double>;
        constexpr double rest_factor = PhysConst::c2_v<double> / PhysConst::q_e_v<double>;
        constexpr double electron = 510998.95069;
        double const projectile = projectile_mass * rest_factor;
        double const neutral = neutral_mass * rest_factor;
        double const ion = neutral - electron + binding;
        double const incident = frame.m_energy;
        auto interval = angularInterval(secondary, binding, incident, projectile, neutral);
        if (!interval.m_valid) {
            double const endpoint =
                molecularMaximumSecondaryEnergy(incident, projectile, neutral, binding);
            // Input/table conversion can round a molecular endpoint upward.
            // Repair only a few input ulps, never a finite portion of the tail.
            if (secondary >= endpoint &&
                secondary - endpoint <= 8 * relative_precision * amrex::max(incident, secondary)) {
                secondary = endpoint * (1 - 32 * std::numeric_limits<double>::epsilon());
                interval = angularInterval(secondary, binding, incident, projectile, neutral);
            }
        }
        if (!interval.m_valid || ion <= 0) {
            return {};
        }
        double const cosine = interval.m_lower + uniform * (interval.m_upper - interval.m_lower);
        auto const direction = directionFromPolarAngle(frame.m_direction, frame.m_transverse1,
                                                       frame.m_transverse2, cosine, azimuth);
        double const p = std::sqrt(incident * (incident + 2 * projectile));
        double const pe = std::sqrt(secondary * (secondary + 2 * electron));
        auto const electron_pc = scale(direction, pe);
        auto const remainder_pc = subtract(scale(frame.m_direction, p), electron_pc);
        double const remainder_p = std::sqrt(dot(remainder_pc, remainder_pc));
        double const tau = incident - secondary - binding;
        double const total = projectile + ion + tau;
        // Factored invariant excess avoids subtraction of molecular rest masses.
        double const excess =
            pe > 0 ? 2 * p * pe * (cosine - interval.m_minimum)
                   : 2 * (incident * (neutral - electron) -
                          binding * (projectile + neutral - electron) - binding * binding / 2);
        if (excess < 0 || remainder_p <= 0) {
            return {};
        }
        double const invariant = (projectile + ion) * (projectile + ion) + excess;
        double const root = std::sqrt(invariant);
        double const star = std::sqrt(excess * (excess + 4 * projectile * ion) / (4 * invariant));
        // Rationalize the small difference between the boosted ion momentum
        // and its COM momentum. The free-electron limit then has zero ion recoil.
        double const defect = (incident + projectile + electron) * secondary +
                              (incident + projectile - secondary) * binding -
                              binding * binding / 2 - p * pe * cosine;
        double const moment = ion * total - defect;
        double const denominator = moment * remainder_p + total * root * star;
        if (denominator <= 0) {
            return {};
        }
        double const ion_p = defect * (2 * ion * total - defect) / denominator;
        auto const ion_pc = scale(remainder_pc, ion_p / remainder_p);
        auto const projectile_pc = subtract(remainder_pc, ion_pc);
        Products products;
        products.m_electron = boostToLab(scale(electron_pc, c / electron), frame.m_velocity,
                                         frame.m_gamma, frame.m_boost);
        products.m_projectile = boostToLab(scale(projectile_pc, c / projectile), frame.m_velocity,
                                           frame.m_gamma, frame.m_boost);
        products.m_ion =
            boostToLab(scale(ion_pc, c / ion), frame.m_velocity, frame.m_gamma, frame.m_boost);
        products.m_ion_rest_energy = ion;
        products.m_valid = true;
        return products;
    }
} // namespace ProtonImpactIonization

#endif // WARPX_PARTICLES_COLLISION_PROTONIMPACTIONIZATION_KINEMATICS_H_
