/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_KINEMATICS_H_
#define WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_KINEMATICS_H_

#include "BackgroundMCCKinematics.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Random.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace BackgroundMCCElasticKinematics
{
    /** Largest internal-energy increase accessible at a neutral-rest-frame angle.
     *
     * Energies and masses (mc^2) are in eV. Maximizing the outgoing molecular
     * invariant mass over the outgoing electron energy gives this bound. It is
     * a kinematic restriction, not the spectator model's rotational-rainbow scale.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    maximumInternalEnergyChange (double energy, double electron, double mass,
                                 double cosine) noexcept
    {
        double const total = energy + electron + mass;
        double const momentum2 = energy * (energy + 2 * electron);
        double const projected2 = cosine > 0 ? momentum2 * cosine * cosine : 0;
        double const root = std::sqrt(total * total - projected2);
        double const shift =
            (energy * (2 * (electron + mass) + energy) - projected2) / (root + electron + mass);
        double const difference = 2 * (mass * energy - electron * shift);
        return difference / (std::sqrt(mass * mass + difference) + mass);
    }

    /** Signed excitation at the laboratory angle sampled from the elastic DCS.
     *
     * The outgoing electron mass-shell quadratic uses its larger-energy branch
     * when both forward solutions exist in the narrow recoil threshold region.
     * Exact accessibility is checked before solving; no angle is projected.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
    computeInternalEnergyChangeLab (
        amrex::ParticleReal incident_ux, amrex::ParticleReal incident_uy,
        amrex::ParticleReal incident_uz, amrex::ParticleReal neutral_vx,
        amrex::ParticleReal neutral_vy, amrex::ParticleReal neutral_vz, double electron_mass,
        double initial_mass, double change, double cosine, amrex::RandomEngine const& engine,
        amrex::ParticleReal& electron_ux, amrex::ParticleReal& electron_uy,
        amrex::ParticleReal& electron_uz, amrex::ParticleReal& recoil_ux,
        amrex::ParticleReal& recoil_uy, amrex::ParticleReal& recoil_uz) noexcept
    {
        using namespace BackgroundMCCKinematics;
        constexpr double c = PhysConst::c_v<double>;
        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double qe = PhysConst::q_e_v<double>;
        Vector3 const velocity{neutral_vx, neutral_vy, neutral_vz};
        double const gamma = 1 / std::sqrt(1 - dot(velocity, velocity) / c2);
        double const boost = gamma * gamma / ((gamma + 1) * c2);
        auto const incoming =
            boostToRest({incident_ux, incident_uy, incident_uz}, velocity, gamma, boost);
        double const energy = kineticEnergyFromProperVelocity(incoming, electron_mass);
        double const m = electron_mass * c2 / qe, mass = initial_mass * c2 / qe;
        if (mass + change <= 0 || change > maximumInternalEnergyChange(energy, m, mass, cosine)) {
            return false;
        }
        double const p = std::sqrt(energy * (energy + 2 * m));
        double const total = energy + m + mass;
        double const projection = p * cosine;
        double const a = total * total - projection * projection;
        double const d = mass * (energy - change) - m * energy - change * change / 2;
        double const discriminant = d * (d + 2 * total * m) + projection * projection * m * m;
        if (discriminant < 0) {
            return false;
        }
        double const center = total * d + projection * projection * m;
        double const splitting = projection * std::sqrt(discriminant);
        double const outgoing =
            projection < 0 && d >= 0 ? d * d / (center - splitting) : (center + splitting) / a;
        if (!(outgoing >= 0)) {
            return false;
        }
        double const norm = std::sqrt(dot(incoming, incoming));
        Vector3 const axis = norm > 0 ? scale(incoming, 1 / norm) : Vector3{0, 0, 1};
        Vector3 transverse1, transverse2;
        transverseDirections(axis, transverse1, transverse2);
        auto const direction =
            directionFromPolarAngle(axis, transverse1, transverse2, cosine,
                                    2 * static_cast<double>(MathConst::pi) * amrex::Random(engine));
        auto const electron_pc = scale(direction, std::sqrt(outgoing * (outgoing + 2 * m)));
        auto const recoil_pc = subtract(scale(axis, p), electron_pc);
        auto const electron_lab = boostToLab(scale(electron_pc, c / m), velocity, gamma, boost);
        auto const recoil_lab =
            boostToLab(scale(recoil_pc, c / (mass + change)), velocity, gamma, boost);
        electron_ux = static_cast<amrex::ParticleReal>(electron_lab.x);
        electron_uy = static_cast<amrex::ParticleReal>(electron_lab.y);
        electron_uz = static_cast<amrex::ParticleReal>(electron_lab.z);
        recoil_ux = static_cast<amrex::ParticleReal>(recoil_lab.x);
        recoil_uy = static_cast<amrex::ParticleReal>(recoil_lab.y);
        recoil_uz = static_cast<amrex::ParticleReal>(recoil_lab.z);
        return true;
    }

/** Scatter at a COM angle with a signed change of molecular internal energy.
 *
 * Positive changes excite the target; negative changes release its stored
 * rotational energy. Factored phase space avoids rest-energy cancellation.
 * Unlike the legacy laboratory-angle excitation closure, this transformation
 * needs no angular projection. Both products conserve four-momentum.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
computeInternalEnergyChange (amrex::ParticleReal incident_ux, amrex::ParticleReal incident_uy,
                             amrex::ParticleReal incident_uz, amrex::ParticleReal neutral_vx,
                             amrex::ParticleReal neutral_vy, amrex::ParticleReal neutral_vz,
                             double electron_mass, double initial_mass, double change,
                             double cosine, amrex::RandomEngine const& engine,
                             amrex::ParticleReal& electron_ux, amrex::ParticleReal& electron_uy,
                             amrex::ParticleReal& electron_uz, amrex::ParticleReal& recoil_ux,
                             amrex::ParticleReal& recoil_uy,
                             amrex::ParticleReal& recoil_uz) noexcept
{
    using namespace BackgroundMCCKinematics;
    constexpr double c = PhysConst::c_v<double>;
    constexpr double c2 = PhysConst::c2_v<double>;
    constexpr double qe = PhysConst::q_e_v<double>;
    Vector3 const velocity{neutral_vx, neutral_vy, neutral_vz};
    double const beta2 = dot(velocity, velocity) / c2;
    if (!(beta2 < 1)) {
        return false;
    }
    double const gamma = 1 / std::sqrt(1 - beta2);
    double const boost_factor = gamma * gamma / ((gamma + 1) * c2);
    auto const incoming =
        boostToRest({incident_ux, incident_uy, incident_uz}, velocity, gamma, boost_factor);
    double const energy = kineticEnergyFromProperVelocity(incoming, electron_mass);
    double const m = electron_mass * c2 / qe;
    double const mass = initial_mass * c2 / qe;
    double const final_mass = initial_mass + change * qe / c2;
    double const threshold = change * (1 + m / mass) + change * change / (2 * mass);
    if (!(final_mass > 0) || energy < threshold) {
        return false;
    }
    double const s = (m + mass) * (m + mass) + 2 * mass * energy;
    double const phase1 = 2 * mass * (energy - threshold);
    double const phase2 =
        4 * mass * m + 2 * mass * energy - 2 * change * (mass - m) - change * change;
    if (phase1 < 0 || phase2 < 0) {
        return false;
    }
    double const pc = std::sqrt(phase1 * phase2 / (4 * s));
    double const incoming_norm = std::sqrt(dot(incoming, incoming));
    Vector3 direction{0, 0, 1};
    if (incoming_norm > 0) {
        direction = scale(incoming, 1 / incoming_norm);
    }
    Vector3 transverse1, transverse2;
    transverseDirections(direction, transverse1, transverse2);
    auto const outgoing_direction = directionFromPolarAngle(
        direction, transverse1, transverse2, amrex::max(-1.0, amrex::min(cosine, 1.0)),
        2 * static_cast<double>(MathConst::pi) * amrex::Random(engine));
    auto const electron_com = scale(outgoing_direction, c * pc / m);
    auto const recoil_com = scale(electron_com, -electron_mass / final_mass);
    auto const com_velocity = scale(incoming, m / (energy + m + mass));
    double const com_gamma = (energy + m + mass) / std::sqrt(s);
    double const com_boost = com_gamma * com_gamma / ((com_gamma + 1) * c2);
    auto const electron_rest = boostToLab(electron_com, com_velocity, com_gamma, com_boost);
    auto const recoil_rest = boostToLab(recoil_com, com_velocity, com_gamma, com_boost);
    auto const electron_lab = boostToLab(electron_rest, velocity, gamma, boost_factor);
    auto const recoil_lab = boostToLab(recoil_rest, velocity, gamma, boost_factor);
    electron_ux = static_cast<amrex::ParticleReal>(electron_lab.x);
    electron_uy = static_cast<amrex::ParticleReal>(electron_lab.y);
    electron_uz = static_cast<amrex::ParticleReal>(electron_lab.z);
    recoil_ux = static_cast<amrex::ParticleReal>(recoil_lab.x);
    recoil_uy = static_cast<amrex::ParticleReal>(recoil_lab.y);
    recoil_uz = static_cast<amrex::ParticleReal>(recoil_lab.z);
    return true;
}

    /** Apply an elastic electron-neutral collision in the neutral rest frame.
     *
     * The sampled polar angle is the electron angle represented by the
     * laboratory-frame IAA/elmolcs DCS for a stationary target. Relativistic
     * two-body kinematics set the outgoing electron energy, and the neutral
     * receives the exact recoil momentum. Both products are then boosted back
     * to the simulation frame.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute (amrex::ParticleReal const incident_ux, amrex::ParticleReal const incident_uy,
             amrex::ParticleReal const incident_uz, amrex::ParticleReal const neutral_vx,
             amrex::ParticleReal const neutral_vy, amrex::ParticleReal const neutral_vz,
             double const electron_mass, double const neutral_mass,
             double const scattering_cosine, amrex::RandomEngine const& engine,
             amrex::ParticleReal& electron_ux, amrex::ParticleReal& electron_uy,
             amrex::ParticleReal& electron_uz, amrex::ParticleReal& recoil_ux,
             amrex::ParticleReal& recoil_uy, amrex::ParticleReal& recoil_uz) noexcept
    {
        using BackgroundMCCKinematics::Vector3;
        using BackgroundMCCKinematics::boostToLab;
        using BackgroundMCCKinematics::boostToRest;
        using BackgroundMCCKinematics::directionFromPolarAngle;
        using BackgroundMCCKinematics::dot;
        using BackgroundMCCKinematics::kineticEnergyFromProperVelocity;
        using BackgroundMCCKinematics::properSpeedFromKineticEnergy;
        using BackgroundMCCKinematics::scale;
        using BackgroundMCCKinematics::subtract;
        using BackgroundMCCKinematics::transverseDirections;

        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double elementary_charge = PhysConst::q_e_v<double>;
        constexpr double two_pi = 2.0 * static_cast<double>(MathConst::pi);

        Vector3 const neutral_velocity{static_cast<double>(neutral_vx),
                                       static_cast<double>(neutral_vy),
                                       static_cast<double>(neutral_vz)};
        Vector3 const incident_lab{static_cast<double>(incident_ux),
                                   static_cast<double>(incident_uy),
                                   static_cast<double>(incident_uz)};

        auto const neutral_speed_sq = dot(neutral_velocity, neutral_velocity);
        auto const neutral_gamma =
            1.0 / std::sqrt(amrex::max(1.0 - neutral_speed_sq * inv_c2,
                                       std::numeric_limits<double>::epsilon()));
        auto const boost_spatial_factor =
            neutral_gamma * neutral_gamma / ((neutral_gamma + 1.0) * c2);
        auto const incident_rest = boostToRest(
            incident_lab, neutral_velocity, neutral_gamma, boost_spatial_factor);

        auto const incident_speed = std::sqrt(dot(incident_rest, incident_rest));
        Vector3 incident_direction{0.0, 0.0, 1.0};
        if (incident_speed > std::numeric_limits<double>::min())
        {
            incident_direction = scale(incident_rest, 1.0 / incident_speed);
        }
        Vector3 transverse_1;
        Vector3 transverse_2;
        transverseDirections(incident_direction, transverse_1, transverse_2);

        auto const cosine = amrex::max(-1.0, amrex::min(scattering_cosine, 1.0));
        auto const sine_sq = amrex::max(1.0 - cosine * cosine, 0.0);
        auto const incident_energy =
            kineticEnergyFromProperVelocity(incident_rest, electron_mass);
        auto const electron_rest_energy = electron_mass * c2 / elementary_charge;
        auto const neutral_rest_energy = neutral_mass * c2 / elementary_charge;

        // This stable form avoids subtracting two nearly equal neutral rest
        // energies for strongly forward scattering from a heavy molecule.
        auto const recoil_root = std::sqrt(amrex::max(
            neutral_rest_energy * neutral_rest_energy -
                electron_rest_energy * electron_rest_energy * sine_sq,
            0.0));
        auto const recoil_factor =
            (incident_energy + electron_rest_energy) * sine_sq +
            neutral_rest_energy * (1.0 - cosine) +
            cosine * electron_rest_energy * electron_rest_energy * sine_sq /
                (neutral_rest_energy + recoil_root);
        auto const incident_pc_sq =
            incident_energy * (incident_energy + 2.0 * electron_rest_energy);
        auto const total_energy =
            incident_energy + electron_rest_energy + neutral_rest_energy;
        auto const denominator =
            total_energy * total_energy - incident_pc_sq * cosine * cosine;
        auto const energy_loss = denominator > 0.0
            ? recoil_factor * incident_pc_sq / denominator
            : incident_energy;
        auto const outgoing_energy =
            amrex::max(incident_energy - energy_loss, 0.0);

        auto const azimuth = two_pi * static_cast<double>(amrex::Random(engine));
        auto const outgoing_direction = directionFromPolarAngle(
            incident_direction, transverse_1, transverse_2, cosine, azimuth);
        auto const electron_rest = scale(
            outgoing_direction,
            properSpeedFromKineticEnergy(outgoing_energy, electron_mass));
        auto const recoil_rest = scale(
            subtract(incident_rest, electron_rest), electron_mass / neutral_mass);

        auto const electron_lab =
            boostToLab(electron_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        auto const recoil_lab =
            boostToLab(recoil_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        electron_ux = static_cast<amrex::ParticleReal>(electron_lab.x);
        electron_uy = static_cast<amrex::ParticleReal>(electron_lab.y);
        electron_uz = static_cast<amrex::ParticleReal>(electron_lab.z);
        recoil_ux = static_cast<amrex::ParticleReal>(recoil_lab.x);
        recoil_uy = static_cast<amrex::ParticleReal>(recoil_lab.y);
        recoil_uz = static_cast<amrex::ParticleReal>(recoil_lab.z);
    }

    /** Apply a discrete-loss electron-neutral collision at a sampled laboratory angle.
     *
     * The configured loss is represented as internal energy of the recoiling target.
     * Solving the two-body energy equation analytically avoids a device-side root solve.
     * Above the recoil-shifted threshold, a sampled angle outside the allowed two-body
     * laboratory range is projected to the nearest kinematically allowed angle.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    computeExcitation (
        amrex::ParticleReal const incident_ux, amrex::ParticleReal const incident_uy,
        amrex::ParticleReal const incident_uz, amrex::ParticleReal const neutral_vx,
        amrex::ParticleReal const neutral_vy, amrex::ParticleReal const neutral_vz,
        double const electron_mass, double const neutral_mass,
        amrex::ParticleReal const excitation_energy,
        double const scattering_cosine, amrex::RandomEngine const& engine,
        amrex::ParticleReal& electron_ux, amrex::ParticleReal& electron_uy,
        amrex::ParticleReal& electron_uz, amrex::ParticleReal& recoil_ux,
        amrex::ParticleReal& recoil_uy, amrex::ParticleReal& recoil_uz) noexcept
    {
        using BackgroundMCCKinematics::Vector3;
        using BackgroundMCCKinematics::boostToLab;
        using BackgroundMCCKinematics::boostToRest;
        using BackgroundMCCKinematics::directionFromPolarAngle;
        using BackgroundMCCKinematics::dot;
        using BackgroundMCCKinematics::kineticEnergyFromProperVelocity;
        using BackgroundMCCKinematics::properSpeedFromKineticEnergy;
        using BackgroundMCCKinematics::scale;
        using BackgroundMCCKinematics::subtract;
        using BackgroundMCCKinematics::transverseDirections;

        constexpr double c2 = PhysConst::c2_v<double>;
        constexpr double inv_c2 = PhysConst::inv_c2_v<double>;
        constexpr double elementary_charge = PhysConst::q_e_v<double>;
        constexpr double two_pi = 2.0 * static_cast<double>(MathConst::pi);

        Vector3 const neutral_velocity{static_cast<double>(neutral_vx),
                                       static_cast<double>(neutral_vy),
                                       static_cast<double>(neutral_vz)};
        Vector3 const incident_lab{static_cast<double>(incident_ux),
                                   static_cast<double>(incident_uy),
                                   static_cast<double>(incident_uz)};

        auto const neutral_speed_sq = dot(neutral_velocity, neutral_velocity);
        auto const neutral_gamma =
            1.0 / std::sqrt(amrex::max(1.0 - neutral_speed_sq * inv_c2,
                                       std::numeric_limits<double>::epsilon()));
        auto const boost_spatial_factor =
            neutral_gamma * neutral_gamma / ((neutral_gamma + 1.0) * c2);
        auto const incident_rest = boostToRest(
            incident_lab, neutral_velocity, neutral_gamma, boost_spatial_factor);

        auto const incident_speed = std::sqrt(dot(incident_rest, incident_rest));
        Vector3 incident_direction{0.0, 0.0, 1.0};
        if (incident_speed > std::numeric_limits<double>::min())
        {
            incident_direction = scale(incident_rest, 1.0 / incident_speed);
        }
        Vector3 transverse_1;
        Vector3 transverse_2;
        transverseDirections(incident_direction, transverse_1, transverse_2);

        auto const incident_kinetic_energy =
            kineticEnergyFromProperVelocity(incident_rest, electron_mass);
        auto const loss = amrex::max(static_cast<double>(excitation_energy), 0.0);
        auto const electron_rest_energy = electron_mass * c2 / elementary_charge;
        auto const neutral_rest_energy = neutral_mass * c2 / elementary_charge;
        auto const incident_total_energy = incident_kinetic_energy + electron_rest_energy;
        auto const incident_pc = std::sqrt(amrex::max(
            incident_kinetic_energy *
                (incident_kinetic_energy + 2.0 * electron_rest_energy),
            0.0));

        // These factored Kallen-function terms remain accurate close to the
        // recoil-shifted excitation threshold, where direct squared energies cancel.
        auto const phase_space_1 = amrex::max(
            2.0 * neutral_rest_energy * (incident_kinetic_energy - loss) -
                2.0 * electron_rest_energy * loss - loss * loss,
            0.0);
        auto const phase_space_2 = amrex::max(
            4.0 * neutral_rest_energy * electron_rest_energy +
                2.0 * neutral_rest_energy * incident_kinetic_energy -
                2.0 * loss * (neutral_rest_energy - electron_rest_energy) - loss * loss,
            0.0);
        auto const kallen_quarter = 0.25 * phase_space_1 * phase_space_2;

        auto cosine = amrex::max(-1.0, amrex::min(scattering_cosine, 1.0));
        auto const angular_scale = incident_pc * incident_pc *
            electron_rest_energy * electron_rest_energy;
        if (angular_scale > 0.0 && kallen_quarter < angular_scale)
        {
            auto const minimum_cosine =
                std::sqrt(amrex::max(1.0 - kallen_quarter / angular_scale, 0.0));
            // Squaring the energy equation also admits a backward cone with
            // negative outgoing momentum. Only the forward cone is physical.
            cosine = amrex::max(cosine, minimum_cosine);
        }

        auto const total_energy = incident_total_energy + neutral_rest_energy;
        auto const coefficient = total_energy * total_energy -
            incident_pc * incident_pc * cosine * cosine;
        // Half of s + m_e^2 c^4 - M_excited^2 c^4, evaluated without
        // subtracting the nearly equal neutral and excited-neutral rest energies.
        auto const linear_term =
            electron_rest_energy * electron_rest_energy +
            neutral_rest_energy * (incident_total_energy - loss) - 0.5 * loss * loss;
        auto const radicand = amrex::max(
            kallen_quarter - angular_scale * (1.0 - cosine * cosine), 0.0);
        // Solve for pc instead of subtracting m_e*c^2 from the outgoing total
        // energy. Rationalize the backward root to avoid cancellation there
        // too.
        auto const projected_pc = incident_pc * cosine;
        auto const root = std::sqrt(radicand);
        auto const linear_minus_rest =
            (neutral_rest_energy - electron_rest_energy) * (incident_kinetic_energy - loss) -
            electron_rest_energy * loss - 0.5 * loss * loss;
        auto const backward_denominator = total_energy * root - projected_pc * linear_term;
        auto const outgoing_pc =
            cosine >= 0.0
                ? (projected_pc * linear_term + total_energy * root) / coefficient
                : (backward_denominator > 0.0
                       ? linear_minus_rest * (linear_term + total_energy * electron_rest_energy) /
                             backward_denominator
                       : 0.0);
        auto const outgoing_pc_sq = outgoing_pc * outgoing_pc;
        auto const outgoing_kinetic_energy =
            outgoing_pc_sq /
            (std::sqrt(electron_rest_energy * electron_rest_energy + outgoing_pc_sq) +
             electron_rest_energy);

        auto const azimuth = two_pi * static_cast<double>(amrex::Random(engine));
        auto const outgoing_direction = directionFromPolarAngle(
            incident_direction, transverse_1, transverse_2, cosine, azimuth);
        auto const electron_rest = scale(
            outgoing_direction,
            properSpeedFromKineticEnergy(outgoing_kinetic_energy, electron_mass));
        auto const excited_mass = neutral_mass + loss * elementary_charge / c2;
        auto const recoil_rest = scale(
            subtract(incident_rest, electron_rest), electron_mass / excited_mass);

        auto const electron_lab =
            boostToLab(electron_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        auto const recoil_lab =
            boostToLab(recoil_rest, neutral_velocity, neutral_gamma, boost_spatial_factor);
        electron_ux = static_cast<amrex::ParticleReal>(electron_lab.x);
        electron_uy = static_cast<amrex::ParticleReal>(electron_lab.y);
        electron_uz = static_cast<amrex::ParticleReal>(electron_lab.z);
        recoil_ux = static_cast<amrex::ParticleReal>(recoil_lab.x);
        recoil_uy = static_cast<amrex::ParticleReal>(recoil_lab.y);
        recoil_uz = static_cast<amrex::ParticleReal>(recoil_lab.z);
    }
} // namespace BackgroundMCCElasticKinematics

#endif // WARPX_PARTICLES_COLLISION_BACKGROUNDMCC_ELASTIC_KINEMATICS_H_
