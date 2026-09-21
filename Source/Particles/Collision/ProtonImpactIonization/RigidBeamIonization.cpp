/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ProtonImpactIonization.H"

#include "Fluids/WarpXFluidContainer.H"
#include "Particles/Collision/BackgroundMCC/BackgroundMCCKinematics.H"
#include "Particles/Collision/IonProductDestination.H"
#include "Particles/Collision/ProtonImpactIonization/IonizationSampling.H"
#include "Particles/Collision/ProtonImpactIonization/ProtonImpactIonizationKinematics.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleCreation/SmartCreate.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/profiler/ProfilerWrapper.H>

#include <AMReX_GpuAssert.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Scan.H>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>

#ifdef WARPX_DIM_RZ
namespace
{
    /** Fixed integer scrambling; no dependence on MPI rank or particle ordering. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE std::uint64_t
    scramble (std::uint64_t value) noexcept
    {
        value += 0x9e3779b97f4a7c15ull;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
        return value ^ (value >> 31);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    radialQuantile (double lo, double hi, double sigma, double uniform) noexcept
    {
        double const span = -std::expm1(-(hi - lo) * (hi + lo) / (2 * sigma * sigma));
        return std::sqrt(lo * lo - 2 * sigma * sigma * std::log1p(-uniform * span));
    }

    /** Locate a positive interval in a cumulative source table. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
    interval (double const* cdf, int count, double value) noexcept
    {
        int lo = 0, hi = count;
        while (lo + 1 < hi) {
            int const mid = lo + (hi - lo) / 2;
            if (cdf[mid] <= value) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return lo;
    }
} // namespace
#endif

void
ProtonImpactIonizationCollision::ProduceFromFluid (amrex::Real cur_time, amrex::Real start_time,
                                                   amrex::Real dt, MultiParticleContainer* mypc)
{
    ABLASTR_PROFILE("ProtonImpactIonizationCollision::ProduceFromFluid");
#ifdef WARPX_DIM_RZ
    if (m_constant_density && m_background_density == 0.0) {
        return;
    }
    auto& warpx = WarpX::GetInstance();
    if (!m_fluid_projectile->getRigidBeam()->active(warpx.Geom(0), start_time, dt)) {
        return;
    }
    auto& electron = mypc->GetParticleContainerFromName(m_product_species[0]);
    IonProductDestination const destination(m_product_species[1], *mypc);
    auto* ion = destination.particles();
    electron.defineAllParticleTiles();
    if (ion) {
        ion->defineAllParticleTiles();
    }
    SmartCreateFactory const electron_factory(electron);
    auto const ion_factory = ion ? std::make_unique<SmartCreateFactory>(*ion) : nullptr;
    auto const create_electron = electron_factory.getSmartCreate();
    auto const create_ion = ion_factory ? ion_factory->getSmartCreate() : create_electron;
    auto const profile = m_fluid_projectile->getRigidBeam()->executor();
    auto const pjg = m_pjg_model->executor();
    auto const sampling = m_sampling_state;
    auto const maximum_transfer =
        pjg.maximumEnergyTransfer(m_fluid_projectile->getRigidBeam()->kineticEnergyEV());
    auto const& geom = warpx.Geom(0);
    auto const domain = geom.Domain();
    int const nr = domain.length(0), nz = domain.length(1);
    int const ir0 = domain.smallEnd(0), iz0 = domain.smallEnd(1);
    double const dr = geom.CellSize(0), dz = geom.CellSize(1);
    double const rmin = geom.ProbLo(0), zmin = geom.ProbLo(1);
    int const resolution = m_source_sampling_points;
    int const quadrature = m_gas_quadrature_points;
    int const bins = quadrature * quadrature * quadrature;
    auto const axial_stride = static_cast<std::size_t>(resolution) + 1;
    auto const gas_stride = static_cast<std::size_t>(bins) + 1;
    auto const constant_gas = m_constant_density;
    auto const gas = m_background_density;
    auto const density_function = m_background_density_func;
    auto const temperature_function = m_background_temperature_func;
    auto const constant_temperature = m_constant_temperature;
    auto const temperature = m_background_temperature;
    auto const ion_mass = destination.getMass();
    auto const neutral_mass = ion_mass + PhysConst::m_e;
    auto const fixed_weight = m_fixed_product_weight;
    auto const cap = m_max_products_per_cell;
    auto const rate = m_source_rate;
    auto const seed = m_sampling_seed;
    bool const fluid_ion = destination.isFluid();
    // Constant inputs were checked at construction. Dynamic parser values
    // must also be checked when device assertions are disabled in Release.
    std::unique_ptr<amrex::Gpu::DeviceScalar<int>> runtime_error;
    if (!constant_gas || !constant_temperature) {
        runtime_error = std::make_unique<amrex::Gpu::DeviceScalar<int>>(0);
    }
    auto* error_pointer = runtime_error ? runtime_error->dataPtr() : nullptr;
    amrex::ignore_unused(error_pointer);

    // The separable air source needs only one axial table per global z cell.
    // Its last entry is the independently integrated physical yield; intermediate
    // entries provide a quiet piecewise-uniform spatial sampler, refined by the
    // source_sampling_points input without changing the physical cell yield.
    m_source_axial.resize(static_cast<std::size_t>(nz) * axial_stride);
    auto* axial = m_source_axial.data();
    amrex::ParallelFor(nz, [=] AMREX_GPU_DEVICE(int j) noexcept {
        double const lo = zmin + j * dz;
        for (int sub = 0; sub <= resolution; ++sub) {
            axial[j * axial_stride + sub] =
                profile.integratedLongitudinal(lo, lo + dz * sub / resolution, start_time, dt);
        }
    });

    auto& remainder = *warpx.m_fields.get(m_remainder_field_name, 0);
    auto& budget = *warpx.m_fields.get(m_budget_field_name, 0);
    auto& counter = *warpx.m_fields.get(m_counter_field_name, 0);
    amrex::MFItInfo info;
    if (amrex::Gpu::notInLaunchRegion()) {
        info.EnableTiling(WarpXParticleContainer::tile_size);
    }
    for (amrex::MFIter mfi = electron.MakeMFIter(0, info); mfi.isValid(); ++mfi) {
        auto const box = mfi.tilebox(amrex::IntVect::TheZeroVector());
        auto const lower = box.smallEnd();
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            box.numPts() <= std::numeric_limits<int>::max(),
            "Rigid-source tiles must contain fewer than INT_MAX cells.");
        int const nx = box.length(0), cells = static_cast<int>(box.numPts());
        m_source_indices.resize(2 * static_cast<std::size_t>(cells));
        m_source_weights.resize(cells);
        if (!constant_gas) {
            m_source_cdf.resize(cells * gas_stride);
        }
        auto* counts = m_source_indices.data();
        auto* offsets = counts + cells;
        auto* weights = m_source_weights.data();
        auto* cdf = m_source_cdf.data();
        auto const rest = remainder.array(mfi);
        auto const totals = budget.array(mfi);
        auto const sequences = counter.array(mfi);
        amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE(int cell) noexcept {
            amrex::ignore_unused(error_pointer);
            int const i = lower[0] + cell % nx, j = lower[1] + cell / nx;
            double const rlo = rmin + (i - ir0) * dr, zlo = zmin + (j - iz0) * dz;
            double source = 0.0;
            if (constant_gas) {
                source = gas * profile.radialIntegral(rlo, rlo + dr) *
                         axial[(j - iz0) * axial_stride + resolution];
            } else {
                // Positive composite quadrature with the exact beam measure in
                // each subcell/time interval. Refinement resolves a varying gas
                // without negative source weights or rejection sampling.
                cdf[cell * gas_stride] = 0.0;
                for (int bin = 0; bin < bins; ++bin) {
                    int const a = bin % quadrature;
                    int const b = (bin / quadrature) % quadrature;
                    int const t = bin / (quadrature * quadrature);
                    double const rb = rlo + a * dr / quadrature, zb = zlo + b * dz / quadrature;
                    double const tb = start_time + t * dt / quadrature;
                    double const background =
                        density_function(rb + 0.5 * dr / quadrature, 0.0,
                                         zb + 0.5 * dz / quadrature, tb + 0.5 * dt / quadrature);
                    if (!(background >= 0.0 && std::isfinite(background))) {
                        AMREX_IF_ON_DEVICE(
                            (amrex::Gpu::Atomic::Max(error_pointer, 1);))
                        AMREX_IF_ON_HOST(
                            (amrex::Abort("Rigid-source neutral density must "
                                          "be finite and non-negative.");))
                        counts[cell] = 0;
                        return;
                    }
                    source += background * profile.radialIntegral(rb, rb + dr / quadrature) *
                              profile.integratedLongitudinal(zb, zb + dz / quadrature, tb,
                                                             dt / quadrature);
                    cdf[cell * gas_stride + bin + 1] = source;
                }
            }
            double const available = rest(i, j, 0) + rate * source;
            AMREX_ALWAYS_ASSERT(available >= 0.0 && std::isfinite(available));
            counts[cell] = 0;
            if (source <= 0.0) {
                return;
            }
            double const expected = available / fixed_weight;
            auto const count =
                static_cast<amrex::Long>(std::min(double(cap), std::floor(expected)));
            auto weight =
                expected >= cap ? static_cast<amrex::ParticleReal>(available / cap) : fixed_weight;
            // Rounding product weights upwards must not make the remainder negative.
            if (count * double(weight) > available) {
                weight = std::nextafter(weight, amrex::ParticleReal{0});
            }
            counts[cell] = count;
            weights[cell] = weight;
            rest(i, j, 0) = static_cast<amrex::Real>(available - count * double(weight));
        });
        auto const total = amrex::Scan::ExclusiveSum(cells, counts, offsets);
        if (!constant_gas) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                runtime_error->dataValue() == 0,
                "Rigid-source neutral density must be finite and "
                "non-negative.");
        }
        if (total == 0) {
            continue;
        }
        auto& electrons = electron.ParticlesAt(0, mfi);
        auto* ions = ion ? &ion->ParticlesAt(0, mfi) : nullptr;
        auto const first_e = electrons.numParticles(), first_i = ions ? ions->numParticles() : 0;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(first_e + total <= std::numeric_limits<int>::max() &&
                                             first_i + total <= std::numeric_limits<int>::max(),
                                         "Rigid-source product tile overflow.");
        electrons.resize(first_e + total);
        if (ions) {
            ions->resize(first_i + total);
        }
        auto e = electrons.getParticleTileData();
        auto ion_data = ions ? ions->getParticleTileData() : decltype(e){};
        amrex::ParallelForRNG(cells, [=] AMREX_GPU_DEVICE(
                                         int cell, amrex::RandomEngine const&
                                                       engine) noexcept {
            amrex::ignore_unused(error_pointer);
            auto electron_data = e;
            auto ions_data = ion_data;
            auto const count = counts[cell];
            if (count == 0) {
                return;
            }
            int const i = lower[0] + cell % nx, j = lower[1] + cell / nx;
            double const rlo = rmin + (i - ir0) * dr, zlo = zmin + (j - iz0) * dz;
            std::uint64_t sequence = 0;
            // Four exact 16-bit digits survive single-precision field checkpoints.
            for (int digit = 0; digit < 4; ++digit) {
                sequence |= static_cast<std::uint64_t>(sequences(i, j, 0, digit)) << (16 * digit);
            }
            auto const phase = scramble(seed + static_cast<std::uint64_t>(j - iz0) * nr + i - ir0);
            double const energy_shift = double(phase >> 11) * 0x1p-53;
            auto const space_shift = static_cast<std::uint32_t>(scramble(phase));
            auto const angle_shift = static_cast<std::uint32_t>(scramble(phase + 1));
            auto const azimuth_shift = static_cast<std::uint32_t>(scramble(phase + 2));
            auto const radial_shift = static_cast<std::uint32_t>(scramble(phase + 3));
            auto const position_shift = static_cast<std::uint32_t>(scramble(phase + 4));
            // Rotate each batch's axial strata through a quiet sequence. A
            // fixed shift repeats the same axial points whenever the cap is
            // reached, leaving persistent density noise across pulse steps.
            double const axial_shift = ProtonImpactIonization::shiftedKronecker<double>(
                static_cast<std::uint32_t>(sequence), space_shift, 0x1f83d9abu);
            double emitted = 0.0, energy = 0.0, binding = 0.0, discarded = 0.0;
            for (amrex::Long product = 0; product < count; ++product) {
                auto const sample_index = static_cast<std::uint32_t>(sequence + product);
                using ProtonImpactIonization::shiftedKronecker;
                double const u_r =
                    shiftedKronecker<double>(sample_index, radial_shift, 0x3c6ef373u);
                double const u_z = (double(product) + axial_shift) / count;
                double r0 = rlo, r1 = std::min(rlo + dr, profile.m_cutoff_r * profile.m_sigma_r);
                double z = 0.0;
                if (constant_gas) {
                    auto const* row = axial + (j - iz0) * axial_stride;
                    double const target =
                        std::min(u_z * row[resolution], std::nextafter(row[resolution], 0.0));
                    int const bin = interval(row, resolution, target);
                    z = zlo +
                        dz / resolution * (bin + (target - row[bin]) / (row[bin + 1] - row[bin]));
                } else {
                    auto const* row = cdf + cell * gas_stride;
                    double const target = std::min(u_z * row[bins], std::nextafter(row[bins], 0.0));
                    int const bin = interval(row, bins, target);
                    int const a = bin % quadrature, b = (bin / quadrature) % quadrature;
                    r0 += a * dr / quadrature;
                    r1 = std::min(r0 + dr / quadrature, profile.m_cutoff_r * profile.m_sigma_r);
                    z = zlo +
                        dz / quadrature * (b + (target - row[bin]) / (row[bin + 1] - row[bin]));
                }
                double const r = radialQuantile(r0, r1, profile.m_sigma_r, u_r);
                auto const stored_r =
                    std::min(static_cast<amrex::ParticleReal>(r),
                             std::nextafter(static_cast<amrex::ParticleReal>(rlo + dr),
                                            static_cast<amrex::ParticleReal>(rlo)));
                auto const stored_z =
                    std::min(static_cast<amrex::ParticleReal>(z),
                             std::nextafter(static_cast<amrex::ParticleReal>(zlo + dz),
                                            static_cast<amrex::ParticleReal>(zlo)));
                int const pe = static_cast<int>(first_e + offsets[cell] + product);
                int const pi = static_cast<int>(first_i + offsets[cell] + product);
                create_electron(electron_data, pe, engine, stored_r, 0.0, stored_z);
                if (!fluid_ion) {
                    create_ion(ions_data, pi, engine, stored_r, 0.0, stored_z);
                }
                auto const theta =
                    2 * MathConst::pi *
                    shiftedKronecker<double>(sample_index, position_shift, 0xa54ff53bu);
                e.m_rdata[PIdx::theta][pe] = static_cast<amrex::ParticleReal>(theta);
                if (!fluid_ion) {
                    ion_data.m_rdata[PIdx::theta][pi] = e.m_rdata[PIdx::theta][pe];
                }
                amrex::ParticleReal secondary, threshold;
                pjg.sample(
                    sampling,
                    ProtonImpactIonization::shiftedRadicalInverse(sample_index, energy_shift),
                    secondary, threshold);
                auto const cosine = ProtonImpactIonization::polarCosine(
                    secondary, threshold, maximum_transfer,
                    shiftedKronecker<amrex::ParticleReal>(sample_index, angle_shift, 0x6a09e667u));
                double const phi =
                    2 * MathConst::pi *
                    shiftedKronecker<double>(sample_index, azimuth_shift, 0xbb67ae85u);
                double const speed = BackgroundMCCKinematics::properSpeedFromKineticEnergy(
                    secondary, PhysConst::m_e_v<double>);
                double const transverse =
                    speed * std::sqrt(std::max(0.0, 1 - double(cosine) * cosine));
                e.m_rdata[PIdx::ux][pe] =
                    static_cast<amrex::ParticleReal>(transverse * std::cos(phi));
                e.m_rdata[PIdx::uy][pe] =
                    static_cast<amrex::ParticleReal>(transverse * std::sin(phi));
                e.m_rdata[PIdx::uz][pe] = static_cast<amrex::ParticleReal>(
                    speed * cosine * (profile.m_velocity > 0 ? 1 : -1));
                double const kelvin =
                    constant_temperature ? temperature : temperature_function(r, 0, z, cur_time);
                if (!(kelvin >= 0.0 && std::isfinite(kelvin))) {
                    AMREX_IF_ON_DEVICE(
                        (amrex::Gpu::Atomic::Max(error_pointer, 2);))
                    AMREX_IF_ON_HOST(
                        (amrex::Abort("Rigid-source neutral temperature must "
                                      "be finite and non-negative.");))
                    return;
                }
                double ion_energy = 0.0;
                for (int dir = 0; dir < 3; ++dir) {
                    double const a = shiftedKronecker<double>(
                        sample_index, static_cast<std::uint32_t>(scramble(phase + 5 + 2 * dir)),
                        0x510e527fu);
                    double const b = shiftedKronecker<double>(
                        sample_index, static_cast<std::uint32_t>(scramble(phase + 6 + 2 * dir)),
                        0x7311c281u);
                    double const velocity =
                        kelvin > 0
                            ? std::sqrt(-2 * PhysConst::kb * kelvin / neutral_mass * std::log(a)) *
                                  std::cos(2 * MathConst::pi * b)
                            : 0.0;
                    if (!fluid_ion) {
                        ion_data.m_rdata[PIdx::ux + dir][pi] =
                            static_cast<amrex::ParticleReal>(velocity);
                    }
                    ion_energy += 0.5 * ion_mass * velocity * velocity;
                }
                auto const weight = weights[cell];
                e.m_rdata[PIdx::w][pe] = weight;
                if (!fluid_ion) {
                    ion_data.m_rdata[PIdx::w][pi] = weight;
                }
                emitted += weight;
                energy += weight * double(secondary) * PhysConst::q_e;
                binding += weight * double(threshold) * PhysConst::q_e;
                if (fluid_ion) {
                    discarded += weight * ion_energy;
                }
            }
            sequence += static_cast<std::uint64_t>(count);
            for (int digit = 0; digit < 4; ++digit) {
                sequences(i, j, 0, digit) =
                    static_cast<amrex::Real>((sequence >> (16 * digit)) & 65535u);
            }
            totals(i, j, 0, 0) += static_cast<amrex::Real>(emitted);
            totals(i, j, 0, 1) += static_cast<amrex::Real>(energy);
            totals(i, j, 0, 2) += static_cast<amrex::Real>(binding);
            totals(i, j, 0, 3) += static_cast<amrex::Real>(discarded);
        });
        if (!constant_temperature) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                runtime_error->dataValue() == 0,
                "Rigid-source neutral temperature must be finite and "
                "non-negative.");
        }
        ParticleCreation::DefaultInitializeRuntimeAttributes(
            electrons, electron, static_cast<int>(first_e), static_cast<int>(first_e + total));
        if (ions) {
            ParticleCreation::DefaultInitializeRuntimeAttributes(
                *ions, *ion, static_cast<int>(first_i), static_cast<int>(first_i + total));
        }
        if (fluid_ion) {
            auto const deposit = destination.deposit(0, mfi);
            amrex::For(total, [=] AMREX_GPU_DEVICE(int p) noexcept {
                auto const index = first_e + p;
                deposit(e.m_rdata[PIdx::r][index], e.m_rdata[PIdx::z][index],
                        e.m_rdata[PIdx::w][index]);
            });
        }
        amrex::Gpu::streamSynchronize();
        setNewParticleIDs(electrons, first_e, total);
        if (ions) {
            setNewParticleIDs(*ions, first_i, total);
        }
    }
    destination.commit(0);
#else
    amrex::ignore_unused(cur_time, start_time, dt, mypc);
#endif
}
