/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ProtonImpactIonization.H"

#include "Particles/Collision/BackgroundMCC/BackgroundMCCKinematics.H"
#include "Particles/Collision/ProtonImpactIonization/ProtonImpactIonizationKinematics.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/ParticleUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/profiler/ProfilerWrapper.H>

#include <AMReX_Array4.H>
#include <AMReX_GpuAssert.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Scan.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <limits>
#include <memory>
#include <string>

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::ParticleReal
    fractionalPart (amrex::ParticleReal const value) noexcept
    {
        using std::floor;
        return value - floor(value);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    normalPair (amrex::ParticleReal const first_uniform, amrex::ParticleReal const second_uniform,
                amrex::ParticleReal& first_normal, amrex::ParticleReal& second_normal) noexcept
    {
        using namespace amrex::literals;
        using std::cos, std::log, std::sin, std::sqrt;

        auto const bounded_first =
            amrex::max(first_uniform, std::numeric_limits<amrex::ParticleReal>::epsilon());
        auto const radius = sqrt(-2.0_prt * log(bounded_first));
        auto const angle =
            2.0_prt * static_cast<amrex::ParticleReal>(MathConst::pi) * second_uniform;
        first_normal = radius * cos(angle);
        second_normal = radius * sin(angle);
    }
} // namespace

ProtonImpactIonizationCollision::ProtonImpactIonizationCollision (
    std::string const& collision_name, MultiParticleContainer const* mypc)
    : CollisionBase(collision_name)
{
    using namespace amrex::literals;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_species_names.size() == 1,
        "Proton-impact ionization must have exactly one projectile species.");

    amrex::ParmParse const pp_collision_name(collision_name);
    pp_collision_name.getarr("product_species", m_product_species);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_product_species.size() == 2,
        "Proton-impact ionization requires product_species = electron ion.");

    auto& projectile = mypc->GetParticleContainerFromName(m_species_names[0]);
    auto& electron = mypc->GetParticleContainerFromName(m_product_species[0]);
    auto& ion = mypc->GetParticleContainerFromName(m_product_species[1]);

    auto const projectile_charge_state = projectile.getCharge() / PhysConst::q_e;
    auto const rounded_charge_state = amrex::Math::round(projectile_charge_state);
    auto const charge_tolerance = 100.0_prt * std::numeric_limits<amrex::ParticleReal>::epsilon();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        projectile.getMass() > PhysConst::m_e && rounded_charge_state >= 1.0_prt &&
            amrex::Math::abs(projectile_charge_state - rounded_charge_state) <= charge_tolerance,
        "The proton-impact projectile must be a positively charged proton or "
        "bare ion.");
    m_projectile_charge_squared = rounded_charge_state * rounded_charge_state;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::Math::abs(electron.getMass() / PhysConst::m_e - 1.0_prt) <= charge_tolerance &&
            amrex::Math::abs(electron.getCharge() / PhysConst::q_e + 1.0_prt) <= charge_tolerance,
        "The first proton-impact product species must be an electron.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        ion.getMass() > PhysConst::m_e &&
            amrex::Math::abs(ion.getCharge() / PhysConst::q_e - 1.0_prt) <= charge_tolerance,
        "The second proton-impact product species must be a singly charged "
        "molecular ion.");

    std::string target_name;
    pp_collision_name.get("ionization_target", target_name);
    m_target = ProtonImpactIonization::PJGModel::parseTarget(target_name);

    auto const neutral_mass =
        (m_target == ProtonImpactIonization::PJGTarget::N2 ? 28.0134_prt : 31.9988_prt) *
        PhysConst::m_u;
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::Math::abs((ion.getMass() + PhysConst::m_e) / neutral_mass - 1.0_prt) < 0.02_prt,
        "The product-ion mass is inconsistent with the selected N2 or O2 "
        "target.");

    utils::parser::getWithParser(pp_collision_name, "fixed_product_weight", m_fixed_product_weight);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(static_cast<double>(m_fixed_product_weight)) &&
            m_fixed_product_weight > 0.0_prt,
        "Proton-impact fixed_product_weight must be finite and positive.");

    utils::parser::queryWithParser(pp_collision_name, "max_products_per_cell",
                                   m_max_products_per_cell);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_max_products_per_cell >= 1,
                                     "Proton-impact max_products_per_cell must be at least one.");

    amrex::ParticleReal background_density;
    if (utils::parser::queryWithParser(pp_collision_name, "background_density",
                                       background_density)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(static_cast<double>(background_density)) &&
                                             background_density >= 0.0_prt,
                                         "Proton-impact background_density must be finite and "
                                         "non-negative.");
        m_background_density_parser =
            utils::parser::makeParser(std::to_string(background_density), {"x", "y", "z", "t"});
    } else {
        std::string background_density_string;
        utils::parser::Store_parserString(pp_collision_name, "background_density(x,y,z,t)",
                                          background_density_string);
        m_background_density_parser =
            utils::parser::makeParser(background_density_string, {"x", "y", "z", "t"});
    }

    amrex::ParticleReal background_temperature;
    if (utils::parser::queryWithParser(pp_collision_name, "background_temperature",
                                       background_temperature)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(background_temperature)) &&
                background_temperature >= 0.0_prt,
            "Proton-impact background_temperature must be finite and "
            "non-negative.");
        m_background_temperature_parser =
            utils::parser::makeParser(std::to_string(background_temperature), {"x", "y", "z", "t"});
    } else {
        std::string background_temperature_string;
        utils::parser::Store_parserString(pp_collision_name, "background_temperature(x,y,z,t)",
                                          background_temperature_string);
        m_background_temperature_parser =
            utils::parser::makeParser(background_temperature_string, {"x", "y", "z", "t"});
    }
    m_background_density_func = m_background_density_parser.compile<4>();
    m_background_temperature_func = m_background_temperature_parser.compile<4>();

    amrex::ParticleReal projectile_energy_min = 1.0e3_prt;
    amrex::ParticleReal projectile_energy_max = 1.0e9_prt;
    utils::parser::queryWithParser(pp_collision_name, "projectile_energy_min",
                                   projectile_energy_min);
    utils::parser::queryWithParser(pp_collision_name, "projectile_energy_max",
                                   projectile_energy_max);
    auto const projectile_rest_energy = projectile.getMass() * PhysConst::c2 / PhysConst::q_e;
    m_pjg_model = std::make_unique<ProtonImpactIonization::PJGModel>(
        m_target, projectile_rest_energy, projectile_energy_min, projectile_energy_max);

    m_remainder_field_name = collision_name + "_product_weight_remainder";
}

void
ProtonImpactIonizationCollision::AllocData ()
{
    auto& warpx = WarpX::GetInstance();
    for (int lev = 0; lev <= warpx.finestLevel(); ++lev) {
        auto const& box_array = warpx.boxArray(lev);
        auto const& distribution_mapping = warpx.DistributionMap(lev);
        warpx.m_fields.alloc_init(m_remainder_field_name, lev, box_array, distribution_mapping, 1,
                                  amrex::IntVect::TheZeroVector(), amrex::Real{0.0},
                                  /*remake=*/true, /*redistribute_on_remake=*/true,
                                  /*checkpoint_restart=*/true);
    }
}

void
ProtonImpactIonizationCollision::doCollisions (amrex::Real const cur_time, amrex::Real const dt,
                                               MultiParticleContainer* mypc)
{
    ABLASTR_PROFILE("ProtonImpactIonizationCollision::doCollisions()");

    using namespace amrex::literals;
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = ParticleTileType::ParticleTileDataType;
    using ParticleBins = amrex::DenseBins<ParticleTileDataType>;
    using index_type = ParticleBins::index_type;
    using SoaDataType = ParticleTileType::ParticleTileDataType;

    auto& projectile = mypc->GetParticleContainerFromName(m_species_names[0]);
    auto& electron = mypc->GetParticleContainerFromName(m_product_species[0]);
    auto& ion = mypc->GetParticleContainerFromName(m_product_species[1]);
    electron.defineAllParticleTiles();
    ion.defineAllParticleTiles();

    SmartCopyFactory const electron_copy_factory(projectile, electron);
    SmartCopyFactory const ion_copy_factory(projectile, ion);
    SmartCopy const copy_electron = electron_copy_factory.getSmartCopy();
    SmartCopy const copy_ion = ion_copy_factory.getSmartCopy();

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceScalar<SmartCopy> device_copy_electron(copy_electron);
    amrex::Gpu::DeviceScalar<SmartCopy> device_copy_ion(copy_ion);
    SmartCopy const* AMREX_RESTRICT copy_electron_pointer = device_copy_electron.dataPtr();
    SmartCopy const* AMREX_RESTRICT copy_ion_pointer = device_copy_ion.dataPtr();
#else
    SmartCopy const* AMREX_RESTRICT copy_electron_pointer = &copy_electron;
    SmartCopy const* AMREX_RESTRICT copy_ion_pointer = &copy_ion;
#endif

    auto const density_function = m_background_density_func;
    auto const temperature_function = m_background_temperature_func;
    auto const pjg = m_pjg_model->executor();
    auto const projectile_mass = projectile.getMass();
    auto const ion_mass = ion.getMass();
    auto const charge_squared = m_projectile_charge_squared;
    auto const fixed_product_weight = m_fixed_product_weight;
    auto const max_products_per_cell = m_max_products_per_cell;

    auto const finest_level = projectile.finestLevel();
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::MFItInfo info;
        if (amrex::Gpu::notInLaunchRegion()) {
            info.EnableTiling(WarpXParticleContainer::tile_size);
        }

#ifdef AMREX_USE_OMP
        info.SetDynamic(true);
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi = projectile.MakeMFIter(lev, info); mfi.isValid(); ++mfi) {
            auto& projectile_tile = projectile.ParticlesAt(lev, mfi);
            auto const projectile_data = projectile_tile.getParticleTileData();
            auto const& geometry = WarpX::GetInstance().Geom(lev);
            auto const bins =
                ParticleUtils::findParticlesInEachCell(geometry, mfi, projectile_tile);

            auto const num_cells = static_cast<int>(bins.numBins());
            auto const* AMREX_RESTRICT cell_offsets = bins.offsetsPtr();
            auto const* AMREX_RESTRICT particle_indices = bins.permutationPtr();
            auto const* AMREX_RESTRICT projectile_idcpu = projectile_data.m_idcpu;
            auto const* AMREX_RESTRICT projectile_weight = projectile_data.m_rdata[PIdx::w];
            auto const* AMREX_RESTRICT projectile_ux = projectile_data.m_rdata[PIdx::ux];
            auto const* AMREX_RESTRICT projectile_uy = projectile_data.m_rdata[PIdx::uy];
            auto const* AMREX_RESTRICT projectile_uz = projectile_data.m_rdata[PIdx::uz];

            // Keep per-cell scratch in two allocations. This path runs once per
            // tile and collision call, so allocation count matters on GPUs.
            amrex::Gpu::DeviceVector<index_type> cell_indices(2 * num_cells, 0);
            amrex::Gpu::DeviceVector<amrex::ParticleReal> cell_reals(3 * num_cells, 0.0_prt);
            auto* AMREX_RESTRICT count_pointer = cell_indices.dataPtr();
            auto* AMREX_RESTRICT offset_pointer = count_pointer + num_cells;
            auto* AMREX_RESTRICT product_weight_pointer = cell_reals.dataPtr();
            auto* AMREX_RESTRICT collision_score_pointer = product_weight_pointer + num_cells;
            auto* AMREX_RESTRICT temperature_pointer = collision_score_pointer + num_cells;

            auto* remainder_multifab =
                WarpX::GetInstance().m_fields.get(m_remainder_field_name, lev);
            auto const remainder_array = remainder_multifab->array(mfi);
            auto const box = mfi.tilebox(amrex::IntVect::TheZeroVector());
            auto const lower = box.smallEnd();
            auto const xyz_min = WarpX::LowerCorner(box, lev, 0.0_rt);
            auto const cell_size = geometry.CellSizeArray();
#if AMREX_SPACEDIM > 1
            auto const length = box.length();
#endif

            amrex::ParallelFor(num_cells, [=] AMREX_GPU_DEVICE(int const cell) noexcept {
                amrex::IntVect grid_index = lower;
                amrex::XDim3 position = {0.0_rt, 0.0_rt, 0.0_rt};
                constexpr auto half = 0.5_rt;
#if defined(WARPX_DIM_1D_Z)
                grid_index[0] += cell;
                position.z = xyz_min.z + (cell + half) * cell_size[0];
#elif defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
                    grid_index[0] += cell;
                    position.x = xyz_min.x + (cell + half) * cell_size[0];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    auto const ix = cell % length[0];
                    auto const iz = cell / length[0];
                    grid_index[0] += ix;
                    grid_index[1] += iz;
                    position.x = xyz_min.x + (ix + half) * cell_size[0];
                    position.z = xyz_min.z + (iz + half) * cell_size[1];
#elif defined(WARPX_DIM_3D)
                    auto const ix = cell % length[0];
                    auto const iy = (cell / length[0]) % length[1];
                    auto const iz = cell / (length[0] * length[1]);
                    grid_index[0] += ix;
                    grid_index[1] += iy;
                    grid_index[2] += iz;
                    position.x = xyz_min.x + (ix + half) * cell_size[0];
                    position.y = xyz_min.y + (iy + half) * cell_size[1];
                    position.z = xyz_min.z + (iz + half) * cell_size[2];
#endif

                auto const density = density_function(position.x, position.y, position.z, cur_time);
                auto const temperature =
                    temperature_function(position.x, position.y, position.z, cur_time);
                AMREX_IF_ON_DEVICE((AMREX_DEVICE_ASSERT(density >= 0.0_prt);
                                    AMREX_DEVICE_ASSERT(temperature >= 0.0_prt);))
                AMREX_IF_ON_HOST((if (density < 0.0_prt || temperature < 0.0_prt) {
                    amrex::Abort("Proton-impact ionization requires "
                                 "non-negative neutral "
                                 "density and temperature.");
                }))
                temperature_pointer[cell] = temperature;

                amrex::ParticleReal score = 0.0_prt;
                for (index_type permutation_index = cell_offsets[cell];
                     permutation_index < cell_offsets[cell + 1]; ++permutation_index) {
                    auto const particle = particle_indices[permutation_index];
                    if (projectile_idcpu[particle] == amrex::ParticleIdCpus::Invalid) {
                        continue;
                    }
                    auto const proper_speed_squared =
                        projectile_ux[particle] * projectile_ux[particle] +
                        projectile_uy[particle] * projectile_uy[particle] +
                        projectile_uz[particle] * projectile_uz[particle];
                    auto const gamma =
                        std::sqrt(1.0_prt + proper_speed_squared * PhysConst::inv_c2);
                    auto const kinetic_energy = projectile_mass * proper_speed_squared /
                                                ((gamma + 1.0_prt) * PhysConst::q_e);
                    auto const speed = std::sqrt(proper_speed_squared) / gamma;
                    score += projectile_weight[particle] * pjg.crossSection(kinetic_energy) * speed;
                }
                collision_score_pointer[cell] = score;

                auto const accumulated_weight =
                    remainder_array(grid_index) + density * charge_squared * score * dt;
                if (accumulated_weight <= 0.0_prt) {
                    remainder_array(grid_index) = 0.0_rt;
                    return;
                }

                auto const expected_products = accumulated_weight / fixed_product_weight;
                if (expected_products >= static_cast<amrex::ParticleReal>(max_products_per_cell)) {
                    count_pointer[cell] = max_products_per_cell;
                    product_weight_pointer[cell] =
                        accumulated_weight /
                        static_cast<amrex::ParticleReal>(max_products_per_cell);
                    remainder_array(grid_index) = 0.0_rt;
                } else {
                    auto const count = static_cast<index_type>(std::floor(expected_products));
                    count_pointer[cell] = count;
                    product_weight_pointer[cell] = fixed_product_weight;
                    remainder_array(grid_index) =
                        accumulated_weight -
                        static_cast<amrex::ParticleReal>(count) * fixed_product_weight;
                }
            });

            auto const total_new =
                amrex::Scan::ExclusiveSum(num_cells, count_pointer, offset_pointer);

            auto& electron_tile = electron.ParticlesAt(lev, mfi);
            auto& ion_tile = ion.ParticlesAt(lev, mfi);
            auto const old_electron_count = electron_tile.numParticles();
            auto const old_ion_count = ion_tile.numParticles();
            electron_tile.resize(old_electron_count + total_new);
            ion_tile.resize(old_ion_count + total_new);
            SoaDataType const electron_data = electron_tile.getParticleTileData();
            SoaDataType const ion_data = ion_tile.getParticleTileData();

#ifdef AMREX_USE_GPU
            amrex::Gpu::DeviceScalar<SoaDataType> device_electron_data(electron_data);
            amrex::Gpu::DeviceScalar<SoaDataType> device_ion_data(ion_data);
            SoaDataType const* AMREX_RESTRICT electron_data_pointer =
                device_electron_data.dataPtr();
            SoaDataType const* AMREX_RESTRICT ion_data_pointer = device_ion_data.dataPtr();
#else
            SoaDataType const* AMREX_RESTRICT electron_data_pointer = &electron_data;
            SoaDataType const* AMREX_RESTRICT ion_data_pointer = &ion_data;
#endif

            auto* AMREX_RESTRICT electron_weight = electron_data.m_rdata[PIdx::w];
            auto* AMREX_RESTRICT electron_ux = electron_data.m_rdata[PIdx::ux];
            auto* AMREX_RESTRICT electron_uy = electron_data.m_rdata[PIdx::uy];
            auto* AMREX_RESTRICT electron_uz = electron_data.m_rdata[PIdx::uz];
            auto* AMREX_RESTRICT ion_weight = ion_data.m_rdata[PIdx::w];
            auto* AMREX_RESTRICT ion_ux = ion_data.m_rdata[PIdx::ux];
            auto* AMREX_RESTRICT ion_uy = ion_data.m_rdata[PIdx::uy];
            auto* AMREX_RESTRICT ion_uz = ion_data.m_rdata[PIdx::uz];

            amrex::ParallelForRNG(num_cells, [=] AMREX_GPU_DEVICE(
                                                 int const cell,
                                                 amrex::RandomEngine const& engine) noexcept {
                auto const product_count = count_pointer[cell];
                auto const total_score = collision_score_pointer[cell];
                if (product_count == 0 || total_score <= 0.0_prt) {
                    return;
                }

                auto const electron_copy = *copy_electron_pointer;
                auto const ion_copy = *copy_ion_pointer;
                auto const first_particle = cell_offsets[cell];
                auto const last_particle = cell_offsets[cell + 1];
                auto const score_spacing =
                    total_score / static_cast<amrex::ParticleReal>(product_count);
                auto score_target = amrex::Random(engine) * score_spacing;
                auto cumulative_score = 0.0_prt;
                index_type permutation_index = first_particle;
                index_type selected_particle = -1;

                auto const energy_shift = amrex::Random(engine);
                auto const angle_shift = amrex::Random(engine);
                auto const azimuth_shift = amrex::Random(engine);
                auto const normal_shift_1 = amrex::Random(engine);
                auto const normal_shift_2 = amrex::Random(engine);
                auto const normal_shift_3 = amrex::Random(engine);
                auto const normal_shift_4 = amrex::Random(engine);

                for (index_type product = 0; product < product_count; ++product) {
                    while (permutation_index < last_particle && cumulative_score <= score_target) {
                        auto const candidate = particle_indices[permutation_index++];
                        if (projectile_idcpu[candidate] == amrex::ParticleIdCpus::Invalid) {
                            continue;
                        }
                        auto const proper_speed_squared =
                            projectile_ux[candidate] * projectile_ux[candidate] +
                            projectile_uy[candidate] * projectile_uy[candidate] +
                            projectile_uz[candidate] * projectile_uz[candidate];
                        auto const gamma =
                            std::sqrt(1.0_prt + proper_speed_squared * PhysConst::inv_c2);
                        auto const kinetic_energy = projectile_mass * proper_speed_squared /
                                                    ((gamma + 1.0_prt) * PhysConst::q_e);
                        auto const speed = std::sqrt(proper_speed_squared) / gamma;
                        auto const particle_score =
                            projectile_weight[candidate] * pjg.crossSection(kinetic_energy) * speed;
                        if (particle_score > 0.0_prt) {
                            selected_particle = candidate;
                            cumulative_score += particle_score;
                        }
                    }
                    AMREX_IF_ON_DEVICE((AMREX_DEVICE_ASSERT(selected_particle >= 0);))
                    AMREX_IF_ON_HOST((if (selected_particle < 0) {
                        amrex::Abort("Proton-impact ionization failed to "
                                     "select a projectile.");
                    }))

                    auto const output_offset = offset_pointer[cell] + product;
                    auto const electron_index = old_electron_count + output_offset;
                    auto const ion_index = old_ion_count + output_offset;
                    electron_copy(*electron_data_pointer, projectile_data, selected_particle,
                                  static_cast<int>(electron_index), engine);
                    ion_copy(*ion_data_pointer, projectile_data, selected_particle,
                             static_cast<int>(ion_index), engine);

                    auto const proper_speed_squared =
                        projectile_ux[selected_particle] * projectile_ux[selected_particle] +
                        projectile_uy[selected_particle] * projectile_uy[selected_particle] +
                        projectile_uz[selected_particle] * projectile_uz[selected_particle];
                    auto const gamma =
                        std::sqrt(1.0_prt + proper_speed_squared * PhysConst::inv_c2);
                    auto const kinetic_energy = projectile_mass * proper_speed_squared /
                                                ((gamma + 1.0_prt) * PhysConst::q_e);

                    // Stratify the SDCS exactly. Distinct irrational rotations
                    // decorrelate angle and thermal sequences without
                    // per-product RNG calls.
                    auto const sequence_index = static_cast<amrex::ParticleReal>(product) + 0.5_prt;
                    auto const energy_quantile =
                        (static_cast<amrex::ParticleReal>(product) + energy_shift) /
                        static_cast<amrex::ParticleReal>(product_count);
                    amrex::ParticleReal secondary_energy;
                    amrex::ParticleReal binding_energy;
                    pjg.sample(kinetic_energy, energy_quantile, secondary_energy, binding_energy);

                    BackgroundMCCKinematics::Vector3 incident_direction{0.0, 0.0, 1.0};
                    auto const proper_speed = std::sqrt(proper_speed_squared);
                    if (proper_speed > 0.0_prt) {
                        incident_direction = {
                            static_cast<double>(projectile_ux[selected_particle] / proper_speed),
                            static_cast<double>(projectile_uy[selected_particle] / proper_speed),
                            static_cast<double>(projectile_uz[selected_particle] / proper_speed)};
                    }
                    BackgroundMCCKinematics::Vector3 transverse_1;
                    BackgroundMCCKinematics::Vector3 transverse_2;
                    BackgroundMCCKinematics::transverseDirections(incident_direction, transverse_1,
                                                                  transverse_2);

                    auto const maximum_transfer = pjg.maximumEnergyTransfer(kinetic_energy);
                    auto const cosine = ProtonImpactIonization::polarCosine(
                        secondary_energy, binding_energy, maximum_transfer,
                        fractionalPart(angle_shift + 0.4142135623730950488_prt * sequence_index));
                    auto const azimuth =
                        2.0_prt * static_cast<amrex::ParticleReal>(MathConst::pi) *
                        fractionalPart(azimuth_shift + 0.7320508075688772935_prt * sequence_index);
                    auto const electron_direction =
                        BackgroundMCCKinematics::directionFromPolarAngle(
                            incident_direction, transverse_1, transverse_2,
                            static_cast<double>(cosine), static_cast<double>(azimuth));
                    auto const secondary_proper_speed =
                        BackgroundMCCKinematics::properSpeedFromKineticEnergy(
                            static_cast<double>(secondary_energy), PhysConst::m_e_v<double>);
                    electron_ux[electron_index] = static_cast<amrex::ParticleReal>(
                        secondary_proper_speed * electron_direction.x);
                    electron_uy[electron_index] = static_cast<amrex::ParticleReal>(
                        secondary_proper_speed * electron_direction.y);
                    electron_uz[electron_index] = static_cast<amrex::ParticleReal>(
                        secondary_proper_speed * electron_direction.z);

                    auto const thermal_speed =
                        std::sqrt(PhysConst::kb * temperature_pointer[cell] / ion_mass);
                    amrex::ParticleReal normal_x;
                    amrex::ParticleReal normal_y;
                    amrex::ParticleReal normal_z;
                    amrex::ParticleReal unused_normal;
                    normalPair(
                        amrex::max(fractionalPart(normal_shift_1 +
                                                  0.2360679774997896964_prt * sequence_index),
                                   std::numeric_limits<amrex::ParticleReal>::epsilon()),
                        fractionalPart(normal_shift_2 + 0.6457513110645905905_prt * sequence_index),
                        normal_x, normal_y);
                    normalPair(
                        amrex::max(fractionalPart(normal_shift_3 +
                                                  0.3166247903553998491_prt * sequence_index),
                                   std::numeric_limits<amrex::ParticleReal>::epsilon()),
                        fractionalPart(normal_shift_4 + 0.4494897427831780982_prt * sequence_index),
                        normal_z, unused_normal);
                    ion_ux[ion_index] = thermal_speed * normal_x;
                    ion_uy[ion_index] = thermal_speed * normal_y;
                    ion_uz[ion_index] = thermal_speed * normal_z;

                    auto const weight = product_weight_pointer[cell];
                    electron_weight[electron_index] = weight;
                    ion_weight[ion_index] = weight;
                    score_target += score_spacing;
                }
            });

            if (total_new > 0) {
                ParticleCreation::DefaultInitializeRuntimeAttributes(
                    electron_tile, electron, static_cast<int>(old_electron_count),
                    static_cast<int>(old_electron_count + total_new));
                ParticleCreation::DefaultInitializeRuntimeAttributes(
                    ion_tile, ion, static_cast<int>(old_ion_count),
                    static_cast<int>(old_ion_count + total_new));
            }
            amrex::Gpu::synchronize();
            setNewParticleIDs(electron_tile, old_electron_count, total_new);
            setNewParticleIDs(ion_tile, old_ion_count, total_new);
        }
    }
}
