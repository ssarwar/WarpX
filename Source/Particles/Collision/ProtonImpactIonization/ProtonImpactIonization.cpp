/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "ProtonImpactIonization.H"

#include "Particles/Collision/BackgroundMCC/BackgroundMCCKinematics.H"
#include "Particles/Collision/IonProductDestination.H"
#include "Particles/Collision/ProtonImpactIonization/IonizationSampling.H"
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
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_REAL.H>
#include <AMReX_Random.H>
#include <AMReX_Scan.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE std::uint32_t
    randomPhase (amrex::RandomEngine const& engine) noexcept
    {
        return (amrex::Random_int(1u << 16, engine) << 16) | amrex::Random_int(1u << 16, engine);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    normalPair (amrex::ParticleReal const first_uniform, amrex::ParticleReal const second_uniform,
                amrex::ParticleReal& first_normal, amrex::ParticleReal& second_normal) noexcept
    {
        using namespace amrex::literals;
        using std::log, std::sqrt;

        // shiftedKronecker supplies open-interval uniforms, including in float.
        auto const radius = sqrt(-2.0_prt * log(first_uniform));
        auto const angle =
            2.0_prt * static_cast<amrex::ParticleReal>(MathConst::pi) * second_uniform;
        auto const [sine, cosine] = amrex::Math::sincos(angle);
        first_normal = radius * cosine;
        second_normal = radius * sine;
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
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_species_names[0] != m_product_species[0] && m_species_names[0] != m_product_species[1],
        "Proton-impact product species must differ from the projectile "
        "species.");

    auto& warpx = WarpX::GetInstance();
    if (warpx.DoFluidSpecies()) {
        m_fluid_projectile = warpx.GetFluidContainer().FindSpecies(m_species_names[0]);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_fluid_projectile || m_fluid_projectile->getRigidBeam(),
        "A fluid proton-impact projectile must use model = rigid_beam.");
    auto* projectile = m_fluid_projectile ? nullptr :
        &mypc->GetParticleContainerFromName(m_species_names[0]);
    auto const projectile_mass = m_fluid_projectile ? m_fluid_projectile->getMass() : projectile->getMass();
    auto const projectile_charge = m_fluid_projectile ? m_fluid_projectile->getCharge() : projectile->getCharge();
    auto& electron = mypc->GetParticleContainerFromName(m_product_species[0]);
    IonProductDestination const ion(m_product_species[1], *mypc);

    auto const projectile_charge_state = projectile_charge / PhysConst::q_e;
    auto const rounded_charge_state = amrex::Math::round(projectile_charge_state);
    auto const charge_tolerance = 100.0_prt * std::numeric_limits<amrex::ParticleReal>::epsilon();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        projectile_mass > PhysConst::m_e && rounded_charge_state >= 1.0_prt &&
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
        m_background_density = background_density;
        m_constant_density = true;
    } else {
        std::string background_density_string;
        utils::parser::Store_parserString(pp_collision_name, "background_density(x,y,z,t)",
                                          background_density_string);
        m_background_density_parser =
            utils::parser::makeParser(background_density_string, {"x", "y", "z", "t"});
        m_background_density_func = m_background_density_parser.compile<4>();
    }

    amrex::ParticleReal background_temperature;
    if (utils::parser::queryWithParser(pp_collision_name, "background_temperature",
                                       background_temperature)) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(background_temperature)) &&
                background_temperature >= 0.0_prt,
            "Proton-impact background_temperature must be finite and "
            "non-negative.");
        m_background_temperature = background_temperature;
        m_constant_temperature = true;
    } else {
        std::string background_temperature_string;
        utils::parser::Store_parserString(pp_collision_name, "background_temperature(x,y,z,t)",
                                          background_temperature_string);
        m_background_temperature_parser =
            utils::parser::makeParser(background_temperature_string, {"x", "y", "z", "t"});
        m_background_temperature_func = m_background_temperature_parser.compile<4>();
    }

    auto const projectile_rest_energy = projectile_mass * PhysConst::c2 / PhysConst::q_e;
    auto const projectile_mass_scale = projectile_mass / PhysConst::m_p;
    amrex::ParticleReal projectile_energy_min = 5.0e3_prt * projectile_mass_scale;
    amrex::ParticleReal projectile_energy_max = 1.0e10_prt * projectile_mass_scale;
    utils::parser::queryWithParser(pp_collision_name, "projectile_energy_min",
                                   projectile_energy_min);
    utils::parser::queryWithParser(pp_collision_name, "projectile_energy_max",
                                   projectile_energy_max);
    auto const mono_energy = m_fluid_projectile ?
        m_fluid_projectile->getRigidBeam()->kineticEnergyEV() : -1.0;
    m_pjg_model = std::make_unique<ProtonImpactIonization::PJGModel>(
        m_target, projectile_rest_energy, projectile_energy_min, projectile_energy_max, mono_energy);
    if (m_fluid_projectile) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_ndt == 1 ||
            m_collision_stepping_mode == CollisionSteppingMode::Subcycle,
            "Rigid-beam ionization supports subcycling, but not supercycling above one.");
        m_sampling_state = m_pjg_model->monoenergeticSamplingState();
        m_source_rate = m_projectile_charge_squared*m_pjg_model->monoenergeticCrossSection()*
            std::abs(m_fluid_projectile->getRigidBeam()->velocity());
        utils::parser::queryWithParser(pp_collision_name, "source_sampling_points", m_source_sampling_points);
        utils::parser::queryWithParser(pp_collision_name, "gas_quadrature_points", m_gas_quadrature_points);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_source_sampling_points > 0 &&
            m_source_sampling_points < std::numeric_limits<int>::max() &&
            m_gas_quadrature_points > 0 && m_gas_quadrature_points <= 16,
            "Source sampling points must be positive; gas quadrature points must lie in [1,16].");
        int seed = 0;
        utils::parser::queryWithParser(pp_collision_name, "sampling_seed", seed);
        m_sampling_seed = static_cast<std::uint64_t>(seed);
        // Stable across processes, restarts, standard libraries and backends.
        for (unsigned char character : collision_name) {
            m_sampling_seed = (m_sampling_seed ^ character)*1099511628211ull;
        }
    }

    m_remainder_field_name = collision_name + "_product_weight_remainder";
    m_budget_field_name = collision_name + "_source_budget";
    m_counter_field_name = collision_name + "_sampling_counter";
    if (m_fluid_projectile || ion.isFluid()) {
        std::ostringstream configuration;
        configuration << std::setprecision(std::numeric_limits<double>::max_digits10);
        configuration << collision_name << ' ' << m_species_names[0] << ' '
            << m_product_species[0] << ' ' << m_product_species[1] << ' '
            << projectile_mass << ' ' << projectile_charge << ' '
            << ion.getMass() << ' ' << static_cast<int>(m_target) << ' '
            << m_ndt << ' ' << static_cast<int>(m_collision_stepping_mode) << ' '
            << m_fixed_product_weight << ' ' << m_max_products_per_cell << ' '
            << projectile_energy_min << ' ' << projectile_energy_max << ' '
            << m_sampling_seed << ' ' << m_source_sampling_points << ' '
            << m_gas_quadrature_points << ' ' << m_constant_density << ' '
            << m_background_density << ' ' << m_constant_temperature << ' '
            << m_background_temperature << '\n';
        if (!m_constant_density) { configuration << m_background_density_parser.expr() << '\n'; }
        if (!m_constant_temperature) { configuration << m_background_temperature_parser.expr() << '\n'; }
        if (m_fluid_projectile) { configuration << "quiet_sampler 2\n"; }
        m_configuration = configuration.str();
    }
}

amrex::Vector<std::string>
ProtonImpactIonizationCollision::CheckpointFields () const
{
    if (m_configuration.empty()) { return {}; }
    if (!m_fluid_projectile) { return {m_remainder_field_name}; }
    return {m_remainder_field_name, m_budget_field_name, m_counter_field_name};
}

void
ProtonImpactIonizationCollision::ValidateRestartState () const
{
    if (m_configuration.empty()) { return; }
    auto const& fields = WarpX::GetInstance().m_fields;
    for (auto const& name : CheckpointFields()) {
        auto const& state = *fields.get(name, 0);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!state.contains_nan() && !state.contains_inf(),
            "Non-finite prescribed-source checkpoint state: "+name);
        for (int comp = 0; comp < state.nComp(); ++comp) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(state.min(comp) >= 0.0,
                "Negative prescribed-source checkpoint state: "+name);
        }
    }
    if (!m_fluid_projectile) { return; }
    auto const& counter = *fields.get(m_counter_field_name, 0);
    for (amrex::MFIter mfi(counter, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const values = counter.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), 4,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                auto const value = values(i,j,k,n);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(value <= 65535.0 && value == std::floor(value),
                    "Invalid prescribed-source checkpoint sampling counter.");
            });
    }
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
        if (m_fluid_projectile) {
            for (auto const& name : {m_budget_field_name, m_counter_field_name}) {
                warpx.m_fields.alloc_init(name, lev, box_array, distribution_mapping, 4,
                    amrex::IntVect::TheZeroVector(), amrex::Real{0.0}, true, true, true);
            }
        }
    }
}

void
ProtonImpactIonizationCollision::doCollisionsInInterval (
    amrex::Real cur_time, amrex::Real start_time, amrex::Real dt, MultiParticleContainer* mypc)
{
    if (m_fluid_projectile) {
        ProduceFromFluid(cur_time, start_time, dt, mypc);
    } else {
        doCollisions(cur_time, dt, mypc);
    }
}

void
ProtonImpactIonizationCollision::doCollisions (amrex::Real const cur_time, amrex::Real const dt,
                                               MultiParticleContainer* mypc)
{
    ABLASTR_PROFILE("ProtonImpactIonizationCollision::doCollisions()");

    using namespace amrex::literals;
    if (m_constant_density && m_background_density == 0.0_prt) {
        return;
    }
    using ParticleTileType = WarpXParticleContainer::ParticleTileType;
    using ParticleTileDataType = ParticleTileType::ParticleTileDataType;
    using ParticleBins = amrex::DenseBins<ParticleTileDataType>;
    using index_type = ParticleBins::index_type;
    using SoaDataType = ParticleTileType::ParticleTileDataType;

    auto& projectile = mypc->GetParticleContainerFromName(m_species_names[0]);
    auto& electron = mypc->GetParticleContainerFromName(m_product_species[0]);
    IonProductDestination const ion_destination(m_product_species[1], *mypc);
    auto* ion = ion_destination.particles();
    bool const fluid_ion = ion_destination.isFluid();
    electron.defineAllParticleTiles();
    if (ion) { ion->defineAllParticleTiles(); }

    SmartCopyFactory const electron_copy_factory(projectile, electron);
    auto const ion_copy_factory = ion ? std::make_unique<SmartCopyFactory>(projectile, *ion) : nullptr;
    SmartCopy const copy_electron = electron_copy_factory.getSmartCopy();
    SmartCopy const copy_ion = ion_copy_factory ? ion_copy_factory->getSmartCopy() : copy_electron;

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
    auto const constant_density = m_constant_density;
    auto const constant_temperature = m_constant_temperature;
    auto const background_density = m_background_density;
    auto const background_temperature = m_background_temperature;
    auto const pjg = m_pjg_model->executor();
    auto const projectile_mass = projectile.getMass();
    // Products inherit the neutral velocity distribution, not an ion Maxwellian
    // at the same temperature. Neglect the binding mass defect in this conversion.
    auto const neutral_mass = ion_destination.getMass() + PhysConst::m_e;
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
            // Empty tiles cannot contribute new weight or select a parent.
            // Leave their checkpointed fractional remainder for the next visit.
            if (projectile_tile.numParticles() == 0) {
                continue;
            }
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
            // Sum counts in 64 bits before checking the tile's int index limit.
            amrex::Gpu::DeviceVector<amrex::Long> cell_indices(
                2 * static_cast<std::size_t>(num_cells) + 1, 0);
            amrex::Gpu::DeviceVector<amrex::ParticleReal> cell_reals(3 * num_cells, 0.0_prt);
            auto* AMREX_RESTRICT count_pointer = cell_indices.dataPtr();
            auto* AMREX_RESTRICT offset_pointer = count_pointer + num_cells;
            auto* AMREX_RESTRICT runtime_error = offset_pointer + num_cells;
            amrex::ignore_unused(runtime_error);
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

            amrex::ParallelFor(num_cells, [=] AMREX_GPU_DEVICE(
                                              int const cell) noexcept {
                // Keep the capture list identical in CUDA host/device
                // compilation.
                amrex::ignore_unused(runtime_error);
                if (cell_offsets[cell] == cell_offsets[cell + 1]) {
                    return;
                }
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

                auto const density = constant_density ? background_density
                                                      : density_function(position.x, position.y,
                                                                         position.z, cur_time);
                auto const temperature =
                    constant_temperature
                        ? background_temperature
                        : temperature_function(position.x, position.y, position.z, cur_time);
                if (!(density >= 0.0_prt && temperature >= 0.0_prt) ||
                    !std::isfinite(density) || !std::isfinite(temperature)) {
                    AMREX_IF_ON_DEVICE((amrex::Gpu::Atomic::Max(
                                            runtime_error, amrex::Long{1});))
                    AMREX_IF_ON_HOST(
                        (amrex::Abort("Proton-impact ionization requires "
                                      "finite, non-negative neutral "
                                      "density and temperature.");))
                    return;
                }
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
#ifdef AMREX_USE_GPU
            if (!constant_density || !constant_temperature) {
                // Constants were validated at construction. Parser values need
                // a host check: device assertions vanish in Release builds.
                amrex::Long error = 0;
                amrex::Gpu::copy(amrex::Gpu::deviceToHost, runtime_error,
                                 runtime_error + 1, &error);
                WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                    error == 0, "Proton-impact ionization requires finite, "
                                "non-negative neutral "
                                "density and temperature.");
            }
#endif
            if (total_new == 0) {
                continue;
            }

            auto& electron_tile = electron.ParticlesAt(lev, mfi);
            auto* ion_tile = ion ? &ion->ParticlesAt(lev, mfi) : nullptr;
            auto const old_electron_count = electron_tile.numParticles();
            auto const old_ion_count = ion_tile ? ion_tile->numParticles() : 0;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                old_electron_count + total_new <= std::numeric_limits<int>::max() &&
                    old_ion_count + total_new <= std::numeric_limits<int>::max(),
                "Proton-impact product tiles must contain fewer than INT_MAX "
                "particles.");
            electron_tile.resize(old_electron_count + total_new);
            if (ion_tile) { ion_tile->resize(old_ion_count + total_new); }
            SoaDataType const electron_data = electron_tile.getParticleTileData();
            SoaDataType const ion_data = ion_tile ? ion_tile->getParticleTileData() : SoaDataType{};

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
                auto const score_shift = amrex::Random(engine);
                auto cumulative_score = 0.0_prt;
                index_type permutation_index = first_particle;
                index_type selected_particle = -1;
                index_type previous_particle = -1;
                amrex::ParticleReal kinetic_energy = 0.0_prt;
                ProtonImpactIonization::NeutralFrame cold_frame;
                auto const thermal_speed =
                    std::sqrt(PhysConst::kb * temperature_pointer[cell] / neutral_mass);
                ProtonImpactIonization::PJGModel::Executor::SamplingState sampling_state;

                // A 53-bit cell shift retains the rare hard-electron tail
                // even in single-precision builds. Only two integer draws
                // per cell are needed, not per emitted electron.
                auto const energy_shift =
                    static_cast<double>(amrex::Random_int(1u << 26, engine)) * 0x1p-26 +
                    static_cast<double>(amrex::Random_int(1u << 27, engine)) * 0x1p-53;
                auto const angle_shift = randomPhase(engine);
                auto const azimuth_shift = randomPhase(engine);
                auto const normal_shift_1 = thermal_speed > 0.0_prt ? randomPhase(engine) : 0u;
                auto const normal_shift_2 = thermal_speed > 0.0_prt ? randomPhase(engine) : 0u;
                auto const normal_shift_3 = thermal_speed > 0.0_prt ? randomPhase(engine) : 0u;
                auto const normal_shift_4 = thermal_speed > 0.0_prt ? randomPhase(engine) : 0u;

                for (index_type product = 0; product < product_count; ++product) {
                    // Avoid the cumulative rounding drift of repeatedly adding
                    // score_spacing, especially with many float products.
                    auto const score_target =
                        (static_cast<double>(product) + score_shift) * score_spacing;
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
                        auto const candidate_energy = projectile_mass * proper_speed_squared /
                                                      ((gamma + 1.0_prt) * PhysConst::q_e);
                        auto const proper_speed = std::sqrt(proper_speed_squared);
                        auto const speed = proper_speed / gamma;
                        auto const particle_score = projectile_weight[candidate] *
                                                    pjg.crossSection(candidate_energy) * speed;
                        if (particle_score > 0.0_prt) {
                            selected_particle = candidate;
                            kinetic_energy = candidate_energy;
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
                    if (!fluid_ion) {
                        ion_copy(*ion_data_pointer, projectile_data, selected_particle,
                                 static_cast<int>(ion_index), engine);
                    }

                    // An independent shift makes every energy quantile uniform
                    // conditional on its selected parent. Ordered quantiles
                    // would correlate energy with ordered parent selection.
                    auto const energy_quantile =
                        ProtonImpactIonization::shiftedRadicalInverse(
                            static_cast<std::uint32_t>(product),
                            energy_shift);
                    // All products selected from one parent share this frame
                    // and sampling state. Rebuild them only when the parent changes.
                    if (selected_particle != previous_particle) {
                        cold_frame = ProtonImpactIonization::neutralFrame(
                            {projectile_ux[selected_particle], projectile_uy[selected_particle],
                             projectile_uz[selected_particle]},
                            {}, projectile_mass);
                        sampling_state = pjg.prepareSampling(kinetic_energy);
                        previous_particle = selected_particle;
                    }
                    using ProtonImpactIonization::shiftedKronecker;
                    amrex::ParticleReal normal_x = 0.0_prt;
                    amrex::ParticleReal normal_y = 0.0_prt;
                    amrex::ParticleReal normal_z = 0.0_prt;
                    if (thermal_speed > 0.0_prt) {
                        amrex::ParticleReal unused_normal;
                        normalPair(
                            shiftedKronecker<amrex::ParticleReal>(product, normal_shift_1,
                                                                 0x3c6ef373u),
                            shiftedKronecker<amrex::ParticleReal>(product, normal_shift_2,
                                                                 0xa54ff53bu),
                            normal_x, normal_y);
                        normalPair(
                            shiftedKronecker<amrex::ParticleReal>(product, normal_shift_3,
                                                                 0x510e527fu),
                            shiftedKronecker<amrex::ParticleReal>(product, normal_shift_4,
                                                                 0x7311c281u),
                            normal_z, unused_normal);
                    }
                    auto frame = cold_frame;
                    auto event_sampling = sampling_state;
                    if (thermal_speed > 0) {
                        frame = ProtonImpactIonization::neutralFrame(
                            {projectile_ux[selected_particle], projectile_uy[selected_particle],
                             projectile_uz[selected_particle]},
                            {thermal_speed * normal_x, thermal_speed * normal_y,
                             thermal_speed * normal_z},
                            projectile_mass);
                        event_sampling =
                            pjg.prepareSampling(static_cast<amrex::ParticleReal>(frame.m_energy));
                    }
                    if (event_sampling.m_energy_index < 0) {
                        AMREX_IF_ON_DEVICE(
                            (amrex::Gpu::Atomic::Max(runtime_error, amrex::Long{2}); return;))
                        AMREX_IF_ON_HOST((amrex::Abort("Neutral-frame proton energy is outside the "
                                                       "calibrated PJG range.");))
                    }
                    amrex::ParticleReal secondary_energy, binding_energy;
                    pjg.sample(event_sampling, energy_quantile, secondary_energy, binding_energy);
                    auto const products = ProtonImpactIonization::compute(
                        frame, secondary_energy, binding_energy, projectile_mass, neutral_mass,
                        shiftedKronecker<double>(product, angle_shift, 0x6a09e667u),
                        2 * static_cast<double>(MathConst::pi) *
                            shiftedKronecker<double>(product, azimuth_shift, 0xbb67ae85u));
                    if (!products.m_valid) {
                        AMREX_IF_ON_DEVICE(
                            (amrex::Gpu::Atomic::Max(runtime_error, amrex::Long{3}); return;))
                        AMREX_IF_ON_HOST(
                            (amrex::Abort(
                                 "Proton-impact spectrum/angle has no on-shell recoil state.");))
                    }
                    electron_ux[electron_index] =
                        static_cast<amrex::ParticleReal>(products.m_electron.x);
                    electron_uy[electron_index] =
                        static_cast<amrex::ParticleReal>(products.m_electron.y);
                    electron_uz[electron_index] =
                        static_cast<amrex::ParticleReal>(products.m_electron.z);
                    if (!fluid_ion) {
                        ion_ux[ion_index] = static_cast<amrex::ParticleReal>(products.m_ion.x);
                        ion_uy[ion_index] = static_cast<amrex::ParticleReal>(products.m_ion.y);
                        ion_uz[ion_index] = static_cast<amrex::ParticleReal>(products.m_ion.z);
                    }

                    auto const weight = product_weight_pointer[cell];
                    electron_weight[electron_index] = weight;
                    if (!fluid_ion) { ion_weight[ion_index] = weight; }
                }
            });

#ifdef AMREX_USE_GPU
            // Device assertions vanish in Release builds. Reuse the tile's
            // error slot and check it before any product deposition can run.
            amrex::Long product_error = 0;
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, runtime_error, runtime_error + 1,
                             &product_error);
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                product_error != 2,
                "Neutral-frame proton energy is outside the calibrated PJG range.");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                product_error != 3, "Proton-impact spectrum/angle has no on-shell recoil state.");
#endif
            if (total_new > 0) {
                ParticleCreation::DefaultInitializeRuntimeAttributes(
                    electron_tile, electron, static_cast<int>(old_electron_count),
                    static_cast<int>(old_electron_count + total_new));
                if (ion_tile) {
                    ParticleCreation::DefaultInitializeRuntimeAttributes(
                        *ion_tile, *ion, static_cast<int>(old_ion_count),
                        static_cast<int>(old_ion_count + total_new));
                }
            }
#ifdef WARPX_DIM_RZ
            if (fluid_ion) {
                auto const deposit = ion_destination.deposit(lev, mfi);
                // Scatter the actual stored electron footprint, after SmartCopy
                // and weight assignment. There is no temporary ion particle.
                amrex::For(total_new, [=] AMREX_GPU_DEVICE(int i) noexcept {
                    auto const p = old_electron_count + i;
                    deposit(electron_data.m_rdata[PIdx::r][p],
                            electron_data.m_rdata[PIdx::z][p], electron_weight[p]);
                });
            }
#endif
            amrex::Gpu::synchronize();
            setNewParticleIDs(electron_tile, old_electron_count, total_new);
            if (ion_tile) { setNewParticleIDs(*ion_tile, old_ion_count, total_new); }
        }
        ion_destination.commit(lev);
    }
}
