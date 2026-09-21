/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "RigidBeam.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

RigidBeam::RigidBeam (std::string const& name, amrex::Real mass, amrex::Real charge)
    : m_charge(charge)
{
    amrex::ParmParse const pp(name);
    auto& p = m_parameters;
    auto const positive = [](double value) { return std::isfinite(value) && value > 0.0; };
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(positive(mass) && positive(std::abs(charge)),
        "A rigid beam requires finite positive mass and finite nonzero charge.");
    bool const energy = utils::parser::queryWithParser(pp, "kinetic_energy", m_energy_ev);
    bool const speed = utils::parser::queryWithParser(pp, "velocity_z", p.m_velocity);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(energy != speed,
        "A rigid beam requires exactly one of kinetic_energy [eV] and velocity_z [m/s].");
    double const rest = mass*PhysConst::c2/PhysConst::q_e;
    if (energy) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(positive(m_energy_ev),
            "Rigid-beam kinetic_energy must be finite and positive.");
        p.m_velocity = PhysConst::c*std::sqrt(m_energy_ev*(m_energy_ev+2*rest))/(m_energy_ev+rest);
    } else {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(positive(std::abs(p.m_velocity)) &&
                                            std::abs(p.m_velocity) < PhysConst::c,
            "Rigid-beam velocity_z must be finite, nonzero and subluminal.");
        auto const beta = p.m_velocity/PhysConst::c;
        auto const gamma = 1.0/std::sqrt(1.0-beta*beta);
        m_energy_ev = rest*gamma*gamma*beta*beta/(gamma+1.0);
    }
    utils::parser::getWithParser(pp, "sigma_r", p.m_sigma_r);
    bool const length = utils::parser::queryWithParser(pp, "sigma_z", p.m_sigma_z);
    double duration = 0.0;
    bool const time_width = utils::parser::queryWithParser(pp, "sigma_t", duration);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(length != time_width,
        "A rigid beam requires exactly one of sigma_z [m] and sigma_t [s].");
    if (time_width) { p.m_sigma_z = std::abs(p.m_velocity)*duration; }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(positive(p.m_sigma_r) && positive(p.m_sigma_z),
        "Rigid-beam sigma_r and longitudinal RMS width must be finite and positive.");
    utils::parser::queryWithParser(pp, "z_reference", p.m_z_reference);
    utils::parser::queryWithParser(pp, "cutoff_r", p.m_cutoff_r);
    utils::parser::queryWithParser(pp, "cutoff_z", p.m_cutoff_z);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(p.m_z_reference) &&
        (positive(p.m_cutoff_r) || p.m_cutoff_r == -1.0) &&
        (positive(p.m_cutoff_z) || p.m_cutoff_z == -1.0),
        "Rigid-beam cutoffs must be finite positive sigma multiples, or -1 for no truncation.");
    if (p.m_cutoff_r == -1.0) { p.m_cutoff_r = std::numeric_limits<double>::infinity(); }
    if (p.m_cutoff_z == -1.0) { p.m_cutoff_z = std::numeric_limits<double>::infinity(); }

    bool const explicit_times = utils::parser::queryArrWithParser(pp, "pulse_times", m_times);
    bool const train = pp.contains("pulse_period") || pp.contains("pulse_count") ||
                       pp.contains("first_pulse_time");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(explicit_times != train,
        "Specify pulse_times or a finite train with first_pulse_time, pulse_period and pulse_count.");
    if (train) {
        double first = 0.0, period = 0.0;
        int count = 0;
        utils::parser::getWithParser(pp, "first_pulse_time", first);
        utils::parser::getWithParser(pp, "pulse_period", period);
        utils::parser::getWithParser(pp, "pulse_count", count);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(first) && positive(period) && count > 0,
            "A rigid-beam pulse train requires a finite start, positive period and count.");
        m_times.resize(count);
        for (int i = 0; i < count; ++i) { m_times[i] = first+i*period; }
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_times.empty(), "A rigid beam requires at least one pulse.");
    m_amplitudes.assign(m_times.size(), 1.0);
    utils::parser::queryArrWithParser(pp, "pulse_amplitudes", m_amplitudes);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(m_times.size() == m_amplitudes.size(),
        "pulse_amplitudes must have one entry per pulse.");
    for (amrex::Long i = 0; i < m_times.size(); ++i) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(m_times[i]) &&
                                            std::isfinite(m_amplitudes[i]) && m_amplitudes[i] >= 0.0,
            "Rigid-beam pulse times must be finite and amplitudes finite and non-negative.");
    }
    amrex::Vector<std::pair<double, double>> pulses;
    for (amrex::Long i = 0; i < m_times.size(); ++i) {
        pulses.emplace_back(m_times[i], m_amplitudes[i]);
    }
    std::stable_sort(pulses.begin(), pulses.end());
    for (amrex::Long i = 0; i < m_times.size(); ++i) {
        m_times[i] = pulses[i].first;
        m_amplitudes[i] = pulses[i].second;
    }
    double current = 0.0, bunch_charge = 0.0;
    bool const by_current = utils::parser::queryWithParser(pp, "peak_current", current);
    bool const by_charge = utils::parser::queryWithParser(pp, "bunch_charge", bunch_charge);
    bool const by_density = utils::parser::queryWithParser(pp, "peak_density", p.m_peak_density);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(int(by_current)+int(by_charge)+int(by_density) == 1,
        "Specify exactly one rigid-beam normalization: peak_current, bunch_charge or peak_density.");
    double const area = p.radialIntegral(0.0, p.m_cutoff_r*p.m_sigma_r);
    if (by_current) { p.m_peak_density = current/(std::abs(charge*p.m_velocity)*area); }
    if (by_charge) {
        p.m_peak_density = bunch_charge/(std::abs(charge)*area*p.m_sigma_z*
            warpx::fluid::gaussianIntegral(-p.m_cutoff_z, p.m_cutoff_z, p.m_cutoff_z));
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(positive(p.m_peak_density),
        "Rigid-beam normalization must specify a finite positive magnitude.");
    p.m_count = static_cast<int>(m_times.size());
    m_times_d.resize(m_times.size());
    m_amplitudes_d.resize(m_times.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_times.begin(), m_times.end(), m_times_d.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, m_amplitudes.begin(), m_amplitudes.end(),
                     m_amplitudes_d.begin());
}

RigidBeam::Executor
RigidBeam::executor () const
{
    auto result = m_parameters;
    result.m_times = m_times_d.data();
    result.m_amplitudes = m_amplitudes_d.data();
    return result;
}

std::string
RigidBeam::configuration () const
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    auto const& p = m_parameters;
    out << p.m_velocity << ' ' << p.m_sigma_r << ' ' << p.m_sigma_z << ' '
        << p.m_peak_density << ' ' << p.m_z_reference << ' ' << p.m_cutoff_r << ' '
        << p.m_cutoff_z << ' ' << m_times.size();
    for (amrex::Long i = 0; i < m_times.size(); ++i) { out << ' ' << m_times[i] << ' ' << m_amplitudes[i]; }
    return out.str();
}

bool
RigidBeam::active (amrex::Geometry const& geom, amrex::Real start, amrex::Real dt,
                   amrex::Real padding) const
{
#ifdef WARPX_DIM_RZ
    auto parameters = m_parameters;
    parameters.m_times = m_times.data();
    int first, last;
    parameters.pulseRange(geom.ProbLo(1), geom.ProbHi(1), start, dt, padding, first, last);
    return first < last && geom.ProbLo(0) < parameters.m_cutoff_r*parameters.m_sigma_r+padding;
#else
    amrex::ignore_unused(geom, start, dt, padding);
    return false;
#endif
}

void
RigidBeam::CacheRadial (amrex::MultiFab const& density, amrex::Geometry const& geom)
{
#ifdef WARPX_DIM_RZ
    if (!m_radial.empty() && m_radial_nodal == density.ixType()[0]) { return; }
    m_radial_nodal = density.ixType()[0];
    m_current_time = std::numeric_limits<amrex::Real>::lowest();
    auto const box = amrex::convert(geom.Domain(), density.ixType());
    m_radial_lo = box.smallEnd(0);
    int const count = box.length(0);
    m_radial.resize(count);
    auto* factors = m_radial.data();
    auto const p = executor();
    auto const dr = geom.CellSize(0);
    auto const rmax = geom.ProbHi(0);
    auto const rmin = geom.ProbLo(0);
    int const nodal = density.ixType()[0];
    int const order = WarpX::nox;
    bool correction = true;
    amrex::ParmParse("boundary").query("verboncoeur_axis_correction", correction);
    double const axis_volume = correction ? 1.0/3.0 : 1.0/4.0;
    amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int i) noexcept {
        double const node = rmin + (i+0.5*(1-nodal))*dr;
        double integral = 0.0;
        for (int mirror = 0; mirror < 2; ++mirror) {
            if (mirror == 1 && (node == 0.0 || rmin != 0.0)) { continue; }
            double const center = mirror == 0 ? node : -node;
            for (int piece = 0; piece <= order; ++piece) {
                double const lo = std::max(rmin, center+(piece-0.5*(order+1))*dr);
                double const hi = std::min(std::min(rmax, p.m_cutoff_r*p.m_sigma_r),
                                            center+(piece+1-0.5*(order+1))*dr);
                if (hi <= lo) { continue; }
                int const subdivisions = 1+static_cast<int>((hi-lo)/p.m_sigma_r);
                double const step = (hi-lo)/subdivisions;
                for (int sub = 0; sub < subdivisions; ++sub) {
                    integral += warpx::fluid::integrate([=] AMREX_GPU_HOST_DEVICE(double r) noexcept {
                        double const x = r/p.m_sigma_r;
                        return r*std::exp(-0.5*x*x)*warpx::fluid::shape(order, (r-center)/dr);
                    }, lo+sub*step, lo+(sub+1)*step);
                }
            }
        }
        factors[i] = node == 0.0 ? 2*integral/(dr*dr*axis_volume) : integral/(dr*node);
    });
#else
    amrex::ignore_unused(density, geom);
#endif
}

void
RigidBeam::UpdateDensity (amrex::MultiFab& density, amrex::Geometry const& geom, amrex::Real time)
{
#ifdef WARPX_DIM_RZ
    CacheRadial(density, geom);
    auto const box = amrex::convert(geom.Domain(), density.ixType());
    m_longitudinal_lo = box.smallEnd(1);
    int const count = box.length(1);
    m_longitudinal.resize(count);
    auto* longitudinal = m_longitudinal.data();
    auto const* radial = m_radial.data();
    auto const p = executor();
    auto const dz = geom.CellSize(1);
    auto const zmin = geom.ProbLo(1);
    int const nodal = density.ixType()[1];
    int const order = WarpX::nox;
    if (nodal != m_axial_nodal) {
        m_density_time = std::numeric_limits<amrex::Real>::lowest();
        m_axial_nodal = nodal;
    }
    if (m_density_time != time) {
        amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int j) noexcept {
            double const node = zmin + (j+0.5*(1-nodal))*dz;
            double value = 0.0;
            int first, last;
            p.pulseRange(node-dz/2, node+dz/2, time, 0.0, order*dz/2, first, last);
            for (int pulse = first; pulse < last; ++pulse) {
                double const center = p.m_z_reference+p.m_velocity*(time-p.m_times[pulse]);
                value += p.m_amplitudes[pulse]*warpx::fluid::projectedIntegral(
                    order-1, node-center-dz/2, node-center+dz/2, dz, p.m_sigma_z, p.m_cutoff_z)/dz;
            }
            longitudinal[j] = p.m_peak_density*value;
        });
        m_density_time = time;
    }
    int const rlo = m_radial_lo, zlo = m_longitudinal_lo;
    for (amrex::MFIter mfi(density, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const a = density.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i,j,k) = static_cast<amrex::Real>(radial[i-rlo]*longitudinal[j-zlo]);
        });
    }
    density.FillBoundary(geom.periodicity());
#else
    amrex::ignore_unused(density, geom, time);
#endif
}

void
RigidBeam::DepositCurrent (amrex::MultiFab& current, amrex::Geometry const& geom,
                           amrex::Real start, amrex::Real dt)
{
#ifdef WARPX_DIM_RZ
    if (!active(geom, start, dt, (WarpX::nox+1)*geom.CellSize(1))) { return; }
    auto const box = amrex::convert(geom.Domain(), current.ixType());
    int const count = box.length(1);
    m_current.resize(count);
    auto* axial = m_current.data();
    auto const* radial = m_radial.data();
    auto const p = executor();
    auto const dz = geom.CellSize(1);
    auto const zmin = geom.ProbLo(1);
    int const nodal = current.ixType()[1];
    int const order = WarpX::nox;
    bool const spectral = WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD;
    auto const charge = m_charge;
    if (m_current_time != start || m_current_dt != dt) {
        amrex::ParallelFor(count, [=] AMREX_GPU_DEVICE(int j) noexcept {
            double const face = zmin+(j+0.5*(1-nodal))*dz;
            double value = 0.0;
            int first, last;
            p.pulseRange(face, face, start, dt, (order+1)*dz/2, first, last);
            for (int pulse = first; pulse < last; ++pulse) {
                double const a = face-p.m_z_reference-p.m_velocity*(start-p.m_times[pulse]);
                double const b = a-p.m_velocity*dt;
                double const sign = p.m_velocity > 0.0 ? 1.0 : -1.0;
                if (dt > 0.0) {
                    value += p.m_amplitudes[pulse]*sign*warpx::fluid::projectedIntegral(
                        spectral ? order : order-1, std::min(a,b), std::max(a,b),
                        dz, p.m_sigma_z, p.m_cutoff_z)/dt;
                } else {
                    value += p.m_amplitudes[pulse]*p.m_velocity*warpx::fluid::projectedIntegral(
                        order-1, a-dz/2, a+dz/2, dz, p.m_sigma_z, p.m_cutoff_z)/dz;
                }
            }
            axial[j] = charge*p.m_peak_density*value;
        });
        m_current_time = start;
        m_current_dt = dt;
    }
    auto const mask = amrex::OwnerMask(current, geom.periodicity());
    int const rlo = m_radial_lo, zlo = box.smallEnd(1);
    for (amrex::MFIter mfi(current, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const a = current.array(mfi);
        auto const owner = mask->const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (owner(i,j,k)) { a(i,j,k) += static_cast<amrex::Real>(radial[i-rlo]*axial[j-zlo]); }
        });
    }
#else
    amrex::ignore_unused(current, geom, start, dt);
#endif
}
