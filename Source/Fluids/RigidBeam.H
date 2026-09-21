/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_RIGID_BEAM_H_
#define WARPX_RIGID_BEAM_H_

#include "RigidBeamProfile.H"

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <limits>
#include <string>

/** Separable, prescribed Gaussian pulses in the laboratory frame. */
class RigidBeam
{
public:
    struct Executor
    {
        double m_velocity = 0.0;
        double m_sigma_r = 0.0;
        double m_sigma_z = 0.0;
        double m_peak_density = 0.0;
        double m_z_reference = 0.0;
        double m_cutoff_r = 8.0;
        double m_cutoff_z = 8.0;
        double const* m_times = nullptr;
        double const* m_amplitudes = nullptr;
        int m_count = 0;

        /** Find pulses whose truncated support intersects a space-time interval. */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
        pulseRange (double lo, double hi, double t, double dt, double padding,
                    int& first, int& last) const noexcept
        {
            double const extent = m_cutoff_z*m_sigma_z+padding;
            double const a = t-(hi+extent-m_z_reference)/m_velocity;
            double const b = t-(lo-extent-m_z_reference)/m_velocity;
            double const lower = std::min(a,b), upper = std::max(a,b)+dt;
            int left = 0, right = m_count;
            while (left < right) {
                int const middle = left+(right-left)/2;
                if (m_times[middle] < lower) { left = middle+1; } else { right = middle; }
            }
            first = left;
            right = m_count;
            while (left < right) {
                int const middle = left+(right-left)/2;
                if (m_times[middle] <= upper) { left = middle+1; } else { right = middle; }
            }
            last = left;
        }

        /** Integral of longitudinal number density over [lo,hi] at time t. */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
        longitudinalIntegral (double lo, double hi, double t) const noexcept
        {
            double sum = 0.0;
            int first, last;
            pulseRange(lo, hi, t, 0.0, 0.0, first, last);
            for (int pulse = first; pulse < last; ++pulse) {
                auto const center = m_z_reference + m_velocity*(t-m_times[pulse]);
                sum += m_amplitudes[pulse]*m_sigma_z*warpx::fluid::gaussianIntegral(
                    (lo-center)/m_sigma_z, (hi-center)/m_sigma_z, m_cutoff_z);
            }
            return m_peak_density*sum;
        }

        /** Space-time integral in z. Subintervals resolve translation by sigma_z. */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
        integratedLongitudinal (double lo, double hi, double t, double dt) const noexcept
        {
            double sum = 0.0;
            auto const profile = *this;
            int first, last;
            pulseRange(lo, hi, t, dt, 0.0, first, last);
            for (int pulse = first; pulse < last; ++pulse) {
                double cuts[6] = {t, t+dt};
                int count = 2;
                for (int boundary = 0; boundary < 2; ++boundary) {
                    double const z = boundary == 0 ? lo : hi;
                    for (int sign = -1; sign <= 1; sign += 2) {
                        double const crossing = m_times[pulse] +
                            (z-m_z_reference+sign*m_cutoff_z*m_sigma_z)/m_velocity;
                        if (crossing > t && crossing < t+dt) { cuts[count++] = crossing; }
                    }
                }
                for (int i = 1; i < count; ++i) {
                    for (int j = i; j > 0 && cuts[j] < cuts[j-1]; --j) {
                        double const value = cuts[j]; cuts[j] = cuts[j-1]; cuts[j-1] = value;
                    }
                }
                for (int segment = 1; segment < count; ++segment) {
                    int const intervals = 1 + static_cast<int>(std::abs(m_velocity)*
                        (cuts[segment]-cuts[segment-1])/m_sigma_z);
                    double const step = (cuts[segment]-cuts[segment-1])/intervals;
                    for (int i = 0; i < intervals; ++i) {
                        sum += m_amplitudes[pulse]*warpx::fluid::integrate(
                            [=] AMREX_GPU_HOST_DEVICE(double time) noexcept {
                                double const center = profile.m_z_reference+
                                    profile.m_velocity*(time-profile.m_times[pulse]);
                                return profile.m_sigma_z*warpx::fluid::gaussianIntegral(
                                    (lo-center)/profile.m_sigma_z, (hi-center)/profile.m_sigma_z,
                                    profile.m_cutoff_z);
                            }, cuts[segment-1]+i*step, cuts[segment-1]+(i+1)*step);
                    }
                }
            }
            return m_peak_density*sum;
        }

        /** Integral of the transverse profile over an annulus, including 2*pi*r. */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
        radialIntegral (double lo, double hi) const noexcept
        {
            lo = lo < 0.0 ? 0.0 : lo;
            hi = hi > m_cutoff_r*m_sigma_r ? m_cutoff_r*m_sigma_r : hi;
            if (hi <= lo) { return 0.0; }
            auto const a = lo/m_sigma_r;
            auto const b = hi/m_sigma_r;
            return 6.2831853071795864769*m_sigma_r*m_sigma_r*std::exp(-0.5*a*a)*
                (-std::expm1(-0.5*(b-a)*(b+a)));
        }
    };

    RigidBeam (std::string const& name, amrex::Real mass, amrex::Real charge);

    [[nodiscard]] Executor executor () const;
    [[nodiscard]] amrex::Real velocity () const { return m_parameters.m_velocity; }
    [[nodiscard]] amrex::Real kineticEnergyEV () const { return m_energy_ev; }
    [[nodiscard]] std::string configuration () const;
    [[nodiscard]] bool active (amrex::Geometry const& geom, amrex::Real start,
                               amrex::Real dt, amrex::Real padding = 0.0) const;

    void UpdateDensity (amrex::MultiFab& density, amrex::Geometry const& geom,
                        amrex::Real time);
    void DepositCurrent (amrex::MultiFab& current, amrex::Geometry const& geom,
                         amrex::Real start, amrex::Real dt);
    void UpdateCurrentDiagnostic (amrex::MultiFab& current, amrex::Geometry const& geom,
                                  amrex::Real time);

    // NVCC requires public enclosing functions for extended device lambdas.
    void CacheRadial (amrex::MultiFab const& density, amrex::Geometry const& geom);

private:
    Executor m_parameters;
    [[maybe_unused]] amrex::Real m_charge;
    amrex::Real m_energy_ev = 0.0;
    amrex::Vector<double> m_times;
    amrex::Vector<double> m_amplitudes;
    amrex::Gpu::DeviceVector<double> m_times_d;
    amrex::Gpu::DeviceVector<double> m_amplitudes_d;
    amrex::Gpu::DeviceVector<double> m_radial;
    amrex::Gpu::DeviceVector<double> m_longitudinal;
    amrex::Gpu::DeviceVector<double> m_current;
    int m_radial_lo = 0;
    int m_longitudinal_lo = 0;
    int m_radial_nodal = -1;
    int m_axial_nodal = -1;
    amrex::Real m_density_time = std::numeric_limits<amrex::Real>::lowest();
    amrex::Real m_current_time = std::numeric_limits<amrex::Real>::lowest();
    amrex::Real m_current_dt = 0.0;
};

#endif
