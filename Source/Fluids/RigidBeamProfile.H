/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_RIGID_BEAM_PROFILE_H_
#define WARPX_RIGID_BEAM_PROFILE_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <algorithm>
#include <cmath>

namespace warpx::fluid
{
    /** Integral of exp(-x^2/2) over [lo,hi], restricted to [-cutoff,cutoff].
     *
     * Positive-tail differences use erfc to retain the small population near
     * the truncation. All profile integrals use double precision, including
     * simulations storing fields and particles in single precision.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    gaussianIntegral (double lo, double hi, double cutoff) noexcept
    {
        lo = lo < -cutoff ? -cutoff : lo;
        hi = hi > cutoff ? cutoff : hi;
        if (hi <= lo) {
            return 0.0;
        }
        constexpr double root_half = 0.70710678118654752440;
        constexpr double root_pi_half = 1.2533141373155002512;
        if (lo >= 0.0) {
            return root_pi_half * (std::erfc(lo * root_half) - std::erfc(hi * root_half));
        }
        if (hi <= 0.0) {
            return root_pi_half * (std::erfc(-hi * root_half) - std::erfc(-lo * root_half));
        }
        return root_pi_half * (std::erf(hi * root_half) - std::erf(lo * root_half));
    }

    /** Compact centered cardinal B-splines, matching particle shape orders 0--4. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    shape (int order, double x) noexcept
    {
        x = std::abs(x);
        if (order == 0) {
            return x < 0.5 ? 1.0 : 0.0;
        }
        if (order == 1) {
            return x < 1.0 ? 1.0 - x : 0.0;
        }
        if (order == 2) {
            if (x < 0.5) {
                return 0.75 - x * x;
            }
            if (x < 1.5) {
                return 0.5 * (1.5 - x) * (1.5 - x);
            }
            return 0.0;
        }
        if (order == 3) {
            if (x < 1.0) {
                return 2.0 / 3.0 - x * x + 0.5 * x * x * x;
            }
            if (x < 2.0) {
                return (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
            }
            return 0.0;
        }
        if (x < 0.5) {
            return 115.0 / 192.0 - 5.0 * x * x / 8.0 + x * x * x * x / 4.0;
        }
        if (x < 1.5) {
            return (55.0 + 20.0 * x - 120.0 * x * x + 80.0 * x * x * x - 16.0 * x * x * x * x) /
                   96.0;
        }
        if (x < 2.5) {
            auto const d = 2.5 - x;
            return d * d * d * d / 24.0;
        }
        return 0.0;
    }

    /** Eight-point Gauss rule on one smooth interval. */
    template <typename F>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    integrate (F const& f, double lo, double hi) noexcept
    {
        constexpr double x[4] = {0.183434642495649805, 0.525532409916328986, 0.796666477413626740,
                                 0.960289856497536232};
        constexpr double w[4] = {0.362683783378361983, 0.313706645877887287, 0.222381034453374471,
                                 0.101228536290376259};
        auto const mid = 0.5 * (lo + hi);
        auto const half = 0.5 * (hi - lo);
        double sum = 0.0;
        for (int i = 0; i < 4; ++i) {
            sum += w[i] * (f(mid - half * x[i]) + f(mid + half * x[i]));
        }
        return half * sum;
    }

    /** Gaussian convolved with a compact particle shape at one coordinate. */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    projectedDensity (int order, double coordinate, double dx, double sigma, double cutoff) noexcept
    {
        double result = 0.0;
        for (int piece = 0; piece <= order; ++piece) {
            double const lo =
                std::max(-0.5 * (order + 1) + piece, (coordinate - cutoff * sigma) / dx);
            double const hi =
                std::min(-0.5 * (order + 1) + piece + 1, (coordinate + cutoff * sigma) / dx);
            if (hi <= lo) {
                continue;
            }
            int const count = 1 + static_cast<int>((hi - lo) * dx / sigma);
            double const step = (hi - lo) / count;
            for (int sub = 0; sub < count; ++sub) {
                result += integrate(
                    [=] AMREX_GPU_HOST_DEVICE(double u) noexcept {
                        double const x = (coordinate - dx * u) / sigma;
                        return shape(order, u) * std::exp(-0.5 * x * x);
                    },
                    lo + sub * step, lo + (sub + 1) * step);
            }
        }
        return result;
    }

    /** Gaussian convolved with a compact particle shape, integrated over an interval.
     *
     * The cumulative shape at an axial current face uses order p-1. Differences
     * of the same cumulative primitive in space give the order-p charge; its
     * differences in time give the step-averaged current. Thus the discrete Yee
     * continuity equation does not depend on quadrature accuracy.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    projectedIntegral (int order, double lo, double hi, double dx, double sigma,
                       double cutoff) noexcept
    {
        double const extent = cutoff * sigma + 0.5 * (order + 1) * dx;
        if (hi <= -extent || lo >= extent || hi <= lo) {
            return 0.0;
        }
        double result = 0.0;
        // Split at spline knots and resolve the Gaussian even for dx > sigma.
        int const subdivisions = static_cast<int>(std::ceil(dx / sigma));
        auto const width = 1.0 / subdivisions;
        for (int piece = 0; piece <= order; ++piece) {
            for (int sub = 0; sub < subdivisions; ++sub) {
                auto const a = -0.5 * (order + 1) + piece + sub * width;
                double cuts[6] = {a, a + width};
                int count = 2;
                double const crossings[4] = {(lo - cutoff * sigma) / dx, (lo + cutoff * sigma) / dx,
                                             (hi - cutoff * sigma) / dx,
                                             (hi + cutoff * sigma) / dx};
                for (auto crossing : crossings) {
                    if (crossing > a && crossing < a + width) {
                        cuts[count++] = crossing;
                    }
                }
                for (int i = 1; i < count; ++i) {
                    for (int j = i; j > 0 && cuts[j] < cuts[j - 1]; --j) {
                        double const value = cuts[j];
                        cuts[j] = cuts[j - 1];
                        cuts[j - 1] = value;
                    }
                }
                for (int i = 1; i < count; ++i) {
                    result += integrate(
                        [=] AMREX_GPU_HOST_DEVICE(double u) noexcept {
                            return shape(order, u) * sigma *
                                   gaussianIntegral((lo - dx * u) / sigma, (hi - dx * u) / sigma,
                                                    cutoff);
                        },
                        cuts[i - 1], cuts[i]);
                }
            }
        }
        return result;
    }
} // namespace warpx::fluid

#endif
