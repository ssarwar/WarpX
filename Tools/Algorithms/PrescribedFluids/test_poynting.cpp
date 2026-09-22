/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Diagnostics/ReducedDiags/FieldPoyntingFlux.H"

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
// Linear fields give exact face interpolation. Inactive coordinates carry
// sentinel variation, so using a nonexistent dimension cannot pass silently.
AMREX_GPU_HOST_DEVICE amrex::Real
electric (int i, int j, int k)
{
    return amrex::Real(2 + i + 3 * j + 5 * k);
}

AMREX_GPU_HOST_DEVICE amrex::Real
magnetic (int i, int j, int k)
{
    return amrex::Real(7 + 2 * i + 11 * j + 17 * k);
}

void
require (amrex::Real actual, amrex::Real expected, char const* message)
{
    if (std::abs(actual - expected) >
        16 * std::numeric_limits<amrex::Real>::epsilon() *
            std::max(amrex::Real(1), std::abs(expected))) {
        amrex::Print() << message << ": " << actual << " != " << expected
                       << '\n';
        throw std::runtime_error(message);
    }
}
} // namespace

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    {
        constexpr int width = 5, cells = width * width * width;
        amrex::Gpu::DeviceVector<amrex::Real> e(cells), b(cells), result(4);
        amrex::Array4<amrex::Real> const ea(e.data(), {0, 0, 0},
                                            {width, width, width}, 1);
        amrex::Array4<amrex::Real> const ba(b.data(), {0, 0, 0},
                                            {width, width, width}, 1);
        amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE(int index) noexcept {
            int const i = index % width, j = (index / width) % width,
                      k = index / (width * width);
            ea(i, j, k) = electric(i, j, k);
            ba(i, j, k) = magnetic(i, j, k);
        });
        auto* values = result.data();
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) noexcept {
            values[0] = PoyntingCellCentered::ExBy(2, 2, 2, 0, ea, ba);
            values[1] = PoyntingCellCentered::EyBx(2, 2, 2, 0, ea, ba);
            values[2] = PoyntingStaggered::ExBy(2, 2, 2, 0, ea, ba);
            values[3] = PoyntingStaggered::EyBx(2, 2, 2, 0, ea, ba);
        });
        std::vector<amrex::Real> host(4);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, result.begin(), result.end(),
                         host.begin());
#if defined(WARPX_DIM_1D_Z)
        amrex::Real const expected = (electric(1, 2, 2) * magnetic(1, 2, 2) +
                                      electric(2, 2, 2) * magnetic(2, 2, 2)) /
                                     2;
#elif defined(WARPX_DIM_RZ)
        amrex::Real const expected = (electric(2, 1, 2) * magnetic(2, 1, 2) +
                                      electric(2, 2, 2) * magnetic(2, 2, 2)) /
                                     2;
#else
        amrex::Real const expected = (electric(2, 2, 1) * magnetic(2, 2, 1) +
                                      electric(2, 2, 2) * magnetic(2, 2, 2)) /
                                     2;
#endif
        require(host[0], expected,
                "Cell-centered ExBy uses the wrong axial dimension");
        require(host[1], expected,
                "Cell-centered EyBx uses the wrong axial dimension");
        // In Yee's axial flux, B is first averaged across the axial face to
        // E's location. Only nodal transverse E then needs a transverse average.
#if defined(WARPX_DIM_1D_Z)
        auto const staggered_exby = electric(2, 2, 2) *
            (magnetic(1, 2, 2) + magnetic(2, 2, 2)) / 2;
        auto const staggered_eybx = staggered_exby;
#elif defined(WARPX_DIM_RZ)
        auto const staggered_exby = electric(2, 2, 2) *
            (magnetic(2, 1, 2) + magnetic(2, 2, 2)) / 2;
        auto const staggered_eybx =
            (electric(2, 2, 2) * (magnetic(2, 1, 2) + magnetic(2, 2, 2)) +
             electric(3, 2, 2) * (magnetic(3, 1, 2) + magnetic(3, 2, 2))) / 4;
#else
        auto const staggered_exby =
            (electric(2, 2, 2) * (magnetic(2, 2, 1) + magnetic(2, 2, 2)) +
             electric(2, 3, 2) * (magnetic(2, 3, 1) + magnetic(2, 3, 2))) / 4;
        auto const staggered_eybx =
            (electric(2, 2, 2) * (magnetic(2, 2, 1) + magnetic(2, 2, 2)) +
             electric(3, 2, 2) * (magnetic(3, 2, 1) + magnetic(3, 2, 2))) / 4;
#endif
        require(host[2], staggered_exby, "Yee ExBy uses the wrong axial dimension");
        require(host[3], staggered_eybx, "Yee EyBx uses the wrong axial dimension");
#ifdef WARPX_DIM_RZ
        amrex::Gpu::DeviceVector<amrex::Real> zero(cells, amrex::Real(0));
        amrex::Array4<amrex::Real const> const empty(zero.data(), {0, 0, 0},
                                                     {width, width, width}, 1);
        amrex::Box const face(amrex::IntVect(2), amrex::IntVect(2));
        auto const p = amrex::lbound(face);
        auto const area = [=] AMREX_GPU_DEVICE(int i, int, int) noexcept {
            return amrex::Real(2 * i);
        };
        auto const flux = Poynting::Kernel<0, PoyntingCellCentered>(
            face, empty, empty, ea, empty, ba, empty, area);
        auto const radial_expected =
            -amrex::Real(2 * p.x) *
            (electric(p.x - 1, p.y, p.z) * magnetic(p.x - 1, p.y, p.z) +
             electric(p.x, p.y, p.z) * magnetic(p.x, p.y, p.z)) /
            2;
        require(flux, radial_expected,
                "Cell-centered RZ flux used Yee staggering");
#endif
        amrex::Print()
            << "PASS: exact cell-centered and Yee axial Poynting interpolation\n";
    }
    amrex::Finalize();
}
