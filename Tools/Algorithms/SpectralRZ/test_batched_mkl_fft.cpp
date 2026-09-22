/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <ablastr/math/fft/BatchedMklFFT.H>

#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Print.H>

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
    void checkTransform (int const length, int const batch)
    {
        using Complex = amrex::GpuComplex<amrex::Real>;
        constexpr int modes = 3;
        int const size = length * batch;
        std::vector<Complex> input(modes * size), result(modes * size);
        for (int m = 0; m < modes; ++m) {
            for (int z = 0; z < length; ++z) {
                for (int r = 0; r < batch; ++r) {
                    // Distinct complex data in every radial sequence and mode.
                    input[m * size + z * batch + r] = Complex{
                        amrex::Real(0.1 * (1 + r) + std::cos(0.7 * z + m)),
                        amrex::Real(0.2 * m + std::sin(0.3 * z - r))};
                }
            }
        }
        amrex::Gpu::DeviceVector<Complex> input_d(input.size()), output_d(input.size());
        ablastr::math::anyfft::BatchedMklFFT plan(length, batch);
        auto const pi = std::acos(-1.0L);
        auto const tolerance = 128 * length * std::numeric_limits<amrex::Real>::epsilon();
        for (int stream = 0; stream < amrex::Gpu::numGpuStreams(); ++stream) {
            amrex::Gpu::Device::setStreamIndex(stream);
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, input.begin(), input.end(),
                                  input_d.begin());
            for (int m = 0; m < modes; ++m) {
                plan.Forward(input_d.data() + m * size, output_d.data() + m * size);
            }
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, output_d.begin(), output_d.end(),
                             result.begin());
            // Direct long-double DFT, independent of oneMKL's FFT algorithm.
            for (int m = 0; m < modes; ++m) {
                for (int k = 0; k < length; ++k) {
                    for (int r = 0; r < batch; ++r) {
                        std::complex<long double> expected{0, 0};
                        for (int z = 0; z < length; ++z) {
                            auto const value = input[m * size + z * batch + r];
                            auto const angle = -2 * pi * k * z / length;
                            expected += std::complex<long double>{value.real(), value.imag()}
                                        * std::polar(1.0L, angle);
                        }
                        auto const value = result[m * size + k * batch + r];
                        auto const actual = std::complex<long double>{value.real(), value.imag()};
                        if (std::abs(actual - expected) > tolerance) {
                            throw std::runtime_error("Batched FFT differs from direct DFT");
                        }
                    }
                }
            }
            for (int m = 0; m < modes; ++m) {
                plan.Backward(output_d.data() + m * size, input_d.data() + m * size);
            }
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, input_d.begin(), input_d.end(),
                             result.begin());
            for (int i = 0; i < modes * size; ++i) {
                auto const error = std::hypot(result[i].real() / length - input[i].real(),
                                              result[i].imag() / length - input[i].imag());
                if (error > tolerance) {
                    throw std::runtime_error("Batched FFT round trip has incorrect normalization");
                }
            }
        }
        amrex::Gpu::Device::resetStreamIndex();
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        for (int length : {7, 12}) {
            for (int batch : {1, 5}) {
                checkTransform(length, batch);
            }
        }
        amrex::Print() << "Interleaved complex FFT, normalization and queue changes pass\n";
    }
    amrex::Finalize();
}
