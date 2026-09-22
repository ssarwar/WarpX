/* Copyright 2026 The WarpX Community
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_BATCHED_MKL_FFT_H_
#define ABLASTR_BATCHED_MKL_FFT_H_

#include <AMReX_Config.H>

#ifdef AMREX_USE_SYCL

#include <AMReX_Gpu.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_REAL.H>

#include <oneapi/mkl/dft.hpp>

#include <complex>
#include <cstdint>
#include <vector>

namespace ablastr::math::anyfft
{
    /** Unnormalized complex one-dimensional FFTs of interleaved sequences.
     *
     * The element at position z in sequence r is stored at r + batch*z.
     * Both directions use the current AMReX in-order queue. A queue change
     * waits for previous transforms before recommitting their descriptor;
     * transforms on the same queue remain asynchronous.
     */
    class BatchedMklFFT
    {
        using Descriptor = oneapi::mkl::dft::descriptor<
#ifdef AMREX_USE_FLOAT
            oneapi::mkl::dft::precision::SINGLE,
#else
            oneapi::mkl::dft::precision::DOUBLE,
#endif
            oneapi::mkl::dft::domain::COMPLEX>;
        using Complex = amrex::GpuComplex<amrex::Real>;

    public:
        BatchedMklFFT (std::int64_t length, std::int64_t batch)
            : m_plan(length), m_stream(amrex::Gpu::gpuStream())
        {
            using oneapi::mkl::dft::config_param;
            using oneapi::mkl::dft::config_value;
            std::vector<std::int64_t> const strides{0, batch};
            m_plan.set_value(config_param::PLACEMENT, config_value::NOT_INPLACE);
            m_plan.set_value(config_param::NUMBER_OF_TRANSFORMS, batch);
            m_plan.set_value(config_param::FWD_STRIDES, strides);
            m_plan.set_value(config_param::BWD_STRIDES, strides);
            m_plan.set_value(config_param::FWD_DISTANCE, std::int64_t{1});
            m_plan.set_value(config_param::BWD_DISTANCE, std::int64_t{1});
            m_plan.commit(amrex::Gpu::Device::streamQueue());
        }

        ~BatchedMklFFT ()
        {
            // The descriptor can own workspace used by queued transforms.
            amrex::Gpu::streamSynchronize(m_stream);
        }

        BatchedMklFFT (BatchedMklFFT const&) = delete;
        BatchedMklFFT& operator= (BatchedMklFFT const&) = delete;
        BatchedMklFFT (BatchedMklFFT&&) = delete;
        BatchedMklFFT& operator= (BatchedMklFFT&&) = delete;

        void Forward (Complex* input, Complex* output)
        {
            SetStream();
            oneapi::mkl::dft::compute_forward(
                m_plan, reinterpret_cast<std::complex<amrex::Real>*>(input),
                reinterpret_cast<std::complex<amrex::Real>*>(output));
        }

        void Backward (Complex* input, Complex* output)
        {
            SetStream();
            oneapi::mkl::dft::compute_backward(
                m_plan, reinterpret_cast<std::complex<amrex::Real>*>(input),
                reinterpret_cast<std::complex<amrex::Real>*>(output));
        }

    private:
        void SetStream ()
        {
            if (!(m_stream == amrex::Gpu::gpuStream())) {
                amrex::Gpu::streamSynchronize(m_stream);
                m_plan.commit(amrex::Gpu::Device::streamQueue());
                m_stream = amrex::Gpu::gpuStream();
            }
        }

        Descriptor m_plan;
        amrex::gpuStream_t m_stream;
    };
}

#endif // AMREX_USE_SYCL
#endif // ABLASTR_BATCHED_MKL_FFT_H_
