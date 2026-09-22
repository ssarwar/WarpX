#!/usr/bin/env bash
# Supplement the stock Perlmutter CUDA build with a oneAPI portability build.
set -eo pipefail
root=$(git rev-parse --show-toplevel)
module load cmake/3.30.2 intel-oneapi-mixed/2025.3
export CC=$(command -v icx)
export CXX=$(command -v icpx)
sycl_include="$(dirname "$CXX")/../include/sycl"
test -f "$sycl_include/sycl.hpp"
# BLAS++ 2024.05.31 includes <sycl.hpp>; oneAPI 2025 installs that header
# below include/sycl. Add the installed include directory without changing
# either dependency's sources.
sycl_flags="-fsycl -isystem $sycl_include -DMKL_ILP64"
export CXXFLAGS="$sycl_flags"
build="$root/build_pm_sycl"
audit="$root/build/reaudit-2026-09-21"
prefix="$audit/software/sycl"
export CMAKE_PREFIX_PATH="$prefix/blaspp:$prefix/lapackpp:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$prefix/blaspp/lib64:$prefix/lapackpp/lib64:${LD_LIBRARY_PATH:-}"

# RZ requires the SYCL BLAS++ backend, not the CUDA installation. These are
# the stock Aurora dependency commands with isolated paths and no cleanup.
CXXFLAGS=-qmkl cmake -S "$audit/sources/blaspp" -B "$audit/dependencies/blaspp-sycl" \
    -DCMAKE_CXX_FLAGS="$sycl_flags -qmkl" \
    -Duse_openmp=OFF -Dgpu_backend=sycl -DCMAKE_CXX_STANDARD=20 \
    -Dblas_int=int64 -DMKL_INTERFACE=ilp64 \
    -DCMAKE_INSTALL_PREFIX="$prefix/blaspp" -DCMAKE_EXE_LINKER_FLAGS=-qmkl
cmake --build "$audit/dependencies/blaspp-sycl" --target install --parallel 16
CXXFLAGS='-DLAPACK_FORTRAN_ADD_ -qmkl' \
    cmake --fresh -S "$audit/sources/lapackpp" -B "$audit/dependencies/lapackpp-sycl" \
        -DCMAKE_CXX_FLAGS="$sycl_flags -DLAPACK_FORTRAN_ADD_ -qmkl" \
    -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -Dgpu_backend=sycl \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX="$prefix/lapackpp" \
    -DCMAKE_EXE_LINKER_FLAGS=-qmkl
cmake --build "$audit/dependencies/lapackpp-sycl" --target install --parallel 16

# Follow the repository's oneAPI CI configuration, using a portable device
# image because Perlmutter has no Intel GPU. This alone is compile acceptance.
cmake -S "$root" -B "$build" -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_CXX_FLAGS="$sycl_flags" \
    -DCMAKE_BUILD_TYPE=Release -DWarpX_DIMS=RZ -DWarpX_COMPUTE=SYCL \
    -DAMReX_SYCL_AOT=OFF -DAMReX_PARALLEL_LINK_JOBS=4 \
    -DWarpX_MPI=OFF -DWarpX_EB=OFF -DWarpX_FFT=ON \
    -DWarpX_PYTHON=OFF -DWarpX_OPENPMD=OFF -DWarpX_QED=OFF \
    -DMKL_INTERFACE=ilp64 -DMKL_THREADING=sequential -DMKL_SYCL_THREADING=sequential
cmake --build "$build" --parallel 16
for entry in 'PrescribedFluids fluid-tests' 'ProtonImpactIonization pjg' \
    'BackgroundMCC mcc-physics' 'SpectralRZ spectral-tests'; do
    read -r source destination <<< "$entry"
    cmake -S "$root/Tools/Algorithms/$source" -B "$build/$destination" \
        -DCMAKE_BUILD_TYPE=Release \
        -DAMReX_DIR="$build/_deps/fetchedamrex-build/lib/cmake/AMReX" \
        -DMKL_THREADING=sequential -DMKL_SYCL_THREADING=sequential
    cmake --build "$build/$destination" --parallel 8
done
