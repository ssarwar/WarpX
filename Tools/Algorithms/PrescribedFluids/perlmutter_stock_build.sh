#!/usr/bin/env bash
# Build the coupled-physics audit with WarpX's documented Perlmutter recipes.
# Run from the validation checkout on Perlmutter; see REAUDIT.md.
set -eo pipefail

root=$(git rev-parse --show-toplevel)
audit="$root/build/reaudit-2026-09-21"
prefix="$audit/software"
phase=${1:-configure}
build_name=${WARPX_AUDIT_BUILD_NAME:-build_pm_gpu_py}
stock_revision=${WARPX_AUDIT_STOCK_REVISION:-66f380f98e44b0cce493a435305009693d0e261e}
mkdir -p "$audit"

# Account and installation paths are the only edits to the stock profile.
sed -e 's/^export proj=""/export proj="m3748_g"/' \
    -e "s@^export SW_DIR=.*@export SW_DIR=\"$prefix\"@" \
    "$root/Tools/machines/perlmutter-nersc/perlmutter_gpu_warpx.profile.example" \
    > "$audit/perlmutter_gpu_warpx.profile"
source "$audit/perlmutter_gpu_warpx.profile"

if [[ "$phase" == prepare ]]; then
    mkdir -p "$prefix" "$audit/sources" "$audit/dependencies"
    # Boost is header-only for these builds; the existing stock version is sufficient.
    if [[ ! -e "$prefix/boost-1.82.0" ]]; then
        ln -s "$PSCRATCH/storage/sw/warpx/perlmutter/gpu/boost-1.82.0" \
            "$prefix/boost-1.82.0"
    fi
    for entry in \
        'c-blosc https://github.com/Blosc/c-blosc.git v1.21.1' \
        'adios2 https://github.com/ornladios/ADIOS2.git v2.10.2' \
        'blaspp https://github.com/icl-utk-edu/blaspp.git v2024.05.31' \
        'lapackpp https://github.com/icl-utk-edu/lapackpp.git v2024.05.31'; do
        read -r name repository revision <<< "$entry"
        if [[ ! -d "$audit/sources/$name/.git" ]]; then
            git clone --branch "$revision" --depth 1 "$repository" "$audit/sources/$name"
        fi
        git -C "$audit/sources/$name" rev-parse HEAD > "$audit/$name-revision.txt"
        git -C "$audit/sources/$name" diff --exit-code
    done

    # These options follow install_gpu_dependencies.sh. Build and install paths
    # are confined to this checkout, so shared user dependencies remain intact.
    cmake -S "$audit/sources/c-blosc" -B "$audit/dependencies/c-blosc" \
        -DBUILD_TESTS=OFF -DBUILD_BENCHMARKS=OFF -DDEACTIVATE_AVX2=OFF \
        -DCMAKE_INSTALL_PREFIX="$prefix/c-blosc-1.21.1"
    cmake --build "$audit/dependencies/c-blosc" --target install --parallel 16
    cmake -S "$audit/sources/adios2" -B "$audit/dependencies/adios2" \
        -DADIOS2_USE_Blosc=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF \
        -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX="$prefix/adios2-2.10.2"
    cmake --build "$audit/dependencies/adios2" --target install --parallel 16
    CXX=$(command -v CC) cmake -S "$audit/sources/blaspp" -B "$audit/dependencies/blaspp" \
        -Duse_openmp=OFF -Dgpu_backend=cuda -DCMAKE_CXX_STANDARD=20 \
        -DCMAKE_INSTALL_PREFIX="$prefix/blaspp-2024.05.31"
    cmake --build "$audit/dependencies/blaspp" --target install --parallel 16
    CXX=$(command -v CC) CXXFLAGS="-DLAPACK_FORTRAN_ADD_" \
        cmake -S "$audit/sources/lapackpp" -B "$audit/dependencies/lapackpp" \
        -DCMAKE_CXX_STANDARD=20 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
        -DCMAKE_INSTALL_PREFIX="$prefix/lapackpp-2024.05.31"
    cmake --build "$audit/dependencies/lapackpp" --target install --parallel 16

    if [[ ! -f "$prefix/venvs/warpx-gpu/bin/activate" ]]; then
        python3 -m venv "$prefix/venvs/warpx-gpu"
    fi
    source "$prefix/venvs/warpx-gpu/bin/activate"
    python3 -m pip install --upgrade pip build packaging wheel 'setuptools[core]' cython
    python3 -m pip install --upgrade numpy pandas scipy matplotlib yt openpmd-api \
        openpmd-viewer cupy-cuda13x pypdf sphinx sphinx-design sphinx-copybutton ruff
    MPICC='cc -target-accel=nvidia80 -shared' \
        python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation \
        --no-binary mpi4py
    python3 -m pip install --upgrade -r "$root/requirements.txt"
    python3 -m pip freeze > "$audit/python-packages.txt"
    module list > "$audit/modules.txt" 2>&1
    exit 0
fi

source "$prefix/venvs/warpx-gpu/bin/activate"
if [[ "$phase" == test-data ]]; then
    mkdir -p "$audit/test-data"
    for entry in \
        'warpx-data https://github.com/BLAST-WarpX/warpx-data.git' \
        'openPMD-example-datasets https://github.com/openPMD/openPMD-example-datasets.git'; do
        read -r name repository <<< "$entry"
        if [[ ! -d "$audit/test-data/$name/.git" ]]; then
            git clone --depth 1 "$repository" "$audit/test-data/$name"
        fi
        git -C "$audit/test-data/$name" rev-parse HEAD > "$audit/$name-revision.txt"
        # Stock inputs resolve these paths relative to their CTest directory.
        for destination in "$root/../$name" "$root/build/$name"; do
            if [[ ! -e "$destination" ]]; then
                ln -s "$audit/test-data/$name" "$destination"
            fi
        done
    done
    exit 0
fi
checkout="$root"
if [[ "$phase" == stock ]]; then
    checkout="$root/build/stock-warpx"
    if [[ ! -d "$checkout" ]]; then
        git worktree add --detach "$checkout" "$stock_revision"
    elif [[ "$(git -C "$checkout" rev-parse HEAD)" != "$stock_revision" ]]; then
        git -C "$checkout" diff --exit-code
        git -C "$checkout" diff --cached --exit-code
        git -C "$checkout" switch --detach "$stock_revision"
    fi
fi
build="$checkout/$build_name"

if [[ "$phase" == float ]]; then
    # Native float particles are needed to test probability narrowing and
    # recoil arithmetic; casting tables in a double build is insufficient.
    float_amrex="$audit/amrex-float-$build_name"
    cmake -S "$root/$build_name/_deps/fetchedamrex-src" -B "$float_amrex" \
        -DCMAKE_BUILD_TYPE=Release -DAMReX_SPACEDIM=1 -DAMReX_GPU_BACKEND=CUDA \
        -DAMReX_BUILD_SHARED_LIBS=ON \
        -DAMReX_CUDA_ARCH=80 -DAMReX_MPI=ON -DAMReX_OMP=OFF \
        -DAMReX_PRECISION=DOUBLE -DAMReX_PARTICLES=ON \
        -DAMReX_PARTICLES_PRECISION=SINGLE -DAMReX_FORTRAN=OFF \
        -DAMReX_EB=OFF -DAMReX_AMRLEVEL=OFF -DAMReX_LINEAR_SOLVERS=OFF
    cmake --build "$float_amrex" --parallel 16
    build="$root/${build_name}_float"
    physics_amrex="$float_amrex/lib/cmake/AMReX"
    phase=physics
else
    physics_amrex="$build/_deps/fetchedamrex-build/lib/cmake/AMReX"
fi

if [[ "$phase" == configure || "$phase" == stock ]]; then
    cmake -S "$checkout" -B "$build" \
        -DCMAKE_BUILD_TYPE=Release -DWarpX_COMPUTE=CUDA -DWarpX_DIMS='1;RZ;3' \
        -DWarpX_FFT=ON -DWarpX_APP=ON -DWarpX_PYTHON=ON -DWarpX_MPI=ON \
        -DWarpX_OPENPMD=ON -DWarpX_QED=OFF -DWarpX_EB=OFF \
        -DWarpX_PRECISION=DOUBLE -DWarpX_PARTICLE_PRECISION=DOUBLE \
        -DWarpX_TESTING=ON -DWarpX_TEST_CLEANUP=OFF \
        -DWarpX_PYTHON_IPO=OFF -DpyAMReX_IPO=OFF \
        -DMPIEXEC_EXECUTABLE="$root/Tools/Algorithms/PrescribedFluids/perlmutter_mpiexec.sh" \
        -DMPIEXEC_NUMPROC_FLAG=-n -DMPIEXEC_PREFLAGS= -DMPIEXEC_POSTFLAGS=
fi

if [[ "$phase" == build || "$phase" == stock ]]; then
    cmake --build "$build" --parallel 16
    # Use the built Python package directly; installing two revisions into the
    # same virtual environment would obscure which library a comparison loaded.
    cp "$build/CMakeCache.txt" "$audit/$(basename "$checkout")-CMakeCache.txt"
    git -C "$checkout" rev-parse HEAD > "$audit/$(basename "$checkout")-revision.txt"
fi

if [[ "$phase" == physics ]]; then
    for entry in 'PrescribedFluids fluid-tests' 'ProtonImpactIonization pjg' \
        'BackgroundMCC mcc-physics'; do
        read -r source destination <<< "$entry"
        cmake -S "$root/Tools/Algorithms/$source" -B "$build/$destination" \
            -DAMReX_DIR="$physics_amrex" \
            -DCMAKE_BUILD_TYPE=Release
        cmake --build "$build/$destination" --parallel 8
    done
fi
