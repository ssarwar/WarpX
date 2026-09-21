#!/usr/bin/env bash
# Build the coupled-physics audit with WarpX's documented Perlmutter recipes.
# Run from the validation checkout on Perlmutter; see REAUDIT.md.
set -eo pipefail

root=$(git rev-parse --show-toplevel)
audit="$root/build/reaudit-2026-09-21"
prefix="$audit/software"
phase=${1:-configure}
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
checkout="$root"
if [[ "$phase" == stock ]]; then
    checkout="$root/build/stock-warpx"
    if [[ ! -d "$checkout" ]]; then
        git worktree add --detach "$checkout" cb5672fae0b9099d2c7631e93ffff3404a553821
    fi
fi
build="$checkout/build_pm_gpu_py"

if [[ "$phase" == configure || "$phase" == stock ]]; then
    cmake -S "$checkout" -B "$build" \
        -DCMAKE_BUILD_TYPE=Release -DWarpX_COMPUTE=CUDA -DWarpX_DIMS='1;RZ;3' \
        -DWarpX_FFT=ON -DWarpX_APP=ON -DWarpX_PYTHON=ON -DWarpX_MPI=ON \
        -DWarpX_OPENPMD=ON -DWarpX_QED=OFF -DWarpX_EB=OFF \
        -DWarpX_PRECISION=DOUBLE -DWarpX_PARTICLE_PRECISION=DOUBLE \
        -DWarpX_TESTING=ON -DWarpX_TEST_CLEANUP=OFF \
        -DWarpX_PYTHON_IPO=OFF -DpyAMReX_IPO=OFF \
        -DMPIEXEC_EXECUTABLE="$(command -v srun)" -DMPIEXEC_NUMPROC_FLAG=-n \
        -DMPIEXEC_PREFLAGS='--cpu-bind=cores;--gpus-per-task=1'
fi

if [[ "$phase" == build || "$phase" == stock ]]; then
    cmake --build "$build" --parallel 16
    # Use the built Python package directly; installing two revisions into the
    # same virtual environment would obscure which library a comparison loaded.
    cp "$build/CMakeCache.txt" "$audit/$(basename "$checkout")-CMakeCache.txt"
    git -C "$checkout" rev-parse HEAD > "$audit/$(basename "$checkout")-revision.txt"
fi
