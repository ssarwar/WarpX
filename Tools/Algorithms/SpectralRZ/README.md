# RZ SYCL transform check

This standalone check compares the production complex FFT with an independent
long-double direct DFT and an unnormalized round trip. It covers odd and even
axial lengths, interleaved radial sequences, three complex azimuthal modes and
changes of AMReX queue. It needs a SYCL device supported by oneMKL.

Use the AMReX package from a SYCL build of WarpX and the same oneAPI environment:

```bash
cmake -S Tools/Algorithms/SpectralRZ -B build/spectral-rz-tests \
    -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx \
    -DCMAKE_CXX_FLAGS=-fsycl -DAMReX_DIR=/path/to/AMReX/lib/cmake/AMReX \
    -DMKL_SYCL_THREADING=sequential -DMKL_THREADING=sequential
cmake --build build/spectral-rz-tests -j 4
ctest --test-dir build/spectral-rz-tests --output-on-failure
```

The existing RZ PSATD Langmuir tests provide end-to-end physics checks of the
Hankel and Fourier transforms. The prescribed-fluid beam and source tests add
current correction, self-fields, collision products and checkpoint coverage.
