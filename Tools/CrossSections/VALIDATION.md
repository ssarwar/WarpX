# Independent-angle rotational validation

The 2026-09-25–26 revision removes angular-accessibility conditioning and recoil
threshold shifts from rotational rates and outcome sampling. The internal
level spacing remains a strict threshold. Sampling uses one alias lookup;
only the cumulative comparison mode stores a cumulative table. A subthreshold
excitation caused by grid rounding becomes an unchanged event without
resampling or increasing another rotational channel's rate.

The angle is drawn independently from the existing elastic DCS. Recoil remains
exact outside the small energy-only continuation band documented in
`Docs/source/theory/mcc_iaa_sources.rst`. Inside that band, the electron retains
`E-loss`; the virtual neutral receives the momentum difference, and its recoil
energy is neglected in the electron update. Tests bound the resulting energy
defect by `2*m_e/M*E`. This selected approximation does not enforce exact
four-momentum conservation in that band or differential detailed balance.

## CPU physics

AppleClang 21, macOS arm64, NOACC; double fields and both double/single particle
precision. The double build includes 1D, 3D and RZ; the single build is 1D.

- 23/23 selected CTest entries passed in each particle precision. They cover
  rotation, anisotropic elastic DCS, RBEQ profiles, relativistic/nonrelativistic
  MCC, the ordinary selector, and invalid inputs.
- After the final recoil error handling and loader changes, 13/13 targeted
  rotation/elastic/invalid-input entries passed again in each precision.
- An independent numerical energy-budget solve checks recoil at retained
  angles, both gases, thresholds, gains, zero energy and energies through
  1 GeV. Separate checks cover the nominal threshold and the bounded
  continuation, without weakening the conservation tolerance elsewhere.
- Deterministic two-dimensional sampling quadrature checks alias/cumulative
  probabilities and verifies that subthreshold excitation does not increase
  any other rotational rate. Canonical loss labels retain double precision.
- Independent Legendre quadrature, Boltzmann tails, state-sum moments and
  heavy-target integral detailed balance are checked against the source-rate
  generator. Midpoint rate and transfer-moment errors remain below 0.2% outside
  separately tested float32 threshold bands.
- The actual supplied elmolcs N2/O2 DCS files passed 16 elastic/excitation cases
  from 0.1 eV through the screened-Rutherford continuation at 1 GeV, with
  independent angular-moment and recoil analysis.
- A CUDA compiler issue exposed a host/device lambda-capture mismatch in the
  existing proton kernel. The capture-only fix passed 21 relevant CPU
  proton, RZ source and rotational-threshold CTest entries.
- A further CUDA-only endpoint discrepancy was repaired within the existing
  input-precision bound, always moving the secondary energy downward. The
  expanded neighboring-endpoint and finite-excess tests passed 20 related CPU
  CTest entries and both float/double device recoil checks.

Deterministic relative-Maxwellian quadrature of the synthetic verification
kernels gives `100*abs(cooling-heating)/(cooling+heating)` below the 0.1% limit:

| Target | 100 K | 300 K | 1000 K |
| --- | ---: | ---: | ---: |
| N2 | 0.00546% | 0.00255% | 0.00130% |
| O2 | 0.00413% | 0.00205% | 0.00122% |

These are numerical verification kernels, not validated production molecular
cross sections or proof of exact differential detailed balance.

## CPU performance

Five seeds, alternating alias/cumulative order, 65,536 particles in each of
14 cases, eight complete timesteps, with a device-fencing scalar reduction.

| Version and sampler | Mean time for eight steps | Particle-steps/s |
| --- | ---: | ---: |
| Previous angular filter, alias | 0.855247 s | 8.582e+06 |
| Independent outcomes, alias | 0.850738 s | 8.628e+06 |
| Independent outcomes, cumulative | 0.842466 s | 8.713e+06 |

The approximately 0.5% before/after alias difference does not establish a
meaningful CPU throughput improvement. The largest alias table decreased from
2.79 to 1.87 MiB in double particle precision, and from 2.32 to 1.40 MiB in
single particle precision. Current maximum table sizes are 1,958,136 and
1,467,332 bytes per host/device copy, respectively.

The driver records estimator variance times computational cost. Five-seed
variance estimates are too noisy to establish a noise advantage. Both samplers
retain the same discrete energy diffusion and pass the same moment checks.

## Perlmutter CUDA validation

Validated on NVIDIA A100-SXM4-40GB devices, driver 580.178.04, CUDA 13.2.78,
GNU 14.3.0 and CMake 3.28.3. The isolated 1D build uses double particle/field
precision, CUDA architecture 80, Python bindings, and no MPI/FFT/OpenPMD/QED.
The baseline is 7dd42a38b with the same capture-only compiler fix as the updated
branch. The rotational implementation is 8271ae895 with its offline-input
workflow in 59ac08393; the later proton endpoint correction does not change
rotational sampling.

- Both baseline and updated rotation/elastic/invalid-input CUDA suites passed
  17/17 CTest entries. An additional device recheck passed 3/3 threshold,
  proton-angle and proton-recoil tests, including the endpoint correction.
- Compute Sanitizer memcheck found zero errors in the rotational threshold
  kernel and a complete MCC run using the prepared verification files.
- Compute Sanitizer racecheck reported zero hazards, errors or warnings in the
  rotational threshold/sampler kernel. This is not a claim of coverage of every
  WarpX kernel or backend.

The complete timestep benchmark uses 524,288 particles in each of 14 cases,
32 steps, five measured seeds and alternating sampler order. The updated run
adds one separate warmup run per sampler, followed by fresh measured ensembles.

| Version and sampler | Mean time for 32 steps | Run-to-run standard deviation | Particle-steps/s |
| --- | ---: | ---: | ---: |
| Previous angular filter, alias | 0.270482 s | 0.001642 s | 8.684e+08 |
| Independent outcomes, alias, warmed | 0.264344 s | 0.002138 s | 8.885e+08 |
| Independent outcomes, cumulative, warmed | 0.265267 s | 0.001239 s | 8.855e+08 |

The baseline followed its CUDA physics tests. The updated timing repeat used
a different node with the same GPU model and driver, in a one-GPU shared
reservation. These conditions do not establish a strong causal claim for the
approximately 2.3% before/after timing difference. Alias and cumulative differ
by only 0.35% in the warmed complete-timestep measurement.

The initial updated run had one first-run outlier (0.461116 s), giving a raw
five-run alias mean of 0.305648 s. That result was not silently discarded.
Instead, the driver gained explicit warmup runs and timing median/standard
deviation reporting, and the entire measurement was repeated. All five
measured runs of the repeat are retained.

Isolated lookup measurements used 1,048,576 concurrent samples per case. Median
aggregate costs were 0.1355 ns/sample for the previous alias sampler,
0.1415 ns/sample for the updated alias sampler and 0.1399 ns/sample for the
updated cumulative sampler. These are throughput-derived costs across the
whole GPU, not single-thread latencies. They do not demonstrate a lookup
speedup; reduced memory and removal of angular conditioning are the established
improvements. Five-seed estimator variance measurements do not establish a
noise advantage between samplers with identical physics.

The successful final device/sanitizer/warmed-benchmark reservation was Slurm
job 58889195. `perlmutter_rotation.sbatch` reproduces the main device suite and
benchmarks; provide the GPU account and Python environment when submitting.
Prepared inputs can be selected with `WARPX_ROTATION_DATA`. HIP and SYCL device
runs, multi-GPU scaling and single-particle-precision CUDA builds were not
performed; single particle precision was validated on the CPU.

## Offline data and source scope

V3 synthetic verification bundles and their integral/DCS inputs are versioned
in `warpx-data/IAA/MCC_cross_sections/IAA/rotation/verification`, with hashes
in `manifest.json`. They match the files used by the CPU tests. The PICMI
input and benchmark only read prepared files. CTest generates its synthetic
fixtures in a separate setup step; production runs do not invoke data tools.

The physical source audits retain their explicit zero-temperature reference
assumption. Additional 300 K audits record negative low-energy unchanged
integral residuals in the current candidate extrapolations. These are not
clipped. Low-energy thermal continuations, omitted-rank moment convergence and
beam-energy coverage remain prerequisites for production rotational bundles.
No physical production bundle or complete downstream cross-section export is
claimed by these verification results. No pull request was opened.
