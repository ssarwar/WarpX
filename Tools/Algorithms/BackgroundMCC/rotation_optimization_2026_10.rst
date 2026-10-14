October 2026 rotational MCC optimization
========================================

The changes reduce the cost of the existing reciprocal rotational model.
They preserve collision acceptance, rates, energy/angle grids, discrete
outcomes, recoil, and the interleaving of coupled collision substeps. Alias
sampling remains the default. No warpx-data file is changed.

Scope and provenance
---------------------

* WarpX branch: ``codex/beam-air-development-sync-2026-10`` in ``ssarwar/WarpX``.
  The comparison baseline is ``305e203ff6befb524e7e5f8277089113da5af0f0``.
  The two latest integration commits, including the development update
  ``b725147c7``, were incorporated before the final build. The fork's
  ``development`` branch was not modified.
* Data: ``ssarwar/warpx-data``, ``codex/elmolcs-elastic-rotation``, commit
  ``61a3947b6ebf365e6b78b63afdc30306791e3898``. Both production 300 K bundles
  and their ordinary elastic tables are used unchanged.
* AMReX: ``b552c85c432defab450688581a1c87aead7b461d``.
* GPU timing, profiling and exact before/after comparisons: Perlmutter
  NVIDIA A100-SXM4-80GB, CUDA 13.2, GCC 13.2.1, driver 580.178.04. Broad
  regression batches also use A100-SXM4-40GB nodes. MPI GPU support and IPC
  remain enabled in performance and regression runs.
* Local: Apple M3 Pro, Clang release build, one OpenMP thread per test,
  native double/double and separate single/single configurations.

``rotation_optimization_2026_10.json`` records the numerical results, source
hashes and test outcomes. Temporary binaries, profiles, native cell dumps,
Slurm scripts and detailed logs are under
``build/validation/rotation-optimization`` in the local worktree and the
dedicated Perlmutter checkout
``/pscratch/sd/s/ssarwar/warpx-beam-air-sync-20261006``.

What costs time
---------------

The rotational path adds joint-distribution searches, discrete outcome
selection and signed internal-energy recoil to ordinary elastic scattering.
An exploratory same-layout comparison measured 13--28% greater complete
MCC time than the ordinary elastic control for thermal, 2.47 eV and mixed
energies. That comparison changes the physical model and is used only to
characterize added computational cost.

The baseline Nsight Systems trace contains 512 MCC launches per PIC step:
four boxes, two gases and 64 subcycles. Summed MCC kernel time is 26.38 ms
for thermal electrons, 37.04 ms at 2.47 eV and 42.21 ms for mixed energies.
The corresponding steady-step trace contains 3186 stream synchronizations
and 271 asynchronous copies. Synchronization API duration includes waiting
for GPU work; it is not an independently additive overhead estimate.

AMReX starts CUDA profiling during initialization. The steady region was
therefore selected after the last explicit ``cudaProfilerStart`` call,
excluding initialization, state restoration and warmup. Full-trace API
percentages would give a misleading account of subcycle overhead.
Nsight Compute could not acquire the driver's profiling resource. There
are no measured cache-hit, bandwidth or occupancy claims in this audit.

Implemented changes
--------------------

1. **Resolve alias alternates once.** The input/preparation alias retains
   its alternate-column format. After validation, initialization replaces
   that column with its outcome identifier in place. Each event then loads
   one alias entry and selects either of its two outcome identifiers.
   The eight-byte entry layout, cutoff and random draw are unchanged; one
   dependent table read is removed without another table-sized allocation.
2. **Narrow the exact energy search.** A small index uses the binary exponent
   and 64 mantissa bins. Its dyadic boundaries are exactly representable.
   Each bin supplies conservative lower/upper intervals for the original
   binary search. Every physical knot and interpolation operation remains.
   Small grids retain the original full search.
3. **Reuse pure-scattering error storage.** A successful call leaves the
   flag zero, so subsequent calls need no allocation or zero upload. The
   host still reads and checks it after every call. Error storage is cleared
   before raising a host assertion. The RNG launch already waits for each
   tile, so the iterator's additional exit fences are omitted; its entry
   synchronization is retained. The iterator configuration is thread-local
   in OpenMP execution.
4. **Remove unreachable kernel branches.** Four compile-time variants
   distinguish reciprocal/ordinary models and scalar/parser backgrounds.
   A reciprocal model cannot coexist with the legacy thermal/spectator
   family in one MCC object; construction enforces that precondition.
   Scalar backgrounds avoid the parser path. The active sampling calls,
   acceptance law and recoil equations remain unchanged.
5. **Accelerate the cumulative reference.** Optional 64-bin per-cell
   bounds narrow the original strict upper-bound search. The stored double
   CDF and its outcome IDs are unchanged. This improves the reference
   sampler, but alias remains faster and smaller.

The empty-particle shortcut applies only to pure scattering. An early
prototype placed it outside that condition; review caught that an empty
MPI rank would then skip required immobile-fluid commits. The final code
keeps those collective operations and has two new empty-rank regressions.
Performance prototypes all had nonempty incident populations; final
three-seed measurements use the corrected implementation.

Why these transformations preserve the model
--------------------------------------------

For an alias entry, both implementations choose exactly the same palette
identifier for a given column and fractional draw. Resolving the alternate
at startup changes the address path only. Primary identifiers are never
modified, so the in-place conversion cannot change a later lookup.

For an energy bin with boundaries :math:`B_l \le E < B_h`, the interval at
:math:`B_l` and one past the interval at :math:`B_h` bracket the same full-grid
answer. Binary-exponent/mantissa binning preserves these boundaries. The
final interpolation still uses the original physical energies and flags,
including square-root threshold intervals.

For a cumulative draw :math:`u \in [k/64,(k+1)/64)`, the first CDF entry
strictly above the lower boundary cannot follow the answer for :math:`u`;
the first entry strictly above the upper boundary cannot precede it.
The final search uses the original comparison, including flat CDF regions.
The existing upper-end draw clamp prevents indexing beyond the sentinel.

The collision-handler schedule is unchanged. Fusing all substeps of N2
before all substeps of O2 would reintroduce the previously corrected
full-PIC-step splitting error. A future fusion must preserve interleaving,
new-particle participation and product-fluid commits. The present changes
add neither a particle pass nor an event-time neutral-state sum.

Measured throughput
--------------------

Both gases are resident on one GPU. The grid has 4096 cells in four boxes.
Each PIC step is 6.4e-13 s and each gas takes 64 collision substeps. Thermal
electrons and neutral translation/rotation are at 300 K. Mixed electron
energies are log-uniform from 0.002 eV to 3 MeV within warps. The quoted
figures are the median of three seed medians, with seeds 19, 2026 and 800.
Each process has four warmup steps followed by 20 synchronized blocks of
eight steps. Execution order alternates between seeds, and every variant
of a workload uses the same physical GPU.

The MCC workloads use the normal PIC-step driver with electron push,
gather and deposition disabled; field-driver work and bookkeeping remain
in their timings. They are not isolated kernel measurements.
Positions and momenta are restored outside each timed block. This controls
changes in the incident distribution as the simulation runs. It does not
require identical seeded particle trajectories across separate GPU runs.
The full-PIC workload includes push, gather, deposition and field advance;
electron density is small enough that self-fields remain negligible.

.. list-table:: Median milliseconds per PIC step
   :header-rows: 1
   :widths: 45 18 18 19

   * - Workload
     - Baseline
     - Optimized
     - Time reduction
   * - 1048576 particles, thermal MCC
     - 35.005
     - 32.038
     - 8.48%
   * - 32768 particles, thermal MCC
     - 14.879
     - 13.080
     - 12.09%
   * - 1048576 particles, mixed energies, full PIC
     - 51.860
     - 47.843
     - 7.75%
   * - Thermal MCC with inactive ionization product channels
     - 43.751
     - 42.158
     - 3.64%
   * - Mixed-energy MCC with parser-defined backgrounds
     - 55.264
     - 43.126
     - 21.96%
   * - Ordinary elastic MCC, mixed energies
     - 39.446
     - 38.511
     - 2.37%

The ionization control has zero cross sections and retains a fixed
population while exercising the product pipeline. The parser control uses
constant-valued expressions. These results do not predict the wall time of
an evolving, fully reactive beam-air simulation. Three seeds also do not
establish a statistical-variance advantage.

An earlier ablation measured each addition independently using 40 timed
blocks, with reversed variant order in a second process. At one million
particles, resolving aliases alone improved resonance/mixed-energy time by
about 1%, with no measurable thermal gain. Adding the energy index brought
the improvement to about 3.4%. Reusing error storage and removing redundant
exit fences accounted for most of the remaining pure-scattering gain.
Specialization added about another 1% with scalar backgrounds. Its larger
parser-path benefit is confirmed by the separate three-seed control.

The independent sampler comparison additionally covers 2.47 eV, 50 eV and
2.5 MeV. Across those constant-background workloads, final alias MCC time
fell by approximately 8--9%. Indexed cumulative sampling remains slower
than optimized alias in every measured case. Full timing summaries for
all intermediate variants are retained in the JSON record.

A separate CPU control compares the saved baseline and optimized native
``benchmark_reciprocal_mcc`` executables with one OpenMP thread, 262144
particles, the same four-box grid and three alternating seeds. Each run
times 512 interleaved N2/O2 operator pairs at 1e-14 s after 256 warmup pairs.
Median elapsed time falls from 7.473 to 6.541 s for thermal electrons and
8.653 to 8.074 s for mixed energies. Every paired run improves. The baseline
thermal timings span 6.846--7.857 s, so these workstation measurements should
be treated as a CPU regression check, not a precise prediction of a PIC
application's speedup. Initialization and field evolution are excluded.

Memory and initialization
-------------------------

Combined table storage is 751225792 bytes (716.425 MiB) for aliases and
954196798 bytes (909.993 MiB) for indexed cumulative sampling. The exact
energy index adds 30728 bytes across both gases. The optional cumulative
index adds 29085290 bytes (27.738 MiB); alias sampling pays none of that
allocation. The combined 1 GiB table budget remains enforced.

Exact alias-cell sharing was also screened as a possible larger memory
optimization. Equal-length cells with different decoded moments cannot be
identical alias arrays. Grouping by entry count and all eight decoded
moments leaves at most 11088 potentially shareable bytes in N2 and 3720 in
O2: only 0.00213% of alias storage. Equal moments do not prove identical
arrays, so this is an upper bound. Adding another cell-mapping structure
would not provide a worthwhile storage reduction for these tables.

In the repeated sampler campaign, warm-filesystem median initialization
changed from 6.95 to 7.07 s for alias and from 7.67 to 7.87 s for cumulative
tables. Median process peak host memory was about 0.775 GiB for alias and
1.287 GiB for cumulative. Device allocator granularity and complete
simulation allocations are separate from the table byte count.

Compiler resource output gives 128 registers and 40 stack bytes for the
scalar reciprocal kernel, and 80 registers and 328 stack bytes for its
parser variant, versus 97 registers and 208 stack bytes in the baseline.
The compiler reported no additional local allocation in these entries.
These are static compilation figures, not hardware-counter measurements;
lower register count alone is not a throughput result.

Correctness and regression evidence
------------------------------------

* CPU and CUDA before/after comparisons are byte-identical for all decoded
  cells and seeded 28-energy sample files, for N2 and O2 with each sampler
  compared to its own baseline. Each energy uses 32768 samples. Decoded
  cells include probability and signed-transfer moment checks, preserving
  rare tails as well as normalization.
* Indexed interpolation agrees exactly with full search at 50439 N2 and
  49539 O2 queries. These cover all physical energy knots and dyadic
  boundaries, their adjacent floating-point values, and a 32768-point
  logarithmic sweep per gas. No interpolation tolerance is introduced.
* Indexed cumulative sampling matches full search at six edge queries in
  each of 116626 N2 and 107107 O2 cells. Paired event checks also compare
  original and indexed searches with the same random uniforms.
* Existing independent equations check signed recoil and laboratory
  four-momentum conservation, threshold continuation, and exact screened
  Rutherford transport moments. Thermal, ionization-competition,
  attachment, subcycling and restart cases exercise the complete operator.
* A delayed invalid-density case succeeds with zero density on its first
  call and must detect a negative density on its next call. This guards
  reuse of the runtime error flag. Two MPI regressions exercise ionization
  and attachment into immobile products with an empty incident rank.
* The full CPU suite has 948 stages: 936 pass and the same 12 previously
  recorded failures remain. Their scientific output matches the prior
  run exactly after removing yt timestamps. All 51 selected all-single
  cases pass, as do 3 double and 3 single analytical MCC tests, 6 radial
  stages, 5 proton-kernel tests and 4 fluid-kernel tests. The 380 checksum
  stages retain the prior 200 passes and 180 differences; no benchmark
  JSON was changed.
* The full CUDA suite has 945 stages: 926 pass and 19 fail, with no skipped
  stages. Eighteen failures occur in the earlier failure set. The earlier
  Vay comparison passes in this run; that does not establish a numerical
  improvement in an unchanged test. The additionally observed CKC restart
  comparison is discussed below. All newly added regressions pass on CUDA.
* The expanded all-single CUDA selection has 74 stages: 73 pass and the
  coarse helium-discharge analysis fails at 7.05% RMS versus 6.50%. Earlier
  single-precision campaigns explicitly excluded that application test;
  it is additional coverage, not a reproduced before/after comparison.
  Its standard input deliberately uses coarse spatial, particle and time
  resolution. No tolerance or reference profile is changed.
* All 3 double and 3 single CUDA analytical MCC tests pass, as do all
  5 proton and 4 fluid analytical GPU tests. Compute
  Sanitizer reports zero errors for both gases with alias and cumulative
  sampling, and for repeated complete rotational-MCC calls. The CUDA
  checksum campaign has 338 passes and 42 differences out of 380 stages;
  differences are retained and ignored under the repository policy.

The original ``test_3d_acceleration_restart.analysis`` failure has relative
``jy`` error 1.1181680e-12 against its unchanged 1e-12 limit. Its input has
no MCC collisions. Five independent repeat pairs use the saved unchanged
development build at ``63fa5a0c1`` and the optimized build, the same inputs,
and the same identity-matched restart analysis. All five controls on each
build also exceed the bound, in ``plasma_p.particle_momentum_y``: reported
errors range from 1.169e-12 to 1.579e-12 on development and 1.026e-12 to
1.410e-12 on the optimized build. This reproduces a precision-sensitive
restart comparison without these optimizations. The original failure is
retained; particle identities are checked exactly and the analysis threshold
is never relaxed.

The single-precision discharge control repeats the coarse input and then
halves only its timestep, doubling the step count and diagnostic averaging
step count to preserve physical duration and averaging time. Both use two
MPI ranks, the same grid and particle count, unchanged WarpX/default and
NumPy seed settings, and byte-identical helium data. The coarse repeat
gives 6.39% RMS and the half-timestep case gives 5.67%; both pass the
unchanged 6.50% limit. These observations demonstrate sensitivity to the
realization and resolution of this deliberately coarse application test;
they are not a replacement for a broader convergence study. Its original
7.05% failure remains in the reported suite result.

CPU/CUDA exactness checks establish preservation of the supplied model;
they do not eliminate its documented experimental uncertainties or the
loss of meV transfers in single-precision MeV particle momenta. Physical
source fitting was not repeated because the model and data are unchanged.
HIP and SYCL toolchains/hardware were unavailable for this campaign.

Reproduction
-------------

Select the baseline or optimized Python package through ``PYTHONPATH``
and run the same driver and seed sequence on the same GPU. For example,
from a build-directory working folder, with absolute paths supplied:

.. code-block:: bash

   mpiexec -n 1 python /path/to/WarpX/Tools/Algorithms/BackgroundMCC/benchmark_rotational_subcycling.py \
       --data /path/to/warpx-data/MCC_cross_sections --output thermal.json \
       --spectrum thermal --particles 1048576 --subcycles 64 \
       --dt 6.4e-13 --warmup 4 --blocks 20 --steps 8 --seed 2026

Use ``--pic --spectrum broad`` for the complete PIC control,
``--inactive-ionization`` for product-buffer overhead,
``--parser-background`` for parser-defined gas inputs,
``--ordinary-elastic`` for the separate physical model, and ``--cumulative``
for the alternative sampler. Startup, timed blocks, configuration and peak
host memory are recorded separately. Linux peak RSS units are KiB; macOS
reports bytes.

Build and run the production C++ sampler with ``lookup_check=1`` and
``cell_output=<path>``; repeat with ``cumulative=1``. Compare each result
to the same sampler built from the baseline. Run the normal CTest suite
without changing analysis thresholds. The Perlmutter regression is divided
into disjoint batches that keep fixtures, dependencies and working
directories together; dry-run test-name sets are checked against the full
945-stage inventory. Each batch has separate CTest driver/log storage and
the original executable commands, inputs and output directories.
