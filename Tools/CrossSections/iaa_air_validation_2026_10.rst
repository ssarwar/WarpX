IAA electron--air validation, October 2026
=========================================

These checks use elmolcs archive SHA256
``b864381086120fff37eb8c658dfcb0f150b80969cf9add626a38d89bab0a3c38``
and the supplied July 2026 MuroranIT attachment export. The implementation is
based on WarpX ``8c0268b56`` and the elastic/rotation data at ``61a3947``.
All synthetic fixtures, dense representations and diagnostic output are kept
in WarpX build directories.

Source and numerical checks
---------------------------

The maximum observed off-grid relative errors, with the declared 1e-6 peak
floor, were 0.0344% for N2 and 0.0382% for O2, below the 0.05% target.
The production rate files contain 53,786 N2 points and 144,378 O2 points across
all non-RBEQ channels and alternatives. Source knots are retained; redundant
knots introduced by other channels are omitted from individual files.

The independent archive-total comparison over 5--1000 eV gives differences
up to 4.93% for N2 and 2.00% for O2. This comparison includes the completed
channels and the selected thesis RBEQ total. These physical-model differences
are reported, not removed by changing the elastic residual. Excitation rate
and first/second loss-moment interpolation errors remain below 0.03% on the
independent audit grid.

All 26 angular profiles pass independent Gauss-quadrature moment and CDF
checks. Analytical checks cover isotropic and d-pi distributions, normalized
mixtures, finite-step angular relaxation, and tiny relativistic deflections.
Native-grid tests check preserved zero prefixes and tails, unresolved-channel
rejection, raw m5 values below binary32 range, and a channel with probability
1e-10 next to a large channel.

The actual MCC loop is tested against finite-step binomial cooling ladders,
their Poisson limit, exponential attachment survival, competing reaction
branching and third-body density scaling. Full-data cases cover N2, O2 and
dry-air cooling, an imposed-field swarm, and two homogeneous 800 MeV proton
pulses with coupled electron collisions. This pulse test is a collision
regression, not validation of the experiment's spatial beam evolution.

CPU builds with both double- and single-precision particles pass these tests.
Selected existing attachment, selector, RBEQ, many-channel excitation and
rotational-equilibrium regressions also pass. No existing assertion tolerances
were widened. The native reader was corrected to preserve physical zero
intervals and threshold constraints under binary32 energy rounding.

Performance measurement
-----------------------

The recorded workload used macOS arm64, AppleClang 21, Release/NOACC,
double-precision particles, 8,192 initial electrons, 16 cells, dry air,
64 timed steps and three independent seeded ensembles. Each representation
had a separate warmup run; each trial excluded four warmup timesteps.
The dense reference inserts redundant knots into the same interpolants.

============================== ============ ============
Measurement                    Compact      Dense
============================== ============ ============
Mean timed steps (seconds)      0.4024       0.4356
Timing standard deviation      0.0048       0.0027
Selector arrays (bytes/copy)    13,663,752   66,317,592
Angular arrays (bytes/copy)     25,921,480   99,717,592
============================== ============ ============

Compact tables reduced timed-step cost by approximately 7.6% in this CPU
workload. The optimized selector remained enabled. Energy observables agree
within ensemble uncertainty. Three ensembles do not establish a reliable
variance advantage; variance-times-cost is recorded rather than claimed to
improve. Memory figures exclude the unchanged reciprocal rotational families.

CUDA, HIP and SYCL compilers/devices were unavailable on this host. The tests
use host/device kernels and the benchmark supports GPU particle arrays, but
accelerator correctness and performance remain to be measured on those backends.

Physical-source limitations
---------------------------

Endpoint completion can contribute up to 48.85% of an accepted N2 angular
profile and 46.26% of an O2 profile. Alternative endpoint completion changes
individual angular moments by up to about 0.14 and 0.084, respectively.
These physical-model uncertainties exceed the numerical sampling errors.

The MuroranIT m5 peak is preserved exactly. It differs by 0.0177724% from the
nominal normalization in the paper; this is a source discrepancy, not a unit
conversion adjustment. Two possible decimal transcription outliers in the
archived Shyn--Sweeney O2 b-state DCS are preserved and flagged for confirmation
against the original paper. The author's current public source retains the
same values (blob ``b3500bfeb03ab802ecdcd9219244a46b162ccfa1``).

Reproduction
------------

Configure with ``WarpX_IAA_AIR_TEST_DATA`` set to the external
``warpx-data/MCC_cross_sections`` directory, then run CTest with ``-R iaa``.
Use ``WarpX_PARTICLE_PRECISION=SINGLE`` in a separate worktree-local build
directory for the particle-precision checks. Set ``WarpX_COMPUTE`` to CUDA,
HIP or SYCL on the corresponding accelerator host.

.. code-block:: bash

   python Tools/CrossSections/benchmark_iaa_air.py \
       --data-dir /path/to/warpx-data/MCC_cross_sections \
       --output build/iaa-benchmark --particles 8192 --steps 64 --repeats 3
   python Tools/CrossSections/audit_iaa_air.py \
       --data-dir /path/to/warpx-data/MCC_cross_sections \
       --archive /path/to/elmolcs-main.tar.gz --output build/iaa-audit.json
