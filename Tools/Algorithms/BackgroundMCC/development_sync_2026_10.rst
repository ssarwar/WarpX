October 2026 development integration and MCC audit
==================================================

Branch and scope
----------------

The integration branch is ``codex/beam-air-development-sync-2026-10`` in
``ssarwar/WarpX``. It starts at the fork's ``development`` commit
``c904424bbe8f6871320ec96da72af40aaef3dc98`` and merges the complete feature
history through ``97faee560e9b8973f5270aacef0b58b4b33ad062``. The merge is
``0df76c697cc1d6777ef440d61223eb9b221d8efd``. Neither the fork's development
branch nor the original feature branch is changed.

After the fork's development branch advanced again, merge
``ac0125555f55a9c0d016d7a3b09b39d07d2ccea5`` incorporates both new commits:
``8293885bf`` (weekly dependencies) and ``63fa5a0c1`` (the ``Direction``
constructor recursion fix). This second merge has no textual conflicts.
The AMReX pin is now ``b552c85c432defab450688581a1c87aead7b461d`` and the
PICMI pin is ``a2fc467f3125d57ea0183562e69f414b84abe675``. Generated checksum
updates are taken directly from development, without local regeneration.

The PICMI Git revision is metadata distinct from WarpX's actual package
requirement, which remains ``picmistandard==0.34.0``. A direct smoke check
of that Git revision exposes an upstream ``GriddedLayout`` regression:
both the singular and plural parameter spellings raise
``KeyError: 'n_macroparticle_per_cell'`` in its deprecated-argument checker.
The checker inspects a parameter removed from the constructor signature.
Full regression runs use the released package required by ``Python/setup.py``;
its constructor passes. The broken optional Git revision is not patched
silently or substituted into the supported runtime.

The data audit uses ``ssarwar/warpx-data`` branch
``codex/elmolcs-elastic-rotation`` at ``61a3947``. All 94 tracked files
remain unchanged. In particular, the elementary rotational tables are not
added as extra processes on top of the inclusive reciprocal family.

Why the merge conflicted
------------------------

* ``ScatteringProcess.cpp``: development renamed ``TWOPRODUCT_REACTION`` to
  ``TWO_PRODUCT_REACTION`` and standardized the input spelling. The feature
  branch added attachment products and prefix-matched ionization/attachment
  channels in the same parser and product classification. The result keeps
  the new spelling and all feature processes.
* ``CollisionHandler.cpp``: development added ``start_step`` and the
  supercycle offset relative to it. The feature branch introduced physical
  source intervals and distinct explicit/implicit substep endpoints. The
  result keeps the offset and ``doCollisionsInInterval``.
* ``SemiImplicitEM.cpp``: development moved reconstruction of the implicit
  field guess into initialization and clarified the stored field times. The
  feature inserted immobile-charge preparation beside that code. Both the
  new field-update algorithm and ``PrepareImmobileCharge`` are retained.
* ``Docs/source/usage/parameters.rst``: three overlapping blocks describe
  the renamed reaction, the additional electron processes, and their energy
  and angular options. The merged text matches the combined implementation.
* One automatic merge also needed repair: ``BinaryCollisionUtils.cpp`` kept
  an obsolete ``TWOPRODUCT_REACTION`` reference in a feature-added condition.
  Git did not flag that line; leaving it would have prevented compilation.

Corrections beyond the textual merge
------------------------------------

Coupled subcycling
~~~~~~~~~~~~~~~~~~

Previously, the handler advanced all substeps of one collision object before
starting the next. Increasing each object's subcycle count could not remove
this full-PIC-step operator-splitting error. The handler now interleaves calls
by physical substep endpoint, retains input order for equal endpoints, and
supports different subcycle counts. Integer fraction comparisons avoid
rounding-dependent ordering. A reusable heap holds one entry per collision
object; the change adds no particle passes or collision kernel launches.

An independent regression uses hot electrons that can cool once and then
attach in a second collision object. For equal rates ``nu``, the exact
populations at ``x = nu*t`` are

.. math::

   H=e^{-x},\qquad C=xe^{-x},\qquad A=1-(1+x)e^{-x}.

At ``x=2``, the continuous survivor fraction is 40.6006%. The old ordering
produced 25.3296% in the recorded run; its split-operator limit is 25.2355%.
With 64 interleaved substeps, the recorded forward-order result is 39.9963%,
consistent with sampling error and the exact finite-substep prediction

.. math::

   H+C=e^{-x}\left[1+N(1-e^{-x/N})\right].

Both operator orders, unequal 64/96 subcycles, delayed starts, and explicit
and semi-implicit evolution are covered. Strongly coupled particle sources
also need adequate subcycling; fine electron subcycling does not resolve a
source that remains coarse in time. Collision subcycling alone does not
establish convergence of the electromagnetic/PIC timestep.

Start times and restart guards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The custom proton-impact PICMI interface now forwards ``start_step``.
Non-default start times are recorded in prescribed-source checkpoint
metadata. Default-zero metadata remain compatible with existing checkpoints.
The guard test rejects a changed source start time before advancing.

MCC probability and angular precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The majorant attempt probability is computed in double precision once per
collision call and shared across tiles. This removes the repeated reciprocal
``expm1`` evaluation per tile and avoids narrowing its host early-return
check. Ordinary MCC attempt and acceptance draws also use the existing
native double-uniform helper, so fine subcycling in all-single builds is
not restricted to a binary32 random grid. Ordinary rate tables retain their
configured particle precision.

IAA angular draws use the same double-uniform helper. The new portable
``test_elastic_angles.cpp`` verifies random-number resolution and compares
screened-Rutherford transport moments with their exact integral:

.. math::

   \langle1-\cos\theta\rangle
     =2\eta\left[(1+\eta)\log(1+1/\eta)-1\right].

Logarithmic quadrature in ``1-u`` resolves the large-angle tail. The test
covers N2/O2 screening radii at 10 keV, 2.5 MeV and 1 GeV in double and single
precision. At 1 GeV, binary32 angular uniforms cannot sample the backward
tail; numerical tabulation refinement alone cannot repair that limitation.

Validation fixture repairs
~~~~~~~~~~~~~~~~~~~~~~~~~~

* A stale negative test attempted to replace an obsolete header token that
  no longer existed. It now verifies the original identifier and explicitly
  writes an unsupported one. The rejection assertion is unchanged.
* A manufactured constant-field restart test divided the error of an exactly
  zero centroid by ``1e-30``. An MPI reduction difference of about
  ``1.1e-18 m`` consequently failed. Translating the periodic domain and its
  probes gives a nonzero, analytically known centroid. The uniform Maxwell
  solution, exact surface/volume integrals, and original ``2e-13`` restart
  criterion are unchanged; an additional exact-centroid assertion is checked.
* A GPU particle-scrape check attempted a NumPy conversion of a CuPy buffer.
  The final repair, ``30c2d90f3``, evaluates ``(arr > 40).all()`` with the
  array's own backend. This preserves the original particle counts and
  ``stepScraped > 40`` condition, and avoids copying the entire buffer.
  The local test and both GPU simulation/analysis stages pass.
* RZ-only OpenMP test properties are set only when that test exists. This
  permits 1D-only OpenMP configurations.

Data verification
-----------------

Both low-energy source references were regenerated and independently refined
before comparison with the unchanged readable production data. The existing
acceptance criteria were retained: 0.2% for rate/transfer moments and 0.1% for
relative equilibrium heating/cooling imbalance.

.. list-table:: Maximum numerical errors
   :header-rows: 1

   * - Check
     - N2
     - O2
   * - Packed low-energy row moments
     - 0.09145%
     - 0.06738%
   * - Intermediate-energy moments
     - 0.11988%
     - 0.09952%
   * - High-energy moments
     - 0.10227%
     - 0.08416%
   * - Angle-weighted equilibrium imbalance
     - 0.00089%
     - 0.00650%

The implicit reverse-normalization checks agree to approximately
``4.3e-16`` (N2) and ``2.3e-16`` (O2). Source-anchor and independent
high-momentum phase-average checks also pass. These are numerical checks
against the adopted model, not measurements of molecular-model uncertainty.

Local regression results
------------------------

The main configuration builds 1D, 2D, 3D and RZ with MPI, OpenMP, Python,
FFT, embedded boundaries, QED and openPMD. AppleClang 21 and the repository's
updated AMReX commit ``b552c85c4`` are used, with pyAMReX/PICSAR 26.10 and the
required PICMI 0.34.0 release. The installed openPMD 0.17.1
provides both HDF5 and ADIOS2. ``FI_PROVIDER=tcp`` avoids a demonstrated
MPICH sockets-provider hang in ``MPI_Finalize``. ``MPLBACKEND=Agg`` avoids
interactive plotting waits in two existing analyses.

All 942 configured non-checksum CTest stages were rerun on integration
commit ``ac0125555``: 930 pass and 12 physics analyses fail. A fresh,
unchanged development build at ``63fa5a0c1`` with the same compiler and
dependency versions reproduces 11 of those:

* both lepton differential-luminosity comparisons;
* three focusing/rotated Gaussian-beam comparisons;
* the diffraction field-probe comparison;
* RZ Ohm-mode amplitudes and the 1D Ohm ion-beam growth-rate reference;
* the time-averaged PSATD acceleration restart comparison;
* effective-potential electrostatic density and spacecraft-charging fits.

The coarse helium MCC comparison is retained as a failure: 6.82% RMS density
error versus the original 6.50% bound (unchanged development gives 3.84%).
Independent refinements preserve physical run/averaging duration and the
original reference and assertion. Halving the timestep gives 3.43%, fourfold
particles/cell gives 4.42%, and joint refinement gives 3.13%; each passes.
A refined result does not replace the failed original coarse case.

All 380 main-configuration checksum stages were also run against completed
outputs: 200 pass and 180 differ. These platform-dependent failures are
ignored as instructed; no local benchmark reset or manual JSON edit was made.

Additional checks pass: 64 Python proton/source tests, three Python fluid-audit
tests, independent Gaussian-field quadrature, 12 standalone double-precision
C++ physics tests, three standalone all-single MCC tests, and 51 all-single
MCC/rotation stages. The 21 scheduling tests and 17 repaired-fixture tests
pass as part of the full run. Separate RCYLINDER and RSPHERE builds also
compile with the updated dependency; all six registered radial-geometry
stages pass, including both checksum comparisons.

A separate production-data elastic test uses the actual N2/O2 differential
and total cross sections in 22 cases from 0.001 eV to 1 GeV, with 65,536
electrons per case. Event probabilities, angular moments, inverse-CDF
statistics and finite-mass recoil conservation pass the existing analyses.
The optional legacy elmolcs rotation tests require legacy bundles absent
from this data branch; the supplied reciprocal and spectator families are
covered by the configured tests.

CPU timing baseline
-------------------

A warmed two-thread run at ``ac0125555`` with the updated dependency uses
262,144 electrons, both gases resident, 16 PIC
steps and three repetitions in alternating sampler order. The broad case
mixes 0.002 eV--3 MeV electrons. Times below are medians per step; they measure
complete PIC evolution and the separate complete MCC operator, respectively.

.. list-table:: CPU milliseconds per step
   :header-rows: 1

   * - Distribution
     - Sampler
     - PIC
     - MCC
   * - thermal
     - alias
     - 17.171
     - 6.950
   * - thermal
     - cumulative
     - 17.177
     - 6.963
   * - broad
     - alias
     - 20.077
     - 8.428
   * - broad
     - cumulative
     - 20.240
     - 9.864

The complete PIC medians are within 1% between samplers. The isolated broad
MCC median favors alias sampling by about 15%, while the thermal MCC medians
are nearly equal. These are three-repeat CPU measurements; they do not
establish a GPU speedup or isolate the small host-probability cache change.

Perlmutter results
------------------

The verified campaign uses A100-SXM4-80GB GPUs, GCC 13.2.1 and CUDA 13.2. It first
checks distinct device UUIDs, transfers a 4096-element integer buffer directly
between the GPUs through MPI, and runs the analytical coupled-MCC regression.
All three preflight checks pass. Four geometries (1D, 2D, 3D, RZ) build in
double precision; a separate 1D build uses single precision for fields and
particles. The C++ physics code tested is ``ac0125555``; the subsequent
``30c2d90f3`` change only repairs the particle-scrape test.

All 939 configured non-checksum stages were run. After the private pandas/NumPy
environment repair and the scraper rerun, 920 pass and 19 remain failed.
Seventeen reproduce on unchanged development with the same CUDA dependencies
and runtime. These include six stages associated with the decomposition guards
below and eleven other diagnostic, energy, reflection-count or restart
comparisons. The two remaining coarse stochastic discrepancies are examined
separately below. The exact failed names and comparison results are included
in the companion JSON. No assertion tolerance was relaxed.

All 380 checksum stages were rerun against the final outputs: 337 pass and
43 fail. These are ignored as requested; some lack complete output because
their simulation stops at a decomposition guard. No checksum was reset.

All targeted new MCC, prescribed-beam, source, immobile-ion and restart
regressions pass. Additional GPU checks pass: 12 standalone double-precision
physics tests, three standalone single-precision MCC tests, all 72 configured
single-precision MCC/rotation stages, and the 22-case production-data elastic
test. CUDA Compute Sanitizer reports zero memory errors for both N2/O2
reciprocal samplers and for the single-precision coupled excitation/attachment
test with unequal 64/96 subcycles.

The launcher exposes the allocated GPU set to a Slurm step and assigns one
entry of that set to each local MPI rank through ``CUDA_VISIBLE_DEVICES``.
This avoids inherited controller masks and per-task device cgroups that can
block Cray MPI GPU peer access, an interaction described in the `HPE MPI
guide <https://h41374.www4.hpe.com/docs/25.09/mpt/mpich9/intro_mpi.html>`_.
GPU-aware MPI and IPC remain enabled in the verified runs. Remaining checks
also use the repository's full-node launcher on a dedicated four-GPU node.
The initial attempt
with incorrect binding and subsequent IPC errors is retained as diagnostic
evidence, not used as proof of cross-GPU validation.

Two standard implicit JFNK cases and one refined Langmuir case abort the
original two-GPU run on a high performance warning: their 25 and 29 boxes
mean averages of 12 and 14 per GPU. Supplemental
four-GPU runs preserve the exact mesh, input parameters and analysis assertions.
All three cases pass simulation and analysis on both integration and unchanged
development (12 supplemental stages). The original two-GPU guard failures
remain recorded.

The original Vay collision energy test gives a maximum relative energy error
of ``6.0148969e-5``, slightly above its unchanged ``6e-5`` bound. Three further
runs of the same input on integration give ``[5.4483, 5.6253, 5.1709]e-5``;
unchanged development gives ``[5.2480, 5.2892, 5.7525]e-5``. All six pass.
Halving the timestep while preserving physical duration and the diagnostic
averaging window gives ``1.3662e-5`` on integration and ``1.3082e-5`` on
development, both passing. The overlapping repeats and approximately fourfold
error reduction support marginal stochastic/numerical variation; they do not
erase the first failed run.

The coarse DSMC discharge comparison gives 7.54% RMS error on integration
versus the unchanged 6.50% bound (development gives 5.09%). Independent
integration refinements retain the same reference and assertion: timestep
4.93%, particles 4.60%, joint 3.67%. All pass. This remains a coarse-run
failure, with convergence evidence rather than a changed tolerance.

GPU subcycle timing
~~~~~~~~~~~~~~~~~~~

The benchmark uses 1,048,576 electrons, both molecular gases resident,
4096 cells, eight measured PIC steps, 64 collision subcycles per gas per PIC
step, and ``dt_PIC = 6.4e-13 s``. Three repetitions follow a warmup, with
alternating sampler order. Each timed interval includes a completed particle
reduction to synchronize device work. These are complete PIC-step times on
one A100 80 GB GPU, not launch-only timings.

.. list-table:: GPU milliseconds per PIC step
   :header-rows: 1

   * - Distribution
     - Alias
     - Cumulative
   * - thermal
     - 36.656
     - 36.752
   * - broad (0.002 eV--3 MeV)
     - 52.114
     - 56.579

Alias sampling takes about 7.9% less elapsed time in the broad-energy case;
thermal times are nearly equal. The corresponding alias throughput is
3.66 billion (thermal) and 2.58 billion (broad) particle/gas/substep visits
per second. These measurements compare samplers in the final implementation;
they do not separately measure the effect of caching the host probability.

Reproduction and limits
-----------------------

Detailed local logs, CMake caches, JUnit results, input/source hashes and
reference outputs are retained under ``build/validation``. A compact result
and revision inventory is in `development_sync_2026_10.json
<development_sync_2026_10.json>`_. The unchanged
source comparison is under ``build/stock-source``. Local build directories are
inside this worktree. The Perlmutter validation uses a separate checkout, with
its build directories contained there, and
build-local CUDA 13 BLAS++/LAPACK++ libraries; it does not replace site or user
installations. The validation launcher assigns a distinct allocated device
to each rank, with the `NERSC affinity guidance
<https://docs.nersc.gov/jobs/affinity/>`_ and the Cray MPI peer-access
constraints above.

Local evidence is grouped under ``latest-development-local``,
``latest-discharge-refinement`` and ``latest-cpu-performance`` within
``build/validation``. Retrieved GPU evidence is under the
``perlmutter-gpu-59472340``, ``perlmutter-single-59476427``, four-rank-control,
scraper, Vay-control and DSMC-refinement directories. The single-precision
configuration is focused MCC coverage; HIP and SYCL runtimes were not tested.

The model retains its documented fixed neutral bath, prescribed-beam and
immobile-ion approximations, inferred rotational source continuations, and
small nominal-threshold recoil continuation. The paper's 34 eV mean energy
per ion pair is not substituted for molecular ionization thresholds.
The beam-example JSON still requires physical inelastic/attachment input
files; this audit does not invent missing cross sections or claim that the
experimental BPM response has been validated.
