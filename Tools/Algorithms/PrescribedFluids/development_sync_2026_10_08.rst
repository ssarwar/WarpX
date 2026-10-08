Fork development integration, 2026-10-08
========================================

Source provenance
-----------------

The new branch is ``codex/beam-air-development-sync-2026-10-08``.  Merge commit
``9863aadbc3766896b50e4452135cb85fc5f35104`` has these two parents:

* Fork development: ``321fb9bf257f36dc4f9b8dedb3ae4be1c0685267``.
* Previous integration: ``566bf52d5637033cd8f427d353d7e1922f8094fa``.

The merge base is ``0f72512ad0a51bc4995fe7165b856704577d0859``.  The first five new
development commits add implicit retries with smaller timesteps, MLMG solver
controls, and corrections to radial fluid vacuum edges, imaginary RZ implicit
electric-field gathering, and coarse-patch spectral divergence cleaning.
The previous integration's six recent development bug fixes are also retained.
An additional fork update arrived during validation:
``a5b7340f0f80feeface1e5ceec7acdec035b19fa`` corrects Darwin reduced diagnostics.
It was merged in ``2e3360f51`` after the RZ storage correction ``6ee8d5c5d``.
Neither source branch was changed.  The collision source and cross-section
tools are identical to the previous integration.  The data repository remains
at ``61a3947b6ebf365e6b78b63afdc30306791e3898``.

Conflicts and corrections
-------------------------

``FieldPoyntingFlux.H`` and ``FieldPoyntingFlux.cpp``
    Development adds an accepted-substep duration to the diagnostic call.
    The feature branch already tracks sampled/integrated PIC step indices to
    avoid integrating initial output, checkpoint output, and end-of-step
    duplicates.  Keeping the old equality guard would discard all but the
    first accepted substep.  The merge integrates each accepted substep with
    its own duration while retaining the ordinary-output and restart guards.
    Rejected nonlinear solves do not invoke the diagnostic.  The existing
    checkpoint format is preserved.

``SemiImplicitEM.cpp``
    Development splits ``OneStep`` into setup, solve, reset, and finish hooks.
    The feature branch's immobile-charge preparation is preserved in setup;
    the common base solver now owns the timestep.  Review also found that the
    new semi-implicit retry path omitted the particle timestep scale already
    present in the theta solver.  It now passes ``1 / m_nsubsteps`` to
    ``PreRHSOp``, so field and particle updates use the same accepted duration.
    The newly introduced magnetic rollback buffer also hard-coded one
    component.  Its allocation now matches the magnetic field, preserving all
    packed real/imaginary RZ modes and avoiding an out-of-bounds restore copy.

``RelativisticExplicitES.cpp`` (automatic merge)
    The rigid-beam self-field overload retained the removed scalar Poisson
    options and the old ``computePhi`` signature.  It now copies ``MLMGOptions``
    and applies the existing per-fluid tolerance, iteration, and verbosity
    overrides.  New global bottom-solver/coarsening controls are retained.

``FieldProbe.cpp`` (second merge)
    The feature branch adds ``<limits>`` for checkpoint metadata precision;
    the new Darwin diagnostic adds ``<memory>`` for temporary total fields.
    Both includes are retained.  The automatically merged body retains
    ``E_total = E_electrostatic - dA/dt`` and the feature branch's elapsed-time
    integration and restart handling.

Analytical regression
---------------------

Two new MPI regressions force the theta and semi-implicit solvers to retry a
uniform cold-plasma step.  Each accepts eight substeps.  Neutral tracers check
the exact relativistic free-streaming distance ``Delta z = v_z dt``.
For uniform magnetic field ``B_y``, the integrated Lorentz equation gives
``integral(E_x dt) = (m/q) Delta u_x + B_y Delta z``.  Multiplying by
``B_y / mu_0`` independently predicts the integrated face Poynting flux and
checks that every accepted interval, and no initial/duplicate interval,
contributes.  The local relative discrepancy is approximately ``9e-10``.
Restoring the missing-scale behavior makes the semi-implicit regression fail
to converge even at eight substeps; the corrected version passes.

The RZ storage regression initializes the production semi-implicit solver
with three azimuthal modes, then checks exact magnetic save/restore for all
three vector directions and all five packed components, including ghost cells.
This is a storage test, not a claim that full multimode implicit dynamics is
validated: the existing ``WarpXSolverVec`` copy/vector algebra still operates
on one component.  Exploratory full multimode evolution exposed that older
limitation; it is left outside this merge's scope.  No tolerances were changed
to accommodate it.  Both the storage regression and the analytical 1D
regressions pass locally.  Their append-only Picard diagnostic is removed
before each run so repeated CTest executions are independent.

Validation
----------

Local validation is complete on source revision ``2e3360f51``.  All configured
1D, 2D, 3D, RZ, RCYLINDER, and RSPHERE libraries/executables rebuilt successfully,
as did the separate all-single-precision 1D build.  The primary build enables
MPI, OpenMP, Python, FFT, openPMD, embedded boundaries, and QED.  Dependency pins
are unchanged; AMReX remains ``b552c85c432defab450688581a1c87aead7b461d``.

.. list-table:: Local results
   :header-rows: 1

   * - Suite
     - Passed / stages
     - Remaining failures
   * - Final complete CPU regression
     - 941 / 953
     - The same 12 previously recorded failures; no new failures or skips
   * - Implicit/fluid rerun after the RZ storage correction
     - 170 / 170
     - None
   * - New substep and storage regressions
     - 3 / 3
     - None; also pass in the final complete run
   * - All-single-precision MCC/rotation
     - 51 / 51
     - None
   * - Independent MCC kernels, double and single
     - 3 / 3 in each precision
     - None
   * - Independent proton and prescribed-fluid kernels
     - 5 / 5 and 4 / 4
     - None
   * - Radial geometry physics and checksums
     - 4 / 4 and 2 / 2
     - None
   * - Main checksum stages
     - 200 / 381
     - 181 informational checksum differences

The final complete run includes the late Darwin update.  The failure names
and source hashes are retained in ``development_sync_2026_10_08.json``.
Recurring failures were recorded without reopening their investigations.
Checksum differences are informational under the repository policy.  No
existing assertion tolerances or checksum references were changed.

GPU validation is pending an external resource hold.  At the recorded
2026-10-08 15:19 Pacific scheduler snapshot, NERSC had an active reservation
covering all ``scratch`` licenses.  The preliminary debug build reported
``PENDING/Licenses``; the final regular build ``59574018`` reports
``PENDING/Priority`` and preparation job ``59574019`` waits on its successful
completion.  No new CUDA compilation or GPU-test result is claimed.

The queued pipeline verifies the dedicated checkout is clean, imports the
exact source bundle, builds all four CUDA geometries and single precision,
then prepares four disjoint regression batches with shared-directory and
dependency grouping.  It subsequently runs checksums, single-precision MCC,
independent kernel checks, N2/O2 sampler memory checks, and a CUDA memory check
of RZ rollback storage.  The original preliminary jobs ``59570210`` and
``59570598`` were canceled before execution when the source was updated.
The future GPU batch IDs will be recorded in ``gpu-job-ids.txt``.

Raw evidence and launch scripts are retained under
``build/validation/development-sync-2026-10-08`` in this worktree and in the
dedicated Perlmutter checkout
``/pscratch/sd/s/ssarwar/warpx-beam-air-sync-20261006``.  The GPU build/test
pipeline remains queued; completion and review of those results are the
remaining validation work.
