# Perlmutter coupled-physics re-audit

This record covers the re-audit requested after the upstream merge at
`97078e27b`. The stock comparison is the merged upstream revision
`cb5672fae0b9099d2c7631e93ffff3404a553821`, not the older feature-branch base.
All new compilation, physics tests, regression tests and timing measurements
run on Perlmutter. Results remain pending until their logs are recorded below.

## Previously committed records

- [Fluid validation](VALIDATION.md), its compact JSON measurements and figures
  record the implementation assumptions, independent references, convergence,
  CPU/A100 measurements and known coverage limits through `3ac32229a`.
- [Proton model theory](../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst)
  and the [PJG research archive](../ProtonImpactIonization/Research/README.md)
  preserve the equations, calibration, input provenance, superseded approaches
  and uncertainty of the 800 MeV extrapolation.
- [Collision theory](../../../Docs/source/theory/multiphysics/collisions.rst)
  and [MCC checks](../BackgroundMCC/README.md) document the event selector,
  relativistic recoil, RBEQ energy sharing, IAA angular models and attachment.
- [Application examples](../../../Examples/Physics_applications/proton_beam_air/README.md)
  and the parameter reference document text/PICMI configuration and limits.

Raw build trees, checkpoints and most runtime logs were not committed. Compact
results and reproducible drivers were committed. This re-audit will retain
machine-readable acceptance summaries and provenance in Git; large raw outputs
remain in the authorized Perlmutter validation directory.

## Acceptance matrix

| Area | Independent check and regression coverage | Status |
| --- | --- | --- |
| Stock build | Current stock Perlmutter profile, documented CMake and dependency recipes, matching branch/stock configurations | Pending |
| Rigid beam | Gaussian normalization, exact space-time yield, continuity, self-fields, axis/walls, pulse overlap, shape orders | Pending |
| Immobile ions | Identical-event frozen-particle footprints, signed charge transfer, nonnegativity, no repeated synchronization/filtering | Pending |
| Proton impact | PJG N2/O2 spectra/totals, IAA angular closure, bare-ion scaling, tail precision, caps/remainders | Pending |
| Electron MCC | Elastic/excitation IAA DCS, RBEQ SDCS and IAA ionization angles, two-/three-body attachment, recoil and null events | Pending |
| Coupled operation | All channels together, both particle/fluid representations, subcycling and all four solver configurations | Pending |
| Diagnostics/restart | Native/plotfile/openPMD/reduced outputs, new upstream per-species particle counts, immediate restored state and continuation | Pending |
| GPU ownership | One/two/four GPUs on one node, eight GPUs on two nodes, CUDA-aware MPI and changed restart decomposition | Pending |
| Performance/noise | Synchronized strong/weak scaling, matched electron work, particle references and repeat-seed statistics | Pending |
| Stock regressions | Collision, fluid, implicit, diagnostics and restart; reproduce failures on the merged upstream revision | Pending |

The previous coupled timing driver included ionization and attachment but did
not include elastic or excitation. The new combined fixture must include both.
In the input API, the electron SDCS model is named `RBEQ`; `IAA` names the
ionization angular closure and the tabulated elastic/excitation DCS option.
Synthetic tables establish numerical correctness, not experimental air-chemistry
accuracy. No assertion tolerance will be relaxed to accommodate a failure.

## Build environment

The saved user profile refers to older Python/HDF5 modules. A fresh copy of the
current repository's `perlmutter_gpu_warpx.profile.example` is used with account
`m3748_g`. Existing BLAS++ links CUDA 12 and existing ADIOS2 links an unavailable
HDF5 library, while the stock profile now loads CUDA 13.2/HDF5 1.14.3.9. These
dependencies must be rebuilt with the stock recipes into an isolated validation
prefix. The user's shared installation and home profile are preserved.

The relevant optional features are MPI, Python, FFT and openPMD/HDF5/ADIOS2.
The validation build covers 1D, RZ and 3D, with double fields/particles. QED and
embedded boundaries are outside this simulation's scope. GPU-aware MPI and
rank-to-GPU placement follow the stock Perlmutter batch example. The earlier
single-GPU timings do not establish full-node or multi-node performance.

## Review findings and results

Pending. Record every fixed defect, rejected hypothesis, unresolved limitation,
test revision, compiler/dependency version, job ID and measured result here.
