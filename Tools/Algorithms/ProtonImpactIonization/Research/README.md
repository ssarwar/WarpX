# PJG research archive: September 2026 calibration

This archive preserves the recoverable calculations, figures, numerical inputs,
fit outputs and research notes behind the N2/O2 proton-impact source. It was
assembled on 2026-09-14 on `codex/proton-impact-ionization-development-sync`,
starting from `967f48efa1a1f473916c9019118901d75b97b228`. Archiving does not
change the SDCS, coefficients, angular closure or collision implementation.

Start with the [complete theory chapter](../../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst)
for the equations, printed parameters, 1977 corrections, derivations, final
parameters, fit objective, data comparisons and physical limitations.
The [production reference guide](../README.md) describes the current code and
tests. The notes below explain how that model was selected, including failed
and superseded approaches. They are not alternative runtime models.

## Finding a result

| Question | Record |
| --- | --- |
| What are the final parameters and all candidate-reduction results? | [Original final fit](final/fit-results.json) and [final formulation note](history/reports/FINAL_PJG.md) |
| What does the frozen production model predict? | [Frozen evaluation](final/frozen-validation.json) and [figure gallery](FIGURES.md) |
| Was the optimizer rerun independently? | [Production-tool refit](final/production-refit.json); not replacement coefficients |
| Which data were measurements, evaluations or recommendations? | [Dataset guide](DATA_SOURCES.md), [detailed source audit](history/reports/DATASET_AUDIT.md), [numerical inputs](inputs/) |
| Why were the original and first corrected models insufficient? | [Initial limit/consistency audit](history/reports/PJG-limit-and-consistency-audit-2026-09-06.md) |
| Where is the algebra of the hard-continuum correction? | [PJG-preserving repair](history/reports/PJG_REPAIR.md), then [matched model](history/reports/MATCHED_PJG.md) |
| What happened to an explicit GOS reconstruction? | [Reconstruction assessment](history/reports/RECONSTRUCTION.md) and [exploratory results](history/results/exploratory/) |
| What were the optical/stopping tradeoffs? | [Optical/stopping audit](history/reports/OPTICAL_STOPPING_AUDIT.md), [source-refresh results](history/results/source-refresh/) |
| Where are means, quantiles, tails, straggling and convergence? | [Property audit](history/reports/PROPERTIES_AUDIT.md), final and source-refresh JSON files |
| Where are the original scripts and saved PIC/table data? | [Reproduction instructions](REPRODUCING.md), [validation inventory](VALIDATION.md), [checksum manifest](manifest.json) |

## Model lineage and decisions

1. **Printed PJG, including the 1977 erratum.** Throughout the final comparison,
   "original" means this baseline, not PJG with the erratum omitted. In Eq. (16),
   the raised plus means addition, not positive-part clipping. The original
   cutoff and signed SDCS must be retained when diagnosing the printed model.
   The current offline evaluator is [original_pjg.py](../original_pjg.py).

2. **Kinematic/algebraic repair and total-only refit.** Exact free-electron
   maximum transfer and the corrected Bhabha remainder were installed, with
   the empirical cutoff offset delta removed. J-only and J/K adjustments could
   improve totals but did not repair the coefficient of the full hard tail:
   the empirical Lorentzian itself already contains a T^-2 contribution.
   Adding a complete hard cross section to it would double-count that term.
   This stage is retained in Git at `665031685` and in the initial audit.
   The printed low-energy N2 Fig. 6 inconsistency was not resolved by restoring
   the printed errors. Its underlying source/transcription cause remains
   unproven; no new dataset or correction was invented to explain it.

3. **Limit and numerical audit.** The September 6 study found high-energy
   negative printed spectra, inadequate hard normalization, single-precision
   cancellation in the free endpoint and a mixed-parent/secondary-quantile
   correlation. These findings are historical, not defects asserted against
   the final implementation. The sampler and kinematics fixes are retained
   in the branch history (`4b4a52128`, `cddd09e89`). The independent invariant
   and Dirac-trace checks do not rely only on transcribing Bhabha's paper.

4. **Exploratory hard/dipole and GOS models.** The `rebuild-*`, `gos-*`,
   `difference-fit-*` and `minimal-*` results record attempts to repair the
   spectrum with few adjustable coefficients. A positive continuation in
   momentum transfer is not a uniquely determined molecular GOS or a complete
   PWBA calculation. Fits to Rudd recommendations are not fits to hundreds
   of independent measurements. These exploratory candidates were not
   promoted. Some require the separately obtained GLOW data to rerun.

5. **Positive PJG-preserving decomposition.** A common-center difference of
   narrow and broad Lorentzians gives a positive soft term falling as T^-4.
   A separately normalized broad term retains the complete T^-2 hard kernel.
   The repair note derives the coefficient by equating the inverse-square
   terms, rather than fitting away the Bhabha normalization. The relativistic
   dipole logarithm is fixed in the fast-projectile limit; its molecular
   coefficient remains an empirical optical approximation. This is not an
   exact pointwise optical spectrum or full finite-q response.

6. **Molecular support and matched candidate.** The exact free-collision Tmax
   stays fixed. It is not imposed as a hard cutoff on bound-electron emission.
   The smooth above-free tail is bounded by molecular energy conservation.
   Low-velocity distortion fades out of the hard term at high velocity.
   Earlier matched fits used Rudd recommendations and a provisional O2 optical
   proxy. They remain in the matched and optical/stopping notes, not in
   production. Improving an optical residual alone could overfill the NIST
   stopping budget; those trials were retained as evidence against accepting
   them solely on that residual.

7. **Source-aware refresh and property sensitivity.** Uploaded papers allowed
   measured N2/O2 information to replace or weaken recommendation-only
   constraints. Rudd (1983) Table V is an authors' fit, not raw measurements;
   Cheng (1989) spectra have a shared, adjusted absolute normalization. NIFS
   absorption and Leiden ionization are distinct optical observables. Mahla
   and Mehnen (2025) were audited but not adopted as the absolute anchor.
   The property study records stopping decomposition, means, second moments,
   quantiles, bound-tail fractions, optical moments, angular-support caveats,
   and sensitivity to the N2 mean-energy weight.

8. **Final selection and production port.** The joint fit selects six N2
   parameters and five O2 parameters, with O2's edge scale fixed to one.
   Delta, the separate low-energy width ratio and independent logarithm-scale
   adjustments are not free fit parameters. The total remains the SDCS
   integral. The final objective balances measured spectral shape, totals,
   N2 means and optical information; NIST stopping is checked independently,
   not fitted to equality. The frozen reference entered at `b2da4ebe8`, the
   fit reproduction at `42c8172e3`, the production tables at `35fcb37ff`, and
   the full theory documentation at `3afc063a9`.

9. **Implementation validation.** Later commits test the device table executor,
   low-noise source, cold backgrounds, mixed energies and budgets. The
   forward-angle point mass was removed by conditioning the closure's interval
   on physical cosines (`d940c19c6`); that does not establish agreement with
   backward experimental DDCS. Checkpoint budget preservation is tested at
   `27cb037d8`. Saved CPU results do not establish CUDA/HIP/SYCL performance.

## What is final, and what is not?

`final/fit-results.json` preserves the original final research run, including
all six/five-parameter candidate results. Its selected parameter dictionaries
match the production [PARAMETERS](../pjg_model.py) exactly. Its status string
still says "not yet installed in production": that was true when the file
was generated. It is preserved verbatim, not a statement of present status.

`final/frozen-validation.json` evaluates those same frozen parameters with the
production reproduction tool. `final/production-refit.json` is a later
optimization run with very small numerical differences in coefficients. The
latter demonstrates repeatability but does not redefine the final calibration.
The tests check this distinction and the input-hash chain.

The final model is a useful compromise, not an exact fit to all evidence.
For example, its Rudd-1983 total ratios span 0.834--1.096 for N2 and
0.879--1.023 for O2; the largest effective-pair/electronic-PSTAR ratios on
the saved grid are about 0.9385 and 0.9653. N2's predicted mean secondary
energy at 5 keV is about 6.94 eV versus the measured 5.42 eV. The full
residuals, grids and qualifications are preserved, including disagreements.

The observable is inclusive electron yield represented by effective pairs,
not exclusive single molecular ionization. The beam is ballistic and the
ions are neutral-thermal. A molecularly allowed electron does not make this
prescribed source event-by-event energy/momentum conserving. All secondary
energies are represented kinetically. Measured total constraints cover
5--4000 keV; the 5 keV--10 GeV numerical audit is not experimental validation
of relativistic molecular spectra. Higher-order radiative QED, projectile
structure, arbitrary-density dielectric effects and a complete finite-q GOS
were not fitted or validated by this work.

## Scope and preservation policy

The manifest identifies each payload by size, SHA-256, recovery location and
research stage. It also enumerates and hashes each compressed tar member.
Historical scripts are isolated in archives so ordinary test discovery and
formatters cannot silently change them or run superseded kernels. Human-
readable reports add a historical banner; the original bytes remain inside
the source archives. The initial audit's links are made portable in its
readable copy. Existing published final figures are copied into this archive
so later documentation edits cannot overwrite this calibration's evidence.

All recoverable saved, authored PJG results found in the worktree's `tmp/`,
`output/`, relevant `build/` diagnostics and the preproduction backup are
represented. Duplicate result/figure copies are identified by aliases.
Full papers, paper-page screenshots/OCR, restricted provisional source tables,
browser state, credentials, caches, executables, unrelated MCC research and
full checkpoint payloads are not republished. The saved diagnostic arrays and
input decks are included. See [source provenance](DATA_SOURCES.md) and
[validation provenance](VALIDATION.md) for details. Runs that were overwritten
or never saved cannot be reconstructed as original output; in particular,
the five-run performance summary is not five preserved raw timing records.
