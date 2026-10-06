# Reproducing the saved PJG work

Run commands from the repository root. Python needs NumPy, SciPy and
Matplotlib for fitting/plots; the archive checksum verifier itself uses
only the standard library. The historical host environment used Python
3.12, NumPy 2.5.1 and SciPy 1.18.0. Last-bit optimizer or image-byte identity
is not promised across library, renderer or platform versions.

## Verify and redraw the frozen final model offline

```sh
python Tools/Algorithms/ProtonImpactIonization/verify_research_archive.py
python -m unittest discover -s Tools/Algorithms/ProtonImpactIonization -v
python Tools/Algorithms/ProtonImpactIonization/fit_pjg.py \
  --nifs Tools/Algorithms/ProtonImpactIonization/Research/inputs/nifs-optical.json \
  --o2-leiden Tools/Algorithms/ProtonImpactIonization/Research/inputs/leiden-o2.txt \
  --output Tools/Algorithms/ProtonImpactIonization/Research/final/frozen-validation.json \
  --figures build/pjg-archive-redrawn --read-results
```

With `--read-results`, the output argument is an **input** result file and is
not overwritten. The figures go into `build/`. These instructions do not
need an original PDF, network access, a WarpX build or GLOW/LXCat inputs.

To recompute all frozen-model properties and the alternate-start check:

```sh
python Tools/Algorithms/ProtonImpactIonization/fit_pjg.py \
  --nifs Tools/Algorithms/ProtonImpactIonization/Research/inputs/nifs-optical.json \
  --o2-leiden Tools/Algorithms/ProtonImpactIonization/Research/inputs/leiden-o2.txt \
  --output build/pjg-archive-recomputed.json \
  --figures build/pjg-archive-recomputed-figures --frozen
```

Omit `--frozen` to repeat the six/five-parameter candidate refits. Neither
mode automatically rewrites production coefficients. Compare the selected
parameter names, SDCS/totals, ratios, moments and source hashes before
considering any numerical differences. The original final record includes
`initial_fits`, the SHA-256 of
`history/results/source-refresh/source-refits.json`; the current tool starts
from the frozen calibration and therefore does not have that provenance key.

## Historical source snapshots

The archives preserve original source bytes and original relative paths.
They contain regular files only, with per-member digests in `manifest.json`.
Use a **new** staging directory, separate from production, and verify first:

```sh
python Tools/Algorithms/ProtonImpactIonization/verify_research_archive.py
mkdir build/pjg-history-replay
tar -xzf Tools/Algorithms/ProtonImpactIonization/Research/history/preproduction-source.tar.gz \
  -C build/pjg-history-replay
tar -xzf Tools/Algorithms/ProtonImpactIonization/Research/history/exploratory-source.tar.gz \
  -C build/pjg-history-replay
python -m unittest discover \
  -s build/pjg-history-replay/Tools/Algorithms/ProtonImpactIonization -v
```

Do not extract over the live `Tools/Algorithms/ProtonImpactIonization` tree.
The first archive is the cleaned preproduction backup, without `__pycache__`
or macOS AppleDouble resource-fork sidecars.
It includes `finalize_pjg.py`, `pjg_matched.py`, `pjg_repair.py`, all fit/audit
drivers, data readers, tests, old CMake/C++ tests and the original research
notes. The second contains the earlier `tmp/pjg_*.py`/C++ diagnostics, input
deck, build helper and the original September 6 report. The manifest is also
the exact member inventory.

| Historical driver/module | Purpose and saved output |
| --- | --- |
| `pjg_limit_audit.py`, `pjg_consistency_audit.py`, `pjg_audit_extras.py` | Original/corrected/total-only-refit limits, analytic moments, lookup and sampling audit; `history/results/exploratory/pjg_*results.json` |
| `pjg_rebuild_experiment.py`, `pjg_gos_experiment.py` | Explicit dipole/hard and momentum-continuation trials; `rebuild-*`, `gos-*` |
| `pjg_difference_fit.py`, `pjg_minimal_experiment.py` | Parameter-count, subtraction-power, cutoff and Rudd-weight trials; `difference-fit-*`, `minimal-*` |
| `pjg_repair.py`, `fit_pjg_repair.py` | Positive common-center algebra and early nonrelativistic bound-tail fits |
| `pjg_matched.py`, `fit_pjg_matched.py`, `audit_pjg_matched.py` | Relativistic matched candidate, multistart checks and early matched figures |
| `stopping_pjg_matched.py`, `study_optical_stopping.py`, `fit_pjg_optical.py` | Independent stopping integration and optical-weight trials |
| `optical_dataset_audit.py`, `study_source_refresh.py`, `audit_source_refresh.py` | NIFS/Leiden/Mahla comparisons, source-aware refits and convergence/digitization sensitivity |
| `audit_pjg_properties.py`, `study_mean_energy.py` | Full property grids and mean-energy-weight sensitivity |
| `finalize_pjg.py` | Original final objective, candidate reductions, properties and figures |
| `pjg_table_audit.cpp`, `build_pjg_table_audit.py`, `pjg_table_design.py` | Superseded production table diagnostic and design study; see raw-data layout below |

Use the CLI help and original reports for historical options. Supply the
archived optical paths explicitly (absolute paths are convenient in a staging
directory). `finalize_pjg.py --initial-fits` can use the archived source-refit
JSON to reconstruct the original initialization. Some old scripts expect
their `tmp/` working-directory layout. Do not mistake the scratch
`pjg_reference.py`'s `correction="printed"` option alone for the fully printed
baseline: the later limit driver also selects the printed endpoint. Use the
current `original_pjg.py` for the authoritative baseline comparison.

The scratch C++ build helper contains the original host compiler paths and
links the **old** production PJG object; it is provenance, not a portable
build command for the new implementation. Rebuilding that diagnostic requires
the code at `665031685` and a matching build inside that checkout. The saved
raw table allows inspection without recompiling. Older GLOW-dependent fits
need separately acquired tables from the cited revision; final reproduction
does not. No claim is made that every historical executable can run unchanged
against today's production API.

## Saved numerical formats

Final JSON is organized as `targets.N2` and `targets.O2`, with `selected`
naming the accepted candidate dictionary. It includes full coefficients,
objective norm, source-specific ratio summaries, individual measured checks,
alternative starts, distributions, convergence, hard-tail and stopping rows.
The readable reports and current driver specify the grids and weights.

Unless a key or column says otherwise, research energies are eV, integrated
cross sections cm2 per molecule, SDCS cm2/eV, and ratios dimensionless. The
`parameters.amplitude` is K/K_printed, not K in cm2. Stopping comparisons
use MeV cm2/g. The source-refresh `properties.json` includes a `units` object;
it describes a pre-final fit, not the selected calibration. Rudd-ratio arrays
in exploratory files use recommendations unless explicitly labeled measured.
Preserve signed legacy spectra when integrating; logarithmic plotting cannot
display zero or negative values.

`validation/legacy-numerical-data.tar.gz` contains
`tmp/pjg_table_audit.bin`, which is numerical data, **not an executable**.
Read it with NumPy dtype `<f8`, shape `(148, 2049, 10)` (little-endian
float64 from the original host). There are 74 incident energies per gas.
Columns are target index (0=N2, 1=O2), incident energy eV, quantile, host total
m2, table total m2, host free Tmax eV, executor free Tmax eV, sampled secondary
energy eV, effective binding eV and host SDCS m2/eV. The generating C++ and
independent Python audit are in the exploratory source archive. These are
outputs of the superseded total-only refit, not production calibration tables.

The same archive retains initial NPZ diagnostics, used-input files and
OpenPMD HDF5 smoke outputs. The separate production archive retains final
source NPZ diagnostics and used-input files. Their original relative paths
and hashes identify the associated test. Load NPZ with `allow_pickle=False`;
use the corresponding analysis script in `Examples/Tests/proton_impact_ionization`
to interpret named fields. Do not treat particle samples as experimental data.
