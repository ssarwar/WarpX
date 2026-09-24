# Electron collision tables

The RBEQ exporter shares its host source definition with the WarpX sampler.
Build it inside the WarpX worktree:

```sh
cmake -S Tools/CrossSections -B build/cross-sections
cmake --build build/cross-sections
build/cross-sections/export_rbeq N2 iaa_thesis_2023 N2_ionization.txt
python Examples/Tests/collision/analysis_rbeq_sources.py build/cross-sections/export_rbeq
```

Use `O2` for oxygen and `elmolcs_b8643810` for the archived fitted model.
An optional final `raw` argument exports signed source totals for comparison;
these are not production MCC tables. Production tables sum positive partials,
span the outer-shell threshold to 1 GeV and target 0.02% linear-interpolation
error (with an absolute floor of 1e-8 of the peak). Their metadata are checked
by WarpX. The independent check requires NumPy and SciPy.

See `Docs/source/theory/mcc_iaa_sources.rst` for the source distinctions,
near-threshold continuation, and rotational production gates.
