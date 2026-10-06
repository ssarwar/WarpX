.. _rotation-data-format:

Production files and initialization format
==========================================

What belongs in each repository
-------------------------------

WarpX contains model equations, source readers, exporters, independent tests,
initialization and GPU sampling code. warpx-data contains numerical source
inputs and production collision distributions. Temporary reference
checkpoints, synthetic parser fixtures, validation outputs, timing records,
commit identifiers and benchmark results do not belong in warpx-data.

The parent ``IAA`` directory holds readable ordinary cross sections and
elastic DCS tables. ``reciprocal_sources`` holds source-specific numerical
constraints, including their units and interpretation. Each
``reciprocal_hybrid_300K`` directory represents one combined collision family,
not a collection of additional processes to add to the elastic rate.

Why the original files were large
---------------------------------

V6 stored final device arrays directly. The N2 bundle contained 55,163,303
8-byte alias entries; O2 contained 31,779,555. Many neighboring angular cells
repeat closely related distributions over the same thermally recoupled
rotational transitions. Aliases therefore occupied 441,306,424 bytes for N2
and 254,236,440 bytes for O2. These are sampling acceleration arrays, not
that many measured cross sections. Including grids, outcomes and search
indices, the two V6 payloads totaled 751,195,072 bytes (716.4 MiB).

The 15 N2 and nine O2 ``.bin`` files were arbitrary slices of a concatenated
payload, each limited to 32 MiB. A file boundary had no physical meaning.
Splitting kept individual files below Git's large-file limit; it did not
reduce total data size. Binary storage made startup simple but made source
review and interpretation unnecessarily difficult.

Readable V7 probability inputs
------------------------------

V7 stores decimal numerical tables and a compact factorization of the
conditional probability matrices. Initialization expands that representation
and builds the device aliases. The numerical factorization is an encoding
of a previously validated probability distribution, not a new molecular
model, and its matrix rank has no relationship to rotational tensor rank
:math:`L`. Every discrete outcome and its exact support are retained.

The index ``thermal_rotation.rot`` has the form::

   WARPX_THERMAL_ROTATION_V7
   N2 reciprocal_hybrid probabilities
   temperature mass maximum_energy high_energy rutherford_energy separation screening_radius
   array_count
   array_name scalar_type element_count relative_filename
   ...
   probabilities alias alias_entry_count probabilities.txt

The third line contains numbers, in K, kg, eV, eV, eV, Bohr radii and Bohr
radii respectively. Filenames refer to the index's directory. The final
probability line is additional to ``array_count``. There is no implicit source
lookup, downloading, temperature change or extrapolation in this index.

Each named array is a text file with ``#`` comments. ``values N`` is followed
by N decimal values; ``repeat N value`` repeats one value; ``copy N start``
reuses N already written values beginning at zero-based index ``start``.
The copy range must end before the current output position. These operations
are a transparent way to avoid repeating identical angular rows and zeros.
They change no values. The declared scalar types are ``f64``, ``f32``,
``u32`` and ``u16``. Binary64/32 values use sufficient decimal digits for
exact round trips; integer ranges are checked.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Array
     - Meaning and units
   * - ``energies``, ``rates``
     - Electron energy in eV and combined :math:`K` in m3/s; both start at zero energy.
   * - ``coordinates``
     - Interval flag: 0 is linear; 1 is square-root interpolation. The final flag is unused.
   * - ``angular_offsets``
     - Start/end positions of each energy row in the angular arrays, with a final sentinel.
   * - ``angular_u``, ``deflection``
     - Angular CDF probability and :math:`1-\cos\theta`, dimensionless.
   * - ``changing``
     - Probability of nonzero rotational energy change at that angular node.
   * - ``conditional_offsets``, ``conditional_u``
     - Row offsets and angular quantiles for the spectrum conditional on changing rotation.
   * - ``conditional_cells``
     - Zero-based probability-cell identifiers at those conditional quantiles.
   * - ``cell_offsets``
     - Start/end of each cell's expanded sampling entries, including the final sentinel.
   * - ``outcomes``
     - Flattened pairs :math:`(\Delta,\epsilon_i)` in eV. Outcome zero is unchanged.
   * - ``high_edges``, ``high_cells``
     - Momentum-transfer bin edges in dimensionless :math:`z` and their cell identifiers.

Angular rows above 10 keV may be empty because their inverse CDF is analytic.
Conditional rows above 1 keV may be empty because the momentum-transfer bank
supplies the discrete spectrum. Empty rows, sentinels and supports are checked
before any device use.

Probability reconstruction
--------------------------

``probabilities.txt`` starts with::

   WARPX_PROBABILITY_FACTORS_V1 number_of_cells number_of_alias_entries

Each consecutive block contains at most 256 cells::

   block first_cell cell_count outcome_count factor_rank
   outcome_id interval_count first end ... left_factor_0 ...
   ... one line for each outcome ...
   ... one right-factor row for each cell ...

Support intervals are zero-based, half-open cell ranges within the block.
They preserve exactly which outcomes can occur, including closed excitation
channels and the zero-change outcome. If :math:`L_{ik}` and :math:`R_{jk}`
are the two stored factor matrices and :math:`S_{ij}` is the support mask,
initialization reconstructs

.. math::
   :label: rotation-text-reconstruction

   q_{ij}=S_{ij}\max\!\left(\sum_k L_{ik}R_{jk},0\right),\qquad
   p_{ij}=q_{ij}/\sum_l q_{lj}.

Signed factor coefficients are allowed; they are not signed physical
cross sections. The small projection in this numerical encoding is distinct
from clipping a negative physical unchanged residual, which the source model
rejects. The reader rejects nonfinite values, inconsistent supports,
excessive negative reconstruction and normalization defects. It then builds
aliases with the standard small/large-column algorithm: scale probabilities
by cell size, pair a column below one with a column above one, and transfer
the excess until all columns represent equal-area sampling slots.

The offline encoder uses a weighted singular-value decomposition, selecting
its rank against reconstructed, decimal-rounded probabilities. It checks
every cell's total variation and unchanged/upward/downward probabilities,
positive and negative first moments, second moment and fourth moment.
The maximum L1 probability error is :math:`2\times10^{-7}` and the maximum
relative error in the listed nonzero moments is :math:`2\times10^{-6}`.
These are additional numerical encoding budgets, much smaller than the
0.2% total discretization target. A zero exact moment remains zero because
support is stored explicitly. Singular values alone do not determine the
accepted rank. The runtime cannot independently certify source accuracy;
that remains an export/validation responsibility.

The same readable input can prepare aliases or a cumulative reference through
``rotation_sampling``. The cumulative construction decodes the float32 alias
probabilities and accumulates them with extended host precision before
normalizing the entire prefix; this avoids a last-entry rounding inversion.
No second production data set is needed for that choice. Test tools also
reconstruct probabilities directly as an independent cumulative reference.

Initialization, memory and backward compatibility
-------------------------------------------------

Initialization additionally constructs 257 quantile-search bounds per energy
row. These bracket the full binary search without coarsening its grid. The
runtime validates all dimensions, monotonicity, probability bounds, canonical
rotational levels, parity, molecular masses, range coverage, temperature and
threshold accessibility. It checks the combined source rate against the
ordinary elastic table above the cold join. No prepared file is written by
the simulation, and no Python or external numerical library is required to
read V7. A restart prepares the immutable arrays again; repeated collision
objects within a process share them.

The 8-byte device alias layout remains float32 cutoff, uint16 alternate
index and uint16 outcome identifier. Energies and outcome labels remain
double precision. The per-process cache enforces a combined 1 GiB table
budget. One MPI rank per GPU shares these tables within that rank; separate
ranks on one GPU each allocate their own copies. Choosing cumulative storage
uses additional memory and is primarily a reference option.

V6 remains supported with its original little-endian binary64 host
requirement. Its directory lists ``name type count byte_offset`` followed by
relative binary-part filenames and byte counts. Offsets are eight-byte aligned
within the concatenated payload. Optional V6 quantile indices are accepted;
without them the full search is used. Unmarked ordinary cross-section files
retain their legacy endpoint behavior. V7 does not silently change any legacy
rotation option or source continuation.
