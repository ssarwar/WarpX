.. _rotation-algorithm:

MCC selection, joint sampling and recoil
========================================

Preparation, initialization and device work
-------------------------------------------

The expensive physical construction remains outside the simulation: selecting
and fitting the molecular sources, summing the equilibrium rotor states,
solving the reciprocal inclusive constraint, refining energy/angle grids and
checking source accuracy. Production data are temperature-specific results
of that calculation. They do not require elmolcs, Python or a nonlinear solver
on a compute node.

Readable bundles reconstruct their tabulated probabilities and build alias
or cumulative sampling arrays once during initialization. This work is
linear algebra and table preparation, not a new fit to the physics sources.
Exact quantile-search bounds are also constructed at initialization. The
resulting device arrays are immutable, shared by collision objects using
the same canonical path, temperature and sampling choice, and bulk host
staging memory is released after upload. Old prepared binary bundles remain
readable. :ref:`rotation-data-format` specifies both representations.

At each collision there is no neutral-state sum, partition function, Bessel
function, table construction or transition-rejection loop. The new family
uses the existing MCC particle pass, selector, null collisions and recoil
infrastructure. Rotational states are virtual outcome labels, not particles.

Combining rates without double counting
---------------------------------------

The ordinary elastic table is first checked against the bundle's total
above the cold join. Its contribution is then replaced exactly once by
:math:`K_{\rm rot}`. With other processes :math:`r`,

.. math::

   \nu(E)=n_n\left[K_{\rm rot}(E)+v(E)\sum_r\sigma_r(E)\right].

A null-collision majorant bounds this sum over the supported collision-energy
interval. For a bounded ordinary table, each interval contribution is bounded
using the largest cross section and speed on that interval. The reciprocal
rate is a positive mixture of endpoint rates and cannot exceed their maximum.
The union of process and reciprocal knots supplies the interval partition.
There is no fictitious extension to :math:`c\sigma_{\rm last}` beyond the
supported domain. Unmarked legacy constant-extrapolation tables retain that
bound where appropriate.

The target translational velocity is sampled from the prescribed Maxwellian.
The existing electron MCC rate evaluation uses its inexpensive relative
proper-velocity approximation; see :ref:`multiphysics-collisions-mcc` for its
stationary-target exactness and moving-target error. The signed recoil update
below boosts explicitly to and from the chosen neutral rest frame. These two
steps should not be described as one exact relativistic rate integration over
a moving gas.

A preselected particle receives at most one accepted channel per collision
substep. Timestep/subcycle convergence remains necessary when the optical
depth is large. A selected unchanged rotational outcome is a real scattering
event, not a null collision; it retains elastic deflection and recoil.

Energy interpolation as a mixture of measures
---------------------------------------------

Let rows :math:`l,h` bracket :math:`E`, with interpolating fraction :math:`t`.
For most intervals :math:`t=(E-E_l)/(E_h-E_l)`. In an interval immediately
above zero or a canonical opening the stored flag selects
:math:`t=\sqrt{(E-E_l)/(E_h-E_l)}`. This resolves the threshold law without
forcing thousands of nearly identical linear knots.

.. math::
   :label: rotation-row-mixture

   K=(1-t)K_l+tK_h,\qquad
   w_l=\frac{(1-t)K_l}{K},\qquad w_h=\frac{tK_h}{K}.

A row is drawn using :math:`w_l,w_h`. The same chosen row supplies the
angle, the probability of changing rotation and the conditional discrete
outcome. These are mixtures of joint rate measures; choosing angle and
outcome from independently selected energy rows would destroy their
correlation. Energy-loss labels are never interpolated.

If a uniform :math:`u` selects the lower row when :math:`u<w_l`, the
conditional uniform is :math:`u/w_l`; otherwise it is
:math:`(u-w_l)/w_h`. Reusing the uniform in this way is exact because each
conditional interval maps uniformly onto [0,1). The implementation handles
zero-rate endpoints explicitly. The reader verifies that an upper row's
entire excitation support is already open at the lower edge of its interval.
Thus row mixing cannot select a subthreshold excitation. This is an input
consistency check, not an event-time angular veto.

For grids larger than 32 nodes, an exact search index narrows the energy
bracket before the original binary search. ``frexp`` identifies the binary
exponent and one of 64 mantissa bins; initialization stores interval bounds
at the corresponding power-of-two-scaled boundaries. This retains every
physical knot and the same interpolation operations. It introduces no
energy rounding or new interpolation approximation.

Angle and discrete rotational outcome
-------------------------------------

For a chosen low-energy row, the angular quantile locates two neighboring
inverse-CDF nodes. Their interpolated :math:`d=1-\cos\theta` gives the angle.
The same nodes interpolate the changing probability

.. math::

   \rho(E,u)=\frac{\sum_{i\ne f} b_iD_{if}(E,\theta(u))}
                   {D_{\rm inc}(E,\theta(u))}.

With another independent uniform :math:`r`, the outcome is unchanged if
:math:`r\ge\rho`. Otherwise :math:`r/\rho` is a uniform for the changing
spectrum. At neighboring conditional-quantile nodes, the spectrum is
interpolated as a mixture of their distributions. That mixture draw is again
rescaled onto its conditional interval before drawing a discrete outcome.
All zeros and discrete level spacings retain their meanings.

For an alias cell with :math:`n` entries and a uniform :math:`r`, compute
:math:`x=nr`, :math:`j=\lfloor x\rfloor`, :math:`f=x-j`. Select entry
:math:`j` if :math:`f<a_j`, otherwise its alternate. Initialization resolves
the alternate column to its outcome identifier, so either branch directly
addresses the palette of :math:`(\Delta,\epsilon_i)` pairs. This removes a
dependent alias-table read without changing a cutoff, outcome or random draw.
The startup alias construction and conversion cost :math:`O(n)` and are
shared across every subsequent draw from that cell.

A cumulative reference instead performs a strict upper-bound search on the
same decoded cumulative probabilities. An optional 64-bin quantile index
narrows that search, using two uint16 row-relative bounds. The original CDF
and comparison are unchanged, including flat intervals and boundary draws.
The index accelerates the reference option; alias remains the default.

Above 1 keV, the sampled deflection and actual energy determine
:math:`z=R P(E)\sqrt{d/2}/(\alpha m_ec^2)` and the momentum-transfer cell.
No energy-row interpolation of rotational labels is used there. Above
10 keV the angle comes from the analytic inverse in
:eq:`rotation-rutherford`, also using actual energy. The two high-energy
switches therefore have different roles.

Deflection and random precision
-------------------------------

Store and carry :math:`d` in double precision. Computing it as
:math:`1-\cos\theta` after rounding a strongly forward cosine to one can
lose the entire deflection. The transverse component is evaluated from
:math:`\sin\theta=\sqrt{d(2-d)}`. Canonical energies, gaps and sensitive
kinematic intermediates also use double precision; alias cutoffs use float32
and local/selected indices use checked uint16 values.

The reciprocal collision and sampling draws use adequately resolved double
uniforms even in all-single builds. CUDA, HIP and SYCL use their backend's
double-uniform generation; the CPU construction combines random integer
bits. Drawing a binary32 uniform and merely casting it to double would not
recover rare-tail resolution. This does not cure single-precision storage of
particle momentum: meV transfers can still round away for MeV electrons.
Double particle precision is recommended for that combination.

Signed internal energy and recoil
---------------------------------

Use energy units :math:`m=m_ec^2`, :math:`M=M_i c^2`, where the virtual
initial molecule has rest mass including :math:`\epsilon_i`. Its final
rest energy is :math:`M+\Delta`. In the initial molecule's rest frame,
energy and momentum conservation imply

.. math::
   :label: rotation-recoil-equation

   E+m+M=E'+m+
      \sqrt{(M+\Delta)^2+P(E)^2+P(E')^2
                         -2P(E)P(E')\cos\theta}.

The desired root has nonnegative :math:`E'` and the unsquared equation's
energy sign. Define

.. math::

   T&=E+m+M,\quad b=P(E)\cos\theta,\quad a=T^2-b^2,\\
   h&=M(E-\Delta)-mE-\tfrac12\Delta^2,\\
   C&=Th+b^2m,\quad S=b\sqrt{h(h+2Tm)+b^2m^2}.

Squaring once gives
:math:`aE'^2-2CE'+h^2=0`, so the physical solution is
:math:`E'=(C+S)/a`. When :math:`b<0` and :math:`h\ge0`, the equivalent
:math:`E'=h^2/(C-S)` avoids subtracting nearly equal terms. The virtual
molecular momentum is the vector difference between the incoming and outgoing
electron momenta. Both products are boosted back to the simulation frame;
the prescribed neutral bath discards the virtual recoil afterward.
De-excitation uses negative :math:`\Delta` in this same solve, not a clamp
of the loss to zero.

The source thresholds use the nominal internal gap, so a narrow excitation
band does not admit the exact sampled angle at finite target mass. The
agreed continuation is applied only when
:math:`0\le E-\Delta\le2(m_e/M_i)E`. There :math:`E'=E-\Delta` and the
virtual neutral still receives the momentum difference; its recoil energy
is not subtracted a second time. This preserves the drawn nominal-threshold
outcome without a new angle veto. Energy is conserved outside this band;
the small recoil energy defect inside it is tested explicitly. The
nonrelativistic heavy-target recoil scale is :math:`2(m_e/M_i)E`; it is not
an exact arbitrary-energy upper bound for every conceivable transition.
For this model's finite rotational gaps the band occurs at low energies.
Unchanged outcomes use the existing stable elastic two-body formula and
retain recoil even when their internal-energy change is zero.

Where the equations are implemented
-----------------------------------

The source/preparation files below are under
``Tools/CrossSections/reciprocal_rotation``; runtime classes are under
``Source/Particles/Collision/BackgroundMCC``. This division makes it possible
to change a numerical encoding without changing the molecular source model.

.. list-table::
   :header-rows: 1
   :class: rotational-scattering-table
   :widths: 40 60

   * - File or function
     - Responsibility
   * - ``rotation_reference.py``
     - Constants, canonical levels, populations, tail bound and momentum factors.
   * - ``elmolcs_source.Source``
     - Elementary source selection, reduced threshold interpolation and matched Born total.
   * - ``prepare_jung.py``
     - Reproduction of inferred elementary fractions from the readable thermal-branch digitization.
   * - ``hybrid_reference.Priors``
     - N2 Read/Jung and Gote priors, O2 Born components, and energy-dependent bridges.
   * - ``o2_completion.Completion``
     - Positive relative-entropy fit to integral and momentum-transfer constraints.
   * - ``hybrid_reference.Hybrid.solve``
     - Backward energy sweep and implicit reciprocal inclusive normalization.
   * - ``Hybrid.interp`` and ``future_weights``
     - Shared primitive interpolation and reverse evaluation at shifted energy.
   * - ``rotor_basis.py`` and ``high_distribution.py``
     - Bound-rotor recoupling, momentum-transfer quadrature and phase-averaged tails.
   * - ``adaptive_low_grid.py`` and ``export.Exporter``
     - Sampling-grid refinement, angular CDFs, discrete supports and probability cells.
   * - ``text_bundle.py``
     - Readable numerical encoding and independently checked probability reconstruction.
   * - ``BackgroundMCCReciprocalText.H``
     - Host text parsing, probability reconstruction and one-time alias preparation.
   * - ``BackgroundMCCReciprocalRotation.cpp``
     - Validation, immutable storage sharing, range checks and device upload.
   * - ``BackgroundMCCReciprocalRotation.H``
     - Device rate interpolation, row mixture, angle and discrete outcome sampling.
   * - ``BackgroundMCCCollision.cpp``
     - Combined-family rate replacement, bounded majorant and integration with the MCC selector.
   * - ``BackgroundMCCElasticKinematics.H``
     - Signed internal-energy update, two-body recoil and nominal-threshold continuation.

Source calculations are not called by the device executor. The independent
validation entry points are mapped in :ref:`rotation-validation`.

.. _rotation-data-format:

Production files and initialization format
------------------------------------------

What belongs in each repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WarpX contains model equations, source readers, exporters, independent tests,
initialization and GPU sampling code. warpx-data contains numerical source
inputs and production collision distributions. Temporary reference
checkpoints, synthetic parser fixtures, validation outputs, timing records,
commit identifiers and benchmark results do not belong in warpx-data.

Each gas has an ``IAA`` directory with physical families in separate
subdirectories. ``elastic`` holds the residual integral and angular DCS, and
``rotational/state_to_state`` holds the elementary rotational source tables.
``rotational/sources`` holds source-specific numerical
constraints, including their units and interpretation. Each
``rotational/elastic_rotation_300K`` directory represents one combined
collision family,
not a collection of additional processes to add to the elastic rate.

Why the original files were large
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Binary inputs stored final device arrays directly. The N2 bundle contained 55,163,303
8-byte alias entries; O2 contained 31,779,555. Many neighboring angular cells
repeat closely related distributions over the same thermally recoupled
rotational transitions. Aliases therefore occupied 441,306,424 bytes for N2
and 254,236,440 bytes for O2. These are sampling acceleration arrays, not
that many measured cross sections. Including grids, outcomes and search
indices, the two binary payloads totaled 751,195,072 bytes (716.4 MiB).

The 15 N2 and nine O2 ``.bin`` files were arbitrary slices of a concatenated
payload, each limited to 32 MiB. A file boundary had no physical meaning.
Splitting kept individual files below Git's large-file limit; it did not
reduce total data size. Binary storage made startup simple but made source
review and interpretation unnecessarily difficult.

Readable probability inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The production data use decimal numerical tables and a compact factorization of the
conditional probability matrices. Initialization expands that representation
and builds the device aliases. The numerical factorization is an encoding
of a previously validated probability distribution, not a new molecular
model, and its matrix rank has no relationship to rotational tensor rank
:math:`L`. Every discrete outcome and its exact support are retained.

The index ``thermal_rotation.rot`` has the form::

   WARPX_RECIPROCAL_ROTATION
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
   :class: rotational-scattering-table
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
~~~~~~~~~~~~~~~~~~~~~~~~~~

``probabilities.txt`` starts with::

   WARPX_ROTATIONAL_PROBABILITIES number_of_cells number_of_alias_entries

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initialization additionally constructs 257 quantile-search bounds per energy
row. These bracket the full binary search without coarsening its grid. The
runtime also builds a small energy-bracketing index with 64 mantissa bins
per binary exponent. Cumulative sampling can add 65 uint16 search bounds
per probability cell, when the prepared table's memory budget permits it.
These indices change search ranges only; all final searches and physical
interpolation grids remain unchanged.

The runtime validates all dimensions, monotonicity, probability bounds, canonical
rotational levels, parity, molecular masses, range coverage, temperature and
threshold accessibility. It checks the combined source rate against the
ordinary elastic table above the cold join. No prepared file is written by
the simulation, and no Python or external numerical library is required to
read the production data. A restart prepares the immutable arrays again; repeated collision
objects within a process share them.

The 8-byte input/preparation alias layout remains float32 cutoff, uint16
alternate column and uint16 primary outcome identifier. After validation,
initialization replaces the alternate column in place with that column's
outcome identifier. Both uint16 device fields therefore address the outcome
palette directly. Prepared binary inputs retain their original column-index
format and remain readable. Energies and outcome labels remain
double precision. The per-process cache enforces a combined 1 GiB table
budget. One MPI rank per GPU shares these tables within that rank; separate
ranks on one GPU each allocate their own copies. Choosing cumulative storage
uses additional memory and is primarily a reference option.

Binary input remains supported with its little-endian binary64 host
requirement. Its directory lists ``name type count byte_offset`` followed by
relative binary-part filenames and byte counts. Offsets are eight-byte aligned
within the concatenated payload. Optional binary quantile indices are accepted;
without them the full angular/conditional search is used. The energy and
cumulative-cell indices are runtime additions and require no file changes.
Unmarked ordinary cross-section files retain their legacy endpoint behavior.
The readable representation does not change any other rotation option or
source continuation.
