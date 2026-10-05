Prepared rotational distributions, version 6
===========================================

``thermal_rotation.rot`` is a UTF-8 index. Payloads are little-endian arrays
in ``tables-NNN.bin`` files, each at most 32 MiB. There is no runtime scientific
generation, compression step, or dependency on elmolcs. The index begins with::

    WARPX_THERMAL_ROTATION_V6
    target reciprocal_hybrid representation
    temperature mass maximum_energy high_energy rutherford_energy separation screening_radius
    array_count part_count payload_bytes

``target`` is N2 or O2. ``representation`` is ``alias`` or ``cumulative``;
production files contain aliases. Units are K, kg, eV, eV, eV, Bohr radii, and
Bohr radii on the third line. Each array is described by::

    name type element_count byte_offset

The types are ``f64``, ``f32``, ``u32``, ``u16``, and ``alias``. Array offsets
refer to the concatenation of all payload parts and are aligned to eight
bytes. The directory is followed by ``part_count`` lines giving a relative
filename and its byte count. Every part except the last is exactly 32 MiB.
An alias entry occupies eight bytes: float32 cutoff, uint16 alternate local
index, uint16 outcome index. No commit identifiers or verification records
are stored in the production index.

Scientific arrays
-----------------

* ``energies`` and ``rates`` (f64) store collision energy in eV and K in m3/s.
  Energy starts at zero. ``coordinates`` (u32) marks each interval as linear
  (0) or square-root (1); the final entry is unused.
* ``angular_offsets`` (u32), ``angular_u`` (f64), ``deflection`` (f64), and
  ``changing`` (f32) describe the elastic inverse CDF and the probability of
  a rotational change on that row. Deflection is 1-cos(theta), not cos(theta).
* ``conditional_offsets`` (u32), ``conditional_u`` (f64), and
  ``conditional_cells`` (u32) describe the distributions conditioned on a
  change. Interpolation chooses one neighboring distribution probabilistically.
* ``cell_offsets`` (u32) index the sampling entries. Cells are nonempty and
  contain at most 65535 entries. ``outcomes`` (f64) contains pairs of signed
  electron energy loss and ground-relative initial neutral internal energy,
  both in eV. Outcome zero denotes unchanged rotation.
* ``high_edges`` (f64) and ``high_cells`` (u32) describe the momentum-transfer
  bank. Its dimensionless coordinate is k R sin(theta/2).
* Alias bundles contain ``aliases``. Offline cumulative reference bundles
  instead contain ``cdf`` (f64) and ``outcome_ids`` (u16).

The angular and conditional offset arrays contain one sentinel after the last
energy row. Cell offsets likewise have a final sentinel. Angular rows above
the Rutherford energy may be empty because that normalized angular law has
an analytic inverse. Conditional rows above the high-energy transition may
be empty because the momentum-transfer bank supplies their outcomes.

Below the high-energy transition, the same rate-weighted energy row supplies
both angle and outcome. An unchanged outcome is still a scattering event.
Canonical thresholds remain in double precision; energy labels are never
interpolated. The cold row contains the finite superelastic rate, so dividing
it by the zero incident speed to manufacture an ordinary cross section is
not meaningful.

Range and continuation contract
-------------------------------

The index's maximum energy bounds every collision query. Preparation may
evaluate the reciprocal forward primitive beyond this endpoint to obtain
reverse rates; that buffer is not a runtime extension. Out-of-range queries
fail. An eight-ParticleReal-epsilon relative tolerance accommodates endpoint
roundoff and is not an extrapolation interval.

New ordinary source exports declare ``energy_min_eV``, ``energy_max_eV``, and
``outside_energy_range = error`` in comments. Background MCC intersects these
ranges at initialization, checks them before either selector path, and bounds
the collision frequency only within the resulting interval. Known zeros below
an excitation threshold are part of its declared range. Unmarked inputs retain
the existing endpoint-clamping convention.
