# RZ proton ionization and electron MCC in air

This example targets the APIs on `codex/proton-impact-ionization-development-sync`,
starting from commit `c9b6cab3b`. It illustrates an axisymmetric 800 MeV proton
bunch in a prescribed N2/O2 background. Its mesh, time step, particle weights,
domain and run length are starting choices, not a converged LANSCE BPM model.

## What the two operators calculate

For each molecule `A = N2` or `O2`, the proton operator represents

```text
p + A -> p + e + A+ (effective electron-production pair)
```

The calibrated Porter–Jackman–Green-type model supplies the differential
electron spectrum and its integral, so **no proton cross-section file is
needed**. The implementation tabulates these on the CPU at initialization,
copies the tables to device memory, and samples them by interpolation at run
time. It includes a soft molecular contribution, a relativistic hard-electron
tail, and an approximate angular model. It can emit electrons above the free
electron binary-collision edge because the initial electron is bound.

The expected physical pair weight added to one cell in a collision interval is

```text
Delta W = n_A * dt_collision * Z_projectile^2
          * sum_over_projectiles(w_p * sigma_A(E_p) * v_p)
```

There is no extra cell-volume factor: the projectile weights already count
physical particles. In RZ they represent particles in annular volumes.
Fractional weight is retained in a checkpointed per-cell field. Pairs normally
have `fixed_product_weight`; if `max_products_per_cell` limits their number,
their weights increase to retain the entire source. This is a weighted source
with stratified sampling, rather than a Poisson count of individual collisions.
Parents are selected in proportion to `w_p*sigma*v`; products inherit their
positions. Ions get the neutral thermal velocity distribution.

The source changes neither proton momentum nor proton weight and does not
deplete the gas. Its observable is inclusive electron production, represented
by singly charged effective molecular ions. It does not resolve molecular
fragmentation, multiple charge states, electron capture or ion recoil.
The numerical proton interval is 5 keV–10 GeV, but the total-yield calibration
uses 5–4000 keV data. The measured spectral calibration extends to 1 MeV for
N2 and 150 keV for O2. **The 800 MeV molecular spectrum is an extrapolation.**
See the [full model derivation](../../../Docs/source/theory/multiphysics/proton_impact_ionization.rst).

The created electrons then collide with the prescribed gas through
`picmi.MCCCollisions`. For a sampled neutral velocity, the total rate is
`nu(E) = n_A*v_rel*sum_j sigma_j(E)`. A null-collision bound first rejects
unlikely attempts. The accepted one-event probability is
`-expm1(-nu*dt_collision)`, and the chosen channel has probability
`sigma_j/sum_j sigma_j`. The implementation uses an exact union of all
piecewise-linear energy grids for efficient channel selection, with a bounded
memory fallback for unusually large channel sets.

| Electron process | Result in this branch | Supplied data |
| --- | --- | --- |
| Elastic | Changes direction, with finite target recoil | Integral elastic cross section; DCS if using `IAA` |
| Excitation | Removes the channel's discrete energy loss and scatters | One cross section and threshold per channel; optional channel DCS for `IAA` |
| Ionization | Changes the incident electron and creates another electron plus a positive ion | Total ionization cross section and product species |
| Dissociative attachment | Removes the electron and creates O− | O2 attachment cross section in m² |
| Three-body attachment | Removes the electron and creates O2− | Raw cross section in m⁵ plus third-body density in m⁻³ |

`energy_sharing_model="RBEQ"` uses built-in N2/O2 shell weights and conditional
secondary-energy distributions. **It does not replace the supplied ionization
rate table.** `scattering_angle_model="IAA"` selects the corresponding
ionization angular closure and needs no ionization DCS file. Recoil kinematics
checks energy and momentum; inadmissible draws become null events. The equal
sharing default is a less detailed alternative. All electron products inherit
the incident macroparticle's weight, including a proton-source weight enlarged
by its per-cell cap.

For IAA elastic scattering the DCS supplies the angular shape, while the
ordinary integral table supplies the event rate. The built-in screened
Rutherford angular continuation starts at 10 keV for recognized N2/O2 elmolcs
metadata. The neutral gas is a classical Maxwellian; collision-rate lookup
uses the documented approximate relative proper velocity. Excitation and
attachment do not evolve a molecular chemistry model: excited neutrals and
dissociation fragments are not tracked, and attached ions inherit neutral
velocities rather than a resolved dissociation recoil distribution.

## Files to supply

Use the [Python input](inputs_rz_proton_beam_air_picmi.py) with a copy of
[cross_sections.example.json](cross_sections.example.json). The JSON is only
a convenience used by this example; the same dictionaries can be written
directly in your Python input. All paths are resolved relative to the JSON.

Supply these evaluated data:

1. **Electron elastic integral cross sections for N2 and O2.** Use true integral
   elastic cross sections with an explicit angular DCS. A momentum-transfer
   cross section is not the same collision frequency. An effective or total
   cross section may already include inelastic contributions; adding those
   channels again would double count them.
2. **Electron excitation cross sections**, one file per rotational,
   vibrational and electronic loss you include. The manifest shows only the
   first vibrational excitation to illustrate syntax; it is not a complete
   air set. Add each retained channel under a unique `excitation_*` name,
   with its actual threshold. The illustrative thresholds are 0.291 eV for
   N2 and 0.196 eV for O2; use values matching your selected dataset.
3. **Electron ionization totals for N2 and O2.** The RBEQ threshold parameter
   must match 15.58 eV for N2 and 12.07 eV for O2, within 0.05 eV. Use one
   consistent molecular total with the built-in RBEQ shell sampler; do not
   attach that sampler indiscriminately to shell-resolved or fragment-specific
   partial rates, or add both a total and its partials.
4. **O2 dissociative attachment and three-body attachment tables**, as separate
   channels and distinct O−/O2− product species. The example's
   `third_body_density: null` must be replaced by a positive number before
   running. Determine the collider or mixture from the source data: use its
   partial density for O2 or N2, or a justified effective density for air.
   The example deliberately does not assume that the total gas density is
   the right collider density. Recalculate it if pressure or temperature
   changes. With collider-specific tables, use separate `attachment_*`
   entries, each with its own density. A rate coefficient in m⁶/s is not a
   cross section in m⁵ and cannot be supplied unchanged.
5. **Elastic DCS files** if selecting `IAA`: elmolcs `DCS.e-N2` and `DCS.e-O2`
   are supported. Each numeric row has energy in eV and equally spaced
   angular values from 0 to 180 degrees (normally 361 values at 0.5-degree
   intervals). Keep the `SPECIES: e / N2` or `O2` metadata. IAA excitation
   needs a suitable excitation DCS; reusing the elastic DCS is an additional
   approximation. The manifest uses isotropic excitation to make that choice
   explicit.

Ordinary cross-section files contain exactly two numeric columns:

```text
# electron energy [eV]     cross section [m^2]
```

Use m⁵ instead for an explicitly declared three-body attachment table. Blank
lines and `#` comments are allowed. Energies must be finite, nonnegative and
strictly increasing; values must be finite and nonnegative. Loss processes
must be zero at and below their threshold, including an explicit zero at the
threshold when needed to keep interpolation from opening the channel early.
Do not pass an entire LXCat multi-process download, a CSV with headers, or an
elmolcs multi-column elastic export directly to this reader. Convert each
selected channel and retain its source, units and threshold in your records.

The existing local resources include
`resources/elmolcs/Data/dcs/{N2,O2}/DCS.e-{N2,O2}` and
`resources/O2_attachment.txt`. The latter is a multi-block LXCat export: its
first attachment block is explicitly m⁵ and its other block is dissociative
attachment in m². Extract the two numerical tables separately. You do not
need the elmolcs Python package, optical oscillator files, proton fit archives,
or an external ionization-energy CDF at simulation run time once the electron
tables have been exported. The local elmolcs `writewarpx` helper should not be
used unchanged for this setup: its elastic export can have headers/multiple
columns, and it predates this branch's attachment and RBEQ/IAA configuration.

The [elmolcs author's package documentation](https://pypi.org/project/elmolcs/0.0.3/)
describes the IAA integral/elastic-angular data and analytical high-energy
extensions. [LXCat's MuroranIT description](https://fr.lxcat.net/Phelps)
identifies the source and validation approach for its electron cross sections;
retain the references for the chosen species and download.

## Energy coverage and timestep

At 800 MeV the free-electron transfer edge is approximately **2.48074 MeV**.
This is not a strict cutoff for the branch's bound-electron source, and fields
can accelerate electrons further. A table ending at 100 eV or 10 keV does
not describe all of these electrons. WarpX holds the terminal cross section
constant outside the tabulated interval. Supply a physically justified,
resolved high-energy continuation and check the electron-energy diagnostic;
do not create a fictitious abrupt zero just to reduce the automatic bound.
Resolve low-energy resonances and attachment as carefully as the energetic
tail. RBEQ's sampling table is initialized using the last supplied ionization
energy, which makes adequate table coverage important for both rate and shape.

Use a small **physical** optical depth per collision interval, for example
`nu_physical*dt_collision <= 0.05–0.1` as a starting point, then halve it and
check observables. The automatic `nu_max` conservatively covers even the
constant high-energy table extension; it can be much larger than rates of
occupied electron states. Its timestep warning is conservative. A loose
bound costs attempts, but the corrected event probability does not require
using that loose bound as the physical convergence scale.

`ndt_subcycle=N` divides the collision interval by N. It does not refine the
field solver or move particles between those substeps. This branch completes
all substeps of one collision object before the next object. Thus increasing
only the N2 and O2 subcycle counts does not remove splitting error between
the gases; decrease the PIC timestep too and check sensitivity to collision
ordering. Likewise, electrons born during a proton source call are processed
by later MCC objects as if present for that collision interval. This source
timing approximation requires timestep convergence. Avoid collision
supercycling for rapidly changing 100 ps bunches.

At atmospheric density, cold-electron collision times can be extremely short.
Choose the actual time step from the supplied rate tables, field CFL,
electron plasma frequency and spatial transport scales. The example uses
0.1 ps only as a setup value. Check Debye-length and grid-heating sensitivity;
do not infer a resolved model from its ability to run.

## RZ and the Roy beam

[Roy et al., AIP Advances 10, 095023 (2020)](https://doi.org/10.1063/5.0021497)
describe an 800 MeV beam in a 1.2 m air gap. Their estimate uses a 0.6 A peak,
a 100 ps (4σ) microbunch, approximately 5 ns bunch spacing, and a 5 mm beam
radius. The example uses a Gaussian longitudinal bunch with σt = 25 ps,
a uniform transverse disk, and a short 10 cm domain. Its Gaussian charge is
`I_peak*sqrt(2*pi)*sigma_t`; it is not the paper's rectangular estimate
`I_peak*(100 ps)`. Use measured bunch charge/profile to choose the normalization.
The 20 ps default run is a short startup example, not the complete bunch train.

Use `n_A = x_A*p/(k_B*T)`. The example takes dry-air fractions 0.78084 N2 and
0.20946 O2 without folding omitted argon into either species. At 101325 Pa
and 293.15 K these give approximately 1.955e25 and 5.244e24 m⁻³. That pressure
is an illustrative sea-level value; use the measured LANSCE gas conditions.
If you instead intend a binary 79/21 mixture, set those fractions explicitly.
Neutral N2/O2 are background fields, so there are no neutral macroparticle
species to load. Electron, positive-ion and negative-ion species start empty.

The input uses `u_z = gamma*v_z`, as required for the PICMI directed velocity
in this implementation; it is approximately 4.675e8 m/s at 800 MeV and can
exceed c because it is proper velocity. Physical beam speed is approximately
0.842c. The example deposits beam charge/current and transports it ballistically.

Set `n_azimuthal_modes=1`. The collision handler temporarily rotates momenta
into the local curvilinear frame and restores them afterward. RZ resolves
axisymmetric plasma response but cannot represent the four stripline
electrodes, their differential signals, or transverse beam offset. A
quantitative electrode/BPM comparison needs a suitable 3D geometry and
electromagnetic boundary/circuit treatment. The example uses illustrative
PEC walls at the outer radius and z ends, compatible with its initial beam
self-field Poisson solve. Move the boundaries farther away and test their
influence before interpreting the result. A longer open-domain run needs
an appropriate beam injection/initial-field treatment: RZ does not support
z-directed PML, and its Yee Silver–Mueller boundaries cannot be combined
with this multigrid self-field initialization in the current branch.

The paper's approximately 34 eV mean energy per ion pair includes the full
energy-degradation cascade. It is not an N2/O2 ionization threshold and must
not replace the 15.58/12.07 eV losses. Do not add an extra stopping-power/W-value
electron source on top of explicit proton ionization plus secondary electron
MCC; that would double count production. Validate total pair yield and energy
deposition against stopping/W-value estimates as integrated observables.

## Running and improving the calculation

Build the branch with Python and RZ support, then use its installed Python
package or the build's `lib/site-packages` directory:

```sh
python inputs_rz_proton_beam_air_picmi.py \
  --cross-sections /path/to/your/cross_sections.json \
  --pressure-pa 101325 --temperature-k 293.15 \
  --dt 1e-13 --steps 200 --mcc-subcycles 1
```

Add `--write-input /path/to/inputs` to write the corresponding native input
instead of evolving. `--product-weight` controls source macroparticle weight;
the cap is set in the Python input. Test convergence in both. RZ cell volume
is `pi*(r_hi^2-r_lo^2)*dz`, so a single constant weight does not give constant
particles per cell across radius. The cap preserves total charge but can
leave a small number of heavy particles carrying the energetic tail.

For GPU runs, retain the compiled collision operators, group channels for a
given target in one MCC object, and keep table preprocessing out of timestep
callbacks. Table interpolation and product creation run on the accelerator.
Profile the actual GPU before selecting more intrusive allocation, scan or
kernel changes: register pressure, sparse event rates, atomic contention and
per-tile synchronization affect the tradeoffs. No CUDA/HIP/SYCL performance
measurement is implied by the CPU validation of this example.

For quantitative long-time air chemistry, additional models/data may be
needed: electron–ion/ion–ion recombination, detachment, ion-neutral transport,
neutral excitation/dissociation populations, thermal rotational superelastic
collisions, gas evolution, humidity and trace species, and wall/circuit
processes. Ordinary excitation in this branch requires a nonnegative loss;
it cannot represent superelastic energy gain by supplying a negative loss.
The two collision operators alone are not a complete air-plasma chemistry
or beam-stopping model.
