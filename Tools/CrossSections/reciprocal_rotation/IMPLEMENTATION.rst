Implementation contract
=======================

Scope
-----

Implement one optional combined vibrationally elastic/rotational family per
gas. Use a new selector, for example ``reciprocal_hybrid``, rather than silently
changing the existing ``elastic_dcs`` or ``iaa_spectator`` meanings. The target
is the current ``codex/rigid-beam-immobile-ions-development-sync`` branch.
No pull request is authorized.

The prescribed neutral rotational bath is fixed. Neutral rotational states
are neither particles nor evolved variables. Initial-state sums and all
angular-momentum recoupling are offline. A virtual initial/final pair in an
outcome palette provides its exact signed energy gap and, when needed,
initial internal energy for recoil kinematics.

Complete source model
---------------------

The following energy intervals refer to electron kinetic energy, not the
neutral temperature.

N2:

* Up to 1 eV: retain the real elmolcs elementary rates and canonical threshold
  factors; use an isotropic rotational DCS as the leading low-energy angular
  approximation. Fill the unchanged part from the inclusive constraint.
* 1--1.25 eV: smooth transition to the resonance construction.
* 1.25--4 eV: Kutz--Meyer energy strengths with the Read tensors and a direct
  background constrained by Jung's vibrationally elastic branches. The
  missing-angle rank-2/rank-4 corrections approach 0.52823198 and 0.84013224;
  they are applied to forward kernels and their reverse partners together.
* 4--10 eV: smooth positive transition to the 10 eV Gote angular constraints.
* 10--200 eV: interpolate the completed Gote transfer-rank distributions;
  use the two-centre prior to extend outside the measured angular range.
* 200--1000 eV: continue the 200 eV correction at equal momentum transfer and
  blend to the full spectator distribution with a smooth step in log energy.
* Above 1000 eV: use the bounded-rotor two-centre spectator continuation.

O2:

* Up to 1 eV: use IAA Eq. 11.21b with the actual transition momentum transfer
  for every N->N+2 excitation, not a rescaled ground-state curve.
* 1--20 eV: smooth positive bridge from the Born construction to the
  moment-constrained 20 eV model. This bridge is an interpolation model,
  not a measured resonance DCS.
* 20--200 eV: preserve Bhattacharyya's 1->3 integral and momentum-transfer
  cross sections. Infer the unresolved higher-rank budget from its totals,
  distribute it using the spectator prior, and fit nonnegative angular
  distributions under the IAA inclusive angular marginal.
* 200--1000 eV and above: use the same equal-momentum-transfer transition
  and bounded spectator continuation as for N2, with O2 parameters.

The O2 rank completion uses the identities
``sigma_11 = a_0 + (2/5) a_2`` and
``sigma_13 = (3/5) a_2 + (4/9) a_4``. The remaining inclusive strength is
``(5/9) a_4 + sum(a_L, L>=6)``. The analogous identities hold for first angular
moments. Higher ranks are not arbitrarily set to zero. The angular fit is the
positive relative-entropy projection of the stated prior onto the integral
and momentum-transfer constraints, with the IAA angular marginal fixed.

State space and endpoints
-------------------------

Use ground-relative rigid-rotor levels ``epsilon_J = B*(J*(J+1)-Jg*(Jg+1))``.
N2 has even:odd nuclear-spin weights 6:3; O2 has only odd nuclear N levels in
this unresolved-spin model. Generate initial populations until both their
omitted mass and transfer moments converge. The population criterion is
1e-10. At 300 K the conservative state limits used here are J=56 and N=64.

The fast-electron continuation is conditioned on the bound, ground-vibrational
rotor manifold. Model levels at or above the supplied dissociation limit are
excluded before normalization, in both directions. The canonical constants
give final limits J=197 for N2 and N=167 for O2. This prevents unlimited
rigid-rotor excitation; it is a declared high-J closure, not validated
non-rigid molecular spectroscopy.

Keep the physically finite zero-energy superelastic rate. Below 1 meV,
retain the reciprocal changing kernels and join a positive unchanged
background to its s-wave limit. Do not force that rate into a finite
zero-energy cross-section entry. The unchanged rate vanishes at rest.
Use the analytic zero-energy limit and sqrt(E) behavior; linear extrapolation
of a positive-energy primitive to zero is incorrect.

Retain the inclusive elastic rate and normalized IAA DCS in its supported
range. Generate a smooth DCS transition from 8 to 10 keV to the existing
screened-Rutherford continuation. The original hard 10 keV switch must not
be interpolated through as though it were continuous. The model is defined
through 1 GeV, with an explicit analytic guard band for reverse evaluations.

Reciprocity and inclusive normalization
--------------------------------------

For an excitation i->f with gap Delta, use a shared positive primitive X:

``D_if(E,theta) = p(E-Delta)/p(E)^2 * X_if(E,theta)``

``D_fi(E-Delta,theta) = (g_i/g_f)/p(E-Delta) * X_if(E,theta)``.

Then ``g_i*p(E)^2*D_if = g_f*p(E-Delta)^2*D_fi``. The same interpolation of
X must be used in both directions. For generic rotational tensors X_if is a
positive Clebsch--Gordan sum; low-energy O2 additionally uses explicit
transition-specific Born components.

At low energy, protect trustworthy changing rates and obtain a nonnegative
unchanged remainder. Elsewhere solve the coupled forward/reverse inclusive
normalization. The reverse term at energy E uses the corrected forward
primitive at E+Delta. Never independently normalize incoming-state rows
after constructing reverse rates. Never clip a negative physical remainder.
A failed positive decomposition blocks the bundle and retains diagnostic
information outside warpx-data.

The reference satisfies differential reciprocity in the stated heavy-target
model. Finite sampling tables approximate it; report their measured errors
rather than claiming machine-exact reciprocity after arbitrary interpolation.

Offline tables and runtime
--------------------------

Use a new versioned bundle containing fixed-temperature, precomputed data:

* target, bath temperature, canonical constants, units, validity limits;
* one aggregate rate coefficient K=v*sigma, including its zero-energy limit;
* inverse angular CDFs and angular coordinates for the low-energy rows;
* the unchanged/change probability and conditional changing-outcome aliases;
* an outcome palette with discrete signed gaps and recoil metadata;
* the compact momentum-transfer alias bank for the high-energy range.

All population sums, source evaluation, special functions, normalization,
CDF inversion, and alias construction occur offline. WarpX initializes by
reading, validating, sharing, and uploading immutable arrays. Production
cross-section data go in warpx-data. Synthetic fixtures, validation outputs,
benchmark records, and commit hashes do not go there.

For low-energy sampling:

1. Find the shared energy interval once. Use a sqrt coordinate in the first
   interval above a physical threshold, and at zero; otherwise use the
   validated interpolation coordinate.
2. Interpolate the joint rates and select the endpoint row with its
   rate-weighted probability.
3. Draw the angle from that row's inverse elastic CDF and retain the row-local
   angular quantile. The original draw used to select a mixture row is not
   automatically the global angular quantile.
4. Use that same row/coordinate to select unchanged or changing rotation and
   a constant-time alias draw of a discrete outcome. Conditional angular
   spectra may be shared or coarsened within the tested error budget.
5. Apply the selected angle and signed-loss recoil kinematics once.

Unchanged rotation is still an elastic scattering event. It must not become
a null collision or skip recoil. Never interpolate discrete energy labels.
Every active threshold must remain an interval boundary so row mixtures
cannot select a forbidden excitation. Distinguish threshold-roundoff handling
from physical rate renormalization.

Above 1 keV, evaluate the momentum-transfer coordinate from the actual energy
and sampled deflection, then use the precomputed q-bin spectrum. Its entire
discrete distribution is stored. Averaging probabilities in a bin preserves
energy diffusion; replacing them by a mean loss does not. Resolve the q tail
with probability phase averaging where appropriate, and audit its angular
as well as energy moments.

A packed alias entry can use a float32 cutoff plus uint16 local-alias and
outcome indices (8 bytes); validate all limits before conversion. Keep energy
labels, kinetic intermediates, and forward deflections in double precision.
Use sufficient uniform-random resolution for rare tails, including in a
single-precision build. Store or carry 1-cos(theta) directly and compute
sin(theta) from sqrt(d*(2-d)) to retain small forward deflections.

Accuracy and validation
-----------------------

The export pipeline must verify rates and first/second transfer moments below
0.2% error, and equilibrium power imbalance below 0.1% of heating plus cooling.
Use separate budgets for source-grid, energy-row, angular, compression, and
storage errors. Include angle-weighted moments and rare-tail absolute bounds.
A grid tested only against interpolation of itself is not an independent
source check; solve intermediate source rows and refine where necessary.

Source tests include the Jung angular/integral constraints, Gote rank coverage,
Bhattacharyya integral and momentum-transfer constraints, canonical thresholds,
degeneracies, reciprocity, and population convergence. The working reference
has explicit zero-temperature and zero-energy checks in addition to the
250/300/350 K bath checks. Numerical precision does not establish experimental
accuracy: low-energy O2 resonance physics, unresolved spin structure, and
high-J centrifugal/vibrational effects remain physical uncertainties.

Integration must add independent particle tests for cold heating, hot cooling,
Maxwellian equilibrium, unequal rotational/translational temperatures,
thresholds, signed recoil conservation, and angle-energy correlations. Re-run
complete MCC and timestep benchmarks with the physical bundles on Perlmutter;
the supplied CPU/CUDA benchmark measures only table lookup. Keep CUDA/HIP/SYCL
portability, shared tables, and existing ordinary-process selection central.

Deliver the exporter/source checks, bundle reader, sampler/kinematics hookup,
then integration tests, benchmarks, and documentation in normal-sized commits.
Push validated work as implementation proceeds. Do not open a PR without an
explicit user instruction.
