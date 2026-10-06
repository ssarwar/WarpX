.. _rotation-sources:

Source evidence and gas-specific angular kernels
================================================

Source hierarchy and reproducibility
------------------------------------

The common source is A. Schmalzried, *Electron Thermal Runaway in Atmospheric
Electrified Gases: A Microscopic Approach* (Universidad de Granada, 2023),
:cite:`rot-Schmalzried2023`, especially Eq. 2.48 (printed p. 60), Sections
11.1--11.2 (including Eqs. 11.10, 11.16, 11.21--11.24 and 11.29--11.35),
the physical discussion in Section 12.1, and the database and comparison
chapters 15--16. Printed page numbers differ
from PDF page indices. The elmolcs source snapshot and its IAA data provide
numerical curves; the specific papers below constrain their interpretation.

Numerical inputs belong to warpx-data under ``MCC_cross_sections/N2/IAA`` and
``MCC_cross_sections/O2/IAA``. The parent directories contain readable
integral and elastic DCS tables. ``reciprocal_sources`` contains the published
or processed angular constraints, with source-specific README files.
``reciprocal_hybrid_300K`` contains the prepared joint distributions.
The WarpX reference implementation is in ``Tools/CrossSections/reciprocal_rotation``.
Changing a model assumption requires regenerating and revalidating the
joint distributions; editing an unrelated source table does not update them.

The numerical ``iaa`` N2 rotational tables with a resonance near 0.8 eV are
not used: their energy dependence disagrees with the elementary ``iaa*``
tables and the cited resonance near 2.3 eV. The latter elementary curves
supply the N2 strengths. The :math:`0\to6` source header threshold is also
inconsistent with :math:`42B`; thresholds always come from
:eq:`rotation-levels`. The O2 numerical :math:`1\to3` Born integral remains
a source cross section, but the low-energy hybrid evaluates the
transition-specific Born formula instead of scaling that one ground-state
curve to arbitrary :math:`N`.

Elastic data provenance and its low-energy limitation
----------------------------------------------------

The imported angular grids are elmolcs ``Data/dcs/N2/DCS.e-N2`` (source
header dated 22 September 2022) and ``Data/dcs/O2/DCS.e-O2`` (12 December
2022). Their native grid has 361 angles at 0.5-degree spacing. We import the
resulting DCS rather than rerunning the package's scattering-potential fits.
The native source headers record these constructions:

* N2: below 30 eV, least-squares fits to the experimental DCS collection;
  below 0.1 eV, a modified effective-range estimate with scattering length
  0.44 Bohr. Between 4 and 8 eV the backward endpoint uses additional
  assistance from Sun et al. (1995), Table VII, with a stated 15% uncertainty.
  At 30--150 eV the source uses angular-momentum-coupled potential calculations;
  at 200--1000 eV it uses the independent-atom model with energy-dependent
  scaling. Above 1 keV the source removes its correlation potential and uses
  Buckingham polarization; above 8 keV its source switches to a Born treatment.
* O2: modified effective-range behavior below 1 eV, using scattering length
  0.30 Bohr, polarizability and permanent quadrupole. Below 10 eV the source
  names Sullivan et al. (1995), Green et al. (1997) and Linert et al. (2004);
  from 15 eV it adds Woste et al. (1995), Trajmar et al. (1971) and Shyn and
  Sharp (1982). Angular-coupled potential calculations enter from 30 eV,
  a scaled independent-atom construction from 200 eV, and the source's Born
  treatment from 8 keV. See thesis Chapters 11 and 13 for the potential and
  fitting details and their experimental bibliography.

These describe the imported grid's provenance, not additional runtime switches.
The production angular continuation replaces the last source regime by the
explicit smooth 8--10 keV Rutherford join in :ref:`rotation-continuations`.
The residual integral's normalization follows thesis Section 11.1.5 and the
three-step prescription on printed p. 566: determine total and other inelastic
cross sections, subtract to obtain the vibrationally elastic residual, and
use the normalized angular shape to obtain transport moments.

Printed p. 567 explicitly identifies limitations below 1 eV: reliance on total
rather than swarm-inferred momentum-transfer cross sections, a basic modified
effective-range angular model, and sparse low-energy source knots (1, 10 and
100 meV). Refining our interpolation does not create missing experimental
information between those knots. Thermal-balance tests check consistency of
the collision operator, not an independent measurement of the true mobility
or diffusion coefficient at room temperature. The supplied header descriptions
are preserved in the warpx-data documentation so that these limitations are
visible without recovering the original archive.

What the papers establish
-------------------------

* H. Kutz and H.-D. Meyer, *Rotational excitation of N2 and Cl2 molecules by
  electron impact in the energy range 0.01--1000 eV: Investigation of
  excitation mechanisms*, Physical Review A **51**, 3819 (1995),
  `doi:10.1103/PhysRevA.51.3819 <https://doi.org/10.1103/PhysRevA.51.3819>`_.
  The elementary N2 curves in elmolcs were digitized from Fig. 7a. The
  rotationally unchanged :math:`0\to0` calculation is a separate
  fixed-initial-state theoretical cross section, not the IAA residual
  minus rotation and not the room-temperature unchanged cross section.
* F. H. Read and D. Andrick, *Angular distributions for the excitation of
  rotational states by resonant electron molecule reactions*, Journal of
  Physics B **4**, 911 (1971). Their angular tensors follow from the angular
  momentum coupling of incident/captured and emitted partial waves of a
  short-lived molecular resonance. The N2 example is the approximately
  2.3 eV resonance. These are resonance shapes, not a universal elastic DCS
  or a quantitatively established extension down to thermal energies.
* K. Jung, Th. Antoni, R. Muller, K.-H. Kochem and H. Ehrhardt,
  *Rotational excitation of N2, CO and H2O by low-energy electron collisions*,
  Journal of Physics B **15**, 3535--3555 (1982). Fig. 5 and Table 1 constrain
  N2 vibrationally elastic branches at 2.22 and 2.47 eV, with a
  :math:`500\pm30` K beam and angles 15--105 degrees. The 10--18 meV
  experimental energy width limits the information in individual lines.
  Reverse branches were constrained in the spectral analysis; they are not
  independent measurements of detailed balance.
* M. A. Morrison, W. Sun, W. A. Isaacs and W. K. Trail,
  *Ultrasimple calculation of very-low-energy momentum-transfer and
  rotational-excitation cross sections: e-N2 scattering*, Physical Review A
  **55**, 2786 (1997),
  `doi:10.1103/PhysRevA.55.2786 <https://doi.org/10.1103/PhysRevA.55.2786>`_.
  Modified effective-range theory and threshold corrections support treating
  the very-low-energy region separately from the spectator approximation.
  The inspected transcription of Eqs. 24--25 and Table III did not reproduce
  Figs. 5--6 when evaluated literally in eV. This is an unresolved reproduction
  problem, not a demonstrated error in the paper. No speculative coefficient
  repair is used; the supplied elmolcs/Itikawa low-energy N2 rates are retained.
* M. Gote and H. Ehrhardt, *Rotational excitation of diatomic molecules at
  intermediate energies: absolute differential state-to-state transition
  cross sections for electron scattering from N2, Cl2, CO and HCl*, Journal
  of Physics B **28**, 3957--3986 (1995). Table 1 supplies N2 transfer-rank
  fractions at 10--200 eV and 10--160 degrees. The data support large
  angular-momentum transfer at intermediate energies; a rank-two-only model
  misses that physics.
* P. K. Bhattacharyya and K. K. Goswami, *Elastic and rotational excitation
  of the oxygen molecule by intermediate-energy electrons*, Physical Review A
  **28**, 713 (1983),
  `doi:10.1103/PhysRevA.28.713 <https://doi.org/10.1103/PhysRevA.28.713>`_.
  Tables II--III, potential B with cutoff radius :math:`2a_0`, give
  :math:`1\to1`, :math:`1\to3` and inclusive integral and momentum-transfer
  cross sections. They are calculations with a Glauber-type eikonal
  amplitude, adiabatic nuclei, polarization and no electron exchange, at
  20--200 eV. They do not uniquely determine every O2 rotational DCS.

The additional N2 low-energy :math:`0\to2` integral values are the Itikawa
(2006) Table 5 entries included by elmolcs. O2's long-range formula is the
quadrupole/anisotropic-polarization approximation of Takayanagi and Itikawa
(1970), Eq. 38, as given in thesis Eq. 11.21b. None of these data establishes
an accurate O2 resonance model throughout 1--20 eV.

Angular-momentum recoupling
---------------------------

The sudden/adiabatic-nuclei approximation freezes the molecule's orientation
during the fast collision. Expanding the orientation dependence in spherical
tensors and integrating products of rotational spherical harmonics produces
Gaunt integrals, or equivalently Clebsch--Gordan coefficients. An elementary
rank :math:`L` contributes to arbitrary :math:`i\to f` through

.. math::
   :label: rotation-coupling

   C_{ifL}=|\langle i0,L0\mid f0\rangle|^2
      =(2f+1)\begin{pmatrix}i&L&f\\0&0&0\end{pmatrix}^{\!2}.

Only :math:`|i-L|\le f\le i+L` and even :math:`i+L+f` contribute.
For even :math:`i+L+f=2s`, the factorial identity used in the reference is

.. math::

   C_{ifL}=(2f+1)
   \left[\frac{s!}{(s-i)!(s-L)!(s-f)!}\right]^2
   \frac{(2s-2i)!(2s-2L)!(2s-2f)!}{(2s+1)!}.

Logarithmic gamma functions evaluate this identity without factorial overflow.
The identity :math:`(2i+1)C_{ifL}=(2f+1)C_{fiL}` underlies reverse recoupling.
For homonuclear gases the allowed transitions preserve nuclear-spin parity,
so the extra nuclear-spin weight cancels within a pair. Importantly,
:math:`i=f` receives contributions from :math:`L>0` as well as :math:`L=0`.
Thus copying :math:`\sigma_{00}` to every :math:`\sigma_{ii}` is incorrect.

If :math:`A_L(E,\theta)` is an elementary reduced angular strength, the
forward primitive before inclusive reconciliation is
:math:`P(E)\sum_L C_{ifL}A_L(E,\theta)`. An additional outgoing/incoming
momentum factor supplies the threshold opening. This extends a sudden model
near threshold; it is not an exact coupled-channel threshold calculation.
N2 source excitation strengths are reduced using
:math:`A_L(E)=\sigma_{0L}(E)P(E)/P(E-\Delta_{0L})`. Positive source knots are
preserved with log--log interpolation in that reduced quantity; its first
value is held between the canonical threshold and the first positive datum.

N2 angular construction
-----------------------

Below 1 eV, ranks 2, 4 and 6 use the real elementary strengths divided by
:math:`4\pi`. The isotropic rotational DCS is the leading low-energy
approximation, consistent with a quadrupole contribution that is much less
forward-peaked than dipole scattering. The unchanged component is solved
from the inclusive constraint, so the *total* DCS need not be isotropic.

For the resonance, the normalized Read shapes used to complete the missing
angular range are, with :math:`\mu=\cos\theta`,

.. math::
   :label: rotation-read-shapes

   h_0&=\frac{5(3\mu^2-1)^2}{16\pi},\\
   h_2&=\frac{5(9\mu^4-9\mu^2+4)}{56\pi},\\
   h_4&=\frac{5(\mu^2+3)^2}{224\pi}.

Each integrates to unity over solid angle. They describe the rank-zero,
rank-two and rank-four angular tensors of the short-lived d-wave resonance.
They alone cannot represent the direct elastic background. The JSON
``n2_jung.json`` therefore contains inferred nonnegative elementary fractions,
obtained by recoupling the digitized thermal branches at seven angles, rather
than labeling those fractions as raw measurements.

Between 15 and 105 degrees the inferred fractions are interpolated linearly
in angle. Outside this range each is continued by the corresponding
:math:`h_L(\theta)/h_L(\theta_{\rm edge})`, then the fractions are normalized.
Let :math:`s(x)=u^2(3-2u)`, :math:`u=\min(1,\max(0,x))`. In the unmeasured
region the rank-two and rank-four fractions are multiplied by

.. math::

   1-(1-a_L)\left[s\!\left(\frac{15^\circ-\theta}{15^\circ}\right)
                      +s\!\left(\frac{\theta-105^\circ}{15^\circ}\right)\right],
   \qquad (a_2,a_4)=(0.52823198,0.84013224).

The removed fraction is assigned to rank zero. These factors are fitted
missing-angle corrections to the absolute branch sums at 2.47 eV and 500 K;
they are not constants from Read's theory. The reference uses the same
corrections at 2.22 eV. The source check compares all individual upward and
downward branches, allowing 5.5% relative disagreement with Jung's inferred
branch fractions. That physical/source agreement is separate from the 0.2%
numerical discretization target.

At each anchor energy the corrected fractions multiply the IAA DCS and are
divided by the elementary Kutz--Meyer strength. These angular multipliers
are linearly blended between 2.22 and 2.47 eV, held outside those anchors,
and multiplied by the energy-dependent elementary strengths. Rank six retains
the isotropic low-energy form. A smooth transition over 1--1.25 eV switches
from the low-energy construction to this resonance model. A second smooth
transition in :math:`\log E` from 4 to 10 eV joins the completed Gote model.

For Gote Table 1, ``-1`` in the source JSON means a reported contribution
below 1%, not a negative cross section. The nominal completion uses 0.5%.
A column is accepted when its lower/upper censored sums bracket 100% within
0.3 percentage points. Otherwise it is interpolated from consistent
neighboring angles; the original values remain in the source JSON. The final
row at each energy is the published rotationally summed DCS, in
:math:`10^{-16}` cm2/sr (:math:`10^{-20}` m2/sr), not a percentage. It is
retained as source data but not used as the production absolute normalization,
which remains the IAA residual. For
example, the printed 200 eV, 50-degree column exceeds 100% in its first three
entries. Fractions are interpolated in angle and :math:`\log E`. The reported
ranks through :math:`L=8` (or 10 where provided) are used in the measured interval; unreported higher
ranks there are zero in this completion. Outside 10--160 degrees a smooth
blend returns to the spectator prior at 0 and 180 degrees. This completion
is an assumption about incomplete angular measurements, not additional data.

O2 angular construction
-----------------------

Below 1 eV, every allowed :math:`N\to N+2` transition uses thesis Eq. 11.21b
with its own final electron momentum. In atomic units,

.. math::
   :label: rotation-o2-born

   \frac{d\sigma_{N,N+2}}{d\Omega}
   =\frac{k_f}{k_i}\frac{6}{5}
     \frac{(N+1)(N+2)}{(2N+1)(2N+3)}
     \left(\frac{Q}{3}+\frac{\pi\alpha_2 q}{32}\right)^2 a_0^2,
   \quad q^2=k_i^2+k_f^2-2k_i k_f\cos\theta.

Here :math:`k_i=\sqrt{2E/E_h}` and
:math:`k_f=\sqrt{2\max(E-\Delta,0)/E_h}`, :math:`Q=-0.29` and
:math:`\alpha_2=4.93`. The bracket is a reduced strength; the production
primitive supplies the phase-space factor using the common relativistic
momentum convention. At sub-eV energies its difference from the atomic-unit
nonrelativistic ratio is negligible. A reverse channel evaluates the same
forward primitive at the shifted energy, not a separately fitted formula.

Between 1 and 20 eV a smooth step in :math:`\log E` joins that Born strength
to the positive moment completion below. This is explicitly an interpolation
model across poorly constrained O2 resonance physics. At 20--200 eV the
Bhattacharyya integral and momentum-transfer data constrain the completion.
Let :math:`a_L=\int A_Ld\Omega` and
:math:`m_L=\int(1-\cos\theta)A_Ld\Omega`. Recoupling from :math:`N=1` gives

.. math::
   :label: rotation-o2-moments

   \sigma_{11}&=a_0+\tfrac25a_2,\\
   \sigma_{13}&=\tfrac35a_2+\tfrac49a_4,\\
   \sigma_{\rm incl}-\sigma_{11}-\sigma_{13}
        &=\tfrac59a_4+\sum_{L\ge6}a_L.

Identical linear relations apply to momentum-transfer integrals. The
unresolved higher-rank residual is distributed in proportion to the spectator
prior integrated under the IAA angular marginal. Then :math:`a_2` and
:math:`m_2` follow from :math:`\sigma_{13}` and its momentum-transfer value.
The adopted total is still IAA's total; the paper's absolute :math:`\sigma_{11}`
is therefore not independently imposed in the final IAA-normalized model.
The source check verifies the paper's :math:`1\to3` strength and first angular
moment, after removing the known finite-threshold phase factor for comparison
with the paper's adiabatic calculation.

To specify a positive angular distribution without inventing further measured
moments, use a relative-entropy projection. For nonzero ranks the prior
:math:`r_L(\theta)` is the spectator fraction, except that the rank-two prior
is isotropic; rank zero has prior one. The resulting fractions have the form

.. math::
   :label: rotation-entropy

   q_L(\mu)=\frac{r_L(\mu)e^{\lambda_L+\eta_L\mu}}
                   {1+\sum_{K>0}r_K(\mu)e^{\lambda_K+\eta_K\mu}},\qquad
   q_0=\frac{1}{1+\sum_{K>0}r_K e^{\lambda_K+\eta_K\mu}}.

The multipliers enforce
:math:`\int Fq_Ld\Omega=a_L/\sigma_{\rm IAA}` and
:math:`\int\mu Fq_Ld\Omega=(a_L-m_L)/\sigma_{\rm IAA}`.
They minimize the convex log-partition dual. The reference solves this
numerically and checks the achieved constraints; positivity alone does not
establish convergence. Ranks whose integral is below :math:`10^{-13}` of
the total are omitted from that optimization. The later reciprocal
normalization protects these O2 changing kernels through 200 eV, adjusting
the unchanged component instead. The completion is a declared inference
from incomplete theory data, not a unique O2 DCS derived by the paper.

.. bibliography::
   :keyprefix: rot-
