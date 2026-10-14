.. _electron-impact-ionization:

Electron-impact ionization
==========================

The IAA models follow Schmalzried's 2023 thesis, especially Section 11.6,
Table 11.12 and Eqs. (2.66), (11.119)--(11.122), :cite:`b-Schmalzried2023`.
Configuration is described in :ref:`usage-collision-data` and
:ref:`usage-beam-air`; independent checks are in :ref:`beam-air-validation`.

Energy sharing
--------------

The user-supplied integral cross-section table always determines the
ionization event rate. By default, WarpX subtracts the configured threshold
``<process>_energy`` and shares the remaining energy equally between the two
outgoing electrons. The per-process ``RBEQ`` energy-sharing model instead uses
the relativistic binary-encounter Bethe, or RBEQ, singly differential cross
section described by :cite:t:`b-Schmalzried2023`::

    mcc.ionization_N2_energy_sharing_model = RBEQ
    mcc.ionization_N2_rbeq_target = N2
    mcc.ionization_N2_energy = 15.58

``RBEQ`` is available for ``N2`` and ``O2``. The default
``rbeq_model = iaa_thesis_2023`` selects one of five
:math:`\mathrm{N}_2` or six :math:`\mathrm{O}_2` subshells from its positive
RBEQ partial cross section. The selected binding energy :math:`B_i`, rather
than only the outer-shell threshold, is then removed from the incident energy.
The lower-energy electron is sampled from the conditional RBEQ distribution on

    .. math::

       0 \leq T_s \leq \frac{T-B_i}{2},

and the other electron receives the remainder, apart from molecular-ion recoil.
The configured process threshold must match the target outer-shell binding
energy within 0.05 eV: 15.58 eV for :math:`\mathrm{N}_2` and 12.07 eV for
:math:`\mathrm{O}_2`.

Some unit-oscillator-strength RBEQ partials become slightly negative immediately
above their subshell thresholds. A negative partial is unphysical and cannot
define a probability, so WarpX clamps it to zero until the analytic expression
becomes positive; it is never reflected with an absolute value. WarpX stores
the analytic zero-crossing coordinate and interpolates the unnormalized,
non-negative partial cross sections. Consequently, interpolation cannot activate
a shell below its crossing, and normalization occurs only after interpolation
at the collision energy.

For several of these shells, the integrated partial becomes positive slightly
before the published SDCS is non-negative over its complete kinematic interval.
A non-monotone cumulative function cannot be sampled as a probability. In that
narrow interval, WarpX uses a uniform conditional energy distribution and
switches to RBEQ once the complete SDCS is non-negative. This continuation is
also used within 0.1 percent of every binding threshold, where direct evaluation
is ill-conditioned. It preserves a finite, symmetric threshold limit without
turning negative values into artificial positive probability.

The partial cross sections and conditional inverse CDFs are precomputed on
logarithmic energy grids during initialization. Subshell selection uses 2,049
energy points, while the much larger inverse-CDF data use 257. The two grids
cover the same range, so the shell coordinate is obtained from the inverse-CDF
coordinate by one fixed scale factor rather than another logarithm or search.
This independently resolves sharp shell onsets without multiplying the
inverse-CDF memory footprint. The inverse CDF uses 513 samples of the symmetric
probability map

    .. math::

       q(x) = \frac{x^4}{x^4+(1-x)^4}, \qquad 0 \leq x \leq 1,

The random energy probability has 53 bits, even with float particles.
Its complement :math:`1-q` is formed before conversion to particle
precision; otherwise rounding :math:`q` to one can collapse an important
part of the relativistic hard tail onto the endpoint. The fourth roots and
table interpolation still use particle precision. Equal sharing requires
no energy or shell random draw.

This map resolves both high-energy probability tails without linearly
interpolating from an interior quantile to the physical endpoint. The tables
cover incident energies through at least 1 GeV. Each accepted event then
requires one logarithm and fixed-size interpolation, with no device allocation
or iterative root solve.

Ionization angles and recoil
----------------------------

The angular model is selected independently from energy sharing. With
``<process>_scattering_angle_model = IAA``, the collision is evaluated in the
sampled neutral rest frame. Let :math:`T_a=T-B_i`, and let :math:`T_p` and
:math:`T_s` be the unrecoiled primary and secondary energies. The primary
electron follows the relativistic binary-encounter relation

    .. math::

       \cos\theta_p =
       \sqrt{\frac{T_p(T_a+2m_ec^2)}{T_a(T_p+2m_ec^2)}}.

The bound/free interpolation for the secondary electron follows
Schmalzried Eq. (2.66) at all incident energies:

    .. math::

       \cos\theta_s =
       \frac{T_s}{T_s+B_i}
       \sqrt{\frac{T_s(T_a+2m_ec^2)}{T_a(T_s+2m_ec^2)}}
       + \frac{B_i}{T_s+B_i}\xi,
       \qquad \xi\sim\mathcal{U}[-1,1].

This convex combination is bounded without clipping its probability
distribution. The nonrelativistic binary factor tends to
:math:`\sqrt{T_s/T_a}`, while the full model becomes isotropic as
:math:`T_s/B_i\to0`. These approximate angles do not themselves enforce
finite ion recoil; the subsequent three-product solve does that.
The removed ``<process>_secondary_angle_model`` input is rejected with a
migration message. There is no secondary-angle selector within the IAA model.

The electron azimuths differ by :math:`\pi`; the product ion receives the
remaining momentum. WarpX solves for the total electron energy after ion recoil
with three initial Newton updates. The analytic derivative includes both
outgoing electron momenta, while a stable :math:`pc` representation avoids
subtracting the ion rest energy. The energy residual is checked explicitly;
if needed, a bounded safeguarded Newton/bisection fallback takes at most
48 further iterations. A draw that cannot produce an energy-conserving state
is a null collision: the incident electron is unchanged and both reserved
product slots are removed before another collision operator runs. This matters
near a shell threshold, where prescribed angles/sharing can be inadmissible,
and at extreme incident energies where three Newton iterations are insufficient.
Kernel tests extend through 100 GeV, including deliberately adverse backward
and isotropic events. The accepted products are Lorentz transformed back to
the simulation frame. This numerical coverage is not a validation of RBEQ
or its angular closure against data at those energies.

The positive-loss threshold guard also applies to ionization; see
:ref:`mcc-finite-recoil`. Shell-specific product kinematics can reject a draw
as described above, independently of the configured outer-shell threshold.

RBEQ parameterizations
----------------------

The original WarpX N2 and O2 shell parameters agree with thesis Table 11.12.
The alternative ``rbeq_model = elmolcs_b8643810`` uses the supplied archive
with SHA256 ``b864381086120fff37eb8c658dfcb0f150b80969cf9add626a38d89bab0a3c38``.
The identifier pins a source snapshot, not an installed package version.
Its fitted elmolcs ``iaa*`` N2 block splits the thesis 18.72 eV shell
into 18.746 and 23.6 eV contributions, changes 37.30 to 37.8 eV, and changes
16.74 to 16.716 eV. Its O2 inner binding is 531 eV instead of 543.8 eV.
The package's default O2 orbital data are a third, eight-contribution model
and must not substitute for either fitted model.

The thesis uses the asymptotic dipole correction in Eq. (11.121).
The package evaluates its logarithmic moment at 1 MeV instead.
Furthermore, ``genRBEQ`` uses ``1+(E+B)/mc2`` in its exchange and squared-binding
terms, whereas ``_RBEQ`` uses ``1+(E+U+B)/mc2``. The latter is the common
kernel underlying the differential and cumulative package helpers and the
selected snapshot model. This small additional helper inconsistency must
be distinguished from the fitted-parameter differences.
The ``lnBm`` antiderivative also differs from numerical integration of its
stated logarithmic moment: its last logarithmic coefficient is 1/4 instead
of 1/2. The snapshot option preserves that formula; it is not a correction
of the underlying moment integral. Independent quadrature checks include
this explicit difference.

The numerical N2 ``iaa`` table is close to the thesis parameterization;
the analytical ``iaa*`` fit is a different model. O2's numerical table is
close to its package fit. Matching production tables must be regenerated.
Negative shell totals are replaced by zero before summing production rates.
Where the differential formula is negative, the conditional sharing is
uniform until the whole differential distribution is nonnegative, including
the narrow continuation within 0.1 percent of a shell's binding energy.
Signed source totals are retained only for reference comparisons.

IAA ionization uses Eq. (2.66) for the secondary at all energies and
Eq. (2.60) for the primary. For 2.5 MeV impact, binding 15.58 eV and equal
sharing, the secondary mean cosine is approximately 0.8804. The ion recoil
is enforced separately through the three-product kinematics solve.
The omission of Coulomb interference is intentional.
