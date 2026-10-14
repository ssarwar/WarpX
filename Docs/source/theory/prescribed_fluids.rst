.. _theory-prescribed-fluids:

Prescribed beams and immobile ions
==================================

The ``rigid_beam`` and ``immobile`` fluid models store charge sources on the
mesh without representing them by macroparticles. They couple to kinetic
electrons through Maxwell's equations and collision product deposition.
Unlike the :ref:`cold relativistic fluid <theory-cold-fluid-model>`, they do
not solve a fluid momentum or advection equation. Their configuration and
solver restrictions are described in :ref:`usage-beam-air` and
:ref:`running-cpp-parameters-fluids`.

Gaussian pulse train
--------------------

For pulses of relative amplitude :math:`a_k\geq0` arriving at
:math:`z_{\rm ref}` at times :math:`t_k`, the prescribed number density is

.. math::

   n_b(r,z,t)=n_0\exp\!\left(-\frac{r^2}{2\sigma_r^2}\right)
     \sum_k a_k\exp\!\left[-\frac{(z-z_{\rm ref}-v_z(t-t_k))^2}
                               {2\sigma_z^2}\right].

Each factor is truncated at its configured cutoff, unless that cutoff is
disabled. The pulse sum allows overlap. The profile is clipped by the physical
domain without renormalization or periodic replication. Its continuum sources
are :math:`\rho_b=q_b n_b` and :math:`J_{b,z}=v_z\rho_b`, with no transverse
current. The velocity is constant and may have either sign.

For radial cutoff :math:`R\sigma_r` and longitudinal cutoff :math:`Z\sigma_z`,
the complete unit-amplitude pulse has

.. math::

   A_R=2\pi\sigma_r^2(1-e^{-R^2/2}),\qquad
   I_0=|q_bv_z|n_0A_R,\qquad
   Q_0=|q_b|n_0A_R\sqrt{2\pi}\sigma_z
                     \operatorname{erf}(Z/\sqrt2).

Thus :math:`Q_0=I_0\sqrt{2\pi}\sigma_t\operatorname{erf}(Z/\sqrt2)` with
:math:`\sigma_t=\sigma_z/|v_z|`. These are complete-profile normalizations,
not the instantaneous population inside a clipped domain. For disabled
cutoffs take the corresponding limit to infinity. ``sigma_r`` is the
Cartesian RMS width of the untruncated Gaussian; the radial RMS is
:math:`\sqrt2\sigma_r`. Truncation changes the actual profile moments.

Projection and discrete continuity
----------------------------------

The source is projected with the same compact particle shape order used by
the kinetic species. Radial projection includes the cylindrical measure and
axis treatment. Longitudinal density and current share a cumulative Gaussian
and shape primitive. Gaussian tails use differences of ``erfc`` when needed
to avoid cancellation; integration and profile parameters use double precision
even when field or particle storage is single precision.

The discrete construction can be written schematically with a cumulative
profile :math:`C` at axial faces and radial amplitude :math:`R_i`:

.. math::

   \rho_{ij}^{n}=\frac{q_bR_i}{\Delta z}
      \left(C_{j+1/2}^{n}-C_{j-1/2}^{n}\right),\qquad
   \overline{J}_{i,j+1/2}=\frac{q_bR_i}{\Delta t}
      \left(C_{j+1/2}^{n}-C_{j+1/2}^{n+1}\right).

The telescoping differences give the Yee continuity identity

.. math::

   \frac{\rho_{ij}^{n+1}-\rho_{ij}^{n}}{\Delta t}
   +\frac{\overline{J}_{i,j+1/2}-\overline{J}_{i,j-1/2}}{\Delta z}=0

independently of quadrature accuracy in the shared primitive. In the
implementation the current primitive uses shape order one lower than charge;
differences recover the charge shape. The identity includes flux through
longitudinal physical faces. Charge/current filtering and guard ownership must
remain consistent. PSATD uses its native staggering and current correction;
the Yee finite-difference identity is not a substitute for that spectral step.

The solver receives the timestep-averaged current. Diagnostic
``fluid_current_<name>z`` instead reports the instantaneous unfiltered profile.
The optional initial self-field solve uses the prescribed charge and velocity
with the existing relativistic Poisson solver. Finite-domain boundary conditions
remain part of that solution; it is not the unbounded Gaussian field.

In a semi-implicit field solve, the prescribed beam contributes no mass-matrix
response because its current is independent of the trial field. Persistent
immobile charge is cached during residual evaluations. Collision sources are
advanced on physical collision intervals, not on nonlinear solver iterations;
rollback retains all RZ magnetic components. These distinctions prevent
repeated source production or an artificial field response.

Ionization from the continuous beam
-----------------------------------

For a monoenergetic beam and gas :math:`g`, define
:math:`K_g=Z_b^2\sigma_g(E_b)|v_z|`, with the total from the
:ref:`calibrated PJG spectrum <multiphysics-collisions-proton-impact-ionization>`.
The physical electron yield in cylindrical cell :math:`c` and collision
interval :math:`[t_0,t_1]` is

.. math::

   W_{g,c}=K_g\int_{t_0}^{t_1}\int_c
                   n_g(r,z,t)n_b(r,z,t)\,2\pi r\,dr\,dz\,dt.

For constant gas this separates into radial and longitudinal space-time
integrals. For parser-defined gas, a positive composite midpoint quadrature
samples the gas within each radial, axial and time subinterval, weighted by
the integrated beam measure. ``gas_quadrature_points`` refines that
approximation. Gas values must be finite and nonnegative. The rate cache uses
the same interpolated PJG total as the runtime source; it is not an additional
independent physical fit.

An axial CDF with ``source_sampling_points`` subintervals places products
within the cell. The radial inverse uses the truncated Gaussian annular
measure. Refining the spatial sampler changes product placement, not the
integrated constant-gas yield. Independently shifted quiet sequences control
position, energy, angle and thermal samples. Their seeds combine the collision
name, global cell index and user seed, with checkpointed sequence counters;
they do not depend on MPI ownership.

Let :math:`A_c=W_{g,c}+R_c`, where :math:`R_c` is the pending physical count.
With fixed macroparticle weight :math:`w_0` and cap :math:`N_{\max}`, emit

.. math::

   N_c=\min\!\left(\left\lfloor A_c/w_0\right\rfloor,N_{\max}\right).

Below the cap, each product pair has weight :math:`w_0`; at the cap its weight
is :math:`A_c/N_{\max}`. In finite precision, a small nonnegative remainder
can remain after rounding the stored weight downward. The updated remainder
is always the un-emitted physical count. The cap limits work and particle
growth, not physical yield. Pending production creates neither electrons nor
ion charge. Small-yield cells can delay emission; convergence requires reducing
product weight as well as refining mesh and timestep.

Counts are scanned before allocating disjoint output ranges. No per-event
allocation or inter-cell scatter is needed for source counting. A host-visible
total is still needed to resize particle containers. The energy spectrum and
recoil model are the same as for a kinetic projectile source. Quiet sampling
reduces numerical source noise but does not reproduce physical Poisson counting
fluctuations or remove noise from subsequent electron MCC collisions.

Immobile ion charge and conservation
------------------------------------

An immobile fluid stores a persistent nonnegative number density, initially
zero or supplied by a constant/parser profile. It has zero current, momentum
and kinetic energy. Accepted proton/electron ionization and attachment events
deposit the corresponding emitted or removed electron's position, weight and
shape into the positive- or negative-ion density. Rejected events add nothing.
New increments are synchronized once before accumulation; persistent boundary
support is retained separately so redistribution and repeated commits do not
double count guard-cell contributions.

Consequently an ionization pair deposits zero net new charge and attachment
transfers existing electron charge to the negative-ion fluid. This is a
charge-footprint statement, stronger than equality of total counts. Electron
finite-mass collision kinematics are still evaluated. Ion momentum, transport
and recoil energy are omitted from the immobile representation.

The prescribed beam and neutral gas are external reservoirs: they are not
depleted or heated. ``PrescribedSourceBudget`` records pending and emitted
counts, electron and binding energies, and discarded ion kinetic energy.
It is not a closed-system energy balance. Checkpoints preserve those budgets,
ion densities and sequence counters. They validate immutable beam/source
configuration and reconstruct analytic caches without replaying events.

Source files and independent checks
-----------------------------------

``Source/Fluids/RigidBeamProfile.H`` implements profile integrals and shape
projection; ``RigidBeam.cpp`` configures and deposits the beam.
``ProtonImpactIonization/RigidBeamIonization.cpp`` integrates and samples its
ionization source. ``IonProductDestination.H`` and
``Source/Fluids/ImmobileDensityDeposit.H`` route product charge to kinetic or
fluid destinations. ``MultiFluidContainer.cpp`` owns accumulation, field
coupling and lifecycle operations. See :ref:`beam-air-runtime` for shared
collision scheduling and :ref:`beam-air-validation` for analytical Gaussian,
continuity, charge-footprint, restart and GPU comparisons.
