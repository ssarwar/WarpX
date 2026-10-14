.. _usage-beam-air:

Proton beams and electron collisions in air
===========================================

This guide combines proton-impact ionization, electron collisions in N2/O2,
and kinetic or prescribed species. Start here for configuration; use the
linked theory pages for equations and assumptions. All input defaults and
accepted values remain in :ref:`running-cpp-parameters`.

Choose the species representation
---------------------------------

.. list-table:: Representations and their physical meaning
   :header-rows: 1
   :widths: 24 38 38
   :class: rotational-scattering-table

   * - Component
     - Representation
     - Consequence
   * - Beam
     - Kinetic species, or ``FluidSpecies(model="rigid_beam")``
     - A kinetic beam can move and respond to fields according to its particle
       settings. The fluid beam follows a prescribed Gaussian pulse train.
       The proton-ionization operator does not deplete either representation.
   * - Electrons
     - Kinetic ``Species``
     - Field acceleration and stochastic collisions evolve the distribution.
   * - Positive and negative ions
     - Kinetic species, or ``FluidSpecies(model="immobile")``
     - Immobile ions retain charge and the electron's deposited charge shape;
       they omit ion transport and kinetic energy.
   * - Neutral N2 and O2
     - Collision background densities and temperatures
     - No neutral particles, internal-state populations, depletion, fragments,
       recombination or detachment are evolved by these models.

Prescribed fluids require a single-level laboratory-frame RZ grid, one azimuthal
mode, no moving window or embedded boundaries, and native solver staggering.
They support explicit Yee, explicit PSATD with current correction, and
``semi_implicit_em`` with or without mass matrices. Runtime load balancing is
supported for the explicit solvers; disable it for the semi-implicit solver.
Restart can change the MPI decomposition in each supported configuration.
See :ref:`theory-prescribed-fluids` for projection, source integration and
conservation properties, and :ref:`running-cpp-parameters-fluids` for inputs.

Set the beam and gas
--------------------

A fluid beam specifies exactly one energy/velocity, transverse width,
longitudinal width, normalization and pulse-time convention. For example,
this PICMI configuration fragment describes a finite 800 MeV train:

.. code-block:: python

    from pywarpx import picmi

    beam = picmi.FluidSpecies(
        name="beam", model="rigid_beam", particle_type="proton",
        kinetic_energy=800e6,
        sigma_r=1.0e-3, sigma_t=25.0e-12,
        peak_current=0.6, z_reference=0.03,
        first_pulse_time=0.0, pulse_period=5.0e-9, pulse_count=2,
        initialize_self_fields=True,
    )

Here ``sigma_r`` is an illustrative Cartesian RMS width, not a conversion of
the experimental 5 mm radius. For an untruncated Gaussian,
:math:`r_{\rm rms}=\sqrt{2}\sigma_r` and
:math:`\sigma_z=|v_z|\sigma_t`. A 100 ps four-sigma duration corresponds to
25 ps RMS; a quoted duration must have a defined convention before conversion.
The profile translates on the laboratory z axis, is clipped by the domain,
and neither wraps nor renormalizes at a boundary. Use nonperiodic longitudinal
boundaries for entering and exiting pulses. Initial self fields are computed
once and are restored from the checkpoint on restart.

Register fluids with ``sim.add_fluid_species(beam)`` and no particle layout.
Use ``sim.add_species`` for electrons and kinetic ions. The low-level
``pywarpx.new_fluid_species`` interface uses the text-input parameter names;
:ref:`usage-picmi` documents both Python interfaces.

Set a separate background density for each gas. At pressure :math:`p` and gas
temperature :math:`T_g`, :math:`n_g=x_g p/(k_BT_g)`. The production 1D tests
use a binary 79%/21% N2/O2 mixture; the RZ examples use 78.084%/20.946% of total
air density and omit the remaining species. State which mixture is used.
Translational temperature is independent of the rotational bath: the supplied
production rotational data require a 300 K rotational temperature.

Configure the collision operators
---------------------------------

For each gas, configure one ``ProtonImpactIonizationCollisions`` source from the
beam to kinetic electrons and a singly charged positive-ion destination.
``ionization_target`` selects N2 or O2. The PJG total cross section is computed
from the built-in calibrated spectrum; no external proton rate file is needed.
The source conserves produced weight using a per-cell remainder and an optional
product-count cap. ``fixed_product_weight`` controls emission delay and sampling
noise; ``max_products_per_cell`` bounds work by increasing product weights.
The 800 MeV molecular spectrum is a continuation beyond measured calibration.
See :ref:`multiphysics-collisions-proton-impact-ionization`.

For electrons, use one ``MCCCollisions`` object per gas. Obtain its complete
process dictionary with ``load_iaa_air_processes`` from the external data checkout.
The following fragment assumes ``electrons``, ``O2plus``, ``Ominus`` and
``O2minus`` have already been created and registered:

.. code-block:: python

    from pywarpx.iaa_air import load_iaa_air_processes

    processes = load_iaa_air_processes(
        "/path/to/warpx-data/MCC_cross_sections", "O2",
        ion_species=O2plus,
        attachment_species={"dissociative": Ominus, "three_body": O2minus},
        third_body_density=n_O2,
    )
    electron_O2 = picmi.MCCCollisions(
        name="electron_O2", species=electrons,
        background_density=n_O2, background_temperature=300.0,
        background_mass=31.9988 * 1.66053906660e-27,
        scattering_processes=processes, ndt_subcycle=mcc_subcycles,
    )

Use the analogous N2 call with its positive-ion species and neutral mass
28.0134 atomic mass units; attachment arguments are not needed. Add both MCC
objects and both proton sources to ``warpx_collisions`` on the simulation.
Product charges must be exactly :math:`+e` for ionization and :math:`-e` for
attachment. Set neutral masses explicitly, especially for dissociative attachment.

The helper selects completed IAA excitation, reciprocal elastic/rotation,
thesis RBEQ ionization, IAA dissociative attachment, and MuroranIT three-body
attachment. :ref:`usage-collision-data` documents source alternatives and files.
Do not add separate rotations or total/aggregate scattering tables to that
dictionary. The O2 three-body stabilizer density is :math:`n_{O_2}`, not total
air density. Keep the RBEQ total and secondary-energy model paired.

Refine ``ndt_subcycle`` for all strongly coupled source and collision objects.
Substeps are interleaved in time and tie-broken by input order. Each MCC substep
allows one accepted event per electron; increasing the majorant does not resolve
multiple physical events. Refine the PIC timestep independently to resolve
fields, particle motion and short beam pulses. Vary product weight, source
position resolution, gas quadrature and mesh spacing independently before a
joint refinement. See :ref:`mcc-collision-timestep` and
:ref:`beam-air-validation` for the corresponding checks.

Examples and diagnostics
------------------------

The maintained examples are in ``Examples/Physics_applications/proton_beam_air``:

* ``inputs_rz_fluid_beam_air_picmi.py``: prescribed beam and immobile ions.
* ``inputs_rz_proton_beam_air_picmi.py``: kinetic beam and ion species.
* ``inputs_rz_fluid_beam_air``: native input using the same fluid interface.

Their README explains command-line arguments, restart and diagnostic choices.
The optional ``cross_sections.example.json`` is an example-specific dictionary,
not the production gas manifest. It illustrates a subset of channels. Use the
loader above to configure the full production inventory; do not pass a gas
``manifest.json`` as ``--cross-sections`` to those examples. Running without the
optional electron tables does not simulate complete electron--air chemistry.
The 1D production integration driver can be run with an appropriate 1D Python
build from a scratch output directory:

.. code-block:: bash

    python /path/to/WarpX/Examples/Tests/collision/inputs_test_1d_iaa_air_picmi.py \
        --data-dir /path/to/warpx-data/MCC_cross_sections --target air --mode pulse

That homogeneous pulse regression is not a simulation of experimental beam
propagation or a validation of measured yield.

.. list-table:: Outputs to distinguish
   :header-rows: 1
   :widths: 42 58
   :class: rotational-scattering-table

   * - Output
     - Meaning
   * - ``rho``, ``rho_<species>``
     - Deposited charge, including fluids and the usual charge filtering.
   * - ``fluid_density_<name>``
     - Native fluid number density in m-3; output interpolation can change
       centering relative to the stored solver field.
   * - ``fluid_current_<name>z``
     - Instantaneous unfiltered beam current in A/m2. Maxwell uses a
       timestep-averaged current.
   * - ``ParticleNumber``, ``ParticleCharge``, ``ParticleEnergy``, ``ParticleMomentum``
     - Physical totals include fluids. Fluid macroparticle counts are zero;
       immobile kinetic energy and momentum are zero.
   * - ``PrescribedSourceBudget``
     - Pending and emitted counts, emitted-electron and binding energies,
       and discarded ion kinetic energy. Pending counts carry no charge.

Checkpointed densities, fractional yields, sampling counters and source budgets
continue across restart. Immutable source physics, shape and solver settings
must match. Integrated ``FieldProbe`` and ``FieldPoyntingFlux`` diagnostics also
restore their accumulated state; a restart output does not integrate the saved
interval again. The parameter list describes legacy-checkpoint restrictions.
Particle phase-space, temperature and emittance diagnostics require kinetic
species. ``BeamRelevant`` reports dimensionless Twiss alpha with the
:math:`mc\epsilon` normalization. A prescribed beam is an external reservoir;
its diagnostics are not a closed-system conservation statement.

Before quantitative use, consult :ref:`beam-air-validation` for independent
references, precision limits, current accelerator coverage and remaining gaps.
