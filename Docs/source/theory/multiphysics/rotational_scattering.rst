.. _theory-rotational-scattering:

Electron--molecule elastic and rotational scattering
====================================================

``rotation_model = reciprocal_hybrid`` treats a collision of an electron with
N2 or O2 as one vibrationally elastic event. The event can leave the molecule's
rotation unchanged, excite it, or de-excite it. Its scattering angle and its
discrete rotational energy change are sampled jointly. Ordinary elastic recoil
is also applied. Vibrational and electronic excitation, ionization, attachment,
proton-impact collisions and the beam source remain separate physical models.

This is a source-constrained hybrid model, not a complete ab initio molecular
scattering calculation. Its main purpose is transport and thermal energy
exchange in air near room temperature, including a high-energy continuation
for electrons produced by a relativistic beam. The published production data
have a fixed rotational temperature of 300 K and a supported collision-energy
domain of 0--1 GeV. Experimental support does not span this entire domain.
The distinction between measurements, theoretical inputs, fitted closures and
numerical accuracy is part of the specification below.

.. toctree::
   :maxdepth: 2

   rotational_scattering/physics
   rotational_scattering/sources
   rotational_scattering/reciprocity
   rotational_scattering/continuations
   rotational_scattering/algorithm
   rotational_scattering/data_format
   rotational_scattering/validation

Reading and using the model
---------------------------

Start with :ref:`rotation-physics` for the physical quantities and with
:ref:`rotation-sources` for the evidence supporting each energy range.
:ref:`rotation-reciprocity` derives the reconciliation of the IAA inclusive
elastic data with state-resolved rotational data. :ref:`rotation-continuations`
states what is assumed outside the measured ranges. :ref:`rotation-algorithm`
connects those equations to MCC selection, sampling and recoil.
:ref:`rotation-data-format` specifies the production files and initialization.
:ref:`rotation-validation` gives reproducible checks and measured limitations.

Configure one elastic process per gas with its ordinary ``elastic.txt`` as
``cross_section``, ``scattering_angle_model = IAA``, and
``rotation_model = reciprocal_hybrid``. Set ``rotation_file`` to
``MCC_cross_sections/<gas>/IAA/reciprocal_hybrid_300K/thermal_rotation.rot``
in warpx-data, ``rotational_temperature = 300``, and normally
``rotation_sampling = alias``. The combined bundle contains its angular law;
there is no separate elastic DCS input for this option. The ordinary elastic
file checks source consistency. Its rate is replaced, rather than added, by
the bundle's rate. Do not also configure the elementary rotational files as
ordinary MCC processes. Both proton-beam examples resolve the bundle relative
to their cross-section manifest.

The neutral translational temperature may differ from the rotational
temperature. Omission of ``rotational_temperature`` uses a constant
translational temperature, which must match the bundle temperature. Spatially
or temporally varying rotational baths, evolving neutral state populations,
spin-resolved O2, isotope mixtures and vibrationally excited rotors are outside
this model. Additional fixed temperatures require separately prepared data.

The two older options retain their meanings. ``elastic_dcs`` samples a
rotational outcome independently of the elastic angle from V3/V4 rate data;
its integral balance does not guarantee differential balance.
``iaa_spectator`` is the earlier N2 kinetic version of the thesis mean-loss
prescription and does not impose the new coupled balance constraint.
See :ref:`mcc-iaa-sources` before using those compatibility models.
