.. _multiphysics-collisions:

Collisions
==========

WarpX includes several different models to capture collisional processes
including collisions between kinetic particles (Coulomb collisions, DSMC,
nuclear fusion) as well as collisions between kinetic particles and a fixed
(i.e. non-evolving) background species (MCC, background stopping).

.. _multiphysics-collisions-mcc:

Background Monte Carlo Collisions (MCC)
---------------------------------------

Several types of collisions between simulation particles and a neutral
background gas are supported including elastic scattering, back scattering,
charge exchange, excitation collisions and impact ionization.

The so-called null collision strategy is used in order to minimize the
computational burden of the MCC module. This strategy is standard in PIC-MCC and
a detailed description can be found elsewhere, for example in :cite:t:`b-Birdsall1991`.
In short the maximum collision probability is found over a sensible range of
energies and is used to pre-select the appropriate number of macroparticles for
collision consideration. Only these pre-selected particles are then individually
considered for a collision based on their energy and the cross-sections of all
the different collisional processes included.

The MCC implementation assumes that the background neutral particles are
**thermal** and move at non-relativistic velocities in the laboratory frame.
For each simulation particle considered for a collision, a neutral ordinary
velocity :math:`\boldsymbol{V}_n` is sampled from the user-specified classical
Maxwellian distribution. WarpX particle momentum components are normalized
momenta, or proper velocities,
:math:`\boldsymbol{u}=\gamma\boldsymbol{v}`.

For an electron projectile, only scalar relative kinematics are needed for the
collision draw. Expanding consistently for a non-relativistic neutral gives

    .. math::

       \gamma_{\mathrm{rel}}-1 \simeq
       (\gamma_e-1)
       - \frac{\boldsymbol{u}_e\cdot\boldsymbol{V}_n}{c^2}
       + \gamma_e\frac{V_n^2}{2c^2}.

The terms shown are important at very low electron speed: in the
non-relativistic limit they reduce to
:math:`|\boldsymbol{v}_e-\boldsymbol{V}_n|^2/(2c^2)`. The cross-section lookup
energy and collision-rate speed are evaluated as

    .. math::

       E_{\mathrm{lookup}} = (\gamma_{\mathrm{rel}}-1)m_e c^2,
       \qquad
       g_{\mathrm{coll}} \simeq
       \frac{c\sqrt{(\gamma_{\mathrm{rel}}-1)(\gamma_{\mathrm{rel}}+1)}}
            {\gamma_e}.

The omitted neutral-Lorentz-factor correction in the flux prefactor is of
relative order :math:`V_n^2/c^2`. In the stationary-neutral limit,
:math:`g_{\mathrm{coll}}` is the ordinary electron speed, not the proper speed
:math:`|\boldsymbol{u}_e|`.

The lookup energy above is the electron kinetic energy in the sampled neutral
rest frame. The total center-of-momentum kinetic energy is

    .. math::

       E_{\mathrm{COM}} =
       \left[\sqrt{m_e^2+M^2+2m_eM\gamma_{\mathrm{rel}}}
       -(m_e+M)\right]c^2.

For electron collisions with atomic or molecular neutrals, the relative
difference between these two energies is of order :math:`m_e/M` until the
electron becomes extremely relativistic. They can therefore be treated as the
same cross-section argument over the intended range. WarpX evaluates
:math:`E_{\mathrm{lookup}}` directly: using
``ParticleUtils::getCollisionEnergy()`` would require reconstructing the
relative proper speed and would add two square roots to each collision
candidate.

Thus, the frequency for process :math:`i` is
:math:`\nu_i=n_n\sigma_i(E_{\mathrm{lookup}})g_{\mathrm{coll}}`. A full
three-vector Lorentz transformation is not performed.

For impact ionization, neutral thermal motion is neglected in the collision
draw because typical ionization thresholds are far above neutral thermal
energies. This avoids three random-number draws per ionization candidate. For
non-electron projectiles, the existing center-of-momentum collision-energy
convention from ``ParticleUtils::getCollisionEnergy()`` is retained.

The configured ``background_mass`` always denotes the neutral target mass and
is kept separate from the mass of any ionization product species. If it is
omitted for an ionizing electron-neutral collision, the neutral mass is inferred
from the positive-ion product mass plus one electron mass.

Once a particle is selected for a specific collision process, that process determines how the particle is scattered as outlined below.

.. _multiphysics-collisions-pulseddecay:

Pulsed Decay
------------

This collision module can be used to have a parent species decay into two product
species with a user-defined decay rate. Mathematically, it solves the following
rate equations on a cell-by-cell basis:

    .. math::

       \begin{aligned}
        \frac{dn_1}{dt} &= -\nu(t)n_1, \\
        \frac{dn_A}{dt} &= +\nu(t)n_1 = \frac{dn_B}{dt}, \\
       \end{aligned}

where :math:`n_1` is the parent species density, :math:`n_A` and :math:`n_B` are the product species densities,
and :math:`\nu(x,y,z,t)` is the user-specified decay rate.

This can be used, for example, to represent ionization of a parent species by an externally applied laser pulse.

.. _multiphysics-collisions-dsmc:

Direct Simulation Monte Carlo (DSMC)
------------------------------------

The algorithm by which binary collisions are treated is outlined below. The
description assumes collisions between different species.

1. Particles from both species are sorted by grid-cells.
2. The order of the particles in each cell is shuffled.
3. Within each cell, particles are paired to form collision partners. Particles
   of the species with fewer members in a given cell is split in half so that
   each particle has exactly one partner of the other species.
4. Each collision pair is considered for a collision using the same logic as in
   the MCC description above.
5. Particles that are chosen for collision are scattered according to the
   selected collision process.

Scattering processes
--------------------

Charge exchange
^^^^^^^^^^^^^^^

This process can occur when an ion and neutral (of the same species) collide
and results in the exchange of an electron. The ion velocity is simply replaced
with the neutral velocity and vice-versa.

Elastic scattering
^^^^^^^^^^^^^^^^^^

The ``elastic`` option uses isotropic scattering, i.e., with a differential
cross section that is independent of angle.
This scattering process as well as the ones below that relate to it, are all
performed in the center-of-momentum (COM) frame. Designating the COM velocity of
the particle as :math:`\boldsymbol{u}_c` and its labframe velocity as :math:`\boldsymbol{u}_l`,
the transformation from lab frame to COM frame is done with a general Lorentz
boost (see function ``ParticleUtils::doLorentzTransform()``):

    .. math::
            \begin{bmatrix}
                \gamma_c c \\
                u_{cx} \\
                u_{cy} \\
                u_{cz}
            \end{bmatrix}
         = \begin{bmatrix}
                \gamma & -\gamma\beta_x & -\gamma\beta_y & -\gamma\beta_z \\
                -\gamma\beta_x & 1+(\gamma-1)\frac{\beta_x^2}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} \\
                -\gamma\beta_y & (\gamma-1)\frac{\beta_x\beta_y}{\beta^2} & 1 +(\gamma-1)\frac{\beta_y^2}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} \\
                -\gamma\beta_z & (\gamma-1)\frac{\beta_x\beta_z}{\beta^2} & (\gamma-1)\frac{\beta_y\beta_z}{\beta^2} & 1+(\gamma-1)\frac{\beta_z^2}{\beta^2} \\
            \end{bmatrix} \begin{bmatrix}
                \gamma_l c \\
                u_{lx} \\
                u_{ly} \\
                u_{lz}
            \end{bmatrix}

where :math:`\gamma` is the Lorentz factor of the relative speed between the lab frame and the COM frame, :math:`\beta_i = v^{COM}_i/c` is the i'th component of the relative velocity between the lab frame and the COM frame with

    .. math::

        \boldsymbol{v}^{COM} = \frac{m \boldsymbol{u}_c}{\gamma_u m + M}

The particle velocity in the COM frame is then isotropically scattered using the function ``ParticleUtils::RandomizeVelocity()``. After the direction of the velocity vector has been appropriately changed, it is transformed back to the lab frame with the reversed Lorentz transform as was done above followed by the reverse Galilean transformation using the starting neutral velocity.

Back scattering
^^^^^^^^^^^^^^^

The process is the same as for elastic scattering above except the scattering angle is fixed at :math:`\pi`, meaning the particle velocity in the COM frame is updated to :math:`-\boldsymbol{u}_c`.

Excitation
^^^^^^^^^^

The process is also the same as for elastic scattering except the excitation energy cost is subtracted from the particle energy. This is done by reducing the velocity before a scattering angle is chosen.

Forward scattering
^^^^^^^^^^^^^^^^^^

This process operates in two ways:

1. If an excitation energy cost is provided, the energy cost is subtracted from the particle energy and no scattering is performed.
2. If an excitation energy cost is not provided, the particle is not scattered and the velocity is unchanged (corresponding to a scattering angle of :math:`0` in the elastic scattering process above).

See :cite:t:`b-Janssen2016` for a recommended use of this process.

Benchmarks
----------

See the :ref:`MCC example <examples-capacitive-discharge>` for a benchmark of the MCC
implementation against literature results.

Particle cooling due to elastic collisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is straight forward to determine the energy a projectile loses during an elastic collision with another body, as a function of scattering angle, through energy and momentum conservation.
See for example :cite:t:`b-Lim2007` for a derivation. The result is that given a projectile with mass :math:`m`, a target with mass :math:`M`, a scattering angle :math:`\theta`, and collision energy :math:`E`, the post collision energy of the projectile is given by

    .. math::

       \begin{aligned}
       E_{final} = E - &[(E + mc^2)\sin^2\theta + Mc^2 - \cos(\theta)\sqrt{M^2c^4 - m^2c^4\sin^2\theta}] \\
       &\times \frac{E(E+2mc^2)}{(E+mc^2+Mc^2)^2 - E(E+2mc^2)\cos^2\theta}
       \end{aligned}

The impact of incorporating relativistic effects in the MCC routine can be seen in the plots below where high energy collisions are considered with both a classical and relativistic implementation of MCC. It is observed that the classical version of MCC reproduces the classical limit of the above equation but especially for ions, this result differs substantially from the fully relativistic result.

.. figure:: https://user-images.githubusercontent.com/40245517/170900079-74e505a5-2790-44f5-ac84-5847eda954e6.png
   :alt: Comparison of classical and relativistic MCC collision results
   :width: 96%

   Classical v. relativistic MCC.

.. bibliography::
    :keyprefix: b-