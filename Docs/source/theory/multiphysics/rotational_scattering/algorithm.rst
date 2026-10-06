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
:math:`j` if :math:`f<a_j`, otherwise its stored alternate. Entry indices
then address a palette of :math:`(\Delta,\epsilon_i)` pairs. A cumulative
reference instead searches the cumulative probabilities. Alias lookup has
constant event cost and one local alternate lookup; cumulative lookup has a
logarithmic search. The startup alias construction costs :math:`O(n)` and
is shared across every subsequent draw from that cell.

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
