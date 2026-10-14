.. _beam-air-runtime:

Collision selection and particle lifecycle
==========================================

The physics and acceptance law are specified in
:ref:`multiphysics-collisions-mcc`. This page describes the shared runtime
algorithm and its cost. Physical rate and angular data remain in warpx-data;
preparation is documented in :ref:`beam-air-data-preparation`.

Many-channel selection and DCS reuse
------------------------------------

The cumulative selector covers **all** processes listed in one Background MCC
object: elastic, excitation, ionization, attachment, charge exchange and
two-product reactions. It is not an excitation-only table. The original input
order defines process indices :math:`i=0,\ldots,P-1`. During initialization,
WarpX constructs the sorted union

    .. math::

       \mathcal{E} = \operatorname{unique}\!\left(
       \bigcup_{i=0}^{P-1}\mathcal{E}_i\right)
       = \{E_0,\ldots,E_{U-1}\}

of the knots from every integral cross-section grid, then stores the prefix
rows

    .. math::

       C_j(E_k) = \sum_{i=0}^{j}\sigma_i(E_k),
       \qquad j=0,\ldots,P-1.

Each :math:`\sigma_i` is evaluated with the same piecewise-linear interpolation
as an ordinary ``ScatteringProcess`` lookup. Legacy endpoint clamping is retained;
declared production energy bounds are checked separately. Every input
breakpoint is present in :math:`\mathcal{E}`, so all :math:`\sigma_i` and all
:math:`C_j` are linear between adjacent union knots. Linear interpolation of a
prefix row is therefore exact relative to the stored cross-section
representation, apart from floating-point roundoff; resampling onto the union
does not smooth a threshold or otherwise approximate a channel.

For each particle that passes the global null-collision preselection, WarpX
samples one neutral velocity and calculates one collision energy. Process
selection then consists of:

#. one bisection of :math:`\mathcal{E}` to obtain the bracketing index
   :math:`k` and interpolation fraction;
#. one interpolation of the last prefix row
   :math:`C_{P-1}(E)=\sigma_{\mathrm{tot}}(E)`;
#. one uniform draw for the physical-event acceptance above, converted to
   cross-section units only if accepted; and
#. one binary search over :math:`j` for the first interpolated prefix with
   :math:`C_j(E)` above that draw.

When every process was supplied on a common grid, that grid is also
the union and the cross-section energy interval is found only once per
candidate particle. The same remains true when the input grids differ, because
they were combined at initialization. The selector cost is
:math:`O(\log U+\log P)` rather than :math:`P` separate energy bisections plus a
linear process scan. For example, 64 processes require at most six prefix
comparisons after the energy search.

"One energy search" refers specifically to the integral cross sections used
for event acceptance and process choice. A selected IAA angular DCS has its own
independent energy grid and performs one additional DCS bisection, but only for
the chosen process. RBEQ energy sharing likewise uses its own fixed logarithmic
coordinate. WarpX never searches the angular or RBEQ tables for processes that
were not selected.

The prefix table contains cumulative cross-section weights only; it does not
contain or average energy losses. The binary search returns the original
process index :math:`j`.
WarpX then reads that process's unchanged discrete ``<process>_energy``, type,
angular model, ionization model and product-species group. An excitation event
therefore subtracts exactly the loss configured for its selected channel. An
ionization or attachment event records the same process index for the grouped
particle-creation pass, so cumulative selection cannot disconnect a product
from its channel.

The storage in one host or device copy is

    .. math::

       U\left[\mathrm{sizeof}(\mathtt{ParticleReal})
                +P\,\mathrm{sizeof}(\mathtt{double})\right],

where the union energy grid uses particle precision and the :math:`P` prefix
rows use double precision. Selection draws and cumulative sums in the fallback
also use double precision, preserving rare channels alongside large rates.
WarpX enables the selector only when this total is at most
64 MiB. Consequently,

    .. math::

       U_{\max} = \left\lfloor
       \frac{64\,\mathrm{MiB}}
            {\mathrm{sizeof}(\mathtt{ParticleReal})+
              P\,\mathrm{sizeof}(\mathtt{double})}
       \right\rfloor.

For 64 processes this permits 129,055 union points with double-precision
particles or 130,055 with single-precision particles. For 128 processes the
limits are 65,027 and 65,280 points, respectively. Typical cross-section sets
are far smaller: if 64 channels each use the same 100-point grid, then
:math:`U=100`, not 6,400. Only completely distinct knots can make :math:`U`
approach the sum of the individual grid sizes. A GPU build retains one host
copy and one device copy, each subject to the 64 MiB limit. Original process
grids, angular DCS tables and RBEQ tables are separate allocations.

If the limit would be exceeded, WarpX does not truncate or coarsen any grid. It
uses the exact legacy fallback, which evaluates each process cross section and
scans the cumulative probabilities at runtime. The first-step information
message reports the active selection path, union point count and bytes per
copy, making an accidental fallback visible in production logs.

Prefix rows are stored process-major. Particles in a GPU warp begin each
process binary search at the same middle prefix, while nearby particle energies
usually address nearby entries within that row. All selector arrays are fixed,
contiguous device data; selection performs no allocation, virtual dispatch or
iterative solve in the particle kernel.

Elastic and excitation channels that name the same DCS file share one
precomputed inverse-CDF table and one device allocation. This avoids duplicating
initialization work and device memory when many excitation channels use the
same target DCS. Use the same path string for channels that should share a
table.

The elmolcs data are distributed separately from WarpX and are not copied into
the BSD-licensed source tree. Users must obtain the tables separately and comply
with their license. The construction and intended range of the IAA database are
documented in :cite:t:`b-Schmalzried2023`.

Product creation and attachment removal
---------------------------------------

Ionization and attachment preserve the original selected process index in a
compact product-event record. The record also contains the source-particle
index, sampled neutral velocity, collision energy and destination-group offset.
Consequently, the creation kernel still reads the selected channel's discrete
energy loss, ionization model and angular model; grouping never averages channel
physics. Carrying the already calculated collision energy also avoids repeating
its square root for every accepted ionization event.

On a GPU, the selection kernel reserves one slot in the compact event queue and
one offset in the destination group for each product-changing event. These are
device-local atomic increments. At the recommended
:math:`\nu_{\max}\Delta t_{\mathrm{coll}}\lesssim 0.1`, fewer than about ten
percent of the source particles can pass even the null-collision preselection,
and product-changing events are a subset of those candidates. The atomics are
therefore sparse in the converged operating regime. On a CPU, selection first
writes independent per-source records; one host pass both compacts those
records and computes the group offsets, without unsafe updates inside an
``amrex::ParallelFor``.

The host receives only the small vector of per-group event counts. Each source
or destination particle tile is then resized once. All destination groups and
all secondary electrons are populated by one kernel over the number of compact
events, not by a kernel over every source particle and not by one source scan
per process or product group. Contiguous particle-ID ranges are reserved before
that kernel and assigned within it, avoiding separate ID kernels. Immutable
smart-copy and product-group metadata are constructed once and reused across
collision substeps.

An attachment marks its incident electron invalid. WarpX compacts only the
source tile that contained an attachment and only when that tile actually had
an attachment event. It does not rescan every tile merely because an attachment
channel is configured. Each MPI rank and GPU allocates and compacts its local
particle tiles. Immobile-ion destinations additionally synchronize deposited
mesh increments before accumulation; ranks with no local products must still
participate in that commit. A device-to-host count transfer remains necessary
before resizing a particle tile, because the new allocation size is a host-side
container operation.

Combined rotation and implementation map
----------------------------------------

A ``reciprocal_hybrid`` elastic family replaces its ordinary elastic rate.
It competes in rate-coefficient units with the other processes; its finite
zero-energy superelastic rate cannot be represented as an ordinary finite
cross section. The ordinary union selector remains the exact representation
of the ordinary channels. :ref:`rotation-algorithm` specifies family selection,
joint sampling and signed recoil.

The main implementation files under ``Source/Particles/Collision`` are:

* ``ScatteringProcess.cpp``: integral-file parsing, native grids and range metadata.
* ``BackgroundMCC/BackgroundMCCProcessSelector.*``: host prefix preparation,
  device lookup and the 64 MiB limit.
* ``BackgroundMCC/BackgroundMCCCollision.cpp``: majorants, acceptance, selected
  scattering, runtime checks and product dispatch.
* ``BackgroundMCC/BackgroundMCCElasticScattering.*``: shared angular profiles,
  descriptor reading, host inverse CDFs and device mixture sampling.
* ``BackgroundMCC/BackgroundMCCParticleCreation.H`` and
  ``BackgroundMCCProducts.H``: grouped product creation and removal.
* ``CollisionHandler.cpp``: the reusable endpoint-ordered subcycle schedule;
  :ref:`mcc-collision-timestep` defines its temporal semantics.
* ``IonProductDestination.H``: kinetic versus immobile-fluid ion destinations.

Source fitting and angular quadrature stay offline or at initialization.
Pure scattering reuses scratch/error state and avoids redundant iterator-exit
synchronization; required error checks, RNG synchronization and fluid commits
remain. Product creation still requires host-visible counts before resize.
Batching or deferring that synchronization is a possible optimization, not
current behavior. It must preserve pointer lifetimes and the visibility of
new or attached electrons to subsequent coupled substeps. Benchmark complete
MCC and PIC cost, memory and variance as described in :ref:`beam-air-validation`.
