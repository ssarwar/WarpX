MCC acceptance: accuracy and GPU cost
====================================

For the present one-candidate-per-substep implementation, retain the
exponential acceptance ratio for transient beam/plasma simulations with
converged collision subcycling. It avoids an extra majorant-dependent
suppression of low physical collision rates. The measured extra cost is
small in the workload below. This is not a claim of universal accuracy
dominance over every observable or collision kernel.

For the collision integrator itself, full Poisson null-collision sampling
with the usual rate-ratio acceptance is the reference algorithm. Changing
only the acceptance in the current code does not implement that algorithm.
Its runtime has not been measured here.

Three distinct algorithms
-------------------------

Let h be the collision substep, M a constant conservative majorant, and
nu the local collision frequency. Both existing implementations first
attempt a candidate with probability P_M = 1-exp(-M h).

* Development accepts this single possible candidate with nu/M. Thus its
  unconditional real-event probability is P_M nu/M.
* Integration accepts with (1-exp(-nu h))/P_M. Its unconditional first-event
  probability is 1-exp(-nu h), independent of the looseness of M, for a
  frozen local rate. This remains an at-most-one-event approximation.
* Full Poisson null collisions sample successive candidate times at rate M
  (or the corresponding Poisson candidate count for a frozen collision
  operator), accept each candidate using the current nu/M, and update the
  particle state after every accepted event. Candidate times must be used
  when explicit time dependence or remaining evolution of newborn products
  matters. Coupled gases/operators require chronological merging of their
  candidates or sampling a combined majorant; advancing a whole interval of
  one gas before another retains operator-splitting error. This construction
  retains multiple collisions within the interval.

The last construction follows from the Poisson expansion of the collision
propagator; see `Longo (2008), section 3(b), equations 13–18
<https://arxiv.org/pdf/0805.3105>`_. It removes collision-time truncation for
the prescribed collision model under the stated assumptions, not the
separate PIC field/push splitting, cross-section-model, or sampling errors.

Why first-event accuracy is not full operator accuracy
-----------------------------------------------------

For a fixed linear collision generator Q with total event rates nu_i, the
old one-candidate update is I+cQ, where c=(1-exp(-Mh))/M. It preserves a
stationary distribution of Q, although its transient evolution is approximate.
The new update instead replaces each state's nu_i*h by 1-exp(-nu_i*h).
For a state model with stationary weights pi_i and positive rates nu_i, its
stationary weights are proportional to::

    pi_i * nu_i / (1 - exp(-nu_i*h))

They generally differ from pi_i at finite h. This is a simplified
fixed-rate state-model statement; the real thermal-target and multi-operator
problem needs its own convergence checks.

An exact reversible two-state example makes the distinction explicit.
For A->B rate a=1, B->A rate b=4, M=4, and h=0.01, the exact equilibrium is
(A,B)=(0.8,0.2). The exponential one-event rule gives
(0.7975992863,0.2024007137); the old one-candidate rule preserves (0.8,0.2).
Neither one-candidate rule reproduces the complete transient transition
matrix. The exact A->B probability is::

    a/(a+b) * (1-exp(-(a+b)*h))

It is 0.0097541151, versus 0.0098026402 for the old rule and 0.0099501663
for the exponential one-event rule. Summing the Poisson powers of
I+Q/M reproduces the exact transition to 1e-14 in independent checks at
h=0.1, 0.01 and 0.001. The accompanying JSON retains the examples.

Consequently, the current exponential rule is exact for frozen-rate
first-event probabilities, not for arbitrary repeated excitation,
de-excitation, recoil, or angular relaxation over a finite interval.
This qualification is particularly relevant to thermal rotational detailed
balance. Collision subcycling convergence is needed for those observables.

GPU measurement
---------------

Perlmutter job 59501123 completed successfully in 2:46 on an A100-SXM4-80GB.
The comparison uses the production integration extension and the diagnostic
extension differing only in the acceptance expression. Both are double
precision and use the same compiler/dependencies. The source is unchanged.

Each case has 1,048,576 electrons, both N2/O2 reciprocal elastic/rotational
bundles resident, 4096 cells, four boxes, one MPI rank, and 64 collision
subcycles per gas per PIC step. The PIC timestep is 6.4e-13 s. The thermal
spectrum is 300 K; the broad spectrum spans 0.002 eV to 3 MeV.

Positions are frozen and particle push, gather and deposition are disabled.
These are collision-dominated timesteps, not an isolated arithmetic kernel
or a full production PIC simulation. Identical initial momenta are restored
outside every timed block, followed by device synchronization. Initial
moments match exactly across variants. Each process warms up for four
steps, then measures twenty eight-step blocks. Process order is new, old,
old, new, yielding forty samples per rule and spectrum.

.. list-table:: Median milliseconds per collision-dominated PIC timestep
   :header-rows: 1

   * - Spectrum
     - Old acceptance
     - Exponential acceptance
     - Extra cost of exponential rule
     - Exploratory 95% interval
   * - thermal
     - 34.095
     - 34.697
     - 1.77%
     - [1.69%, 1.83%]
   * - broad
     - 49.690
     - 50.360
     - 1.35%
     - [1.33%, 1.41%]

The interval uses a hierarchical bootstrap over two process launches and
five-sample circular blocks within a launch. It is an exploratory timing
interval with limited independent launch replication. Both process-level
comparisons show the same ordering. All final particle counts equal their
initial counts, and the particle energies remain finite and nonnegative.

The old arithmetic is cheaper, but part of its runtime advantage can also
come from processing fewer accepted real collisions. These are equal-step
measurements; they do not establish which integrator is fastest at equal
physical accuracy. A first harness attempt failed before yielding valid
timings and is excluded. The corrected synchronization was also smoke-tested
locally before this GPU run.

Accuracy can dominate subcycling cost
------------------------------------

In this benchmark the N2 majorant is approximately 6.660857636e12 /s, giving
M*h=0.0666086 at 64 subcycles. For nu*h much smaller than one, the old rule
multiplies the physical linearized event probability by::

    (1-exp(-M*h))/(M*h) = 0.967423

This is about 3.26% suppression. Without tightening the majorant, reducing
this particular low-rate bias below 1% requires at least 212 subcycles per
PIC step; below 0.1% requires at least 2131. These are bounds on this specific
error component, not recommended universal subcycle counts or claims that
the current 64-subcycle simulation has converged every observable.

Practical choice
----------------

Keep the current exponential ratio while checking timestep/subcycle
convergence of the observables of interest, including equilibrium and
transport moments. Its measured 1–2% timestep cost is small relative to the
possible cost of compensating for a loose-majorant bias. The earlier noisy
DSMC discharge RMS comparison cannot determine which acceptance is more
accurate universally.

For a future redesign, evaluate a GPU implementation of full Poisson null
collisions with state-dependent rate updates and correct treatment of
newborn particles. It is a candidate for both improved collision-time
accuracy and reduced repeated scans when many fixed subcycles are needed;
no speedup for that unimplemented design is claimed here. Exact algebraic
rearrangements or cached factors can also optimize the current exponential
rule without replacing its probability law by a linear approximation.

``acceptance_comparison_2026_10.json`` contains all 160 timing samples,
process medians, library hashes, commands, analytic checks and statistical
metadata. Scripts and logs remain under
``build/validation/acceptance-performance``. No production collision formula,
input parameter default, test threshold, or cross-section data was changed.
