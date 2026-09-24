# Validation record

Validated locally on 2026-09-24 with AppleClang 21, a macOS arm64 CPU, and the
NOACC backend. The double-particle build includes 1D, 3D and RZ; the
single-particle build is 1D. Both use double field precision. No CUDA, HIP or
SYCL compiler/hardware run was available. Portable device tests and benchmark
drivers are included; CPU results are not GPU validation.

The updated implementation base is
`26b09d96e32953b47164c9cfc0ae91be863e112f` on
`codex/rigid-beam-immobile-ions-development-sync`. Matching RBEQ rates and
rotational source-audit reports are on `IAA` in `ssarwar/warpx-data`.

## Physics and integration

- The selected double-particle suite passed 40/40 CTest entries; the
  single-particle suite passed 20/20. These cover both RBEQ profiles,
  relativistic ionization, invalid inputs, rotational aliases and their
  cumulative reference, complete MCC integration, proton particle/fluid
  sources, recoil, restart, and bounded performance checks. After the final
  configured-mass consistency change, all ten rotation entries passed again
  in each precision.
- Eq. (2.66) is the sole IAA electron-impact secondary-angle prescription.
  The removed selector is rejected by `BackwardCompatibility()`. Independent
  float/double checks cover quantiles, moments, the free-projectile mass
  shell, isotropy, endpoints, and the nonrelativistic limit.
- Proton tests compare the heavy-projectile angular interval with an
  independent invariant phase-space calculation. Recoil conservation covers
  5 keV through 10 GeV, moving neutrals, molecular-spectrum tails, and rounded
  endpoint inputs. Actual PJG table samples are also checked. The 3D particle
  and RZ prescribed-fluid source paths pass their integration tests.
- RBEQ totals agree with independent dimensional differential quadrature
  within the 0.03% comparison tolerance. The four rate tables use 560–643
  knots, with a 0.02% interpolation target and grids checked in float32.
- Independent Legendre quadrature checks rotational angular-momentum
  coefficients. Finite-mass integral detailed balance, canonical thresholds,
  spin weights, and population/moment convergence are checked separately
  from the sampler.
- The accessibility test compares the allowed signed energy changes with
  numerical minimization of the final energy budget. It checks the retained
  scattering angle, recoil conservation, and conditioned energy-row/outcome
  probabilities against explicit quadrature in both precisions.
- An anisotropic elastic-DCS integration test checks analytic angular
  moments, absence of an added angle–rotation correlation away from the
  accessibility boundary, and isotropic emission at zero relative momentum.
- Adaptive grids check rates and absolute first/second transfer moments
  below 0.2% outside separately identified float32 threshold bands.
  Insignificant moment tails use a 1e-4-of-peak floor. The unchanged outcome
  and discrete signed losses are retained.
- Malformed, truncated or negative bundles, insufficient population
  coverage, invalid temperatures, inconsistent inclusive rates, missing
  elastic DCS inputs, out-of-range electrons, and explicitly mismatched RBEQ
  metadata are rejected.

Deterministic relative-Maxwellian quadrature gives the following rotational
power imbalance, `abs(cooling-heating)/(cooling+heating)`, including the
electron/neutral reduced mass. These are analytic verification kernels, not
validated N2/O2 production data. The requirement is below 0.1%.

| Target | 100 K | 300 K | 1000 K |
| --- | ---: | ---: | ---: |
| N2 | 0.00620% | 0.00298% | 0.00154% |
| O2 | 0.00468% | 0.00237% | 0.00140% |

Integral detailed balance does not establish exact differential detailed
balance when both directions reuse an energy-dependent elastic DCS. Exact
accessibility conditioning can alter individual integral rates in narrow
recoil-threshold bands. The tests do not remove these documented model
approximations.

## CPU performance

A five-seed benchmark alternated sampler order and used 65,536 particles in
each of 14 cases, eight full timesteps, and a scalar reduction to fence device
work. Both paths use the same integral rates, elastic DCS, accessibility,
discrete losses, recoil, and MCC acceptance.

| Sampler | Mean total time for eight steps | Particle-steps/s |
| --- | ---: | ---: |
| alias | 0.846736 s | 8.669e+06 |
| cumulative | 0.850236 s | 8.633e+06 |

The cumulative/alias time ratio was 1.00413. This small difference does not
establish a meaningful throughput advantage. The largest verification table
uses 2.79 MiB in double particle precision and 2.32 MiB in single particle
precision per host/device copy, compared with 59.7/44.8 MiB for the former
angular-bin construction. Maximum measured per-table initialization was
about 4.8 ms in double precision and 4.6 ms in single precision.

The benchmark writes per-seed observables and variance times computational
cost for every case. Five-seed variance estimates fluctuate substantially;
these runs do not establish a statistically reliable noise advantage. Both
samplers pass the same first/second-moment checks and retain physical energy
diffusion.

Reproduce in the matching WarpX Python environment:

```sh
python Tools/CrossSections/benchmark_rotation.py --output build/rotation-benchmark --particles 65536 --steps 8 --repeats 5
```

Increase repeats for variance studies. Run this and the portable CTest
kernels separately on CUDA, HIP and SYCL before drawing GPU performance
conclusions.

## Production source gates

The pinned-archive audit verifies the imported reader and integral data
against the supplied archive. Both unchanged integral residuals are
nonnegative over the tested 0–20 eV range; their reported minimum is zero at
zero energy. The negative angular residuals from the former spectator/Born
decomposition do not apply because this model uses only the elastic DCS for
angle sampling.

The tested zero-temperature reference population remains an explicit audit
assumption. N2 uses corrected elementary integral data and sudden scaling.
O2 uses the angular integral of the Born approximation for rates only; using
the elastic DCS does not extend that rate approximation's sub-eV validity.

No production rotational bundle is exported. Validated low-energy/resonant
integral-rate continuations, convergence in omitted elementary ranks, and
high-energy coverage or a bound on omitted effects are still required. The
beam profile uses source-matched RBEQ rates and Eq. (2.66), but does not enable
unvalidated rotational data.

Full vibrational/electronic/attachment export remains downstream of the
source checks, including the O2 longest-band threshold and SR join
corrections documented in the source audit. No pull request was opened.
