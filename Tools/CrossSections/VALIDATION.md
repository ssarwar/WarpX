# Validation record

Validated locally on 2026-09-24 with AppleClang 21, a macOS arm64 CPU, the NOACC backend, and a 1D WarpX build. Field precision was double; particle precision was tested in both double and single. No CUDA, HIP, or SYCL hardware/compiler run was available. Portable device tests and benchmark drivers are included; CPU results are not GPU validation.

The implementation base is `4cd2a28fd639ba73ac917f4181ac9bb51d524731` on `codex/rigid-beam-immobile-ions-development-sync`. Matching RBEQ rates and rotational source-gate reports are on the `IAA` branch of `ssarwar/warpx-data`.

## Physics and integration

- The selected double-precision suite passed 19/19 CTest entries; the corresponding single-particle-precision suite passed 14/14. These include both RBEQ models, relativistic ionization, invalid inputs, rotational aliases, the explicit cumulative reference, and complete MCC integration. Ordinary-selector and loose-majorant regressions also passed in double precision.
- Additional sampler checks cover the zero-energy limit, 1e-12 and 1e-10 eV, moving neutrals, and signed recoil at 10 keV, 2.5 MeV, and 1 GeV. Energy and momentum conservation are normalized to the full available energy/momentum, including released rotational energy.
- RBEQ totals agree with independent dimensional differential quadrature within the 0.03% comparison tolerance. The four rate tables use 560–643 knots, with a 0.02% interpolation target and grids checked in float32.
- Independent Legendre quadrature checks the angular-momentum coefficients. Finite-mass detailed balance, canonical thresholds, spin weights, and population/moment convergence are tested separately from the sampler.
- Adaptive rate grids check rates and absolute first/second transfer moments below 0.2% outside the separately identified float32 threshold bands. Insignificant moment tails use a 1e-4-of-peak floor. Sampling checks retain the unchanged atom, signed changes, angular moments, and angle–energy correlation.
- Malformed/truncated/negative rotational bundles, insufficient population coverage, invalid temperatures, inconsistent inclusive rates, out-of-range electrons, and explicitly mismatched RBEQ metadata are rejected.

Deterministic Maxwellian quadrature gives the following rotational power imbalance, expressed as `abs(cooling-heating)/(cooling+heating)`. These are analytic verification kernels, not validated N2/O2 production data. The requirement is below 0.1%.

| Target | 100 K | 300 K | 1000 K |
| --- | ---: | ---: | ---: |
| N2 | 0.00692% | 0.00340% | 0.00178% |
| O2 | 0.00522% | 0.00269% | 0.00157% |

## CPU performance

A five-seed benchmark alternated sampler order and used 65,536 particles in each of 14 cases, eight full timesteps, and a scalar reduction to fence device work. Both paths use the same rates, angular bins, discrete losses, recoil, and MCC acceptance.

| Sampler | Mean total time for eight steps | Particle-steps/s |
| --- | ---: | ---: |
| alias | 0.706554 s | 1.04e+07 |
| cumulative | 0.717635 s | 1.02e+07 |

The full-timestep cumulative/alias time ratio was 1.0157; the improvement is modest because lookup is only part of the timestep. Isolated sampler measurements showed approximately 1.5–1.8× speedups in these CPU runs. The largest verification table occupied about 59.7 MiB in double precision and 44.8 MiB in single precision per host/device copy.

The benchmark writes per-seed observables and variance times computational cost for every case. Five-seed variance estimates fluctuate substantially; these runs do not establish a statistically reliable noise advantage. Both samplers pass the same first/second-moment checks and retain physical energy diffusion.

Reproduce with `benchmark_rotation.py --output build/rotation-benchmark --particles 65536 --steps 8 --repeats 5` in the matching WarpX Python environment. Increase repeats for variance studies. Run this and the portable CTest kernels separately on CUDA, HIP and SYCL before drawing GPU performance conclusions.

## Production source gates

The pinned-archive audit reconstructs the DCS using its named angle/energy axes and degrees, verifies the actual imported reader/data files against the archive, and normalizes the DCS separately from the inclusive integral. The tested zero-temperature reference is an explicit modeling assumption. Both candidate decompositions have negative unchanged angular-bin rates:

| Target/model | Energy of reported minimum | Minimum unchanged bin rate |
| --- | ---: | ---: |
| N2 / iaa_sudden_spectator | 2.42206 eV | -4.32958e-15 m³/s |
| O2 / iaa_born | 20 eV | -1.13087e-15 m³/s |

No negative residual is clipped, and no production rotational bundle is exported. A consistent angular decomposition, low-energy/resonant closure, convergence in elementary rotational rank, and high-energy continuation or an omitted-effect bound are still required. The new beam profile selects the source-matched RBEQ model and Eq. 2.66; it does not enable these rejected rotational candidates.

Full vibrational/electronic/attachment export remains downstream of these source gates, including the O2 longest-band threshold and SR join corrections documented in the source audit. No pull request was opened.
