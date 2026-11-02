# Current-rate128 photon spectrum and carbon containment ✨

For the selected **A-H3/B1-H2/B2-H2/C13-H4 median-carbon128 quartet**, the broad D = 1.10 photon scenario gives an upper hot-receiver leakage fraction of **4.4173% at 1000 kg/m²** and **0.83850% at 2000 kg/m²**, under the localized-cascade/full-column assumptions below. These are conditional bounds on the specified photon components. They are usable inputs to a receiver-mass versus capture calculation; **99% capture is not imposed**.

The current generalized cycle supplies `H + U = 26.61166509 MeV/net He`. Its actual delayed parent inventories are N13 **0.09112580**, O14 **0.90898811**, and O15 **1.00011514 decays/net He**. The primary photon energy is **4.13199936 MeV/He**; mean positron kinetic energy is **1.48315986 MeV/He**; saved grey prompt energy is **0.37413600 MeV/He**. The current inventory differs from the frozen selected64 packet, especially in its N13/O14 split.

## Column table

All entries are percentages of the same current `H + U` boundary. The nominal-line column retains the measured rest-frame nuclear and 511 keV lines, assigns all mean beta kinetic energy to photons bounded by 5.14304 MeV, and assigns all saved grey prompt energy to 20 MeV laboratory photons. The broad column combines annihilation and beta energy below `1.10 × 5.14304 MeV`, boosts nuclear-line energies and source weights by 1.10, and uses the same 20 MeV prompt reserve. Thus “nominal” describes the delayed line treatment; the grey spectrum is still a conservative scenario.

| Carbon column (kg/m²) | Nominal-line leakage bound (%) | Broad D1.10 bound (%) | Broad plus carried 6.4 ppm (%) |
|---:|---:|---:|---:|
| 500 | 8.796559 | 10.356093 | 10.356733 |
| 750 | 5.532467 | 6.751790 | 6.752430 |
| 1000 | 3.495356 | 4.417263 | 4.417903 |
| 1250 | 2.219292 | 2.900106 | 2.900746 |
| 1500 | 1.416694 | 1.910801 | 1.911441 |
| 2000 | 0.587704 | 0.838498 | 0.839138 |
| 3000 | 0.109774 | 0.168751 | 0.169391 |
| 5000 | 0.005462 | 0.008208 | 0.008848 |

At 1000 kg/m² the broad scenario corresponds to at least **95.58274% hot capture for the included components**, or **95.58210%** after the optional carried reserve. Other losses and the geometry qualification remain separate. The 6.4 ppm reserve is an inherited full-loss sensitivity from the older positive-flow hot-gamma screen. It has **not** been recomputed for the current128 trajectories or updated rates and may overlap the saved grey prompt account. It is exposed as an optional parameter rather than added invisibly.

## Actual endpoint Doppler check

| Stage | Maximum saved speed / c | Maximum forward bulk D | Maximum endpoint T (MK) | D1.10 speed margin / thermal sigma |
|---|---:|---:|---:|---:|
| A | 0.060735 | 1.062696 | 115.26 | 37.86 |
| B1 | 0.041698 | 1.042605 | 38.95 | 101.29 |
| B2 | 0.058350 | 1.060156 | 171.75 | 33.17 |
| C13 | 0.087095 | 1.091242 | 105.17 | 9.16 |

The fastest saved endpoint is C13 at **26.1105 million m/s**. D1.10 corresponds to speed 0.095023c. The thermal diagnostic assigns every radioactive parent the global maximum endpoint temperature and minimum mass number 13, then uses the triangle-inequality speed bound. Its largest Maxwell-tail probability is **4.31 × 10⁻¹⁸**. The old selected64 21.3-sigma statement does not apply to this family. This is a conditional Maxwellian endpoint screen; later pressure acceleration, external fields or a non-Maxwellian energetic population require their own velocity envelope.

The broad source uplift consumes **0.561516 MeV/He** from the already counted **19.59080 MeV/He** bulk kinetic budget. It adds no nuclear energy. The evaluated gamma-plus-beta inventory exceeds the old approximate delayed-heat ledger by **0.00105962 MeV/He** (0.00398% of `H+U`); retaining that tiny surplus in the spectral screen is conservative and does not revise the plant’s total energy. Neutral-atom branch data are used. Completely suppressing electron capture could add about **0.00233793 MeV/He** of annihilation energy before attenuation; charge-dependent branching is a separate small sensitivity, not a solved plasma-decay network.

## Method and assumptions

The decay data are the same primary evaluations documented in the earlier packet, checked again for this reduction: [LNHB N13](https://www.lnhb.fr/nuclides/N-13_tables.pdf), [LNHB O15](https://www.lnhb.fr/nuclides/O-15_tables.pdf), and the explicitly dated [ENSDF July 1990 O14 evaluation](https://www.nndc.bnl.gov/ensnds/14/N/ec_decay.pdf). They supply beta branches, mean positron energies and the O14 gamma cascade. Current hot isotope outputs are weighted by the current generalized cycle’s actual stage attempts. All other modeled unstable isotopes receive a full atomic-Q leakage reserve; this adds only 7.50 × 10⁻¹⁰ MeV/He.

The [NIST carbon table](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html) supplies mass energy-absorption coefficients. The conversion is `1 cm²/g = 0.1 m²/kg`. The algorithm takes the minimum coefficient along the explicit piecewise log-linear curve from 1 keV to each maximum photon energy, rounds that minimum down to 10⁻⁵ m²/kg, and uses `exp(−k_min × Sigma)` for remaining photon energy. The 511 keV and 2.312593 MeV lines both use **0.00203 m²/kg**, reflecting the lower coefficient at degraded energies rather than simply the coefficient at incidence.

NIST’s [coefficient definition](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html) accounts for the radiative fraction of secondary charged-particle slowing. Applying it as a local cascade kernel assumes charged secondaries deposit and regenerate radiation near their interaction point. Explicit electron/positron ranges, boundary escape and displaced bremsstrahlung are not solved here. No new Monte Carlo calculation is claimed.

The energy-floor bound also assumes no significant hard-photon upscatter, absorption of photons below the tabulated 1 keV limit, a closed cavity returning inward-directed photons, and the stated carbon column between the source and external escape. Decay material carried behind the credited shield needs a separate source-depth calculation. Ports, gaps and intercepted energy in cold apparatus cannot be silently credited to the hot receiver. These conditions should travel with any use of the table.

## Hot capture, external escape and cold equipment

Let `L` be the table’s leakage fraction and `f_cold` the fraction of that missed hot energy intercepted by cold equipment. For the chosen bound scenario, the partition of `H+U` is:

- Hot receiver: `1 − L`.
- Cold equipment: `f_cold × L`.
- External escape: `(1 − f_cold) × L`.

`f_cold = 0` and `f_cold = 1` bracket external escape and complete cold interception. Energy intercepted cold contributes radiator load without hot-engine work; it must not be counted simultaneously as escaped energy. For a thermal conversion efficiency `eta`, total rejection from these channels is `[1 − eta(1 − L)](H+U)` in the all-cold-interception case, before separately applying driver and auxiliary recycling according to the plant ledger. This packet exposes the partition; the root-owned capture/column tradeoff retains responsibility for plant mass, rate and cold-sink re-evaluation.

## Callable API and reproduction

`current_photon128.py` is import-safe and uses only the standard library for the table API:

```python
result = evaluate(
    spectrum,
    sigma_kg_m2=1000,
    mode="broad_D1p10",
    hot_bypass_ppm=6.4,
    cold_intercept_fraction=0.0,
)
```

Load `spectrum.json` as `spectrum`. Modes are `nominal_lines`, `broad_rest`, `broad_actual_bulk`, and `broad_D1p10`. Return fields include `leakage_upper_fraction_thermal`, individual source components, `hot_capture_lower_fraction_for_included_components`, and `partition_at_this_bound`. `containment.json` records all requested columns and the optional carry separately.

Raw spectrum regeneration uses `--research-root`, `--inputs` pointing to `mixture-design/current-rate128-reduced-inputs.json`, `--plant` pointing to its `current-rate128-sigma5000.json`, `--sources sources.json`, and a new `--out-dir`. Only source hashes are read from the original records; inventory/profile data are already in the reduced-input archive. `verify.py` passes **257/257** arithmetic, provenance, NIST-knot/floor, Doppler, spectral ordering and partition checks. The source data, callable, outputs and manifest are separate from the frozen photon-containment packet.
