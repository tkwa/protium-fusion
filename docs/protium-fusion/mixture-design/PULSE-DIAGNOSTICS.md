# ✨ Saved-pulse pressure, expansion, radiation and isotope diagnostics

This packet reduces 19 existing completed traces: the four original 64-cell mixtures, four optimized 64-cell mixtures, two explicit source-profile responses, and nine available 128-cell comparisons. The frozen manifest identifies every case. It introduces no hydrodynamic rerun or changed nuclear physics. The mixture optimum and source choices remain conditional on the existing spherical model and its nuclear data. Updated nuclear-rate families can be reduced separately without replacing this packet.

The most useful findings are that C13 hydrogen enrichment changes both the reaction routing and pressure history; the mapped C13 source develops a small central hotspot despite its lower average pressure; and B1 has a substantially larger escaped-photon share than the three denser stages. Trapped radiation provides a large fraction of the modeled pressure and internal-energy reservoir. Energy already present in that reservoir is counted once.

## Definitions and reconstruction

`pulse_diagnostics.py` reads each trace's actual composition, cell masses, boundary radii, velocities, temperature normalization and energy normalization. Spherical cell volume is `4π(r_out³−r_in³)/3`; density is mass/volume. Pressure and internal energy are reconstructed with the original `FastEOS`, including the electron/positron table, ideal ions and LTE radiation. The isotope clip and baryon renormalization match the hydro EOS interface. These definitions are implemented in [fast_eos.py](../fast_eos.py) and [radial_completion.py](../b2-completion/radial_completion.py).

Mass-mean pressure is `Σm_i P_i / M`; volume-mean pressure is `ΣV_i P_i / V`. Their maxima generally occur at different times. Volume-mean density is `M/V`, whereas mass-mean density is `Σm_i ρ_i/M`. Half- and ninety-percent-mass radii interpolate enclosed volume within each constant-density cell. Kinetic energy uses the actual staggered boundary-node masses from the solver. The reported RMS velocity is `sqrt(2K/M)` under that lumped-mass convention.

Pressure tables report thermodynamic EOS pressure. The reducer separately reconstructs the artificial viscosity used in the existing trajectory, including its volume-compression condition. That numerical shock stress enters the original momentum and work equations; its heating is already included in the trajectory. At peak mass-mean thermodynamic pressure it contributes an additional 5.3–8.9% of that pressure in the optimized/source 64-cell cases. This is one reason spatial refinement remains relevant.

All pressure, density and radius extrema here are **maxima or minima over saved samples**. Ordinary mixture traces have 400 samples, and the two source-response traces have 120. The metadata temperature maximum uses every accepted solver state; it is not a continuous-time bound. Across these cases the saved-sample temperature peak is at least 99.44% of the metadata peak, which is a useful sampling check but does not establish similarly small errors in sharp pressure features.

## Compression and pressure

The first eight rows compare mixtures at matched prepared pulse energy. The last two use the actual mapped source profiles, including the controlled preheat chosen for A. Pressure and density peaks are separate time extrema; they should not be combined into an invented uniform thermodynamic state.

| Case | Peak mass-mean P (Pa) | Peak mass-mean density (kg/m³) | Minimum half-mass radius (m) | Deposited nuclear / prepared energy |
|---|---:|---:|---:|---:|
| A original64 | 4.454e+20 | 4.842e+07 | 26.717 | 8.19634 |
| B1 original64 | 4.512e+18 | 7.056e+05 | 24.588 | 6.38197 |
| B2 original64 | 9.502e+20 | 6.674e+07 | 24.575 | 10.55503 |
| C13 original64 | 4.577e+20 | 2.869e+07 | 30.951 | 15.53898 |
| A optimized64 H3 | 4.615e+20 | 3.894e+07 | 28.822 | 10.28882 |
| B1 optimized64 H2 | 4.521e+18 | 5.97e+05 | 26.196 | 7.54562 |
| B2 optimized64 H2p5 | 9.65e+20 | 6.636e+07 | 25.039 | 10.79986 |
| C13 optimized64 H4 | 6.255e+20 | 2.84e+07 | 31.533 | 23.65567 |
| A source64 | 4.63e+20 | 3.911e+07 | 28.781 | 10.25252 |
| C13 source64 | 4.923e+20 | 2.46e+07 | 32.860 | 22.88617 |

Hydrogen enrichment raises peak mass-mean pressure by about 3.6% in A, 0.19% in B1, 1.6% in B2 and 36.6% in C13 at 64 cells. It lowers peak mass-mean density in A and B1, slightly lowers it in the other two, and increases their minimum half-mass radii. Thus the gain improvement is not simply greater compression of the same mixture.

The mapped C13 H4 source reduces peak mass-mean pressure from `6.255e20` to `4.923e20 Pa`, while its highest sampled cell pressure rises from `2.481e21` to `7.018e21 Pa`. The latter occurs at 20.483 µs in the innermost cell: radius 0.393 m, density `1.696e8 kg/m³`, temperature 1.983 GK, and enclosed mass fraction `3.815e−6` (3.815 ppm). At that time 0.9476% of pulse mass exceeds half this maximum pressure. The optimized-profile maximum instead occurs at 5.27–5.96% enclosed mass. This distinction matters for local EOS, reaction and multidimensional checks even though the full-cycle performance change of the source family is small.

For the chosen source-response family, the independently maximized volume-mean pressures are `7.50e19, 7.12e17, 1.38e20, 9.04e19 Pa` for A, B1, B2 and C13. Their maximum cell pressures are `1.74e21, 1.86e19, 3.89e21, 7.02e21 Pa`. The large difference between these averages reflects strongly nonuniform shells.

## Radiation and pressure work

At peak mass-mean thermodynamic pressure:

| Selected case | Ion share of P | Electron/positron share | LTE radiation share | Centre-to-surface grey optical depth |
|---|---:|---:|---:|---:|
| A source64 | 12.75% | 39.33% | 47.92% | 1.65e+06 |
| B1 optimized64 H2 | 6.24% | 18.75% | 75.00% | 1.91e+04 |
| B2 optimized64 H2p5 | 7.48% | 39.53% | 52.99% | 2.18e+06 |
| C13 source64 | 8.12% | 25.89% | 65.99% | 1.41e+06 |

The radiation energy reservoir is `ΣaT⁴V`, and its pressure is `aT⁴/3`. Both are already part of the total EOS. At the time of maximum integrated radiation energy, radiation accounts for about 66.5%, 86.4%, 68.4% and 80.5% of total internal energy in those four selected cases. Electron/positron internal energy includes degeneracy and pair-creation energy; interpreting all internal energy as random kinetic energy of particles would be incorrect.

The simulation uses constant-grey opacity `κ=0.001 m²/kg`, flux-limited intercell transport and the outer boundary flux

`F_out = c a T_out⁴ / [4 + 3 κ ρ_out Δr_out / 2]`.

The reducer reproduces that boundary luminosity. Integrating its sampled values agrees with the solver's cumulative escaped radiation to 0.33% or better across this packet. The cumulative state is the preferred energy ledger. Optical depth and reservoir/luminosity histories are diagnostics of this assumed grey transport law; they do not measure frequency-dependent absorption or establish LTE.

| Selected case | Escaped photons / prepared energy | Escaped photons / (prepared + deposited nuclear) | Peak trapped radiation / prepared energy |
|---|---:|---:|---:|
| A source64 | 0.969% | 0.0861% | 6.086 |
| B1 optimized64 H2 | 50.684% | 5.9310% | 6.309 |
| B2 optimized64 H2p5 | 0.835% | 0.0708% | 6.796 |
| C13 source64 | 3.145% | 0.1317% | 15.586 |

B1's 50.7% of prepared energy corresponds to 5.93% of prepared-plus-nuclear energy, because its deposited nuclear energy is 7.55 times the prepared input. The original B1 mixture gives 5.79%. The corresponding selected A, B2 and C13 photon shares are about 0.086%, 0.071% and 0.132%. These are photons that leave the simulated plasma during the archived interval; a surrounding collector can receive them. Radiation still trapped at the final state is tabulated separately and can evolve during later expansion.

**Radiation production remains an imposed LTE closure.** The model does not compute a frequency spectrum, photon-number production, separate electron/ion temperatures, line absorption, non-LTE populations, or a frequency-dependent thermalization time. This packet therefore adds no independent bremsstrahlung power. A classical optically thin gross-emission formula would count emission before reabsorption and would not be a second loss term. The proper next extension is a frequency-dependent transport/thermalization calculation on these exact mixtures and trajectories, followed by coupled hydro only if it changes the thermal history appreciably.

## Expansion and final energy partition

The initially prepared outer radii are about 171 m (A/B2/C13) and 191 m (B1). Minimum surface radii are about 75.1, 66.0, 65.5 and 85.9 m for the selected family; minimum half-mass radii are much smaller because the outer layers are dilute. At the end of the actual recorded durations:

| Selected case | Duration (µs) | Surface radius (m) | Surface speed (Mm/s) | Equivalent RMS speed (Mm/s) | K / input | Internal / input |
|---|---:|---:|---:|---:|---:|---:|
| A source64 | 68.395 | 613.46 | 15.651 | 8.135 | 10.7127 | 0.5301 |
| B1 optimized64 H2 | 82.278 | 445.68 | 10.287 | 6.553 | 7.6003 | 0.4385 |
| B2 optimized64 H2p5 | 63.063 | 591.39 | 18.317 | 9.012 | 11.2710 | 0.5205 |
| C13 source64 | 70.004 | 998.56 | 22.563 | 11.561 | 22.7088 | 1.1459 |

The figures show coherent kinetic energy rising as the hot, radiation-supported pulse expands and does work. The archived final state remains a plasma near the origin; it is not the eventual arrival state at the billion-metre receiver. Expansion, recoil, recombination and charge-state assumptions for that later flight belong to the existing debris/collector packets.

## Actual isotope routing and minor channels

The summary stores initial, final and sampled peak inventories for every one of the 27 isotopes. The history file retains their time series. Values below are final hot-pulse nuclei per initial catalyst nucleus, before the long cold decay and recovery interval:

| Selected case | Main remaining/output nuclei per initial catalyst | Alpha particles per initial catalyst |
|---|---|---:|
| A source64 | n14: 0.235367, o15: 0.764633 | 6.34911e-08 |
| B1 optimized64 H2 | n15: 0.243803, c12: 0.755968, n13: 0.00013787, o16: 9.10255e-05 | 0.756106 |
| B2 optimized64 H2p5 | c12: 0.0860873, n13: 0.0340468, o14: 0.879866 | 2.47168e-11 |
| C13 source64 | c13: 0.0764637, n14: 0.252914, o15: 0.670622 | 9.46167e-08 |

The major extra route in the H-rich C13 pulse is continued proton capture on newly formed N14 to O15. Original C13 H1 ends with about 0.817 N14 and 0.0316 O15 per initial C13; optimized H4 ends with 0.222 N14 and 0.713 O15; the mapped H4 source gives 0.253 N14 and 0.671 O15. Consequently a scalar desired-product yield or nuclear gain cannot substitute for the full cycle's material routing. The existing process integration uses these endpoints.

In this network, the clearly resolved minor B1 routes leave about `1.38e−4` N13 and `9.10e−5` O16 per initial N15 in the optimized mixture. The modeled F-and-heavier inventories and prompt neutron inventories are much smaller. Their tiny tabulated values are numerical diagnostics, not certified engineering production limits: rate uncertainties, omitted channels, small negative trace inventories and cancellation can dominate them. These pulses also exclude later receiver/material interactions.

The optional reaction-flow reducer uses the archived ReacLib coefficients to estimate gross turnovers and forward-minus-reverse flows on the saved samples. Its signed version reproduces the solver's polynomial convention; a clipped version exposes nominal nonnegative gross flows. The resulting quadrature is less accurate than stored inventories. For example, the sparse C13 source trace has a **7.1% nuclear-energy quadrature discrepancy** and a maximum isotope closure discrepancy of **0.111 nuclei per initial catalyst**. Those sampled turnovers identify mechanisms but cannot support precise branching fractions. The ordinary optimized C13 trace has 1.58% energy discrepancy; most other 64-cell traces are within 0.3%. Gross consumption of clipped numerical trace neutrons is particularly unreliable. Exact small side-channel bounds would require integrated reaction counters or denser accepted-state output with the revised nuclear data.

## Numerical evidence and scope

All 19 traces completed their intended physical durations (the interval is based on four reference sound times, with mixture runs preserving physical duration); no guarded A continuation is included. The reducer passes **266/266** checks: positive cell volume and finite EOS pressure; exact endpoint mass, kinetic/internal/radiated/weak-energy and time reconstruction; baryon/positivity diagnostics; monotone cumulative outward radiation; the full sampled energy ledger; and nonnegative artificial stress that vanishes in expanding cells. The largest sampled total-energy residual is `6.42e−6` of prepared input. An independent read-only review also matched the spherical-volume, EOS unit conversion, nodal kinetic-energy, artificial-pressure and outer-flux formulas to the production solver. These checks establish implementation consistency within the model; the physical closures retain their stated limitations.

The 128-cell comparisons remain distinct scenarios. They change peak average pressures and nuclear gains, and in B1 alter both localized temperature peaks and photon escape. For example, B1 H2's escaped share rises from 5.93% at 64 cells to 6.58% at 128. Original B2's peak mass-mean pressure rises from `9.50e20` to `1.23e21 Pa`. These are discretization effects, not a statistical confidence interval or proof of convergence. The 19-case manifest intentionally excludes a C13 H4 128-cell result that was not present when this manifest was selected; a later available case should be added in a new manifest/output directory. New nuclear-rate and Coulomb-EOS variants likewise require explicit provenance and the corresponding physics when reconstructing pressure and flows.

## Artifacts and reproduction

- [Summary, isotope inventories, extrema and provenance](pulse-diagnostics/summary.json)
- [Compact comparison table](pulse-diagnostics/comparison.csv)
- [Compressed time histories](pulse-diagnostics/histories.json.gz)
- [Frozen 19-case manifest](pulse-diagnostics/manifest.json)
- [Pressure, radius and luminosity figure](pulse-diagnostics/pulse-dynamics.png), also PDF
- [Kinetic/internal/radiation energy figure](pulse-diagnostics/pulse-energy-partitions.png), also PDF

Each case records hashes of the metadata, full raw trace, FastEOS, electron table, nuclear module and numerical nuclear coefficient arrays, plus the pinned pynucastro version. The summary records the reducer hash. This packet's raw traces remain in the manifest's existing local locations; it does not duplicate or alter them.

Run from the website checkout, with the restored pulse runtime:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/tmp/protium-python OPENBLAS_NUM_THREADS=1 NUMBA_CACHE_DIR=/home/tkwa/src/protium-fusion/output/pulse-diagnostics-numba MPLCONFIGDIR=/tmp/protium-mpl /tmp/protium-python-interpreters/cpython-3.10.20-linux-x86_64-gnu/bin/python3.10 docs/protium-fusion/mixture-design/pulse_diagnostics.py --manifest docs/protium-fusion/mixture-design/pulse-diagnostics/manifest.json
```

A new nuclear-rate family can use `--manifest` and a separate `--output`. Set an entry's `reaction_flows` to `false` unless its exact original network arrays are available to this reducer. A `rate_update` record explicitly rejects baseline sampled reaction flows even if its unchanged source-module hash matches. Radiation-off records have zero reconstructed surface luminosity. Pressure reconstruction also requires the correct EOS; this version explicitly supports the archived uncorrected FastEOS family, with Coulomb-EOS cases excluded.
