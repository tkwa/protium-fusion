# Open-loop receiver and cadence bounds ✨

The corrected 64-zone pulse family supports a conditional collector screen around a radius of one billion metres. The limiting surface pulse is C13. The model requires at least **1.022 billion metres** for an allowed temperature rise of 1,000 K, or **1.412 billion metres** for 500 K. These are results within an illustrative carbon, surface-deposition model. Photon penetration can relax the skin-temperature limit; material response and isotope recovery still require a designed receiver. The B2 refinement now changes yield by about 5–6% at 128 zones, so these 64-zone receiver and process numbers remain provisional until the chosen refined pulse is propagated.

[collector.py](collector.py) and [collector.json](collector.json) provide the interface consumed by the station model. [collector-inputs.json](collector-inputs.json) preserves the corrected source records, initial profiles and cumulative prompt-radiation samples, with hashes of the original records and raw traces. It is sufficient to reproduce this screen without rerunning hydrodynamics or retaining the multi-megabyte traces.

## Pulse loading and radius

At a reference radius of 10⁹ m, with illustrative carbon conductivity 100 W/(m K), density 2,000 kg/m³ and heat capacity 1,000 J/(kg K):

| Stage | Ballistic kinetic + internal rise | Prompt rise | Delayed decay rise | Sum of separate maxima |
|---|---:|---:|---:|---:|
| B2 / C12 | 257.03 K | 450.30 K | 89.74 K | 797.07 K |
| C13 | 318.76 K | 727.81 K | 1.12 K | 1,047.69 K |
| A / N14 | 187.88 K | 314.86 K | 17.82 K | 520.56 K |
| B1 / N15 | 1.26 K | 125.41 K | <0.00001 K | 126.67 K |

The semi-infinite conduction solution includes no surface cooling. Separate maxima are added conservatively within that linear model. Prompt photons are placed at the surface; real absorption at depth may substantially reduce the short flash's skin-temperature peak. The grey pulse calculation supplies no escaping photon spectrum with which to qualify that improvement. The beta-decay term places the heat uniformly over the receiver area; decay location and charged-particle collection remain apparatus requirements.

Ballistic flight times scale with radius. Its heat flux scales as R⁻³, so its transient temperature scales as R⁻²⋅⁵. Prompt and fixed-lifetime decay heat scale as R⁻². For compatibility with the station interface, `prompt_peak_rise_K` includes both prompt and delayed-decay maxima; the two components are also recorded individually. Absorber mass and radiating area scale as R². Prepared pulse storage remains fixed, while the source's required peak power increases approximately as R⁻¹.

The complete prepared source peak is **3.1228 × 10²³ W** at R = 10⁹ m and 60% driver efficiency. This includes the exact initial EOS internal energy of each shell, assigned to the same launch-time bins as its kinetic energy. The separate kinetic-only peak is 3.0406 × 10²³ W. This is a required equipment rating; it does not establish a realizable route to the prepared flow. Switching, focusing, internal preparation and the warm-cavity transit must all fit the installed equipment and power allowances.

## Timing and the slow core

Timing includes the incoming flight before focus and the outgoing flight after focus. The interpolation assumes constant mass per velocity interval within each resolved shell. Exact quantiles are taken within that distribution; kinetic-energy quantiles use its cubic cumulative energy.

At R = 10⁹ m, retaining 99.9% of each leg's mass gives combined occupied times of 3,961 s for B2, 4,362 s for C13, 4,357 s for A and 4,856 s for B1. Weighted by the stationary process's shot counts, 80% cavity occupancy permits **4.1605 × 10³² He/s**. The remaining 0.1% slow mass requires protected staging or active cleanup. Discarding it exceeds the process's approximately 0.0472% allowable CNO loss per stage.

The `full_resolved` sensitivity follows every mass-weighted shell representative and permits **4.8971 × 10³¹ He/s** at the same occupancy. It is about 8.5 times slower. The central shell contains 3.815 × 10⁻⁶ of the pulse mass and its continuous velocity profile reaches zero. Consequently even this finite representative-shell schedule requires an active slow-core pickup. Literal 100% ballistic clearance takes no finite time. Neither schedule certifies collision-free recurrence of the unmodified continuous pulse.

The later B2 eight-sound-time sensitivity is retained separately. Its remaining internal-energy fraction falls below 1%, with about 1.9% mass-weighted velocity growth after the four-time baseline. This bounds one part of the ballistic approximation; source mesh refinement remains a separate issue.

## Thickness, heat extraction and material return

At a density of 2,000 kg/m³, areal masses 1,000, 5,000 and 10,000 kg/m² mean thicknesses 0.5, 2.5 and 5 m. A sufficient semi-infinite criterion, L ≥ 4√(αt), permits about 313, 7,813 and 31,250 s respectively at conductivity 100 W/(m K); these times halve at conductivity 200. At 5,000 kg/m², the largest pulse's uniform bulk rise is only about 11.8 K. The much larger skin transient, the delay before heat removal and the actual cooling-channel geometry therefore matter separately. A whole slab maintained by conduction alone also needs a gradient qL/k, so cooling hardware must be distributed near the receiving layers.

[NIST's graphite data](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html) give mass attenuation coefficients of 0.001959 and 0.001575 m²/kg at 10 and 20 MeV. A 1,000 kg/m² slab leaves 14.1–20.7% of those incident photons uncollided; 5,000 kg/m² reduces that to 0.0056–0.038%. These are interception checks. Scattering, secondary escape, ports and the actual spectrum determine total energy retention. The current family screens 5,000 and 10,000 kg/m²; a thinner collector needs an explicit transport and supplemental shielding solution.

At 5 AU, an emissivity-0.9 outer sphere can reject internally generated heat at 591.6, 2,079.1 or 6,600.3 W/m² at 330, 450 or 600 K after full solar absorption. A 40% engine operated at 80% of its internal Carnot efficiency needs finite hot and cold approaches. At a 1,000 K hot reservoir and 450 K radiator the optimal internal temperatures are 952.28 and 476.14 K, requiring total exchanger conductance **0.04391 W/K per watt of hot input**. The JSON includes the full 1,000/1,500/2,000 K × 330/450/600 K table. The 1,000/600 K combination is infeasible under this efficiency assumption. These equations establish heat-flow requirements; the station still allocates exchanger, manifold and pump mass.

Warm walls also illuminate incoming cold pellets. A 1 cm radius, 2,000 kg/m³ sphere with unit absorptivity absorbs about 25 MJ/kg during 3,000 s inside a 1,000 K black cavity, or about 129 MJ/kg at 1,500 K. Reflective coatings, shielded launch paths, staged preparation or a different receiver topology must preserve the intended incoming material state.

The receiver must recover implanted CNO as well as collect heat. Each large pulse deposits only about 9 × 10⁻⁷ kg/m² of fuel on a 5,000 kg/m² substrate. Stripping bulk substrate to retrieve sparse catalyst can therefore greatly exceed the nominal fuel-processing flow. Closed substrate reuse, isotope separation and sputtered-material recapture belong in the process and auxiliary budgets.

At the reference mass999 cadence, fresh hydrogen flow is 2.789 × 10⁶ kg/s. Reserving 10% or 20% of the overall 0.007 carbon/hydrogen allowance permits only **4.90 × 10⁻⁹ or 9.81 × 10⁻⁹ kg/m²/year of net wall-carbon loss**. Refurbishing a 5,000 kg/m² wall at 10⁻⁴/year circulates 0.5 kg/m²/year; its substrate must return at approximately 99.9999990% or 99.9999980% if losses consume those allowances. Gross retired material can be internally recycled. Initial construction material is capital stock. Carbon lost from it during operation needs its own allowance rather than silently entering the fuel balance as free makeup carbon.

## Reproduction

Run from this folder with the research Python dependencies available:

```sh
PYTHONPATH=/tmp/protium-python PYTHONDONTWRITEBYTECODE=1 python collector.py
PYTHONPATH=/tmp/protium-python PYTHONDONTWRITEBYTECODE=1 python collector-check.py
```

`--refresh-inputs` reconstructs the compact inputs from raw pulse traces and should only be used when deliberately replacing those inputs. The source records must match the canonical hashes in `process.json`. The check script tests the analytic rectangular-flux solution, radius scaling, mass and energy quantiles, source EOS reconstruction, integrated ballistic energy, cycle occupancy and current input/code hashes. These checks do not qualify irradiation lifetime, multi-pulse material recovery, cold-flow preparation or multidimensional pulse stability.

## Finite-depth prompt sensitivity

[collector-depth.py](collector-depth.py) and [collector-depth.json](collector-depth.json) resolve conduction from an exponential volumetric deposition profile, q(t) μ exp(−μx), using the saved prompt waveform. For an insulated semi-infinite solid the exact surface impulse response is μ erfcx(μ√(αt))/(ρc). Analytic checks recover both the surface-flux limit and the shallow-conduction, local heat-capacity limit. The script uses the same [NIST graphite table](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html) to set two illustrative depths: total attenuation and energy absorption. Both assign all prompt energy to the chosen exponential profile. A transport solution would instead follow scattering, secondary particles and escaping energy.

With the more concentrated total-attenuation depth, the C13 prompt temperature rise at R = 10⁹ m becomes 649.5, 298.1, 97.6, 13.8, 1.51 and 0.894 K for monoenergetic photons at 1, 3, 5, 10, 30 and 100 keV. Assigning each complete pulse to those bands gives corresponding minimum radii of 0.986, 0.809, 0.694, 0.645, 0.641 and 0.641 billion metres for a 1,000 K total transient allowance. At the harder bands, material arrival and delayed decay replace the prompt flash as the main limitation.

A mixed case makes the compactness requirement more useful: at R = 0.70 billion metres, at most 13.2% of prompt energy may use the original surface-deposition kernel when the remainder uses the 10 keV total-attenuation depth. At 0.65 billion metres, that allowance falls to 1.63%. These are sufficient bounds obtained by summing separate maxima. The required spectral condition applies to the pulse family; it is not an established emitted spectrum. The grey hydrodynamic model supplies integrated luminosity, so the conservative surface baseline remains unchanged. This sensitivity identifies photon-spectrum and deposition calculations that could justify a smaller collector.

Run the extra sensitivity with:

```sh
PYTHONPATH=/tmp/protium-python PYTHONDONTWRITEBYTECODE=1 python collector-depth.py
```

## Reserved-carbon station input and finite-wall accumulation

[collector-scenarios.py](collector-scenarios.py) produces [collector-scenarios.json](collector-scenarios.json) and [station-inputs.json](station-inputs.json) for the corrected 64-zone process with fresh fuel C:H = 0.0056. The remaining 0.0014 of the total 0.007 allowance is reserved for net structural carbon losses. These variants refresh the per-helium pulse occupancy and surface loading; they preserve the original `collector.json` and `process.json`.

[collector-periodic.py](collector-periodic.py) and [collector-periodic.json](collector-periodic.json) address cumulative heating for R = 1.1 × 10⁹ m, 5,000 kg/m² carbon, 80% mass999 occupancy and a back-face sink maintained at 1,000 K. The mean absorbed receiving flux is 107.695 W/m². Across the 2.5 m slab at conductivity 100 W/(m K), that establishes a mean front-face offset of 2.692 K. The wall retains finite heat capacity and conducts its absorbed heat to the maintained sink; no reset of its temperature is granted.

The transient bound requires at least one day between any two large A, B2 or C13 pulses, while the smaller B1 pulses continue. The large stages need about 78.48 shots/year, using 21.49% of the daily guard capacity. Cold routing and stock buffers must realize that schedule. One day's entire gross fuel flow is 3.87 × 10¹² kg; this is an explicit additional stock scale, not a zero-residence assumption.

The front-temperature impulse response of the slab is

`h(t) = 2α/(kL) × Σ exp[−α((n+½)π/L)²t]`.

It is positive and decreases with elapsed time. For each previous large pulse, assigning its entire heat to the end of its first 20,000 s bounds its remaining early contribution from above. Separate ballistic t⁻⁴ and beta-exponential envelopes cover its later heat. The B1 history is bounded using its minimum complete-cavity spacing of 5,341 s, which permits irregular scheduling around large pulses; its mean interval is 6,774 s. The resulting prior-pulse background is at most 7.091 K.

Even allowing the largest current pulse and a B1 pulse to attain their isolated maxima simultaneously, the front temperature is bounded by **1,965.33 K**, with the back held at 1,000 K. This bound is independent of their relative phase within the stated spacing rules. Removing those rules and summing every stage's maximum would yield a much looser bound. The result remains a constant-property, one-dimensional cooling calculation; nonlinear surface lifetime and implanted-fuel recovery remain apparatus requirements.

The grey hydro boundary uses a bolometric flux proportional to the outer cell's radiation energy divided by a grey optical-depth factor. It does not determine frequency-dependent thermalization or the escaping photon spectrum. Choosing a Planck spectrum at the material temperature or flux-defined effective temperature would add a spectral closure. The finite-depth photon calculation therefore remains a sensitivity rather than a replacement for the surface-deposition baseline.

## Completed B2 refinement variants

The separate 128- and 256-zone B2 variants retain the other three stages at 64 zones and use matched initial prepared energy. Their inputs and results are saved under `collector-inputs-b2-128-energy-matched.json`, `collector-b2-128-energy-matched.json` and the corresponding `b2-256-energy-matched` names. [collector-refined.py](collector-refined.py) rebuilds either variant from the completed raw trace and its process scenario. Original 64-zone files remain unchanged.

| B2 zones | B2 material rise | Prompt rise | Delayed rise | Complete source peak at R = 10⁹ m |
|---|---:|---:|---:|---:|
| 64 | 257.03 K | 450.30 K | 89.74 K | 3.1228 × 10²³ W |
| 128, matched | 299.91 K | 394.09 K | 94.87 K | 3.1430 × 10²³ W |
| 256, matched | 288.11 K | 336.52 K | 97.89 K | 3.2653 × 10²³ W |

C13 continues to set the same all-surface isolated-pulse radius threshold. The 256-zone process's mass999 occupied time is 1.921703 × 10⁻³³ s per helium at R = 10⁹ m. Its complete prepared-source peak is 4.56% higher than the 64-zone value. The representative-shell `full_resolved` timing is inherently resolution-dependent as the central cell approaches zero velocity; it remains a sensitivity rather than a converged all-particle clearance time.

[station-inputs-b2-256-energy-matched.json](station-inputs-b2-256-energy-matched.json) selects the refined B2 process with fuel C:H = 0.0056. [collector-periodic-b2-256-energy-matched.json](collector-periodic-b2-256-energy-matched.json) repeats the finite-wall check: mean receiving heat is 107.533 W/m², mean conduction offset 2.688 K and accumulated background at most 7.590 K. The same large-pulse guard gives a front-temperature bound of **1,965.82 K** at R = 1.1 × 10⁹ m. Its large stages require 76.16 shots/year and a one-day gross-flow buffer would hold 3.80 × 10¹² kg.

The compact finite-depth point requires more care. At R = 0.70 × 10⁹ m and a 13.18% surface-like soft component, the 128-zone B2 pulse exceeds the isolated 1,000 K rise allowance, while the 256-zone pulse nearly meets it. That exact point has neither refinement robustness nor demonstrated repeated-pulse headroom. A larger radius or smaller soft component is needed for a concrete compact operating choice. The conservative all-surface, 1.1 × 10⁹ m case retains thermal margin across these refinements. The B2 nuclear gain still changes by approximately 3% between 128 and 256 zones, so completion of the runs does not establish mesh convergence.

Example refinement reproduction, after the relevant process family and raw trace exist:

```sh
PYTHONPATH=/tmp/protium-python PYTHONDONTWRITEBYTECODE=1 python collector-refined.py --family b2_256_energy_matched
PYTHONPATH=/tmp/protium-python PYTHONDONTWRITEBYTECODE=1 python collector-periodic.py --station-inputs station-inputs-b2-256-energy-matched.json --collector-inputs collector-inputs-b2-256-energy-matched.json --output collector-periodic-b2-256-energy-matched.json
```
