# ✨ C13-rich polymer source at H:C=2

The completed 128-zone calculation supports replacing the methane-rich C13 source with isotope-rich CH2 polymer when source simplicity is valuable. A constant contact-collision entropy reproduces the nominal H:C=2 pulse gain to **0.0056% at unchanged total prepared energy**. This is a source-profile comparison within the existing spherical model. The separate capsule screen and three-dimensional assembly work determine whether physical packets can deliver that profile.

The original methane source family remains unchanged. `stage-map.json` supplies one explicit C13 replacement for a separately named integrated family; it does not silently replace A, B1 or B2. The final-family agent will use the full isotope endpoints and actual receiver reduction, rather than multiply the old plant output by a gain ratio.

## Source state and controlled collision

The new FEOS material is C13H2, with atomic masses13 and1, number ratio1:2, reference density **1011.642857 kg/m³ =944.2×15/14**, reference temperature300 K and chosen polymer bulk modulus2 GPa. The density convention preserves the molecular volume of the old natural-CH2 reference. It is an explicit isotope scaling, not a measurement of the proposed capsule.

The delivered hot density grid, total mass, initial composition and physical integration duration are taken from the exact parent nominal128 C13H2 case. The inner0.1% seed keeps its existing hot state. For the remaining parcels, a single normal collision-energy parameter is selected so that the FEOS shock entropy equals the mass-weighted entropy of the nominal profile. Each parcel is then mapped isentropically to its original hot density. A single radial-velocity factor restores the exact total internal plus kinetic energy.

| Quantity | Result |
|---|---:|
| Contact heat per core mass | 3.46666 GJ/kg |
| Contact share of prepared energy | 0.117677% |
| Shock density | 4047.81 kg/m³ |
| Shock temperature | 467802 K |
| Shock pressure | 9.35109 TPa |
| Normal speed of each opposing stream | 83.2665 km/s |
| Pair half-angle, mass10/50/90% | 1.540 /2.061 /2.877 degrees |
| Pair half-angle, full modeled nonseed range | 1.194–5.263 degrees |
| Mapped/reference temperature, mass10/50/90% | 0.996633 /0.998389 /1.002148 |
| Mapped/reference temperature, full nonseed range | 0.996523–1.032228 |
| Radial velocity scale for exact energy | 1.000058472 |
| Nonseed mass requiring postshock rarefaction | 0 |

These collision angles follow a stated local energy allocation equal to the hot internal plus radial kinetic energy. They are a kinematic construction to test, not an outlet-count requirement or a demonstrated pressure-work trajectory. Controlled transverse energy pays for the prescribed contact entropy; it is not added on top of the original prepared energy. No extra external preheater or temperature clipping is used.

## Independent contact comparison and completed burn

`contact_primary_check.py` evaluates the same shock energy against the archived primary CH2 FPEOS table, with15u per C13H2 formula. It gives density4042.63 kg/m³, compression3.99608 and pressure9.35519 TPa. FEOS compression differs by0.129% and pressure by−0.044%. FPEOS temperature is520748 K, about11.3% above FEOS. The mechanical agreement is useful; temperature/ionization/entropy accuracy remains a material-model uncertainty. The primary table is not used to claim a validated subsequent isentrope. Its initial density differs from the chosen reference by only8×10⁻⁶.

The hot calculation uses the full27-species network, the frozen median-carbon update, FastEOS, opacity0.001 m²/kg, corrected compression-cap viscosity, rtol10⁻⁷ and atol10⁻¹². Its initial isotope abundances are exactly p=2/15 and C13=1/15. Neither structural carbon nor the small Al optical coat is hidden inside those abundances; the capsule packet reports their actual composition corrections separately.

| Matched128 result | Nominal prescribed profile | Polymer contact-mapped profile |
|---|---:|---:|
| Prepared energy | 3.3282663865×10²⁵ J | Same |
| Pulse mass | 1.1309733553×10¹³ kg | Same |
| Physical duration | 70.0035781 microseconds | Same, completed |
| Nuclear/prepared energy | 22.65872140 | 22.65998125 |
| Recorded target product/catalyst | 0.47860762 | 0.47845983 |
| Solver accepted-step peak T | 1.66578091 GK | 1.66597162 GK |
| Final energy residual/prepared | — | −2.48×10⁻⁶ |

Nine source/reference checks and fourteen independent saved-trace diagnostics pass. The checked quantities include mass, energy, duration, composition, rate identity and the corrected viscosity. The trace reducer's extrema are sampled extrema; the temperature above is the solver's accepted-step maximum. `audit.json` preserves the full diagnostics and hashes.

The earlier nominal128 whole-plant comparison found only a0.237% specific-output difference between C13H2 and C13H4 with identical other stages. That comparison remains separate from the integrated source family and the new receiver-thickness choices. The present result shows that the chosen polymer entropy mapping introduces no meaningful extra pulse-gain penalty. It does not select a final plant mass or capture fraction.

## Capsule connection and limits

The [capsule-survival packet](../capsule-survival/README.md) finds that the reference1cm,80%-filled C13 polymer absorbs about227 kJ/kg over the representative458s flight including fully retained field energy. Pre-cooling it makes the300K reference a plausible warm end of flight: the isotope-scaled primary20→300K enthalpy range is251–282 kJ/kg. A300K-start pellet instead needs a melting/expansion allowance. Manufactured absorptivity, finite surface gradients, seams, gas permeability and neutralization energy deposition remain actual qualifications.

A carbon shell carrying all charging pressure plus an illustrative1MPa internal allowance adds about0.201% core mass at a10GPa allowable, yielding H:C=1.9954; at1GPa it adds2.03%, yielding H:C=1.9543. These are explicit nearby compositions, not exactH2 inputs to this run. A solid polymer core carrying part of its own charging tension can reduce the separate shell. The Al coat is another small, unmodeled species contribution.

The contact model uses condensed polymer density. Capsule void collapse, shell fragments and their thermal histories have not been evolved into the hot snapshot. The incoming trajectory still needs the full pressure-work, geometry, neutralization and symmetry calculation. The useful next test is therefore a finite representative multi-stream sector with this cold material and controlled collision speed, compared against the actual same-energy gain and cycle response. No fixed95% yield threshold or particular angular-channel count is inferred from this successful scalar mapping.

## Reproduction and handoff

`run_source.py` takes the exact nominal prefix, a task-owned FEOS cache and a new output prefix. It writes initial profile and input metadata before integration. `package.py` verifies the completed result, copies metadata byte-for-byte into `records/`, and records the raw trace path/hash. `contact_primary_check.py` is a separate primary-table comparison. The median-carbon kernel hash is9e3252c7f70bb4fa6c931c4264f3df583fe2d6fe19b375c2e1f461e646aae865; table hash9166934c7f82868daf7ff95054c78f6078994ab1c0a8e53040a37c1d651f4b99.

The new exact stage map is `component-design/c13-polymer-source/stage-map.json`, relative to `docs/protium-fusion`. Raw output remains at `/home/tkwa/src/protium-fusion/output/c13-polymer-source/C13-H2-source128-carbon-median.*`. The upstream FEOS implementation and primary-data limitations are documented in [launcher.md](../launcher.md), [launcher-feos.json](../launcher-feos.json) and the [primary table source](https://militzer.berkeley.edu/FPEOS/).
