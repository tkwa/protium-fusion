# Compact opaque station: selected 2049 boundary ✨

This is a conditional engineering allocation supporting a possible in-universe history. It is not a prediction from an engineered reactor simulation. Four small prescribed hot inflows completed (A, B1, C13 and the 3 m helium factory). The B2 ledger uses the preserved late solver-failure endpoint; its late product abundance nearly saturates, but neither its final expansion nor its receiver history is validated. The numerical status is recorded in every output packet. Controlled preparation, material recovery and hardware masses remain chosen future engineering achievements.

The useful redesign is to run the helium factory at 3 m reference radius. Its mass and pulse energy are a thousand times smaller than the old 30 m factory, while carbon conversion is ten times lower. Requiring better catalyst recovery pays for the poorer factory yield and removes the rare enormous storage discharge. All factory equipment and energy are counted in this station; there is no allocation to future customers.

## Complete selected point

| Quantity | Selected value |
|---|---:|
| Whole installed project mass, including explicit reserve | 1.00 × 10^24 kg |
| Useful electrical output after auxiliaries, support and replacement | 1.00 × 10^21 W |
| Useful specific output | 1.00 mW/kg |
| Main catcher radius | 1.00 × 10^9 m |
| Electrical-to-prepared driver efficiency | 60% |
| Captured-heat conversion efficiency | 40% |
| Non-neutrino energy capture | 99% |
| Catalyst loss per pass | 10^-5 |
| Outer radiator temperature/emissivity | 330 K / 0.9 |
| Annual replacement fraction / production energy | 10^-4 / 100 MJ/kg |

The linked cycle supplies 24.7222 MeV of non-neutrino nuclear heat per net helium atom, with 3.58968 MeV of delivered preparation. The small factory accounts for 1.28947 MeV of that preparation. Nuclear/prepared gain is 6.887. Radiation escaping the modeled plasma is intercepted by the closed station: its 0.267777 MeV/He is included before applying the station's single 1% energy-escape allowance. This avoids subtracting prompt photons from useful heat while silently depositing them in the radiator load.

With `H` denoting total non-neutrino nuclear heat and `U` preparation, the preauxiliary electrical margin is `0.40 × 0.99 × (H+U) − U/0.60 = 5.22870 MeV/He`. Prepared energy is recovered exactly once.

| Power flow | W |
|---|---:|
| Captured heat | 8.002 × 10^21 |
| Gross electricity | 3.201 × 10^21 |
| Driver electricity | 1.708 × 10^21 |
| Other auxiliaries | 1.493 × 10^20 |
| Active support loss | 2.656 × 10^19 |
| Replacement manufacture | 3.169 × 10^20 |
| Useful export | 1.000 × 10^21 |
| Total rejected heat | 5.977 × 10^21 |

The auxiliary reserve is 10% of preauxiliary electrical surplus. Replacement power is substantial and has not been hidden in that reserve. At 330 K the outside radiator emits 605.2 W/m². Its required area, 9.876 × 10^18 m², occupies 78.6% of the collector's outside area. Distributed engines and local thermal manifolds avoid billion-metre coolant trunks. A hot collector near 1000 K, 330 K radiating surface and finite internal engine temperatures can support the chosen 40% conversion; the independent receiver audit calculates finite exchanger approaches. The balance relies on insulation and controlled thermal bridges, explicitly allocated below.

## Mass and pulse service

| Installed allocation | kg |
|---|---:|
| Absorber, 10,000 kg/m² | 1.257 × 10^23 |
| Driver, guidance and peak-power hardware | 2.500 × 10^23 |
| Isotope separation, decay, recovery and working inventory | 1.000 × 10^23 |
| Structural guides, control and redundancy | 1.000 × 10^23 |
| Local heat exchange, manifolds, insulating bridges | 1.000 × 10^23 |
| Two maximum electrical pulses in storage, at 100 kJ/kg installed | 1.358 × 10^21 |
| All other computed components, including support streams, insulation, generation, conditioning, panels, manufacture and ten-year replacement stock | 3.596 × 10^21 |
| Installed design reserve for unresolved source, processing and service equipment | 3.194 × 10^23 |
| Total | 1.000 × 10^24 |

The large fixed allocations are requirements, not masses derived from a component design. They provide a concrete engineering boundary that the fictional builders must meet. A 10,000 kg/m² catcher can use carbon-rich receiving material and bulk thermal storage; the independent constant-property carbon screen is not a proof for a substituted laminate.

The station circulates about 2.29 × 10^8 kg/s across all stages, consumes 1.18 × 10^7 kg/s of hydrogen, and makes approximately 6.0 × 10^3 kg/s of replacement carbon. Its stage pulse intervals are A 1.59 days; B1 23.9 minutes; B2 2.77 days; C13 39.8 days; helium factory 18.3 hours. The average interval across all stages is 23.0 minutes. B1 has a smaller pulse mass because its prescribed density differs from the other 30 m stages.

The maximum electrical pulse is B2, 6.791 × 10^25 J. A two-pulse buffer stores 1.358 × 10^26 J. The independent initial-profile audit obtains kinetic-source electrical peaks of 2.415 × 10^23 W (A), 2.032 × 10^21 W (B1), 3.041 × 10^23 W (B2), 2.246 × 10^23 W (C13), and 4.418 × 10^23 W (He). These correspond to about 1.8 W/kg of the fixed driver allowance or 325 W/kg of the two-pulse storage. Internal-energy preparation and overlap need additional rating; even summing all five peaks is about 4.9 W/kg of driver allowance. Storage supplies peak power locally; the average driver demand is already charged to generation.

Cadence is an engineering discriminator. Roughly 99% of B1 outgoing kinetic energy arrives in 109–332 seconds. The corresponding source's 1–99% energy flight-time span is roughly 255–851 seconds. Most energy can therefore fit the 23.9-minute slot, with a small timing margin. The exact saved profile has a much longer low-energy, slow core tail. A practical source needs staged delivery or a redesigned central velocity profile, and protection from other shots' outgoing material. Merely firing the exact saved inflow every 24 minutes would not demonstrate a collision-free schedule. The fictional station's commissioning success includes this redesign and its measured gain; this calculation supplies its required throughput, not that trajectory optimization.

## Recovery and receiving surfaces

Completed small-pulse constant-property receiver screens at 10^9 m have upper sums of separate temperature-rise maxima of approximately 249 K (A), 94 K (B1), 607 K (C13), and 16 K (He3). B2 is rejected by the ballistic screen because its endpoint velocity profile permits shell crossing. Its continuing expansion and heat pulse remain unresolved. See `opaque-receiver-audit.json` for complete conditions and limitations.

99.999% catalyst return must include nuclear changes during stopping, implanted inventory, collection ports, sorting losses and unreturned gases. Geometric enclosure alone does not establish it. Appropriate development choices include segmented slowing before solid impact, light-element receiving layers, offset capture ports, recovery of secondary gas and isotope accounting across every stage. Their stopping spectra and transmutation yields are part of the chosen in-universe commissioning achievement. Ordinary heavy-element wall impacts cannot silently be assumed to meet the number.

Cold preparation also needs a thermal and materials design. At a 1000 K cavity, centimetre parcels with low infrared absorptivity survive the shorter flight more readily than the submillimetre example; the older driver screen quantifies the size/field/gradient trade. Conductive coatings enter the isotope ledger. The helium factory needs a preparation architecture whose container material does not consume more catalyst than the factory creates. These are concrete design tasks, not consequences of an unlimited intelligence assumption.

## Sensitivity and history

At catalyst loss 10^-4, the same small factory raises input to approximately 6.86 MeV/He. The 40% thermal point remains positive before auxiliaries, but at 1 mW/kg it would require 4.26 times the available 330 K outside radiator area and a much faster cadence. Increasing thermal efficiency to 45% reduces that to 1.67 times. Thus 99.999% return is a load-bearing process achievement for this selected compact station.

At 45% heat conversion with baseline recovery, margin rises to 6.630 MeV/He, radiator coverage falls to 59.2%, and the average pulse interval increases to 29.2 minutes. A 25% preparation-energy penalty at 40% conversion requires 104% of the baseline radiator area; a 50% driver requires 106%. These small thermal shortfalls could be addressed by a warmer radiator or extra area, while source/cadence changes still need engineering. None is an automatic physical failure.

A coherent selected history is: build the collection and source programme from 2043–45; conduct open pulse/merger work at the remote 120 AU range in 2048; commission the first enclosed specialist station near 5 AU in late 2049. By then the programme has selected smaller repeated pulses, a small local makeup factory, high-return receiving surfaces, efficient field-energy recovery and a revised low-velocity delivery schedule. The station is a strategic hydrogen-power capability and research asset rather than the cheapest source of electricity in a solar-swarmed system.

A 10^24 kg installation is 20.5% of Venus's mass. At 30 km/s, the straight cruise from Venus's orbit to 5 AU takes approximately 0.676 years. Carbon shipped from a works near 20 AU to the same station takes 2.37 years. Construction and sourcing begun in 2043–45 can therefore deliver material during 2047–49. These cruise times exclude launch, braking and nonradial trajectory allowances. Accelerating all installed mass once to 30 km/s costs 4.5 × 10^32 J of kinetic energy; at 50% launch efficiency, 9 × 10^32 J is 30.2 full-swarm days. That total counts both source routes once at the same speed, with their different distances affecting travel time. Receiving works and braking remain required. Fabrication at 100 MJ/kg costs 10^32 J, another 3.36 full-swarm days. Carbon lifting and processing are separately charged below.

Venus cannot supply 10^23 kg of carbon from its atmosphere. The narrative therefore needs a separate carbon source or a validated composite absorber. A concrete permitted choice is nearby outer-planet methane processing, supplying the carbon-rich receiving material while Venus supplies most structural mass. At an assumed 3% recovered-carbon yield, 1.26 × 10^23 kg of carbon requires about 4.2 × 10^24 kg of feedstock. That is a new industrial sourcing campaign, not a minor atmospheric skim. If lifting/processing energy averaged 10^8–10^9 J/kg of feed, its energy cost would be 4.2 × 10^32–4.2 × 10^33 J (14–141 full-swarm days); these are design parameters, not a measured extraction process. The carbon delivery follows the separate 20-to-5 AU route above. Its cruise acceleration is already included in the whole-installed-mass transport total; lifting and feedstock processing are additional.

Even the optimistic 100 MJ/kg construction energy takes about 3,170 years to repay at 1 mW/kg. The illustrative 30 km/s transport kinetic energy adds about 14,300 output-years before efficiency losses. This makes the selected early plant economically credible as a strategic demonstration in the established civilization, while giving a strong reason to optimize mass and direct recovery before widespread deployment.

## Reproduction and corrections

Run `opaque_station.py --root docs/protium-fusion --inputs docs/protium-fusion/reassessment-2049/small-stages --output <path>`, using the saved Python dependency environment. The JSON contains source hashes, solver status, seven sensitivity cases, electricity closure and six assertions. Physics preparation and receiver validation are outside those assertions.

This packet corrects the old `opaque.md` 30 m factory paragraph: the direct 10^9 m carbon-surface screen gives roughly **11,337 K** ballistic-plus-internal rise, not 567 K; its electric input at 60% is **4.045 × 10^28 J**, not 6.79 × 10^28 J. The 3 m redesign avoids both inherited figures. The earlier failed-station narrative and old eleven-Earth collector are superseded design choices, not scientific necessities.
