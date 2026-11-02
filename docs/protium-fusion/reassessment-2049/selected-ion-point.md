# Conditional ion operating point and startup inventory ✨

The newly selected accounting scenario uses **1.75 MeV of extra electrical processing per primary reaction** and closes at **272.37 GW net, 1.34000×10^15 kg, and 0.20326 mW/kg**. It retains 60% heat conversion, a 50 K radiator, 10^14 kg and 100 GW of additional integration reserve, fixed original general auxiliaries, and the full solar radiation budget at 80 AU. A rounded 0.2 mW/kg is therefore a consistent complete-plant operating target. The author may select a late-2040s broad-capture/prepared-acceleration system achieving the required burden as an explicit future engineering result. The arithmetic does not establish that hardware. The simple monochromatic CARA calculation does not provide its evidence: an independent audit found a resonance-capture bandwidth mismatch.

[The script](selected_ion_point.py) and [saved results](selected-ion-point.json) preserve the inputs, all 20 comparison cases, source hashes, complete mass breakdowns, startup inventories and 118 passing arithmetic checks. This small calculation uses the saved source family without a new plasma search or nuclear evaluation.

## Power and mass acceptance range

All rows below include the 10^14 kg / 100 GW reserve, full solar radiation at 80 AU and fixed original auxiliary draw.

| Extra electrical processing | Net output | Accounted plant mass | Net specific power |
|---:|---:|---:|---:|
| 1.0 MeV/primary | 595.67 GW | 1.33875×10^15 kg | 0.44494 mW/kg |
| 1.2 MeV/primary | 509.45 GW | 1.33908×10^15 kg | 0.38045 mW/kg |
| 1.5 MeV/primary | 380.13 GW | 1.33958×10^15 kg | 0.28377 mW/kg |
| **1.75 MeV/primary** | **272.37 GW** | **1.34000×10^15 kg** | **0.20326 mW/kg** |
| 2.0 MeV/primary | 164.60 GW | 1.34041×10^15 kg | 0.12280 mW/kg |

The underlying plant function assigns other auxiliaries as 10% of power remaining after the extra processing charge. Thus raising processing automatically reduces general auxiliaries. Preserving that behavior gives 561.18 GW and 0.41914 mW/kg at 1.2 MeV. Holding the original 172.269 GW general allowance fixed restores 51.664 GW of draw at 1.2 MeV. At the selected 1.75 MeV it restores 75.343 GW. The table holds the allowance fixed in every row. This is a budget-policy distinction, rather than a new physical loss.

At the selected table row, primary throughput is 2.68718×10^24/s. Extra processing consumes 753.434 GW; added integration consumes 100 GW; replacement manufacturing consumes 424.620 GW. The retained allocation also includes the original extraction, atomic-preparation, channel/end-transfer, neutron-capture and refrigeration costs. The complete specific-power figure includes radiators, the new integration reserve, manufacturing equipment and replacement inventory. No heat recovery credit is added for the new processing burden. Its waste heat remains in the total rejection calculation.

The retained 140 T nominal field, 137.15 T source field, 16 GPa hoop stress, 8 GPa compression stress, 1 TPa compression modulus, 95% heater/direct conversion, 99.9% potential/end transfers, 99.5% neutron capture and prescribed coil leakage remain coupled engineering requirements. Reducing export below 1 mW/kg does not itself qualify those components. The added reserve is a selected capacity envelope, rather than a demonstrated bill of equipment.

## The 80 AU radiator closes with an explicit external load

[NASA's irradiance value](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html) gives 1,361 W/m² at 1 AU and 0.212656 W/m² at 80 AU. The calculation assumes outward spherical emitting area, unit bolometric absorption of the incoming luminosity, and projected area equal to one quarter of emitting area. It credits no reduction from the swarm's solar interception: captured light may return as infrared waste. Directional utility beams and nearby large industrial emitters require siting controls in addition to this central-luminosity budget.

At 50 K, emissivity 0.9 and a 3 K background, emission capacity is 0.318954 W/m². Mean absorbed central radiation is 0.053164 W/m², leaving 0.265790 W/m² for plant heat. The completed row therefore has:

- Internal waste heat: **6.65356 TW**.
- Emitting area: **2.50331×10^13 m²**, at the retained 1 kg/m² panel allowance.
- Equivalent spherical radius: **1,411 km**.
- Absorbed external power: **1.33086 TW**.
- Total emission: **7.98442 TW**.

The global radiation and mass ledgers close independently. Relative to the same fixed-auxiliary point under a 3 K-only sky, the external load adds 4.182×10^12 kg and reduces net output by 1.325 GW through replacement work. It changes 0.20489 to 0.20326 mW/kg. A smaller assumed mass per emitting area is unnecessary to obtain the stated result.

An isothermal sphere needs sector-specific loading. The subsolar face can reject only 0.106298 W/m² of internal heat at 50 K, versus 0.318954 W/m² on the antisolar face. Their coolant loads therefore differ by about threefold. A uniform internal heat flux on every panel would not preserve a uniform 50 K surface. The geometry requires coolant distribution, selective surfaces, shades or a temperature distribution; the current arithmetic uses coolant distribution and has not designed its manifolds.

## Finite heat transfer and reserve size

The saved source/circulation equations credit 8.45483 TW entering the warm heat engine and 5.07290 TW of thermal-converter electricity. These include recirculated energy, so the heat input may exceed the 6.92593 TW net nuclear source. The engine itself rejects 3.38193 TW. Other electrical losses and the new burdens bring total internal rejection to the 6.65356 TW above.

At a selected performance of 85% of internal Carnot, a 60% engine between a 300 K source and a 50 K radiator can use internal temperatures 239.99 K and 70.58 K. Minimum hot-plus-cold exchanger conductance in the simple two-exchanger model is 3.052×10^11 W/K. A working-fluid cycle and the mass achieving that conductance remain component-design requirements.

The existing ideal-helium, two-leg cold-loop screen was rerun at the enlarged 1,411 km radius. For 0.1 MPa helium, a 10.292 K temperature span, 1 GPa pressure-wall stress and Darcy friction 0.01:

- Engine heat alone requires 3.335×10^13 kg of pipe wall and helium, and 17.35 GW pumping.
- Routing all 6.654 TW of internal waste through this cold loop requires 5.952×10^13 kg and 31.00 GW. This intentionally counts the loop's own friction twice because that work is already part of the fixed 100 GW integration reserve.

The second screen leaves approximately 4.05×10^13 kg and 69.0 GW of the added reserve for manifolds, heat exchangers, RF hardware, insulation, pumps and other integration. Existing equipment allowances also remain present. The idealized single-pair bore is 4.90 km, which plainly needs a branched layout at the small central source. Dividing the flows raises wetted area and friction; proximity of distributed engines can reduce routing. These screens establish the scale of a feasible accounting reserve, not a finished piping topology. The earlier central 30 K loop penalty cannot be applied unchanged to this 50 K architecture.

## Startup inventories and the recycling clock

The neutral decay stores contain **9.4646×10^6 kg** in steady operation, including 8.9860×10^6 kg of V49, 261,177 kg of Fe55 and 194,011 kg of Cr51. Here “cold” means removed from the MeV ion plasma. These stores belong in warm chemical-processing equipment, with their own decay heat rejected at useful process temperatures; they do not belong on the 4.2 K magnet circuit.

The saved nuclear network uses a V49 mean life of 1.30346 years, corresponding to a 330-day half-life. [DOE's isotope product data](https://www.isotopes.gov/sites/default/files/2021-02/V-49.pdf) give that half-life and Ti49 daughter. A [primary measurement](https://pubmed.ncbi.nlm.nih.gov/10724407/) obtained 328.6±3.3 days, consistent with the scale. The existing plant ledger conservatively assigns all weak-decay ground-state Q and a known-isomer margin to neutrinos, so it does not rely on delayed weak-decay heat for electrical gain. Real electromagnetic decay heat still needs a warm-store heat load calculation.

If constant primary operation starts with no V49, the inventory and return-rate fractions follow 1−exp(−t/τ). They reach 53.57% after one year, 78.44% after two, 89.99% after three and 99.00% after six. Initial Ti49 feed to bridge the V49 return is about 0.21845 kg/s; cumulative bridge use approaches **8.9859×10^6 kg of Ti49**. This is a temporary catalyst inventory, compared with only 11.08 kg/year of steady boundary makeup in the saved reaction network. The bridge becomes the radioactive inventory already allocated in the mass ledger, so its atoms must not be counted twice. Carrier chemicals, holding tanks and handling apparatus remain equipment costs.

A station can therefore run before the V49 return reaches equilibrium by loading a stable Ti49 bridge, or can preload activated V49 from earlier isotope plants. Demonstrating near-steady catalyst return requires an appropriate operating history: about three years gives 90% V49 return, while six gives 99%. Fe55 reaches 90% only after about 9.16 years, but its side stream is below 1% of primary flow and can carry its own much smaller bridge. These are fixed-flow inventory fill laws. A coupled transient network would determine the exact mixture of initial bridge isotopes and hot-processor feeds.

The cold neutron-capture bank is separately **1.91490×10^11 kg** in the homogeneous-flux lower-bound screen, dominated by stable Ti46, Ti48 and Ti47. Its long per-atom capture times do not impose a million-year commissioning wait: preloaded stable targets supply the required capture rate immediately. The inventory is 0.0143% of whole-plant mass, but is isotope-specific. Total absorption is still used as a useful-capture upper bound for most entries; spectral self-shielding, undesired branches, extraction time and neutron competition can increase its required stock or processing power. Tiny uncovered Ca/K capture flows remain in the underlying network and should not be erased by a claim of exact closure.

## Economic meaning and history boundary

The conservative row consumes about 363 tonnes of ordinary hydrogen per year and exports approximately 0.27 terawatts. At the assigned 100 MJ/kg construction energy, its 1.340×10^23 J capital takes about **15,590 output-years** to repay. Annual replacement at the assigned 10^-4 mass fraction is 1.340×10^11 kg. Against the established 0.9-Lsun industrial supply, the construction energy is only 0.389 milliseconds of total captured power; that comparison permits funding but does not supply the remote manufacturing, transport and qualification schedule. Charging the 1.261×10^21 J field at 95% efficiency and 1 PW takes 15.36 days.

Thus economic restraint can coexist with successful physics: a specialist research or autonomy installation is consistent with the canon, while broad power-market adoption would need a different capital or service case. The selected future history requires a revised broad-capture and prepared-acceleration clearing/extraction system, with total extra electrical processing near 1.75 MeV/primary, apparatus inside the integration reserve, qualified materials and neutron return, a coolant topology matching the stated external heat environment, and sufficient catalyst startup stock. The monochromatic CARA bandwidth failure remains a specific discarded implementation; it neither demonstrates the replacement system nor establishes a fundamental prohibition on one.
