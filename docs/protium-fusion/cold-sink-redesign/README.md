# Shaded receiver and improved energy conversion ✨

The working reactor design now uses **60% heat-to-electric conversion**, a **1,000 K hot receiver** and a **200 K radiator**, with **zero direct solar heating** under Thomas's Dyson-shading premise. Infrared heating remains an explicit environmental input. The former 40% conversion and directly illuminated radiator are retained as historical comparisons.

At the existing 1.1×10⁹ m radius and pulse rate, the revised conditional account exports **6.574×10²⁰ W from 1.219×10²³ kg, or about 5.4 mW/kg**. It uses an illustrative **11.0254 W/m² absorbed infrared background**, **1 W/m² heat bypass around the engine**, and the separate **2 MJ/kg routine refurbishment sensitivity**. Radiators occupy **79.9% of the outer sphere**, within a specified 90% available area. The pulse preparation, source storage, 0.007 fresh C:H budget, source count reference and required material return remain as before.

This is a concrete operating requirement supported by energy and mass calculations. The component ratings, insulation and collection performance still need their detailed designs. The [station account](ACCOUNTING.md) and [36 saved cases](cases.json) retain those assumptions explicitly; no new fusion yield or pulse cadence is credited.

## Why 60% is a reasonable working choice

The external Carnot bound rises from 55% at a 1,000/450 K temperature pair to 80% at 1,000/200 K. The [heat-engine study](heat-engine.md) goes beyond that bound: a closed regenerative gas cycle with six intercooled compression groups, six reheated expansion groups, 93% isentropic stage efficiencies, a 95% recuperator and finite pressure and mechanical losses produces **61.4% package electrical output**. Stronger named components produce **66.9%**. These are calculated future machinery specifications.

The working station selects 60% package output, including the converter's internal loss allowance. It also allocates dedicated exchanger mass, including recuperation, while retaining the existing generic services. At the selected 10 W/(K kg) exchanger rating, the requirement is **8.86×10¹⁸ kg**. The exchanger and generator account therefore accompanies the efficiency assumption.

The complete plant exports **40.55% of total nuclear energy**, or **43.87% of its non-neutrino portion**. Driver electricity, ordinary auxiliaries, support and refurbishment explain the difference from the converter's 60%. Prepared pulse energy returns through the receiver and is charged once through the driver ledger.

| Conditional comparison | Net output | Radiator coverage |
|---|---:|---:|
| 40% conversion, 450 K radiator, revised IR and refurbishment | 3.03 mW/kg | 4.0% |
| **60%, 200 K, 1 W/m² bypass** | **5.39 mW/kg** | **79.9%** |
| 60%, 200 K, 10 W/m² bypass | 4.79 mW/kg | 85.9% |
| 65%, 200 K, stronger engine components | 5.99 mW/kg | 74.1% |
| 70%, 200 K, higher-performance requirement | 6.58 mW/kg | 68.3% |

The 65% row has a corresponding constructive cycle with stronger components. The 70% row has finite temperature approaches and an explicit exchanger requirement, but a complete matching cycle has not been solved. These rows form a sensitivity comparison, rather than a probability interval or a global optimum. At the working 60% point, varying routine refurbishment from 2 to 100 MJ/kg changes output from 5.395 to 5.084 mW/kg without shrinking factory capacity.

## Shade, infrared and insulation

The [environment calculation](environment.md) follows the canonical orbiting swarm, whose detailed radiator sky at the reactor is unspecified. An isotropic compact source reradiating 90% of solar luminosity inside 5 AU gives the 11.0254 W/m² absorbed benchmark for a grey spherical receiver. Directional panels, extended emitters and local structures change this value. A warm enclosing bath is a different geometry.

The 200 K radiator emits about **81.65 W/m² per unobstructed emitting face** at emissivity 0.9. The station charges absorbed environmental power over the entire outer sphere, then assigns enough emitting area to reject it along with plant waste heat. With its other working assumptions, 90% radiator coverage permits about **19.25 W/m²** of full-sphere environmental absorption. Self-view and warm obstructions reduce that allowance.

A 1,000 K hot receiver needs a deliberate thermal barrier behind it. The heat-engine study gives ideal multilayer-shield and support-bridge requirements for 1–10 W/m² bypass. Thousands of thin reflective shields can fit within the existing 100 kg/m² insulation allocation in an ideal mass calculation using abundant materials; their optical performance, contacts, penetrations and lifetime determine actual leakage. Bypassed heat loses its engine conversion opportunity and enters the cold load once.

The cold exterior does not determine the temperature of the cavity-facing receiver or protect incoming fuel from that receiver's radiation. Those remain separate parts of the established fuel-handling design.

## Direct recovery and debris state

The frozen whole-cycle ledger puts **73.68% of outgoing energy in bulk kinetic motion**, 4.03% in endpoint internal energy, 1.28% in prompt radiation and 21.01% in delayed decay heat. The bulk kinetic term averages **79.35 W/m²**. Its largest saved pulse reaches **636 kW/m²**, so direct-converter equipment must handle a very different peak from a thermally buffered generator.

The [direct-recovery study](direct-recovery.md) assesses a time-ramped multistage collector ahead of material collection. Recovered electricity removes that same energy from the thermal engine's input. Converter mass, electron handling, charge-state acceptance, output storage and round-trip losses determine the net advantage. At 80% direct recovery, a sufficient output store for ideal independently periodic stage streams holds about **1.01×10²⁷ J**, separately from the existing driver store. The energy-only improvement cannot establish a completed direct-converter operating point.

The [charge-state assessment](charge-state.md) separates the initial hot ionized plasma from its later recombination history. At an assumed 50% threshold-delivery efficiency, giving every arriving atom one positive charge consumes **0.000911 W/m²**, while complete stripping consumes **0.05506 W/m²**, averaged over the cycle. The latter is about **0.0694% of outgoing kinetic power**. Absorption, electron handling, charge-state control and apparatus mass remain separate requirements.

The [debris-temperature calculation](debris-temperature.md) distinguishes local random motion from variation in the coherent expansion. For C13, 98% of kinetic energy arrives between about **59 and 269 seconds** after launch, with a pulse root-mean-square bulk speed near **9,600 km/s**. Its calculated oxygen-15 daughter recoil gives a pulse-averaged ion-temperature equivalent of about **2,800 K** if that random energy remains with the nuclei. Sharing it with the fully stripped electron population instead gives a common temperature near **400 K**. These are conditional velocity-moment calculations: a cold core and a minority of fast daughters can carry the same random energy without forming a common thermal distribution.

For the ion-retained, weakly heated branch, **3,000 K with a subjective 80% interval of 1,000–10,000 K** is a conditional judgment about the characteristic pulse average. It is not a spatial temperature range or an unconditional arrival prediction. The calculation keeps ion/electron sharing and positron heating as separate branches. Immediate local deposition of 1% of the emitted positron kinetic energy, followed by ideal adiabatic cooling and full particle equilibration, gives about **32,000 K**; the corresponding 100% deposition branch gives about **3.2 million K**. Even that branch retains approximately **99.95% of particle kinetic energy in bulk motion** under its expansion assumptions. Later lepton thermalization, charge evolution, radiation and collective fields require a coupled calculation before selecting an unconditional electron or ion temperature.

The **73.68%** whole-cycle energy share includes prompt radiation and delayed decay energy in its denominator. The much larger bulk fraction within the arriving particles answers a different, local kinetic-energy question.

The 60% thermal-engine premise is compatible with either charged or neutral arriving material. The selected fuel-return account still requires low-loss slowing and collection; a neutral stream may require charge control for that function even when direct electrical recovery is omitted.

## Economic implications

The revised working point requires about **185 kg of station per net watt**, compared with 362 kg/W in the original 2.77 mW/kg account. For a specified construction recipe and shipment cost, this nearly doubles output per installed kilogram and approximately halves energy-amortization time. It does not turn the 2 MJ/kg refurbishment sensitivity into an initial construction cost: fresh-feed chemistry, material acquisition and transport keep their own accounts.

With the old comparison assumptions of 30% solar conversion and 1 kg/m² complete collector mass, unobstructed solar power matches the revised specific output near **275 AU**. That is a reference for an exposed collector. At the shaded installation, the relevant Solar-system alternative includes power delivered from the existing swarm. The improved reactor remains most economically interesting for fuel-cycle development, remote or independent power, and long-lived hydrogen use where external energy delivery is costly.

## Evidence and reproduction

The old source packets remain reproducible historical snapshots. The working update is contained here. Run `heat_engine.py`, `environment.py` and `cases.py` in that order; the environment script needs NumPy. The other linked scripts reproduce the direct-recovery and debris-state screens from their named archived inputs. Code and input hashes are recorded with each result.

The station family checks installed mass, electrical power, nuclear energy including neutrinos, separate cold-heat components and radiator/environment exchange. It reproduces the old numerical controls and rejects zero-approach heat-engine limits. The physical acceptance criteria are positive whole-plant output above 1 mW/kg, a realizable emitting-area requirement, finite temperature approaches and the retained fuel/material budgets. Passing these accounting checks does not replace a component or plasma transport calculation.
