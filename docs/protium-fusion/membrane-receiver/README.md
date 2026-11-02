# Membrane and receiver assessment ✨

This checkpoint assesses a barrier on the cavity-facing side of the reactor absorber: pulse center → supported thin film → energy receiver and collection space. Its purpose is to keep slow receiver material out of subsequent incoming fuel. The reference retains a 1.1×10⁹ m receiver radius, 0.007 fresh C:H mass budget and the provisional [1,000 launcher complexes](../launcher-reference/README.md). Geometry, energy accounting and lifetime have distinct requirements.

## Outgoing pulse and stopping

The corrected C13 pulse distributes 7.438×10⁻⁷ kg/m² over the receiver. It carries about 3.584×10¹⁹ nuclei/m², predominantly CNO, with roughly 34.1 MJ/m² kinetic fluence and 48.85 kJ/m² prompt photon fluence. The smooth ballistic reconstruction places most kinetic energy between about 59 and 269 seconds after focus. The exact staggered-node kinetic accounting and a separate mass-conserving continuous reconstruction are retained in [pulse-loading.json](pulse-loading.json). They differ by about 0.14% in kinetic energy; this is recorded rather than hidden by rescaling mass.

The [species-dependent stopping calculation](stopping.py) uses measured electronic-stopping data with identified interpolation/extrapolation, an elastic-transfer model, and the finite-flight isotope mixture. Its intact-silica results for C13 are:

| Film thickness | Incident kinetic energy removed from projectiles | Peak projectile energy-loss flux |
|---|---:|---:|
| 10 nm | 0.20% | 0.742 kW/m² |
| 100 nm | 2.00% | 7.43 kW/m² |
| 1 μm | 19.2% | 75.0 kW/m² |

Projectile energy loss is an upper account for local heat before secondary-particle escape is modeled. The thin film's radiation damage, sputtering, reaction products and gas permeability require separate transport or material data. At 1 μm, about 2.9% of incident mass stops in this silica calculation, requiring collection from or processing of the film as part of fuel return. The stopping screen alone establishes neither a usable lifetime nor the required integrated CNO/H return.

## Thermal vapor and the reason for a barrier

The [thermal plume calculation](plume.py) combines the archived receiver histories with NIST-JANAF graphite-vapor data for C, C₂ and C₃. It keeps the prompt and delayed particle heating at their actual times. For the prescribed uniform wall and fixed material properties, it gives **0.00572 kg/s** average gross thermal evaporation across the whole sphere, versus 3548.6 kg/s in the entire structural-carbon makeup allowance. The C13 prompt surface reaches about 1713 K; the largest slow-heating maximum is about 1387 K. Adding maxima that occur at different times produces a much larger apparent temperature and evaporation rate.

This result addresses thermal evaporation from the modeled graphite wall. Irradiation sputtering, dust, localized heating and malfunctioning collection can inject additional material. A useful receiver design must quantify those sources before choosing a whole-surface barrier for contamination control. The [gas assessment](GAS-ASSESSMENT.md) adds accumulation, flight times, interaction columns and packet-size criteria. It replaces the original factor-three instantaneous column screen with a checked 3.26 envelope where that model applies. Conditional on the modeled thermal source and the explicitly assumed source, capture and interaction uncertainties, its central encounter optical depth is 6.33×10⁻⁸, with a subjective 80% credible interval of 1.14×10⁻⁹–2.33×10⁻⁵, against a selected 10⁻³ screening limit. Dust, sputtering and localized plumes require their own source estimates.

The [thermal and mechanical assessment](THERMAL-MECHANICAL.md) evaluates the film itself. Assigning a Planck spectrum to each saved grey escape increment gives a 1 μm silica film a calculated prompt sensible-temperature rise of roughly 23–107 K across the four stages at full local retention. Assigning all unresolved sub-keV energy to the film adds about 2–20 K. These are spectrum-proxy results with fixed material properties; the unresolved soft component gives a much larger bound for a 10 nm film. The earlier all-1-keV photon case is retained as a separate stress test.

Radiation emitted inward returns through the closed cavity; cooling needs an actual outward heat path. At the C13 75 kW/m² projectile-loss peak, a backed 1 μm silica layer has a calculated through-thickness temperature drop of about 0.058 K, and a contact conductance of 750 W/m²/K would hold the contact rise to 100 K under full heat retention. This points to a cooled catcher surface after electromagnetic braking. A film placed before a downstream braking gap instead requires lateral support and cooling: a solid backing would intercept the transmitted ions. These two arrangements have different energy-conversion and collection requirements.

Supported cassettes can close pressure forces locally through ties between opposing loaded surfaces. Tile edges alone do not remove the global pressure load. Photon spectrum, secondary-electron heat retention, support geometry and material lifetime remain inputs to the selected cooling arrangement.

## Renewal after the manufacturing correction

[renewal_energy.py](renewal_energy.py) applies explicit 1, 2, 3.1, 10 and 100 MJ/kg recycling sensitivities. The 100 MJ/kg value is retained as the former accounting case. The accepted manufacturing premise is a minimum 50% conversion efficiency for specified positive free-energy upgrades, often higher; same-material renewal instead needs its actual processing losses. The old blanket charge does not select a renewal cost.

A 1 μm silica film contains 3.345×10¹⁶ kg. Renewing it after every modeled pulse processes 5.021×10¹² kg/s on average. At 2 MJ/kg this consumes 1.004×10¹⁹ W. With all original costs retained, the energy-only station output falls from 2.766 to 2.683 mW/kg. If routine station refurbishment also takes the separate 2 MJ/kg sensitivity, the corresponding result is 2.994 mW/kg. Extra cooling, support, heat rejection and larger processing equipment still need their allocations.

The retained factory would need 14.04 times its original total material throughput for routine refurbishment plus this film renewal. Its original flow would take about 86,871 seconds to process one sheet, compared with a B1 interval near 6,774 seconds. Higher throughput, local renewal and material buffers need an actual schedule. The model makes those requirements explicit rather than shrinking factory mass when electricity consumption falls.

Gross reprocessing remains inside the station. For silica, even an unreturned fraction of 10⁻⁸ at this renewal rate requires about 50,212 kg/s of fresh silicon and oxygen combined. Those elements require their own supply account. For a carbon film, irreversible carbon escape competes directly with the 0.0014 structural share of the C:H budget. Material removal toward the center would also undermine the barrier's contamination-control function.

## Integration with 1,000 launchers

The count changes the layout of source complexes and their protected openings. Chamber-facing sources can receive slow feed through shielded transfer systems and emit high-speed packets on the cavity side. Behind-receiver sources require fast-packet apertures through the receiver and any film. The [port assessment](../launcher-reference/ports-and-shielding.md) retains both arrangements and calculates aperture losses and closure timing. Repeated packets from one compact site retain nearly the same approach direction; angular assembly remains a separate requirement.

A continuous membrane has not yet earned a smaller radius or faster pulse cadence in the selected account. The calculations remove the former energy-only rejection of renewal, and they substantially reduce the thermal-vapor source compared with a sustained-maximum-temperature estimate. The present design comparison should retain locally shielded source heads and a backed receiver as its reference, with supported membrane cells as an alternative contingent on material injection, heat removal and lifetime results. The power benefit of that alternative must come from a demonstrated change in collection, receiver size or cadence that exceeds its added loads.

## Reproduction and evidence status

`renewal_energy.py` uses only the frozen station JSON and the Python standard library. Its 120 cost rows were independently reconstructed. `station_extension.py` retains the earlier coupled penalty model with the old routine 100 MJ/kg maintenance charge; its default membrane-renewal rows are historical sensitivities. Use the newer renewal-energy table for the corrected cost comparison, respecting its energy-only scope.

The archived `pulse-loading.json`, `stopping-data.json`, receiver input snapshots and their source hashes support the downstream calculations. Regenerating `pulse_loading.py` itself also requires the original full trajectory files named in its provenance. The raw pulse model, a grey photon spectrum proxy, material stopping models and subjective engineering judgments are identified separately in the linked assessments. No new multidimensional implosion or membrane irradiation experiment is represented by these calculations.
