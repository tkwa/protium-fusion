# ✨ Independent plant-account review

The current catcher-interface ledger passes independent review. It contains **47 scenarios: 42 complete plants and five correctly rejected auxiliary-overload cases**. The preceding revision passed 1,267 reconstructed identities. After adding the full inactive cryobank cycle and the passive cold-pressure interface, a further **646 independent delta checks** pass, with maximum normalized residual **2.67×10⁻¹⁵**. All **14 current source/code hash checks** pass. The unchanged baseline reproduces exactly.

The original audit summary remains in `output/catcher-followup/plant-account-independent-audit.json`; the current extension is in `output/catcher-followup/plant-account-independent-delta-audit.json` in the research scratch checkout. Reviewed fingerprints below refer to the current version.

## Final cryobank and pressure-interface extension

The inactive bank waits at the actual radiator temperature, then cools to 5 K over the specified 100 seconds. Its account includes the complete sensible-cooling work, radiation absorbed during cooldown, and electrical warming. The analytic exponential-cooldown integral independently agrees with the numerical quadrature. The integration adds only `new_full_inactive_cycle − old_85K_sensible_cycle` to the prior entry; the old sensible term is therefore replaced without duplication. The active-bank radiation charge remains separate.

At the warm reference's 200.18862 K sink, the increment is **0.00909747 W/m²**. It consumes the existing auxiliary allowance, leaving **0.99287024 W/m²** unassigned. Reference net power, total mass and radiator temperature remain unchanged; demand-limited sensitivity rows retain their explicit additional debit.

The passive cold-panel relation also closes: `Pc = Peq + (Tc/Th)*Qhot/Spanel`, with combined hot conductance `Sc/[1 + Sc*sqrt(Tc/Th)/Spanel]`. The retained panel areas support a mean summed cold-interface pressure **0.000593865 Pa** in the zero-equilibrium-pressure case. The revised buffer pressure incorporates that finite sink impedance. Physical conduit area and its radiation/conduction heat are unchanged, so adding the pressure correction's provenance does not warrant another photon or gas-enthalpy debit.

## Energy routing and auxiliary allocation

For captured pulse energy C, prepared source power U, hot bypass B and recovered auxiliary heat R, the engine input is `hot = hot_eligible_fraction*C − B + R`. The electrical account still subtracts the complete ordinary-auxiliary input. Its cold dissipation is `aux − R`, so recovered heat appears once. The saved cases satisfy `0 ≤ R ≤ aux`.

The effective auxiliary debit is `max(0.1*process_margin, already_named + new_named)`. A named pump or repair load that fits inside the retained 10% allowance therefore consumes that allowance rather than being charged a second time. Heat recovery increases engine conversion and may also increase the 10%-of-margin allowance; that is the explicit retained budgeting convention. Above the allowance, actual named demand controls the debit.

The 1 W/m² global and 0.13 W/m² primary-booster entries are **total motor input targets**, including molecular windage. The windage packet is included in provenance, and its mechanical load is not added again. The 0.52 W/m² booster sensitivity is confined to the warm-buffer architecture. Repair sensitivities change electrical work per serviced event while keeping the specified serviced throughput and factory requirement separate.

The actual buffer conduit contributes **0.0141133322 W/m²** of hot-to-cold radiation plus conduction. The reference adds 0.5 W/m² bearing/support bypass and **3.2043 µW/m²** of H₂/He sensible enthalpy, giving **0.5141165365 W/m²** additional bypass. The model removes this heat from engine supply and deposits it in the cold account once. The stipulated fraction reaching the 5 K bank changes refrigerator work; it does not create a second photon-heat source.

Radiator temperature is solved with the selected finite heat-engine approaches. Cryogenic work, inherited temperature-dependent loads, radiator area, exchanger and converter masses, support, spares, replacement capacity and installed reserve all retain their feedback. The component heat sum and global non-neutrino nuclear balance agree independently.

## Mass allocation

The **14.2879806 kg/m²** rework-factory requirement is explicitly included within the 300 kg/m² new-service allocation. That allocation fits inside the existing 678.606 kg/m² unassigned service allowance, leaving **378.606 kg/m²**. It is not added again as an independent plant mass. Larger service scenarios add only the amount exceeding the available allowance, after which the station recalculates mass-dependent loads and reserve.

The warm-buffer reference gives **20.65764 mW/kg** at **200.18862 K** radiator temperature. The direct-cold-port comparison gives **19.27292 mW/kg** at **202.22611 K**. These are conditional accounting results for the unchanged pulse rate and source family. The review establishes bookkeeping closure; component power targets and the service allocation retain their stated engineering scope.

## Reviewed fingerprints

- `calculate.py`: `25db98145b97a17811dfe4badcac951bf7371b71fbd927dfff28351b04918e2e`
- `station_catcher.py`: `6175e57a538622b34959acb40af1a346677fb891e1745c9117b8deb179a1008a`
- `results.json`: `1e4f725953035424880f50c104e1ba11668df56e682d15daa49f9abe7e624197`
