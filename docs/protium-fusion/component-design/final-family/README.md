# Delivered B1 and polymer C13 family integration ✨

This packet composes actual 128-zone source and burn records. The two explicit families differ only in C13 feed: H/C13 = 4 for the comparison, and H/C13 = 2 for the polymer-compatible candidate. Both use the delivered B1 profile after the κ = 0.001 m²/kg preparation calculation. A and B2 retain the frozen source128 profiles.

| Family | Prepared at hot burn, MeV/He | Prepared at precursor, MeV/He | Preparation radiation, MeV/He | Fresh C:H mass ratio |
|---|---:|---:|---:|---:|
| Delivered B1, C13H4 | 1.927311030 | 1.927325295 | 0.0000142655 | 0.0003129094 |
| Delivered B1, C13H2 | 1.936702645 | 1.936716910 | 0.0000142655 | 0.0003129181 |

These material balances assume eventual return of the modeled reusable hydrogen and CNO streams. Temporary retention and recycled extraction rejects do not create a fresh-material debit. Oxygen held outside the selected catalyst-return states still contributes to the fresh-carbon balance; this is a selected process boundary, not evidence that oxygen escaped the station. The available 0.007 C:H budget need not be consumed.

The delivered B1 precursor contains 3.196326437256525×10²³ J, of which 7.303449010114459×10¹⁸ J leaves as preparation radiation. Its hot burn starts with 3.196253402764988×10²³ J. The combined ledger charges the precursor once and preserves both radiation stages. The original FEOS/FastEOS reference difference remains the separate model-handoff qualification in [assembly-burn](../assembly-burn/README.md); it is not silently assigned heater work or recovered energy.

B1's incoming source rating is recomputed from its actual common FastEOS precursor profile. At the 10⁹ m reference radius, the source peak is 2.59427×10²¹ electric W under the existing 60% preparation convention, compared with 2.16023×10²¹ W inferred from the delivered hot state. The mass99.9% incoming interval is 3320.43 s rather than 3313.73 s. Outgoing velocity, composition, receiver heating and timing come from the completed delivered burn. Other stages retain their own actual source profiles. The main-window convention preserves the separately handled slow tail rather than treating it as matter lost from the station.

## Interface

Each family directory contains `integration.json`, `spectrum.json`, `hot-reduced-inputs.json`, `hot-spectrum-interface.json`, and `independent-checks.json`.

- `hot_cycle` is the unmodified actual-composition hot-burn material and energy ledger.
- `precursor_cycle` increases `U_prepared_MeV_per_He` by the actual preparation input increment. Nuclear, isotope, outgoing kinetic/internal, burn-radiation and delayed-decay fields remain unchanged. The energy-partition residual includes the separate preparation radiation.
- `preparation_escaped_MeV_per_He` and `preparation_energy_increment_MeV_per_He` are separate named channels. Their tiny difference is the recorded numerical handoff residual.
- `hot_selected` and `stage_inputs` retain physical hot-burn records and exact inputs.
- `plant_selected` is a documented allocation view. `source_prepared_energy_J` gives the actual precursor input where present. The legacy plant calculator's `record.initial_energy_J` field is overridden **only in this view for storage sizing**. It must not be mistaken for a changed hot-burn state or used to reconstruct the material cycle. The receiver's `source_rating` and incoming timing use the precursor; `hot_start_source_rating` and `hot_start_timing_at_radius` retain the comparison values.
- `preparations` contains the actual precursor profile, composition, source rating, incoming timing, saved preparation-radiation waveform and source hashes.
- `spectrum` is embedded, and `spectrum_path` identifies the adjacent standalone copy. It is normalized to **H + U_hot**, with preparation metadata included. Preparation radiation is excluded from its hot-source weights.

For a downstream hot transport loss fraction `L_hot`, useful hot-boundary supply before apparatus bypass is `(1 − L_hot) × (H + U_hot)`. The driver draws `U_precursor / eta_driver`. Preparation radiation is assigned separately to cold equipment or external escape. An additional apparatus miss fraction such as 0.005 must name its denominator; a total-source convention would subtract `0.005 × (H + U_precursor)` from useful hot supply once and add it to the selected cold/external channel. Changing source denominators must not create an additional thermal credit for prepared energy.

## Comparison plant points

The supplied comparison points retain the older 99% hot-boundary energy-capture convention, 60% heat conversion, 60% preparation efficiency, radius 1.1×10⁹ m and the existing radiator/dispatch constraints. Preparation radiation receives no recovery credit and is assigned externally for these comparison points. The effective composite capture is `0.99 × (H + U_hot)/(H + U_precursor)`; its definition is stored explicitly.

| Family | 5000 kg/m² receiver | 2000 kg/m² receiver |
|---|---:|---:|
| Delivered B1, C13H4 | 6.6070 mW/kg | 12.4472 mW/kg |
| Delivered B1, C13H2 | 6.5971 mW/kg | 12.4286 mW/kg |

These are integration controls. The parent coupled-photon/thermal calculation selects receiver thickness, photon capture and the preparation-radiation cold/external partition explicitly. The polymer case is about 0.15% lower in these control points. The C13H2 source mapping itself completes at Q = 22.65998; its physical material and capsule qualifications are in [c13-polymer-source](../c13-polymer-source/README.md).

## Verification and reproduction

Each family passes 35 integration checks and 28 independent pulse-energy/material/partition checks. The latter reconstruct prepared input from exact pulses per helium, verify precursor and hot partitions, require unchanged isotope and nuclear fields, and independently check driver, hot capture, external escape and plant mass. The hot spectra use the actual completed isotope outputs.

`integrate_family.py` accepts an explicit four-stage map. A map entry can contain a `preparation` record/initial/trace triplet; no stage is discovered implicitly. `verify_family.py` verifies a saved family without rerunning hydrodynamics. Full source traces remain at their hashed research paths; the reduced inputs, precursor waveform and physical endpoint records required by downstream material, spectrum and plant calculations are embedded here. Old source128, assembly and radiation-transport packets are preserved.

The delivered precursor is an auxiliary spherical whole-fluid construction. The other mapped source profiles, cold packet assembly, multidimensional merger, actual opacity and material interfaces retain their separate qualifications. This integration closes the declared numerical and material interfaces without establishing a complete source apparatus.
