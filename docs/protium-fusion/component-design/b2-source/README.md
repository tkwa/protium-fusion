# ✨ B2 source shaping with the updated carbon rates

The first B2 source-map test exposed a consequential distinction between a mean entropy budget and its distribution. The existing N=2000 parallel-sector CH2 mapping completed its full simulation but produced very little fusion. A constant contact-shock entropy, obtained by controlling the normal collision speed, restored productive burning at the same total prepared energy. Both results are preserved.

All cases here use actual C12:H2 abundances, the 64-zone compression-cap solver, FastEOS during burning, grey opacity 0.001 m²/kg, and the same median LUNA carbon-rate revision. Physical duration is 63.0627249 microseconds. The existing CH2 mapping's initial radii, temperatures, velocities and cell masses are bit-identical to the H2 reference profile. Its condensed reference is 944.2 kg/m³ and 300 K, with FEOS's chosen 2 GPa reference bulk modulus.

| Case | T/reference mass 10/50/90% | Deposited nuclear / prepared energy | O14 per initial C12 |
|---|---|---:|---:|
| Updated-rate reference | 1 / 1 / 1 | 10.271886 | 0.799328 |
| Local N2000 mapping | 0.465 / 0.899 / 1.537 | 0.176361 | 0.00138883 |
| Constant-contact-entropy mapping | 0.99705 / 0.99801 / 1.00352 | 10.275337 | 0.799649 |

The first mapping's mass-mean temperature ratio is about 0.964, concealing a large variation between parcels. Its coldest nonseed parcel is at 0.0477 of reference, and the hottest relative parcel at 2.223. A small mean difference therefore cannot certify ignition. The low-gain run completes without a guard, with peak temperature 1.3983 GK and energy residual −1.96e−9 of input. It is a completed unsuccessful **chosen input profile**, not evidence that an adaptive source cannot work.

The saved temperature histories make the distinction clearer: the N2000 case brings at most **0.0824% of mass above 0.5 GK**, equal to the retained seed mass, despite its 1.40 GK local peak. Its maximum mass-mean temperature is 0.317 GK. The constructive case brings **95.39% above 0.5 GK** and **86.73% above 1 GK**, with peak mass-mean temperature 1.045 GK. `temperature_coverage.py` / `temperature-coverage.json` retain the 400-sample diagnostics. They describe spatial coverage and do not isolate a unique causal parcel or replace reaction-rate integration.

## A constructive collision-energy choice

The reference FEOS entropy is nearly uniform: approximately 38.23–38.56 kJ/kg/K. The replacement uses one shock entropy, 38.52495 kJ/kg/K, chosen from the active reference mass. Each parcel is then compressed isentropically from that contact state to its existing prescribed hot density. The resulting temperature is calculated; no simulated temperature is clipped, and no additional external heater is inserted.

The contact shock deposits **3.45103 GJ/kg**, only **0.095715% of the whole prepared energy**. Equal streams need normal/transverse speeds **±83.079 km/s**. The FEOS state is **3761.2 kg/m³, 0.439707 MK and 8.7013e12 Pa**. The global radial-speed correction restoring the original prepared energy is only **1.00005438**. The inner approximately 0.1% seed retains its separate reference temperature.

An illustrative symmetric-pair geometry with local kinetic allocation equal to the eventual hot internal plus radial energy gives mass 10/50/90% half-angles **1.38°/1.86°/2.63°**. These are collision geometry scales. They do not specify a minimum number of physical outlets or random pointing tolerances. Relative particle directions near contact, first-contact locations, pressure work and the actual cold-to-hot trajectory remain to be solved in three dimensions.

The constructive run completes at peak temperature **1.6176 GK** with energy residual **1.21e−6** of prepared input. Its gain is within 0.034% of the updated reference. This is strong evidence that the entropy *program* can be changed constructively within the 1D model; it does not qualify the program's physical delivery.

## Primary EOS comparison

`contact_primary_check.py` solves the primary CH2 FPEOS Hugoniot at the same deposited shock energy. It gives density **3771.0 kg/m³**, temperature **0.486349 MK** and pressure **8.6938e12 Pa**. Relative to this interpolation, FEOS temperature is 9.59% lower and compression 0.259% lower. The close pressure/compression agreement supports the energy/geometry scale, while the temperature difference and different entropy models remain a preparation uncertainty. The FPEOS check does not independently integrate the full contact-to-hot isentrope.

## Files and evidence boundary

- `run_mapped.py`: applies the existing N2000 map with exact-profile provenance and the median carbon kernel.
- `B2-H2-N2000-carbon-median.{json,initial.json,inputs.json}`: preserved completed low-gain control; raw trace remains under `/home/tkwa/src/protium-fusion/output/b2-source/`.
- `run_constant_collision.py`: constructs the alternative common shock entropy and records all energy/angle assumptions before the full-network run.
- `contact_primary_check.py` / `contact-primary-check.json`: independent primary CH2 Hugoniot comparison.
- The constructive median-rate result and exact initial profile are archived in the separate four-stage `source-family-carbon-median` packet, whose `stage-map.json` identifies all records and raw traces. Earlier source maps are unchanged.

The next smallest physical discriminator is a finite multi-stream contact patch with the selected collision-energy distribution. It must measure the delivered entropy, density, transverse motion and timing together; replacing a whole source by its mean entropy would miss the failure seen here. Receiver hardware and whole-plant consequences must use the constructive result's actual endpoint isotopes and spectra, not a scalar gain correction. The lower gain of one geometry creates no universal angular-count or near-perfect-yield requirement.
