# ✨ Local guide shielding and 65 K heat

The cold guides intercept a small part of the thin receiver's radiation leakage. Their **pulse** heat capacity is more demanding than their original average interception suggests. Local hot-side shielding around the guide runs, together with an explicit cold buffer or sufficiently fast refrigeration, can address this without putting a heavy shield over the entire station.

This bounded screen reads the actual 500/750/1000 kg/m² parent rows with 0.005 additional missed-hot fraction and all missed-hot energy charged to cold equipment. It does not alter those power results. `guide_gamma.py` and `guide-gamma.json` retain geometry, source hashes, material inputs and angular/shield sensitivities.

## Intercepted average heat

The existing guide topology has six bores per 1 km meridional pitch: 0.006 m of bore per receiver m². At the selected winding specification the outer cryostat radius is 0.432 m. The summed radial projected coverage is therefore `2r×0.006=0.005186`, or **0.519%**. For isotropic outward intensity above a local plane, a tangent cylinder intercepts an effective width πr, giving **0.815%**. This angular result follows by integrating the cylinder's projected area; it is not a bound on an arbitrarily concentrated or grazing field. The saved factors 1, π/2, 3 and 10 expose that sensitivity.

Assigning all intercepted radiation to 65 K gives the following average-only result for the original geometry. Actual 200 K shields may absorb some first.

| Main receiver | Missed hot heat entering cold system | Refrigerator electricity, radial field | Refrigerator electricity, outward diffuse field |
|---|---:|---:|---:|
| 500 kg/m² | 12.126 W/m² | 0.871 W/m² | 1.368 W/m² |
| 750 kg/m² | 8.395 W/m² | 0.603 W/m² | 0.947 W/m² |
| 1000 kg/m² | 5.830 W/m² | 0.419 W/m² | 0.658 W/m² |

These use the existing cooler COP of 0.07222 (15% of Carnot between 65 and 200 K). For an allowed extra refrigerator input P, the effective intercepted-and-transmitted fraction must satisfy `f_eff ≤ COP×P/q_missed`. At P=1 W/m² the limits are 0.596%, 0.860% and 1.239%. The entire 12 W/m² should not be treated as a 65 K load merely because some guides are cold.

## Pulse buffering and a concrete local shield

Average heat alone does not establish quench margin. The B2 pulse carries about 16.68 MJ/m² of eventual non-neutrino decay heat before shielding. The script conservatively converts **all** that delayed heat, including beta kinetic energy, into a ≤20 MeV photon reserve with a 1.10 energy margin; it treats the saved prompt energy similarly. Main-wall attenuation uses the minimum carbon mass energy-absorption coefficient through 20 MeV, 0.001198 m²/kg. The additional 0.005 missed-hot fraction is separately retained as a photon-like acceptance scenario. Hard-particle leakage or an unbounded photon spectrum requires a separate transfer calculation.

A local continuous closed-ring graphite enclosure is placed on the hot side of the cryogenic assembly. For the illustrated 2000 kg/m² radial column it is 1 m thick at density 2000 kg/m³. Its conditional transmitted-energy factor is exp(−0.001198×2000)=0.0911. This inherits the photon packet's assumptions of local charged-secondary deposition, no upscatter and sources outside the enclosure. Joints, openings, secondary ranges and actual spectrum still need transport verification. This factor is applied to the **outer enclosure's intercepted area**, allowing photons scattered toward the guide from a wider collection region; the calculation does not retain the smaller bare-coil cross-section after adding shielding.

The cold buffer is sized to absorb the whole worst pulse with no simultaneous refrigeration credit, for a chosen 65–70 K temperature interval. NIST copper heat capacity gives about 155 J/kg/K at 65 K and 819 J/kg from 65 to 70 K; after a 5% fit margin the calculation uses **778 J/kg**. Copper is a transparent calorimetric reference, not a selected universal buffer recipe. Its inventory cost and conductivity/contact performance remain part of a real design. [NIST copper properties](https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm), [NIST carbon coefficients](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html).

The buffer is modeled as an annulus surrounding the winding, and **its increase of cryostat radius, shield mass and radiation interception is solved self-consistently**. Existing stabilizer heat capacity is not credited as free storage. For outward diffuse radiation:

| Main receiver | Local hot graphite shield, averaged over receiver | 65 K copper-equivalent buffer | Total local addition | Additional refrigerator electricity including enlarged 200→65 K radiative area |
|---|---:|---:|---:|---:|
| 500 kg/m² | 87.3 kg/m² | 37.7 kg/m² | 125.0 kg/m² | 0.583 W/m² |
| 750 kg/m² | 83.3 kg/m² | 27.3 kg/m² | 110.5 kg/m² | 0.401 W/m² |
| 1000 kg/m² | 80.2 kg/m² | 19.9 kg/m² | 100.1 kg/m² | 0.278 W/m² |

The shielded outer radius is approximately 1.56–1.66 m; neighbouring guides are kilometres apart, so this is local shielding rather than a continuous metre-thick outer shell. At angular factor 3 the 500 kg/m² row requires about 181 kg/m² and 1.20 W/m² instead. In contrast, the unshielded diffuse-case full-pulse buffer enlarges to about 416 kg/m² at the 500 kg/m² point, with nearly 5.85 W/m² additional refrigeration. That is why the original average-only estimate should not be promoted into a completed cold-guide design.

The 1000 kg/m² services allowance had roughly 829 kg/m² unallocated in the named-sector screen. The illustrated 100–125 kg/m² local shield/buffer can be assigned from that remainder while retaining the original guide package and frame. It is finite apparatus and material stock; it is not zero additional hardware. The extra refrigerator input must likewise be assigned within the retained ordinary auxiliary allowance or debited once. At a radiator-saturated point, any unallocated extra electrical dissipation changes the allowed rate.

Hot local shielding also needs insulation from the 200 K radiation shield. An effective emissivity of 10⁻⁴ across its inner surface would add approximately 0.12–0.14 W/m² station-wide hot-to-cold leakage in these examples. That is an explicit insulation requirement to budget within the retained 1 W/m² bypass, or add if that allowance is already spent. Returning some intercepted radiation to the hot loop could improve the parent's conservative all-cold energy account, but no such credit is taken here.

The 5 K buffer interval is an acceptance scenario, not a measured conductor quench margin. Spatial energy deposition, heat diffusion into the buffer, conductor critical-current margin, radiation damage to insulation and the actual angular leakage pattern remain to be checked. The useful conclusion is a **local shielded-guide architecture with finite cold storage**, rather than either charging all receiver leakage to 65 K or assuming mean cooling alone protects the magnet.
