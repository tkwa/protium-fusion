# Additional thermal cavity-source robustness ✨

The finite-cap diffusion result describes one ideal, already-mobile channel population. Total gas release also depends on prompt surface interactions, energetic ions that cross the cavity and reimplant, later thermal liberation, damage, seams, process exhaust and source leakage. This extension exposes those additions as independent thermal source rates. It leaves the frozen gas packet unchanged.

For the actual final350 source, a total common H/He thermal cavity fraction of **0.003**, with effective global capture **g=0.01**, passes the selected bulk 300 K comparison for the median and fast 2 cm C13 capsules. A fraction **0.01** exceeds that comparison for the fast capsules. Raising actual global capture or changing delivery/capsule design changes that boundary. These are conditional heat margins, rather than estimates of how much gas the catcher releases.

## Reusable source balance

`calculate.steady_column(R,T,q,m,g)` accepts arrays of species source mass rates, gas-particle masses and effective capture probabilities. For a uniform Maxwellian gas,

`Sigma = (4R/A) sum_i[q_i/(g_i vbar_i)]`, with `vbar_i=sqrt(8kT/(pi m_i))` and `A=4piR²`.

It returns columns, pressures, inventory and residence times for each species. Capture includes physical access and net transfer into a retained return line. Heating, ionization or a wall bounce alone does not remove the gas. The formula also permits a source exhaust composition different from the returned H/He mixture.

At 1000 K, full returned H₂+He with unit capture gives **6.35242×10⁻⁷ kg/m²**. If the same fraction `f` of both full streams becomes thermal cavity gas and both have capture `g`, the column is that value times `f/g`. `results.json` contains all **18 combinations** of total source fraction 0.001, 0.003, 0.01, 0.03, 0.1 and 1 with capture 0.01, 0.03 and 0.1. Atomic-H cases retain the escaped hydrogen mass while changing particle mass and speed. Nitrogen and hydrocarbon chemistry remain in the [separate chemistry screen](../chemistry.json).

The chronological calculation rescales frozen full-source cases at g=0.01 and 0.1. One additional linear inventory response at g=0.03 uses the archived 150 s arrival bins and exact 4,228-event schedule. All cases begin with the corresponding steady reservoir; a poor pump does not receive an empty-start advantage. This prescribed immediate thermal-release history is distinct from the still-unresolved delay and angular distributions of reflected ions or localized plumes.

## Capsule-specific heat margins

The reference is a **2 cm-radius, 80%-filled C13CH₂ capsule**, density **1011.64 kg/m³**, with the archived absorptivity, full-flight IR exposure and surface-field energy. The isotope-scaled 20→300 K enthalpy cases are **251.07 and 281.87 kJ/kg**. They represent material comparisons, not an 80% probability interval. Gas energy below assumes full incident kinetic-energy retention per core mass.

| Incoming mass-speed quantile | Speed | IR + field | Remaining gas heat | Maximum total gas column |
|---|---:|---:|---:|---:|
| Median | 2.207 Mm/s | 130.71 kJ/kg | 120.36–151.16 kJ/kg | (1.067–1.340)×10⁻⁶ kg/m² |
| 99.9% (“fast”) | 4.151 Mm/s | 77.17 kJ/kg | 173.89–204.69 kJ/kg | (4.355–5.126)×10⁻⁷ kg/m² |

The limit follows directly from `e_gas=3 Sigma v²/(8 rho fill a)`. Above the ideal 100 nm channel-routing column at g=0.01, the added-column allowance is **(1.007–1.279)×10⁻⁶ kg/m²** for the median grain and **(3.754–4.526)×10⁻⁷ kg/m²** for the fast grain.

For a common added H/He source fraction at g=0.01, the fast-grain allowance is **0.00591–0.00712 in addition to the ideal routing source**. The corresponding steady total common fraction is **0.00686–0.00807**; the actual chronological peak gives **0.00684–0.00806**. The median-grain total fraction allowance is 0.01679–0.02109. H-only and He-only additions have separate limits in the JSON. No combination of unknown release mechanisms is assigned a probability merely because it fits this allowance.

## Actual-schedule examples

These are the maximum scheduled exposure for each source-speed quantile, with both gas species assigned the tabulated total source fraction and capture. Energies include IR, field and gas.

| Total thermal fraction f | Capture g | Median 2 cm total | Fast 2 cm total |
|---:|---:|---:|---:|
| 0.001 | 0.01 | 137.89 kJ/kg | 102.58 kJ/kg |
| 0.003 | 0.01 | 152.25 kJ/kg | 153.40 kJ/kg |
| 0.010 | 0.01 | 202.51 kJ/kg | 331.25 kJ/kg |
| 0.010 | 0.03 | 154.72 kJ/kg | 162.14 kJ/kg |
| 0.030 | 0.10 | 152.53 kJ/kg | 154.39 kJ/kg |
| 0.100 | 0.10 | 203.45 kJ/kg | 334.57 kJ/kg |

The fast grain supplies the tighter mass-column constraint in these two samples because gas-impact heating grows as v². Slow capsules have a different limitation: the 1st-speed-percentile 2 cm grain already receives **392.57 kJ/kg** from full-flight IR and field energy, and the 0.1st percentile receives **829.30 kJ/kg**. Protected delivery or a qualified contained-melt route remains relevant to those populations.

At the maximum allowed total column, the geometric uniform-drag delay is **12.32–15.47 µs** for the median grain and **2.67–3.15 µs** for the fast grain. Repeatable mean shifts can be precompensated; spatial column variation, grain-size dispersion, local plume crossings, skin ablation and asymmetric ejecta affect actual assembly timing and pointing. The bulk heat comparison does not qualify those effects or other stage capsule materials.

## Lower anisotropy with a separate leakage allowance

`anisotropy-comparison.json` evaluates the independently saved high-resolution **anisotropy 10,000** case at the same 100 nm cap, 20 µm pitch and 150 µm active depth. Its flow-weighted ideal mobile front fractions are **0.0022127935 for H** and **0.0018190368 for He**. A separate plenum leakage fraction of **0.001 of the gas remaining locally captured** gives total fractions **0.0032105807 / 0.0028172178** and a steady column **1.90295×10⁻⁷ kg/m²** at g=0.01.

| Case at g=0.01 | Median 2 cm total | Fast 2 cm total |
|---|---:|---:|
| Anisotropy 10,000 plus 0.001 plenum leakage | 152.18 kJ/kg | 153.16 kJ/kg |
| Also release 0.003 of the remaining retained gas | 173.62 kJ/kg | 229.02 kJ/kg |
| Also release 0.004 of the remaining retained gas | 180.77 kJ/kg | 254.31 kJ/kg |

After the 0.001 plenum term, the fast capsule permits a further conditional release fraction **0.003872–0.005090** across the two material enthalpy cases. The lower end still leaves more than 0.003 for other thermal paths in this particular steady bulk screen. This comparison does not infer a minimum anisotropy. The relevant acceptance variables are the total species-resolved cavity source, physical cleanup, capsule/delivery design and surface response.

The same anisotropy-10,000 plus 0.001-plenum case at **1 cm** gives **287.96 kJ/kg** for the median grain and **289.90 kJ/kg** for the fast grain, exceeding both full-flight solid-300 K enthalpy cases. Its use therefore needs the already identified delivery or capsule changes. The anisotropy comparison uses mean fluxes; it does not claim a new packet-specific transient calculation or resolve energetic reflection and local plume geometry.

## Physical scope and reproduction

An energetic reflected ion may traverse the cavity and stop in another part of the retained receiver. It becomes part of the thermal source only if and when it is liberated as thermal gas. Actual source populations should therefore distinguish immediate thermal release, energetic return/reimplantation and delayed thermal release. All remain station inventory unless material physically escapes.

The global capture values are apparatus requirements. The warm ion-cleanup study still has a field-interference and late-pulse collection interface to close; this packet supplies no capture credit from the ideal front channel. At these dimensions, neutral collisions can matter even at very low pressure. Uniform steady Maxwellian wall balance remains useful, while localized pumps and directed plumes require spatial transport. Higher capture and larger grains also require their own power, mass and delivery changes before a station rating can be revised.

Run `calculate.py` with the existing scientific Python environment. It reads only archived files in the parent gas directory and writes in this new folder. The original hydro, chemistry, source and plant records are unchanged. **29 arithmetic, source-normalization and enthalpy-closure checks pass, plus four leakage-additivity checks for the anisotropy comparison.** `manifest.json` pins this extension and the calculation records all frozen input hashes.
