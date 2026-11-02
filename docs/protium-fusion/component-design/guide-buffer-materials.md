# ✨ An abundant cold-buffer branch

The additional pulse buffer can use a carbon-and-hydrogen solid instead of a station-scale copper inventory. **Measured polyethylene calorimetry supports the required stored heat; heat access during a pulse remains a component-design requirement.** This extension preserves `guide-gamma.json`, its copper benchmark, and all original station accounts. The new `guide_buffer_materials.py` and `guide-buffer-materials.json` carry the inputs, geometry and checks.

## Measured temperature interval

Chang and Bestul measured linear polyethylene SRM 1475 by adiabatic calorimetry. Table 1 gives 6.008, 6.662 and 7.252 J/mol/K at 61.44, 67.41 and 73.12 K, per 14.027 g base mole of CH₂. Piecewise-linear interpolation gives **456–494 J/kg/K at 65–70 K**, and an integrated sensible heat of **2376.7 J/kg**. A chosen 5% engineering reserve leaves **2257.9 J/kg**. This uses the solid's actual low-temperature data. The measured sample's room-temperature density was 954 kg/m³; the calculation selects an effective packed density of 900 kg/m³, with a 600 kg/m³ sensitivity. Those effective densities include volume allowance and are not claimed measurements at 65 K. [Chang and Bestul, NBS (1973), original paper](https://nvlpubs.nist.gov/nistpubs/jres/77A/jresv77An4p395_A1b.pdf).

## Same guide geometry and heat loads

The calculation retains six guide bores per kilometre, the selected winding radius, the 65–70 K storage interval, the existing cooling COP, and a local graphite enclosure with 2000 kg/m² radial column. It recomputes buffer volume, enclosure radius, outer radiation interception, shield mass and refrigerator power together. It retains all the preceding photon-transport qualifications. Existing conductor stabilizer is unchanged; only the additional calorimetric stock is replaced.

For the outward diffuse angular case, with the exact parent rows carrying 0.005 additional missed-hot fraction into cold equipment:

| Main receiver | Polyethylene buffer | Local graphite enclosure | Combined addition | Extra refrigerator electricity |
|---|---:|---:|---:|---:|
| 500 kg/m² | 16.53 kg/m² | 121.29 kg/m² | **137.82 kg/m²** | **0.925 W/m²** |
| 750 kg/m² | 11.46 kg/m² | 109.81 kg/m² | **121.27 kg/m²** | **0.636 W/m²** |
| 1000 kg/m² | 8.06 kg/m² | 100.90 kg/m² | **108.97 kg/m²** | **0.445 W/m²** |

At 500 kg/m², the copper reference was 37.70 kg/m² buffer plus 87.29 kg/m² enclosure, totaling 124.99 kg/m² and 0.583 W/m². Polyethylene lowers buffer mass but enlarges the shielded radius from 1.66 to 2.11 m because its volumetric heat capacity is lower. The additional copper inventory of **5.73×10²⁰ kg** becomes **2.51×10²⁰ kg of C/H polymer**, while local shield stock increases. This is a material-access alternative, with a modest larger mass and power requirement.

The lower-density 600 kg/m³ branch at 500 kg/m² totals 161.56 kg/m² and 1.140 W/m². Allowing only half the measured heat capacity to be accessible during a pulse gives 205.09 kg/m² and 1.346 W/m². Combining half-accessible capacity with the larger angular factor 3 gives 341.54 kg/m² and 2.997 W/m². These are explicit stress scenarios, not probability intervals.

The main branches fit within the previously unassigned approximately 829 kg/m² of services. They require a named allocation there, including additional contact and spreading hardware, rather than a second full services allowance. Extra refrigerator electricity must be assigned within retained auxiliaries or debited once; no revised plant output is claimed here. At effective emissivity 10⁻⁴, the enlarged local hot enclosure leaks another 0.179–0.237 W/m² to the 200 K shield. That requires part of the retained hot-to-cold bypass allocation or an explicit additional charge.

## Qualification and replacement boundary

A metre-scale cold polymer annulus cannot be treated as instantaneously mixed. The acceptance relation `l ≲ sqrt(k t/(ρ cp))` gives a diffusion length of only 4.83 mm in 10 s for an **assumed effective composite conductivity** of 1 W/m/K. Neither that conductivity nor that time is a measured complete-guide response. Actual spectral energy deposition, local spreaders, contacts, conductor temperature and the time-dependent refrigerator response must establish how much capacity is accessible. Distributed direct photon deposition also needs transport calculation. The half-accessible case is a useful accounting sensitivity, not proof of adequate pulse response.

The diffuse main branch receives a mean ionizing-energy equivalent of approximately 72–84 kGy/year if all assigned photon heat is deposited in the polymer. That establishes a radiation-test requirement. It does not establish a failure dose. Gamma irradiation can create persistent radicals in polyethylene; cold irradiation and subsequent thermal annealing require material-specific treatment. [Original 77 K radical study](https://www.sciencedirect.com/science/article/abs/pii/S2352492820325320).

The saved 0.1, 1 and 10 MGy replacement-dose rows are acceptance curves only. For example, a **hypothetical** 1 MGy qualified limit would imply about 12 years and 1.39 kg/m²/year of buffer reprocessing at the 500 kg/m² point. Retired polymer remains recoverable C/H stock; the calculation creates no fresh-carbon loss or purge requirement. Process throughput, cooling and energy must be allocated for whatever regeneration method is chosen, respecting the manufacturing efficiency premise.

All six copper mass comparisons reproduce the frozen reference within 5.2×10⁻¹² relative error. The new storage identities close within 2.7×10⁻¹³. The supported outcome is an abundant-material calorimetric branch with explicit volume and electrical penalties; pulse heat transfer and radiation lifetime remain open component qualifications.
