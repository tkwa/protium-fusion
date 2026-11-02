# ✨ Repeatable receiver sector: first component checkpoint

This checkpoint turns several allowances in the 5.395 mW/kg thermal reference into a candidate physical layout. The frozen pulse family, radius, 1000/200 K thermal boundaries, 60% converter output, shaded site and explicit infrared load remain the input. Updated pulse mixtures can replace those inputs; the present numbers describe the archived family. `sector_account.py` and `sector-account.json` record input hashes, equations and requirements. No article or earlier result is changed.

## Layout and interfaces

Use replaceable receiver cassettes mounted on a hot load-spreading layer, distributed thermal pickup behind them, a compartmented radiation-shield stack, local regenerative engines and cold heat-rejection circuits, and an outboard support-guide deck. The support deck attaches to the hot receiver through slender thermally interrupted supports. The radiator faces open sky between guide lines and service structures. Source outlets penetrate designated gaps; the incoming accelerator hardware continues behind the receiver rather than filling its cavity.

The initial candidate has approximately **1 km engine/maintenance territories**, subdivided into much smaller independently replaceable catcher tiles. The territory is neither one launcher nor one incident particle packet. A 1 km territory receives 105.6 MW of mean useful hot input, contains 5×10⁹ kg of receiver, and has 10⁹ kg of generic services available. Source hardware, driver storage, engines, dedicated exchangers and radiator panels retain their separately counted budgets. No extra 1000 kg/m² is added for the same services.

The baseline catcher remains backed graphite while the recovery team compares implanted-surface extraction and chemically distinct reusable catchers. Catcher choice can change nuclear interactions, extraction depth, thermal contacts and repair flows; it does not automatically require replacing the 2.5 m absorber after each pulse. A global separating membrane is not assumed. Delayed radiation and particles require a temperature-resolved capture account before all of their energy can be assigned to the 1000 K heat reservoir.

My provisional judgment for a useful territory scale is **1 km, with a subjective 80% credible interval of 0.2–3 km**, conditional on distributed gas circuits, repeatable engines and allowing up to roughly 300 kg/m² of the existing services for thermal distribution. This is a preference inferred from the simple transport tradeoff, not a measured lower bound or an optimization result. Tile size, guide spacing and launcher spacing are independent design variables.

## Local heat transport closes much more easily than a global coolant network

The explicit screen uses closed helium headers at 2 MPa. Hot fluid runs from 970 to 990 K, with mean speed 10 m/s; cold fluid runs from 205 to 225 K at 5 m/s. These leave finite temperature differences to the 1000 K receiver, 950 K turbine inlet, 230 K compressor inlet and 200 K radiator of the constructive gas-cycle study. Helium is treated as an ideal monatomic gas. Darcy friction factor 0.015, six additional velocity heads and 80% pump efficiency are declared requirements. Two full-flow pipes of total length twice the territory width conservatively represent supply/return headers. Terminal channels, heat exchangers and isolation hardware remain separate.

Pressure walls use density 8000 kg/m³ and an **assumed 100 MPa allowable stress**. This is a material requirement at the actual temperature and life, not an asserted property of a particular alloy.

| Territory width | Hot/cold header wall and helium | Extra external pump load per receiver area |
|---|---:|---:|
| 100 m | 8.55 kg/m² | 0.656 W/m² |
| 1 km | 85.5 kg/m² | 0.656 W/m² |
| 10 km | 855 kg/m² | 0.656 W/m² |

The increase with territory width comes from longer full-flow headers. A distributed branching network can improve this bound. Very large territories consume the service allowance; kilometre or smaller loops have substantial room for other functions. The extra pump load is approximately 1.5% of the reference net export. It is charged outside the converter package and is not a second charge for the compressor or the cycle's existing internal loss allowance.

Adding these coolant circuits creates additional heat-transfer interfaces. The script explicitly evaluates their logarithmic mean temperature differences rather than retaining the old minimum conductance unchanged. Conservative isothermal hot/cold cycle bounds give **18.11 W/K/m²** total external conductance, versus 2.71 W/K/m² in the original abstract external-exchanger minimum. At the selected 10 W/(K kg) equipment rating, this requires **1.54 kg/m² additional exchanger mass**. Existing recuperator conductance stays in the account. The finite-loop addition therefore modestly changes the mass budget; it does not require abandoning 60% conversion. A component flow map must still show that heat-transfer coefficients and pressure drops can coexist.

Through-thickness mean conduction across 2.5 m of graphite at the inherited k=100 W/(m K) gives approximately **2.64 K** mean gradient. The corresponding diffusion time at the inherited diffusivity is about 1.45 days. This supports slow heat pickup from the existing bulk while the saved transient surface calculation handles the pulse. It does not justify lateral conduction across a kilometre or an instantaneous uniform temperature.

The current material-property inputs are scenarios. As a present-day scale reference, SGL's R6650 graphite datasheet gives 12.5 GPa dynamic modulus and 95 W/(m K) conductivity at 20°C; the latter cannot be silently carried to irradiated 1000 K service. [SGL R6650 primary datasheet](https://www.sglcarbon.com/pdf/SGL-Datasheet-SIGRAFINE-R6650-EN.pdf). The NASA Brayton study supplies the relevant component/pressure-loss methodology, not these future ratings. [NASA converter study](https://ntrs.nasa.gov/api/citations/20050244475/downloads/20050244475.pdf).

## Insulation and openings: separate the real constraints

The inherited 1 W/m² bypass still requires a specified multilayer insulation system. The companion cold-sink analysis supplies finite-emissivity shield examples. The support strength floor itself does **not** force large heat leakage: even four times the combined static plus maximum pulse pressure carried at 10 MPa only needs a solid tie cross-section fraction of 5.6×10⁻⁸. With k=10 W/(m K), a 1 m path leaks approximately 0.00045 W/m². Actual joints, buckling stability, rails, service penetrations and fault loads will enlarge that floor. The conclusion is that poor attachment geometry can spoil insulation, while gravity alone does not require it.

The source team's current 100 m-radius apertures, numbering 10⁵–10⁶, occupy only 2.1×10⁻¹⁰–2.1×10⁻⁹ of the outer area. Even a direct black opening from 1000 to 200 K then contributes only **1.2×10⁻⁵–1.2×10⁻⁴ W/m²** globally. Thus their area is not a major station-wide thermal leak. Individual cold fuel channels still need protection from cavity irradiation; small global leakage does not make their contents cold. The script also preserves larger-opening sensitivities so a changed source geometry can be tested.

A 200 m-diameter outlet can fit between guide lines on an approximately 1 km lattice, with positions selected to avoid all three guide families. The source positions have millions-of-metres spacing, so kilometre-scale local placement freedom should be examined before inventing stream overpasses through every outlet. Outboard accelerator tubes and their supports must be included in the actual sky-view model.

## Surface processing is a more promising design discriminator

The reference factory capacity is 3.86×10¹¹ kg/s, or **0.801 kg of handled material per receiver m² per year**. Its mass-to-flow coefficient is retained at 10⁸ kg of factory per kg/s processed. Lower processing electricity does not make that factory smaller. The coefficient is an inherited equipment allocation, not a physical throughput limit for nanotechnology.

For hypothetical graphite skin handling at 2 MJ/kg effective total cost:

- Treating **7.20 µm after every pulse** consumes 10% of reference net electricity.
- Treating **435 µm after the approximately 78.5 large pulses/year** consumes the same fraction.
- Ten micrometres after every pulse therefore costs approximately 14% of net export. A much thicker whole-skin process on each frequent B1 pulse changes the architecture or requires a lower-energy recovery route.

These are handled depths, not predicted erosion and not fresh imported carbon. Selective release of implanted material, reuse of a chemically distinct thin catcher, treating deeper layers only after large pulses, and higher factory throughput are concrete alternatives. The recovery team's actual spectrum initially places 99.9% of CNO stopping within approximately 5 µm for B1 and 22 µm for C13 in graphite, before straggling/roughness allowances; those preliminary depths motivate an extraction calculation rather than a millimetre replacement assumption. Their source packet is the authority for those stopping results.

The ≥50% useful free-energy conversion premise applies to declared feed/product transformations. Same-material repair can have a small bulk free-energy change. The 2 MJ/kg used here is a process-cost sensitivity; it is neither a mandatory atomization bill nor an automatic prediction from that efficiency premise. A hot-to-cold processing route must account for sensible-heat recovery and reheating within its actual total cost.

## What this checkpoint resolves and what it leaves

Local heat transport, finite exchanger approaches and aperture leakage fit within the existing allocations for a useful family of sectors. They do not currently justify a major reduction in specific output. The consequential next choices are the catcher/extraction route, an actual magnetic guide and outward-fault topology, and a coupled nonradial structural model. The guide topology and common/differential motion distinction are developed in [support.md](support.md). An updated pulse mixture can be inserted without imposing the old source energy or carbon-budget split as a physical constraint.

## Capital carbon and the operating supply ceiling

The 0.007 C:H figure is an **available operational fresh-feed ceiling**. Retired skins, displaced target material and contaminated inventory can remain inside the closed receiver and be recycled. The former 0.0014 structural share is an allocation to revisit, not a requirement to lose that much carbon. Gross processing work, unavailable inventory during treatment and factory capacity still matter even when all mass returns.

Capital stock is a different scale. The reference graphite receiver contains **7.60×10²² kg of carbon**, about 62% of total installed mass. Twenty years of its reference 2.53×10⁶ kg/s fresh hydrogen feed amounts to only **1.60×10¹⁵ kg**. Thus startup receiver carbon divided by twenty-year hydrogen consumption is approximately **4.75×10⁷ kg/kg**. This comparison describes resource amortization; the operating 0.007 ceiling is not imposed on construction.

If that receiver carbon is made by methane decomposition, stoichiometry also supplies approximately **2.55×10²² kg of hydrogen**—over 300 million years of the reference station's hydrogen flow, if all were usable and retained. The actual project can distribute that byproduct elsewhere. Its existence illustrates why initial material acquisition, long-lived stock and ongoing fuel purchases require separate ledgers.

At the documented standard-state methane-to-graphite free-energy change, the carbon chemical step alone requires about **6.43×10²⁹ J at 50% efficiency**, equivalent to about 31 years of reference export. This is one declared recipe step, not a total construction estimate. Extraction, delivery, assembly, temperature preparation and any chemical-energy coproducts need their own consistent boundaries. `capital_stock.py` and `capital-stock.json` preserve these examples; the underlying thermochemistry is sourced in the [manufacturing interpretation](../manufacturing-efficiency/interpretation.md).

## Named allocation at the improved mixture point

[sector-allocations.json](sector-allocations.json) applies the same sector equations to the frozen 5.395 mW/kg point and the separately saved `plant-scan-003` mixture point. It does not change either receiver denominator or cadence. At the improved point the hot loop carries 123.40 W/m²; the cold loop carries 73.49 W/m² including the specified environmental radiation. The named subdivisions of the existing 1,000 kg/m² service allowance are:

| Service | kg per m² of receiver |
| --- | ---: |
| Hot and cold helium headers and gas inventory | 99.116 |
| Added loop-interface heat exchangers, above retained external exchanger allocation | 1.778 |
| Cold cross-braced frame package requirement | 20 |
| Guide conductor, cryostats, trim/control and mechanical package requirement | 50 |
| Named subtotal | 170.894 |
| Remainder for other services | 829.106 |

The last two packages are requirements with partial component calculations, rather than complete bills of materials. The guide conductor itself occupies 8.02 kg/m² at the selected 30 A/mm² engineering current density. A lower-current-density alternative can use a 200 kg/m² guide package while leaving 679 kg/m² of services. Dedicated station converter hardware, the retained internal/external exchangers, the 100 kg/m² insulation allowance and the radiator panels remain in their existing separate allocations.

The external helium pumps draw 0.74368 W/m² and the reference guide's radiative cryocooling draws 0.31946 W/m². Itemizing their 1.06314 W/m² inside the existing 5.86379 W/m² ordinary auxiliary allowance leaves 4.80065 W/m². This preserves the current plant output and radiator constraint. Charging these loads again as incremental consumption is a different sensitivity, requiring a smaller pulse rate or warmer/more extensive radiators at the presently saturated point. Guide AC/control losses, pipe hardware at temperature, and the catcher recovery process remain to be closed against the remainder.

The incoming source interface has meanwhile expanded to 100,000–1,000,000 spatial outlets with a candidate exit radius up to 250 m. Their whole-shell opening fraction is only 1.29×10⁻⁹–1.29×10⁻⁸. Even a black, unbaffled 1000 K view through those openings transfers only 0.000073–0.000731 W/m² averaged across the sphere. This small global load does not protect cold incoming grains from the hot cavity; their individual shielded trajectory remains a separate launcher requirement.

The support package now includes the local guide-control requirement described in [guide-local-control.md](guide-local-control.md). A stabilized global shell alone is insufficient to establish a stabilized moving-stream rail.
