# Manufacturing energy correction ✨

Thomas clarified the manufacturing premise on 30 October 2026: advanced nanotechnology converts lower-energy reactants plus electricity or readily accessible chemical energy into higher-energy products at **at least 50% efficiency**, with higher efficiencies common. This is the governing capability assumption. The earlier uniform 100 MJ/kg charge remains a reproducible sensitivity case; it does not establish a necessary manufacturing cost or a selected construction-energy repayment time.

The 0.007 kg fresh carbon per kg fresh hydrogen operating budget is unchanged. Energy efficiency and material return are separate quantities.

## Energy accounting

Specify the feed, all products and byproducts, temperature, pressure and material concentrations. For a positive useful free-energy increase ΔG, a conversion with efficiency η requires usable input ΔG/η and dissipates ΔG(1/η−1). The accepted floor is η≥0.5. For the same conversion, ordinary processing losses already represented by η are counted once. If supplementary chemical fuel powers the conversion, count its usable free energy alongside electricity and retain its material products.

A transformation that returns the same material to the same bulk chemical state can have very little bulk free-energy change. Its cost follows the actual separation, new surfaces, defects, thermal and mechanical handling, and dissipation. Applying an efficiency ratio to a zero bulk change does not determine these costs. Atomizing every kilogram and discarding the released bonding energy on reassembly is one possible process; its energy bill is not a universal requirement for recycling.

The distinction between reversible separation work and losses is standard thermodynamics; see the [Bureau of Reclamation's thermodynamic separation analysis](https://www.usbr.gov/research/dwpr/reportpdfs/report078.pdf), section 1.3.

## Drexler source check

Drexler's *Nanosystems*, [§13.3.8](https://nanosystems.contact.ms/molecular_sorting_processing_and_assembly/), predicts greater than 99% efficiency for certain mechanochemical energy conversions, using hydrogen oxidation as an example and extending the argument to small organic molecules. This is a theoretical apparatus estimate. Its feed ordering and entropy account belong in the complete process.

[§14.4.8](https://nanosystems.contact.ms/molecular_manufacturing_systems/) lists 1.5 MJ/kg for molecular operations and 0.5 MJ/kg for block assembly, with computation separately allowed. Its acetone/oxygen/diamond example subsequently states 3.1 MJ/kg free-energy dissipation and about 4.8 MJ/kg waste heat; the two totals are retained as stated. A dissipation estimate excludes usable energy retained in a product. The example also consumes chemically energetic feedstock and oxygen. These figures support a much more efficient manufacturing model than the blanket 100 MJ/kg assumption, while preserving a process-specific balance.

## Fresh carbon and returned graphite

For CH₄(g) → C(graphite) + 2H₂(g) at 298.15 K and standard 1 bar states, the useful free-energy increase is 4.227 MJ per kilogram of carbon. The usable input is 8.454 MJ/kg at 50% efficiency or 4.270 MJ/kg at 99%. Both carbon and hydrogen remain useful products. Returned graphite instead needs a refurbishment process whose energy follows its actual damage, contamination and handling. [NIST-JANAF methane data](https://janaf.nist.gov/tables/C-067.html).

The independent [manufacturing audit](interpretation.md) and [reproducible examples](energy_examples.py) distinguish those flows. At the existing station's entire structural-carbon makeup allowance of 3548.6 kg/s, fresh methane conversion at 50% would consume 30.0 GW. Gross receiver refurbishment processes 2.409×10¹¹ kg/s, so its per-kilogram cost has a much larger effect. The examples also report thermal, surface and separation scales, and label factory mass reallocation as a secondary sensitivity.

## Station sensitivity

[sensitivity.py](sensitivity.py) changes the recurring manufacturing energy at fixed station mass and installed hardware. It retains the original factory capacity, radiators, reserve and support allocations. No smaller factory or faster throughput is granted by reducing its electricity consumption. The total refurbished mass remains 3.85×10¹¹ kg/s.

| Effective total refurbishment cost | Refurbishment electricity | Net station specific output |
| --- | ---: | ---: |
| 100 MJ/kg | 38.5 EW | 2.766 mW/kg |
| 10 MJ/kg | 3.85 EW | 3.051 mW/kg |
| 3.1 MJ/kg | 1.19 EW | 3.073 mW/kg |
| 2 MJ/kg | 0.770 EW | 3.076 mW/kg |
| 1 MJ/kg | 0.385 EW | 3.079 mW/kg |

EW denotes 10¹⁸ W. These rows vary **total effective cost**; the Drexler dissipation numbers do not automatically supply that total. At fixed hardware, eliminating the entire old manufacturing-electricity charge raises output to 3.083 mW/kg. Construction repayment is much more sensitive: the saved table gives about 10.3 years at 1 MJ/kg and 104 years at 10 MJ/kg, conditional on those being the full manufacturing energy per installed kilogram. Feed preparation and transport still require their specified accounts.

The former 1,146-year fabrication repayment and the finding that a 1 μm silica membrane renewed after every pulse exhausts the power surplus both depend on the 100 MJ/kg convention. They cannot determine the revised outcome under the accepted capability premise.

## Reproduction and next integration

Run `python3 docs/protium-fusion/manufacturing-efficiency/sensitivity.py`. It uses the frozen `open-loop-redesign/operating-cases.json`, saves input/code hashes, reproduces the 100 MJ/kg case and checks its power balance. Results are in [sensitivity.json](sensitivity.json).

A complete revised selected station requires component feed/product recipes and maintenance routes. The completed [membrane assessment checkpoint](../membrane-receiver/README.md) applies the corrected accounting to renewal and records the remaining cooling, collection and lifetime requirements. The newer [shaded-receiver account](../cold-sink-redesign/README.md) integrates the 2 MJ/kg routine-refurbishment sensitivity with revised heat conversion, infrared surroundings and explicit exchanger mass. Its working point is about 5.4 mW/kg. Initial construction recipes and component-specific maintenance costs still require their own material and free-energy accounts.
