# Manufacturing energy and economic implications ✨

The current plant is an infrastructure-intensive way to obtain energy from hydrogen. Its fuel carries ample energy; the economic burden is assembling, supplying and maintaining a very large receiver and its services. Reducing the receiver from 5,000 to 2,000 kg/m² nearly doubles conditional output per installed kilogram. Efficient manufacturing helps both versions, but construction feedstock and delivery remain consequential at these specific powers.

This calculation uses the **64-zone median-carbon-rate, fixed-reference-profile family**: A-H3/B1-H2/B2-H2/C13-H4. The 5,000 kg/m² case is the controlled mixture comparison; 2,000 kg/m² is a separate component candidate with unresolved capture, surface and apparatus qualification. Source-profile and finer-mesh revisions are being integrated separately. [economics-inputs.json](economics-inputs.json) names and hashes the precise plant records; [economics.json](economics.json) contains the derived quantities. The calculator accepts replacement input records without changing its formulas.

## What the installed mass buys

| Quantity | 5,000 kg/m² receiver | 2,000 kg/m² receiver candidate |
|---|---:|---:|
| Conditional net electricity | 7.994 × 10²⁰ W | 8.006 × 10²⁰ W |
| Complete installed mass | 1.2192 × 10²³ kg | 6.4801 × 10²² kg |
| Specific power | 6.557 mW/kg | 12.355 mW/kg |
| Installed mass per net watt | 152.5 kg/W | 80.94 kg/W |
| Graphite receiver carbon | 7.603 × 10²² kg | 3.041 × 10²² kg |
| Receiver share of installed mass | 62.36% | 46.93% |
| Services share | 12.47% | 23.46% |
| Installed reserve share | 20% | 20% |

The receiver is the largest named mass item. Receiver, services and reserve together account for **94.8% and 90.4%** of the respective totals. This makes thermal/capture qualification of a thinner receiver economically important. It also means that optimizing minor isotope-processing machinery has limited effect on the denominator. A detailed services bill could still change the result substantially; its 1,000 kg/m² allowance is retained in both comparisons.

Each kilogram of freshly consumed hydrogen produces about **74.3 GWh of net electricity** in these ledgers. The 2,000 kg/m² point consumes about **2.99 million kg/s H and 1,073 kg/s fresh C**. Preparing the pulses already draws **2.36 × 10²⁰ W**, about 29.5% of net export. That pulse-preparation cost is an operating energy account distinct from factory construction or cold chemical sorting.

## Manufacturing efficiency and energy payback

Use Thomas's premise directly for a specified positive useful free-energy change: `E_supplied = ΔG / η`, with η ≥ 0.5 and often higher. Supplied chemical free energy and electricity belong in the same balance. The conversion's dissipative loss is already included in this expression. For recycled graphite returning to the same bulk state, the net bulk chemical ΔG can be small; actual purification, repair, thermal cycling and handling determine the process cost.

Drexler's molecular-manufacturing analysis separates chemical feedstock energetics, molecular processing, assembly and reuse. Its illustrative process dissipation supports investigating modest fabrication costs, while the required feedstock and throughput still need to be specified. The existing [manufacturing interpretation](../manufacturing-efficiency/interpretation.md) preserves the distinct 2 and 3.1 MJ/kg dissipation examples. [Nanosystems, chapter 14](https://nanosystems.contact.ms/molecular_manufacturing_systems/).

For actual whole-plant embodied energy `e` joules per installed kilogram and net specific power `q`, simple energy payback is `t = e/q`. Thus every **1 MJ/kg** takes **4.83 years** of export in the 5,000 kg/m² case or **2.56 years** in the 2,000 kg/m² case.

| Actual average embodied energy | Payback, 5,000 kg/m² | Payback, 2,000 kg/m² |
|---|---:|---:|
| 1 MJ/kg | 4.83 years | 2.56 years |
| 2 MJ/kg | 9.67 years | 5.13 years |
| 10 MJ/kg | 48.3 years | 25.6 years |
| 100 MJ/kg | 483 years | 256 years |

These are cost scenarios, not estimates of the actual construction recipe. Twenty years of continuous operation supplies **4.14 or 7.80 MJ per installed kilogram**; a 100-year life supplies **20.7 or 39.0 MJ/kg**. Those are useful acceptance thresholds for any proposed acquisition-and-construction chain. Availability below one lengthens payback in proportion to `1/availability`, before accounting for fixed loads during outages.

For one explicit chemical step, standard-state methane decomposition, CH₄(g) → graphite + 2 H₂(g), JANAF gives ΔG = **50.768 kJ/mol at 298.15 K and 1 bar**. This is **4.227 MJ/kg C** of useful free-energy increase, requiring **8.454 MJ/kg C at 50% efficiency** or **4.270 MJ/kg C at 99%**. [NIST-JANAF methane](https://janaf.nist.gov/tables/C-067.html).

Applied once to the receiver carbon, that step takes **25.5 or 10.2 years** of the respective plant's export at 50% efficiency; at 99%, **12.9 or 5.14 years**. Existing condensed carbon changes the initial chemical endpoint and can avoid that methane upgrade. Extraction, selected product structure, delivery, installation and thermal preparation then determine their own costs. Applying fresh methane chemistry to every kilogram of recirculated graphite would charge repeatedly for carbon already present.

Energy payback does not specify construction time. For illustration, an external **10²⁴ W** industrial power supply could provide the 2,000 kg/m² receiver's 50%-efficiency chemical-step energy in about **three days**; the material throughput and factory deployment would set a much longer schedule. A late deployment date therefore needs an account of the existing industrial base, rather than treating the completed reactor's own output as its only construction supply.

## Construction carbon and the 0.007 C:H ceiling

At complete eventual ordinary-material return, this operating family needs **fresh C:H ≈ 0.00035864 by mass**. Nuclear conversion to retained side products produces that demand. The original 0.0056 fuel / 0.0014 structural division is an allocation within the 0.007 fresh-feed ceiling. Retained skins, residues, sorting tails and material waiting for processing use inventory and work; they consume fresh material only when an identified terminal stream leaves the reusable system.

The startup graphite inventory is a separate resource requirement. With a resource having C:H = 0.007, obtaining the receiver's carbon corresponds to accessing at least **1.086 × 10²⁵ kg H** for the 5,000 kg/m² case or **4.344 × 10²⁴ kg H** for the 2,000 kg/m² case, even with complete carbon recovery. This is the hydrogen associated with the accessed resource; it need not all be transported, stored at the reactor, or burned. A carbon-rich source can substantially reduce the material-access burden.

By comparison, twenty years of reactor operation consumes only **1.89 × 10¹⁵ kg H**. The thinner receiver's startup carbon is roughly **2.30 billion times** the carbon allowed by applying 0.007 solely to that twenty-year hydrogen flow. This does not violate the operating ceiling: it establishes that construction requires a much larger resource programme than fuel supply. If every receiver kilogram came from methane, the resulting **1.02 × 10²² kg H coproduct** would cover about **108 million years** of this reactor's hydrogen consumption. The resource programme must therefore have uses or destinations for its coproducts.

Receiver carbon alone is a lower bound on the selected plant's total carbon stock. The remaining services, drivers, reserves and insulation require an elemental bill before total resource access can be stated. The 0.007 resource composition is Thomas's specified comparison, not a universal composition assigned to every astronomical source.

## Recurring work, retained material and factory throughput

The existing routine-repair account handles 10⁻⁴ of installed mass annually. At **2 MJ/kg**, its electricity is **0.0967%** of reference export for the 5,000 kg/m² point or **0.0513%** for the 2,000 kg/m² point. At 100 MJ/kg, the corresponding fixed-throughput burdens would be **4.83% and 2.56%**. Only the increase over the included repair load would be a new debit. These fixed-ledger comparisons do not resize the thermal system or change its capacity rate.

Returned excess hydrogen is a smaller energy burden: even **1 GJ/kg** for its 3.72 million kg/s flow costs about **0.000464%** of net export. Local gas capture and short residence remain important because hydrogen in the incoming-grain path affects source heating. The energy calculation cannot replace that geometric requirement. Likewise, annual skin processing retains its actual mass throughput and can occupy a named part of the ordinary-auxiliary allowance; it must not be added again if that allowance already covers it. [Material-loop account](material-loop.md), [return-gas account](return-gas.md).

The selected repair-factory mass coefficient is **10⁸ kg machinery per kg/s processed**. The 2,000 kg/m² plant therefore includes about **2.05 × 10¹⁹ kg** of repair machinery for **2.05 × 10¹¹ kg/s** throughput. Building that entire plant in twenty years requires an average **1.03 × 10¹⁴ kg/s** installation rate, **500 times** the installed routine-repair throughput. If the same productivity applied to new construction, the equivalent external construction machinery would weigh about **1.03 × 10²² kg**. Those are explicit throughput scenarios; a different construction process or growing fleet changes them. Manufacturing energy efficiency alone supplies no factory growth rate or completion date.

External delivery can also exceed fabrication energy. A parcel's kinetic energy at 30 km/s is **450 MJ/kg**. That is a mechanical-energy scale, not a compulsory dissipative shipping charge: orbital trajectories, exchanged momentum, braking recovery and where the material is sourced determine the supplied work. Long-lived local infrastructure and local recycling can therefore be economically preferable even when fabrication is highly efficient.

## Alternatives and what the comparison establishes

An exposed solar collector with the declared comparison assumptions of **30% conversion and 1 kg/m² complete collector mass** supplies 16.33 W/kg at 5 AU. Its mass-specific output equals these reactor cases near **250 AU and 182 AU**, respectively. The inverse-square calculation uses the IAU nominal irradiance **1,361 W/m² at 1 AU**. [IAU 2015 Resolution B3](https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf).

Collector mass, illumination and delivery are explicit variables: the parity distance scales as the square root of `conversion × delivery efficiency / collector areal mass`. At the shaded reactor site, power supplied by an existing swarm needs a transmission and receiver account. This simple collector comparison motivates that economic study; it does not by itself price delivered electricity.

The strongest potential economic role is long-lived, independently supplied power where external energy delivery is costly, together with development of a hydrogen fuel cycle. Near a star with accessible power infrastructure, the reactor needs a concrete advantage in delivered service or independence to justify its bulk. No money price is inferred from thermodynamic efficiency. A financial estimate needs the chosen materials and delivery chain, manufacturing capacity, usable lifetime, availability and the opportunity cost of those resources.

## Scope and reproduction

The arithmetic is deterministic conditional on its inputs. No empirical or subjective 80% interval is claimed for the unknown complete construction cost. The scenario table can be compared with any subsequently proposed 80% cost range, with payback endpoints obtained by the same linear conversion. An actual range should specify the acquisition route, finished material recipes, recycling boundary and transport recovery assumptions.

Run `PYTHONDONTWRITEBYTECODE=1 python3 docs/protium-fusion/component-design/economics.py`. Five identities per plant check specific power, receiver area/mass, the complete mass sum, included refurbishment and the energy-payback conversion. Existing research snapshots and plant operating records remain unchanged.
