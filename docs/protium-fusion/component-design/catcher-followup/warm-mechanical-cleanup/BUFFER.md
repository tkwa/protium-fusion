# ✨ Hot foreline buffer and isolated cold outlet

A **0.5 m³ hot buffer per square metre of receiver**, fed through a finite warm molecular/drag booster, smooths the selected plant's gas pulses to **0.679–1.691 Pa** in the H₂/He reference. Its passive cold outlet carries **0.64–1.78 times the continuous mean throughput**. The buffer and first compression stages remain near 1000 K. A small gas conduit then passes through a 300 K intercept to a separate H/N/O recovery branch and the helium sorption train. The calculation retains the cold sorber's working area, radiation, cycling and refrigeration requirements.

This is a concrete conditional alternative to the broad direct cold port. The warm booster needs an assigned inlet speed and compression ratio; its inlet pressure rises with the pulse. The reference does not require the catcher to hold all its gas until a later refurbishment.

The calculation is [buffer.py](buffer.py), its complete numerical record is [buffer.json](buffer.json), and downsampled pressure histories are [buffer-curves.json](buffer-curves.json). It uses the exact delivered-B1/C13H2 source128 inputs and 4,228 dispatched events from [the frozen gas packet](../gas/README.md). It does not change the selected plant or previous pump calculations.

## Finite source and buffer equations

For each species, pressure is governed by

\[
h_s\dot p_s=q(t)-S_e(p_s-p_b/K_0),\qquad
V_b\dot p_b=S_e(p_s-p_b/K_0)-S_c[p_b-p_c\sqrt{T_h/T_c}].
\]

All volumes and speeds are per receiver area. The source plenum has void-equivalent depth **38.25 μm**, excluding the long vias. Channel and booster conductances combine in series. The specified booster inlet speed is **0.001 m³/(s·m²)** per species; with the saved channel geometry, effective speeds are **0.0008625 H₂ and 0.0008166 He**. Their source-plenum relaxation times are **0.0443 and 0.0468 s**. Thus the source response is fast, but it is explicitly finite.

The warm booster's selected zero-flow compression is **10⁶**, with a **10⁵** sensitivity. This is a future apparatus requirement. The high-compression room-temperature pump evidence and its measured operating-temperature limits are recorded in [the benchmark supplement](benchmarks/POWER-TEMPERATURE.md). Source pressure is allowed to reach the independently calculated pulse peaks: **A 4.581 Pa, B1 0.05183 Pa, B2 3.270 Pa and C13 2.341 Pa**. No constant 0.001 Pa inlet or free pressure step is imposed. At these inlet pressures, the detailed last stages must be checked for transitional gas flow and actual heat generation.

The molecular outlet has one physical effective aperture, so H₂ and He speeds scale with their thermal molecular speeds. The selected effective area fraction is **1.73596×10⁻⁹**, giving speeds **1.40656×10⁻⁶ and 9.98150×10⁻⁷ m³/(s·m²)**. The separate cold partial pressures sum to 10⁻⁴ Pa at 300 K. Thermal-transpiration backflow is included.

The linear system is integrated exactly within each constant-source time bin. Sources preserve all H and He in the dispatch, with causal arrival binning. The principal 50 s/1,024-quantile result was checked against 100 s/256 and 25 s/2,048. The finest refinement changes the peak by **0.010%**. Reported extrema are sampled bin endpoints; within-bin extrema are not searched. Exact modal integrals independently check the source, output and inventory balances.

| Hot volume per area | Minimum total pressure | Maximum total pressure | Peak outlet / mean throughput |
|---:|---:|---:|---:|
| 0.1 m³/m² | 0.556 Pa | 5.249 Pa | 5.776 |
| **0.5 m³/m²** | **0.679 Pa** | **1.691 Pa** | **1.776** |
| 1.0 m³/m² | 0.796 Pa | 1.341 Pa | 1.383 |

The operating rows start with an explicitly precharged mean-flow inventory: **1.664×10⁻⁷ kg/m²**, or about **2.53×10¹² kg** station-wide, for the 0.5 m case. Empty startup is saved separately and initially provides negligible cold-outlet pressure. Gas arriving after the finite 370-day record remains source inventory; it is not counted as lost. Lowering booster compression from 10⁶ to 10⁵ increases the small source backpressure floor while leaving the buffer swing nearly unchanged.

## The cold conduit and gas separation

One **1 mm-radius, 0.1 m-long dogleg tube** per existing one-metre-radius global-pump-module area supplies more conductance than required. Its full-open effective aperture is 3.26446×10⁻⁹; a **53.18% transmission setting** produces the reference outlet. The mouth's entire admitted black radiation is charged even while its molecular conductance is throttled:

- 1000→300 K admitted radiation upper screen: **0.0137706 W/m² receiver**.
- Tube-wall conduction for 10 μm wall, conductivity 10 W/(m·K): **0.0003428 W/m²**.
- Combined shared conduit heat: **0.0141133 W/m²**.

At the actual operating pressure minima, the reverse-flux fractions are **0.04196% H₂** and **0.01480% He**. The outlet does not deliver its mean flow at every instant: the buffer explicitly supplies the calculated varying flow. The downstream system must accept its peak, or add downstream storage. The retained He active area is not automatically resized downward.

The global cleanup foreline stays physically separate from the pulsed primary source. A small additional warm booster transfers it into the same buffer. At the specified simultaneous 10⁻¹⁰ Pa H₂ and He cavity partial pressures, the fixed H₂-sized global pump geometry supplies 0.127% of the primary mean PV flow; its addition raises mean buffer pressure by about **0.00136 Pa**. Its required booster speeds are about **8×10⁻⁷ m³/(s·m²)**. There is one shared cold outlet, not a duplicate 0.0141 W/m² debit for each inlet stream.

A **separate chemical separation/recovery stage precedes the helium sorber**. It removes hydrogen- and nitrogen/oxygen-bearing streams and routes them to the existing material loop. Some eventual nitrogen reaches the receiver as oxygen before beta decay. This packet assigns no present separation efficiency or measured volatile fraction. The parent plant ledger carries the separator's positive-work allowance and component allocation. Radioactive hold-up retains the separate shielding/heat-location requirement of the material-processing design.

If all eventual stable nitrogen forms mobile N₂, its additional continuous mean PV is **4.606×10⁻⁷ Pa·m³/(s·m²)**. The same outlet then adds **1.243 Pa** to mean pressure. A deliberately bursty sensitivity releases each event's full nitrogen inventory at its earliest resolved material arrival; the 0.5 m buffer spans **1.580–3.523 Pa** when added to H₂/He. At 10% mobile N, it spans **0.769–1.874 Pa**. These are release-profile sensitivities. The N₂ molecular-speed proxy uses the actual continuous N14/N15 inventory-weighted mass; stage-dependent isotopologues and N/O chemical species remain unresolved. The finite N booster can smooth their impulses and needs its own throughput qualification. Keeping a roughly 2.2 Pa mean in the all-N case avoids enlarging the cold mouth. Keeping a 1 Pa mean instead requires 2.243 times the reference outlet conductance, modestly exceeding the presently specified full-open tube; that alternative needs a larger or shorter conduit.

## Hot storage box and allocated hardware

The gas inventory occupies **a half-metre gap**, which is negligible relative to the station radius. It is built as locally tied cassettes, rather than an unsupported global pressure vessel. One explicit mechanical scale is:

- One-metre cassette span, 0.5 m depth, two 0.5 mm hot dense-composite faces and 0.5 mm shared sidewalls: **2.7 kg/m²** at density 1800 kg/m³.
- 0.1 m face grid with 10 mm-deep, 0.5 mm-thick ribs on both faces: **0.36 kg/m²**.
- 1 mm² tension ties on that grid across the gap: **0.09 kg/m²**.

The total is **3.15 kg/m² before hermetic liner, seals, valves and fittings**. At a chosen **10 Pa** design gauge pressure, a simple supported-strip face calculation gives **0.30 MPa** bending stress and **0.125 mm** deflection between 0.1 m supports for E=10 GPa. Tie stress is **0.10 MPa**. These are low mechanical requirements, not qualified refractory-laminate properties. The ties connect the two faces on the same grid; the face ribs are not assumed to span a free metre. Off-design handling, thermal mismatch and joint design remain to be specified.

Both walls and the buffer gas are hot. The station's cold insulation remains outside the buffer; no 0.5 m bulk conductive link to a cold radiator is added. Small cold conduits, bearing supports and instrumentation have their own heat-path account. Dense liner/seals or a collected guard envelope are necessary: ordinary porous graphite is not granted gas-tightness. For example, diverting 0.1% of the mean primary flow into a guard permits aggregate leak conductance **7.41×10⁻¹⁰ m³/(s·m²)** at the reference peak pressure. A guard's gas is retained and recompressed. Routing such leakage into the central cavity instead would invoke its separate gas-column limit.

The primary booster needs only **1.984×10⁻⁵ pump-face area per receiver area** if it preserves the benchmark capture coefficient at 1000 K. Scaling the benchmark package areal mass gives **0.01165 kg/m²**. Scaling the older measured 48 W/260 mm pump-face idle input gives **0.01794 W/m²**; scaling the official modern nEXT85H 20 W ultimate-pressure input gives **about 0.1273 W/m²**. Neither small-pump figure demonstrates a hot implementation. They are transparent hardware-scale comparisons, separate from the global pump's chosen low-idle acceptance requirement. The primary mean-state ideal compression energy at 50% useful-work efficiency is only **1.63×10⁻⁵ W/m²**; it does not replace motor/bearing/windage accounting.

The 3.15 kg/m² box and small booster fit within the parent's **300 kg/m² total new-service allowance**, alongside roughly 144 kg/m² global pump, 35 kg/m² inlet baffle and 14.3 kg/m² declared rework apparatus. They are named uses of that allowance, not additions on top of it. The remaining allowance covers actual liner/seal/valve/thermal-support hardware; a complete bill of materials and evidence-based probabilistic mass bound are not available. No new whole-layer atomization or mandatory full-component replacement cycle is introduced.

## Controlled release during regeneration

The prompt-flow calculation and a trapped cohort are alternative routes for the same atoms. One mean calendar month of the entire H₂/He source is **3.293 Pa·m³/m² at 1000 K**. Releasing it at once into the 0.5 m buffer produces a **6.587 Pa pressure increment**. This can materially alter local cold-bank loading even though total annual throughput is unchanged. The nitrogen all-mobile month adds a further roughly **2.42 Pa** increment; the chosen 10 Pa mechanical screen is not a universal allowance for arbitrary simultaneous releases.

Use a metered valve from the holding cell or pressure buffer, and stagger release among tiles sharing an actual header and sorber bank. A scenario with total discharge limited to **twice the continuous mean** preserves a definite downstream capacity requirement. If fraction *f* is retained for a month while fraction *(1−f)* arrives promptly, the shortest discharge interval at that cap is *f/(1+f)* months: **7.23 hours at f=0.01, 2.77 days at f=0.1, and 15.22 days at f=1**. Atoms can remain in their isolated traps/holding cells until admitted. A cold-outlet regulating valve can also hold excess gas in the hot buffer while respecting its pressure rating.

Staggering across independent cassettes does not average the load on an isolated local sorber; it works only over the group that shares the header and cold-bank capacity. The parent retains hourly two-bank helium cycling and its active area. Neither the hot buffer nor efficient refurbishment grants unbounded regeneration throughput.
