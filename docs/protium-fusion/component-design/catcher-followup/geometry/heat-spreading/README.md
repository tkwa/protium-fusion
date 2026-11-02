# ✨ Fixed-mass receiver heat-spreading screen

Increasing the backing's effective conductivity substantially lowers its material-arrival temperature peaks while preserving the **350 kg/m² carbon inventory, actual dispatch, constant mean heat pickup and 1000 K minimum pickup temperature**. In the seven saved representative histories, the highest front temperature falls from **1496.0 K at k=50 W/(m·K)** to **1219.7 K at k=500**, and **1198.9 K at k=1000**. The 150 μm k=4 stopping/gas-collection layer and its existing support-post resistance remain in every case.

This is a conditional heat-spreading comparison. **1200 K and 1250 K are chosen screens**, useful for comparing radiation-enhanced sublimation exposure; neither is imposed as a failure temperature or a universal onset of carbon release. Previous thermal packets and the selected plant remain unchanged.

[calculate.py](calculate.py) reconstructs the finite-backing diffusion response using the existing validated kernel and actual source profiles. [results.json](results.json) contains parameters, source hashes, per-profile temperatures and durations, solved conductivity requirements, and links to complete new histories.

| Effective backing k, W/(m·K) | Highest material-phase front, K | Warmest A: seconds above 1200 K | Hottest C13: seconds above 1200 K |
|---:|---:|---:|---:|
| 50 | 1496.03 | 619.5 | 352.0 |
| 100 | 1377.17 | 324.25 | 187.25 |
| 200 | 1293.43 | 179.75 | 108.25 |
| 500 | 1219.74 | 73.5 | 42.25 |
| 1000 | 1198.89 | 0 | 0 |
| 5000 | 1189.65 | 0 | 0 |
| Isothermal backing | 1189.64 | 0 | 0 |

Durations cover the 0.25–4000 s material/decay and resolved prompt-tail window for the actual selected events. The microsecond photon flash remains a separate finite-depth spectral calculation. It is not replaced by a constant high-k assumption or held at its peak for the material pulse.

## What changes in the calculation

The source heat, areal heat capacity **630 kJ/(m²·K)**, and mean pickup **116.519 W/m²** remain fixed. Changing k alters only the nonuniform diffusion kernel and the small steady gradient correction needed to keep pickup at least 1000 K. The whole-schedule stored-energy mode remains the same, with its maximum around **1189.6–1189.8 K**. There is no fixed-temperature back face and no requirement that the engine absorb the instantaneous pulse power.

The front's series resistance is retained:

- 150 μm stopping/collection layer at k=4: **3.75×10⁻⁵ m²·K/W**.
- Existing 50 μm post gap, 25% supporting fraction and post k=50: **4.00×10⁻⁶ m²·K/W**.
- Total: **4.15×10⁻⁵ m²·K/W**.

The k=50 body histories reproduce the archived temperatures to numerical precision. All three C13 events, the selected warmest-background A/B1/B2 events and an ordinary B1 event are retained. This is not an exhaustive fine-time resimulation of every dispatched event, and the solved comparison screens apply to those representative histories. Increasing k reduces the nonuniform temperature contribution monotonically in this family.

## Conductivity and interface requirements

The **1250 K** screen needs effective backing **k=325.14 W/(m·K)**, corresponding to **1858 W/(m²·K)** through 175 mm and a nonuniform fundamental decay time of **34.36 s**. The controlling event is C13.

The **1200 K** screen needs **k=933.64 W/(m·K)**, or **5335 W/(m²·K)** through 175 mm, with a **11.96 s** fundamental time. Warm A becomes controlling. The specified surface layer and posts do not prevent this screen: in the isothermal-backing limit, the tightest allowable combined surface resistance is **1.039×10⁻⁴ m²·K/W**, above the existing 4.15×10⁻⁵. That margin is not permission to omit new contacts.

Primary measurements of compression-annealed pyrolytic graphite show basal-plane conductivity of approximately **423 W/(m·K) at 1000 K, 386 at 1200 K and 368 at 1250 K**, with conductivity falling as temperature rises. These values interpolate the report’s numerical table; they are material measurements rather than assembled-component ratings. The specimen density is about 2240 kg/m³; achieving the model's 2000 kg/m³ with parallel voids gives ideal axial values of approximately **377, 345 and 328 W/(m·K)** at those temperatures. Consequently the k≈325 assembled requirement is a credible materials direction with limited contact/damage margin. k≈934 requires substantially enhanced transport. The high-conductivity direction must run through the backing; the shallow gas-barrier laminate retains its separate basal-plane orientation. [Roth et al., Sandia SAND88-2057](https://www.osti.gov/servlets/purl/6188316).

At the largest incoming flux, a chosen 10 K allowance for any new interface permits total contact resistance only **1.13×10⁻⁵ m²·K/W**. A heat-spreader claim must include that contact, its evaporator/condensing interfaces, and its dynamic transport limit. The companion material-evidence packet distinguishes measured conductivity from assembled performance and irradiation qualifications.

## Phase-change or distributed-spreader boundary

The full incoming material-plus-delayed peak heat fluxes are **0.395 MW/m² A, 0.00206 MW/m² B1, 0.529 MW/m² B2 and 0.881 MW/m² C13**. Their total stage fluences, including prompt energy, are approximately **32.95, 0.188, 45.92 and 55.58 MJ/m²**. These define the evaporator-side requirement for a spreader transporting all incoming heat. Transporting only the nonuniform component can reduce its load, but must be modeled explicitly. The microsecond prompt energy is first buffered near its absorption depth and can enter the spreader later; its instantaneous photon power is not assigned to the vapor channel.

For a declared latent-energy scale of 4 MJ/kg, the C13 full-flux requirement is **0.220 kg/(m²·s)** and **13.89 kg/m² of cumulative phase circulation during the event**. Cumulative circulation is not installed fluid inventory: a closed return path can circulate the same liquid repeatedly. At the unchanged mean heat rate, the corresponding circulation is only about **2.91×10⁻⁵ kg/(m²·s)**. No new instantaneous heat-engine capacity is inferred from the evaporator peak.

A Na/Li heat-pipe network or other two-phase spreader needs its own vapor area, pressure drop, sonic/entrainment margin, liquid-return route, evaporator wetting, sealed inventory and contact conductance. A fine-pore capillary surface and a high-throughput liquid return can be distinct engineered structures. The working fluid stays in its own sealed circuit, transferring heat into the graphite storage body or the existing conversion-fluid pickup; it does not carry catcher gas or nuclear ash through the engine.

Voids and pipe walls require explicit geometry. Retaining 350 kg/m² carbon by increasing thickness does not by itself preserve the old photon-capture distribution: blind or staggered channels, residual shielding columns, and a changed transport calculation are possible ways to close that interface. Pipe/wick/fluid mass must be assigned within an existing service allowance or added to the plant denominator. The present k sweep grants no free apparatus or material conversion.

Heat spreading is therefore a useful numerical route to reducing the duration and severity of hot-surface carbon emission. Its effect on actual emission still depends on the measured or deliberately selected yield law; lowering the temperature does not replace that material uncertainty.
