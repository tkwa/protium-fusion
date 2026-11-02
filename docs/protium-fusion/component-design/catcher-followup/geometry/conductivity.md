# ✨ Conductivity and gas-interface supplement

The 150 µm active layer does not require exactly 4 W/(m K) to carry the material-arrival pulse. On the frozen exact-family temperature histories, 0.5 W/(m K) produces a maximum of 1,687 K with the existing conductive-post geometry. This result describes the slower material/decay heating; the prompt spectral-depth transient remains a separate calculation at the candidate conductivity.

[conductivity.py](conductivity.py) reads the preserved [temperature histories](temperature-histories.json) and [geometry](results.json), then writes [conductivity.json](conductivity.json). It does not modify those inputs. The seven detailed profiles contain all three C13 events and the named A/B1/B2 histories. The archived finite-backing model retains a constant mean pickup; the new skin resistance does not require an engine to follow peak heat flux.

## Active layer and post gap

The added temperature at each saved time is

`ΔT(t) = q_material+decay(t) [d/k_normal + h_gap/(f_post k_post)]`.

Here d=150 µm, h_gap=50 µm, f_post=0.25 and k_post=50 W/(m K). The post contribution is 4×10⁻⁶ m² K/W. A small intact-layer resistance already represented in the backing solution is not subtracted, making this additive screen conservative.

| Active-layer normal k, W/(m K) | Total added resistance, m² K/W | Highest saved material-phase front, K |
| ---: | ---: | ---: |
| 0.2 | 7.54×10⁻⁴ | 2,064.03 |
| 0.5 | 3.04×10⁻⁴ | 1,687.48 |
| 1 | 1.54×10⁻⁴ | 1,572.11 |
| 2 | 7.90×10⁻⁵ | 1,520.16 |
| 4 | 4.15×10⁻⁵ | 1,496.03 |
| 10 | 1.90×10⁻⁵ | 1,482.16 |

C13 event 1 controls every row. The through-layer diffusion scale `ρ cp d²/k` ranges from 0.405 s at k=0.2 to 0.0081 s at k=10, short compared with the resolved material-arrival pulse. This supports the quasi-steady resistance comparison. Temperature-dependent conductivity, interface resistance and changes in heat deposition within a damaged layer remain material qualifications.

If **1,700 K is chosen merely as an illustrative material-arrival screen**, the conductivity requirement follows directly from all saved times:

`k_min = max_t { q(t)d / [1700 − T_backing(t) − q(t)R_post] } = 0.474777 W/(m K)`.

This number is not a measured allowable, a universal operating temperature, or a complete prompt-heating requirement. The source's microsecond photon flash must be recalculated or bounded at the chosen low conductivity before using this screen as a whole-pulse temperature claim. No empirical conductivity credible interval is assigned to the irradiated laminate.

The post-gap term can be reduced by wider or better-conducting posts. At fixed 50 W/(m K), post fractions 0.1/0.25/0.5 give resistances 1×10⁻⁵/4×10⁻⁶/2×10⁻⁶ m² K/W. At the largest saved material flux, the present gap adds 3.52 K at most. The active material's resistance dominates the low-conductivity rows. Void fraction also affects gas volume and conductance, so post area cannot be optimized as an isolated thermal scalar.

## Supported cap

The 100 nm continuous cap spans a 100 nm-wide buried slit. Its in-plane conductivity is a different tensor component from the active layer's normal conductivity. With fixed-temperature edges, the upper center rise under an absorbed flux is

`ΔT_cap = q_abs w²/(8 k_parallel t_cap)`.

At k_parallel=4 W/(m K), its coefficient is 3.125×10⁻⁹ m² K/W. Even assigning the full incident C13 photon flux to the cap gives a 26.04 K span rise. At k_parallel=0.5/1/2/10, the same deliberately generous absorption screen gives 208.36/104.18/52.09/10.42 K. The cap span's diffusion time remains below 0.18 µs for the tabulated values; it is short relative to the approximately 70 µs prompt duration. Actual cap absorption can be lower, but no credit is assumed here.

The unsupported slit projection is 0.5% of area. Most of the cap rests directly on solid graphite. Its through-thickness resistance, t_cap/k_normal, is already inside the 150 µm active-layer resistance and must not be added again. The local span temperature rise is above the adjacent solid's temperature during the same instant; independently timed maxima cannot be combined as an achieved temperature. The spectral-depth response in the adjacent low-k material remains separate.

## Completed gas handoff

The [joint gas-port packet](../gas/joint-ports.json) supplies a concrete conditional interface for the [tile layout](README.md), while the [barrier review](../barrier-review/README.md) supplies its pressure/leakage requirements. The comparison uses a shared warm admission mouth occupying 10⁻⁴ of receiver area, molecular transmission 0.1, and a downstream speed target of 0.001 m³/(s·m² receiver) for each of H2 and He. Its active getter/cryo area fractions are 6.248×10⁻⁵ and 4.867×10⁻⁵.

Continuous extraction keeps the mean source-side H2+He pressure near 0.001477 Pa. The source-rate peak sums in the immediate thermal-release model reach approximately 4.58 Pa for A, 0.0518 Pa for B1, 3.27 Pa for B2 and 2.34 Pa for C13. Concentrated later desorption can have a different peak, so these are source-model values rather than universal pressure bounds. The barrier review's earlier B1 rounded comparison is superseded by the gas packet's final source-rate value.

Gas cooling from 1,000 K to 300 K is only 2.8044×10⁻⁶ W/m². The declared 50%-efficient compression example from 0.001 Pa to 1 bar takes 1.3841×10⁻⁵ W/m². The radiative interface is the larger load:

- Counting all 1,000 K radiation admitted by the shared mouth gives a **5.67037 W/m² hot-to-cold bypass upper bound**.
- A full 200 K view of the specified 5 K He sorber gives 0.00441536 W/m² of cold heat and **0.344398 W/m² refrigerator electricity** at 50% of Carnot performance with a 200 K sink.
- If 10⁻⁴ of admitted hot radiation also reaches 5 K, it adds **0.044229 W/m²** refrigerator electricity.

The hot bypass is removed from engine-eligible heat once. Radiation transferred internally from the 200 K shield to the 5 K sorber is not a second heat source; the refrigeration work is a separate electrical input. Existing insulation, bypass and service allocations must be matched to these actual components before updating the plant. The gas packet's baffle transmission, cold view, loaded capture, pump speed and regeneration isolation remain joint hardware requirements.

A panel-only mass sensitivity of 1/10/100 kg per active m² gives approximately 0.000111/0.001112/0.01112 kg per receiver m² for the two active surfaces. Those numbers exclude manifolds, valves, thermal shields, compressors, electrical connections and regeneration hardware. The existing large service allowance does not establish their construction mass.

Nine new algebra/replay checks pass, including the k=4 replay and the exact chosen-screen inversion. A separate supplement checkpoint preserves the input hashes and leaves the temperature packet unchanged.
