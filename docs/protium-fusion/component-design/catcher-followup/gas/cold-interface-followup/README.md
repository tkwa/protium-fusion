# ✨ Cold-panel pressure, two-bank operation and regional refrigeration

The retained pumping areas support the primary H₂/He flow at a higher cold-interface pressure than the original prescribed 0.0001 Pa. At 300 K the compatible mean is about **0.000594 Pa**, with a **0.001057 Pa** sampled combined peak. Closing this finite pumping speed raises the warm-buffer peak from 1.69060 to **1.69143 Pa** in the separate [passive replay](../../warm-mechanical-cleanup/cold-interface-correction/README.md).

An explicit two-bank sequence adds **0.009070 W/m²** to the old sorber-cycle electrical allowance when its full sensible warmup is supplied electrically. Its instantaneous cold capacity is a larger apparatus requirement: a shared station serving 1 km² of receiver needs **4.42 kW** for the active bank's 200 K radiation and up to **17.81 kW at 5 K** at the selected bank handoff. These figures preserve the existing active areas and original local hot buffer.

## Compatible pressure and conductance

The existing active areas are 6.2483931e−5 m²/m² receiver for the H getter and 4.8666936e−5 for He cryosorption. The H reference is a selected 10 m³/(s·m² active area) at 300 K within the [CERN technical range](https://cds.cern.ch/record/2929324/files/document.pdf). The He reference is the approximately 13 m³/(s·m²) low-load graph reading from [Day et al.](https://www.osti.gov/etdeweb/servlets/purl/20047879), whose speed normalization is 273.15 K. Its loading, temperature and equilibrium-pressure qualifications remain in the [primary evidence](../../helium-evidence/README.md).

At cold-reservoir temperature `Tc`, the table uses `S ∝ sqrt(Tc/Tref)`, conditional on unchanged net capture. The surface temperature of the He adsorbent remains 5 K. The reservoir-equivalent scaling does not determine getter kinetics at 200 K.

| Cold reservoir | H₂ equivalent speed, m³/(s·m² receiver) | He equivalent speed | H₂ mean pressure excess, Pa | He mean pressure excess, Pa |
|---|---:|---:|---:|---:|
| 200 K | 0.000510179 | 0.000541367 | 0.000343191 | 0.000139241 |
| 300 K | 0.000624839 | 0.000663036 | 0.000420321 | 0.000170535 |

These means use the continuous design flow. The actual finite dispatch has a slightly different mean: the passive replay gives total 0.000593865 Pa at 300 K. Pressure entries are excess over the actual surface-equilibrium pressure. The roughly 1e−6 Pa loaded-He reference at 6.5 K supplies an explicit sensitivity; the 5 K operating equilibrium still depends on coverage and gas history. Maintaining the old partial pressures instead would require at least 6.01 times the H area and 5.67 times the He area, before any additional manifold resistance.

For hot throughput `Qh`, the finite-panel boundary is

\[
Q_h=S_c[P_b-\sqrt{T_h/T_c}P_c],\qquad
P_c=P_{eq}+\frac{T_c}{T_h}\frac{Q_h}{S_{panel}},
\]

so

\[
S_{eff,h}=\left[\frac1{S_c}+\frac{\sqrt{T_c/T_h}}{S_{panel}}\right]^{-1},\qquad
P_{floor,h}=\sqrt{T_h/T_c}P_{eq}.
\]

The combined hot speeds are 1.40483248e−6 and 9.97328031e−7 m³/(s·m²), reducing the original conduit speeds by 0.123145% and 0.0823875%. The passive replay includes this change in the buffer time constant. A separate cold manifold with conductance at least ten times the active-panel speed adds at most 10% to the panel-side pressure rise. That manifold is distinct from the already charged hot outlet.

The retained He area receives only 19.04% of the selected low-load capacity scale in one mean hour; the maximum saved He outlet rate raises this to 19.17%. Smaller active area is an optional trade against pressure and loading-dependent speed. The reference retains the original area.

## Explicit two-bank operation

One bank captures at 5 K for 3,600 s. The other remains isolated from the gas feed and direct 1000 K view while it warms and regenerates, then waits near 200 K. The example reserves 300 s for warming/regeneration, 3,200 s for standby, and 100 s for cooldown. The final cooldown is `T(t)=200 exp[−ln(40)t/100]` K. Regenerated gas goes to an isolated product line.

The model explicitly extends the selected constant thermal capacity of 1,000 J/(m² active area·K) from the old 5–85 K comparison to 5–200 K. Actual temperature-dependent heat capacity, adsorbed gas, supports and exchanger material should replace this choice. A finite 300 s approach to nominal 200 K requires a slightly warmer switched source or final electrical trim; exact 200 K is the bookkeeping endpoint. The fully electrical warmup below supplies a concrete conservative alternative to an unspecified thermal switch.

For each bank, absorbed radiation is `a σ(200⁴−T⁴)`, and refrigerator work is the extracted heat times `(200/T−1)/ηCarnot`. Radiation absorbed during warmup contributes to warming; it is not also refrigerated. At the central `ηCarnot=0.3`:

| Electrical term | W/m² receiver |
|---|---:|
| Active bank, 200 K radiation | 0.5739961 |
| Inactive bank, sensible cooldown from 200 K | 0.0244586 |
| Inactive bank, radiation during cooldown | 0.0039042 |
| Inactive bank, 200 K standby | 0 |
| Full sensible electrical warmup | 0.0026361 |
| New inactive-bank total | **0.0309989** |
| Old 85 K sensible-cooldown term | 0.0219291 |
| Additional electrical allowance | **0.0090698** |

With a qualified switched warm source supplying the sensible heat from the existing 200 K waste reservoir, the increment is 0.0064337 W/m². Borrowed sensible heat and its later return cancel in the global cold balance; refrigerator work remains an electrical debit. A permanently cold inactive bank would instead add another 0.5740 W/m² of refrigeration. An 85 K standby bank has a larger radiative holding cost than the chosen near-200 K standby.

The 100 s cooldown gives an instantaneous two-bank electrical maximum of **2.31491 W/m²**, and an end-of-cooldown 5 K load of **0.0178070 W/m²**. These are local or physically connected shared-station capacity requirements. Staggering independent units only in a plant spreadsheet does not reduce each unit's equipment rating.

`bank_cycle.py` provides a pure callable `bank_cycle(sink_K=..., fraction_Carnot=..., electrical_warming=True)` for the plant ledger. Its `additional_to_old_85K_sensible_cycle_W_m2` replaces the missing inactive-bank terms at the actual radiator temperature. Existing active-bank radiation, admitted hot photons, gas cooling, compression and getter regeneration remain separately charged. There is no new hot-to-cold photon source in this cycle correction. Any new hot view, support conduction, adsorption heat or electrical warmup beyond these terms must be included once.

The central efficiency 0.3 and subjective 80% interval 0.15–0.5 apply conditionally to a large industrial helium cycle rejecting near 200 K. [CERN's industrial refrigerator report](https://cds.cern.ch/record/1507611/files/CERN-ATS-2013-004.pdf) gives the reference scale of 230 W/W at 4.5 K and 28.5% of Carnot; the linked report's indexed primary passage is the retained evidence. The six deterministic cases in `results.json` span both 85/200 K standby and those three efficiencies. They are scenario calculations, not measurements of the proposed equipment.

## Optional regional station

A station serving 1 km² of receiver has 62.48 m² of active H getter, 48.67 m² of active He panel, and 97.33 m² of installed He panels across two banks. The active 5 K view load is 4.415 kW. The unphased handoff requires 17.807 kW at the end of cooldown; associated maximum electrical input is 2.315 MW in this thermal cycle. Connected subdivisions and shared cryogenic headers can stagger bank transitions physically. This offers an industrial refrigeration topology while retaining the local 0.5 m³/m², 3.15 kg/m² hot buffer.

The region's mean warm throughput is 1.25235 Pa·m³/s at 1000 K. An illustrative 500 m-long, 1 m-radius warm header has molecular-model conductances 13.576 m³/s for H₂ and 9.634 for He, giving a summed mean pressure drop of 0.10361 Pa. With a selected collision area of 1e−19 m² at 1 Pa, its mean free path is 0.0976 m and diameter-based Knudsen number is 0.0488. This lies in the transitional regime: the molecular result sets a useful comparison, while a multicomponent transitional-flow calculation must establish the actual installed conductance.

A 1 mm wall at 2000 kg/m³ contributes only 0.00628 kg/m² receiver as bare pipe skin; stiffness, supports, bends, joints, insulation and service access require their own mass. Thermal placement matters more than that skin estimate: a bare header with full cold view would radiate about **178 W/m² receiver**. It must lie inside a suitable warm enclosure or have charged insulation. The header volume is only 0.00157 m³ per receiver m² and does not replace the local buffer.

Large warm headers end upstream of the shielded cold intake. The existing 1 mm-radius terminal cold-conduit equivalent is retained, either distributed or as parallel channels at the regional station. No metre-scale open cold aperture is introduced. Regional refrigerator, header and shield mass remain installed-component requirements; the existing generic service allocation alone does not demonstrate their fit.

## Same-cohort material closure

The 100%-incident warm-buffer response screens maximum throughput. Its separate global-cleanup capacity increment cannot simultaneously represent additional physical atoms from that same incident flow. For incident flow `I`, cavity branch `J=f_packet I`, direct local branch `L=(1−f_packet)I`, cavity inventory `G`, buffer inventory `B`, later cavity capture `C`, and product outlet `O`,

\[
\dot G=J-C,\qquad \dot B=L+C-O,\qquad
\dot G+\dot B=I-O.
\]

Thus steady `O=I`, with any external supply, escape or changing stock stated separately. If a solid cohort holds atoms, replace the immediate local release with `dS/dt=L−R` and `dB/dt=R+C−O`. Regeneration and delayed release transfer atoms between those reservoirs.

`lineage.py` reconstructs total incident arrivals on the same 4,096-node quadrature as the frozen routed cavity source, then integrates exponential cavity capture exactly between the saved 150 s source impulses. At `gH=gHe=0.01`, over the 31,968,150 s archived interval:

| Atomic mass ledger | H, kg | He, kg |
|---|---:|---:|
| Incident `∫I` | 1.03966482047e14 | 8.82031604519e13 |
| Cavity injection `∫J` | 1.07093858984e11 | 7.72033618603e10 |
| Later cavity capture `∫C` | 1.06606450046e11 | 7.72121447690e10 |
| Cavity stock change `ΔG` | +4.87408937333e8 | −8.78290841772e6 |

`∫L+∫C+ΔG=∫I` closes within 7e−14 of incident mass. The He return includes a small withdrawal from the initial steady cavity stock. The check also covers the mechanical cleanup's distinct He capture coefficient and a 0.1/0.1 capture comparison. It supplies a same-cohort numerical identity without rerunning the hot buffer on a newly partitioned physical source history. All source fractions retain the ideal mobile-routing scope; prompt reflection, permeability, chemical gases, traps and later release remain separate physical populations.

## Reproduction and provenance

Run `calculate.py`, `regional.py`, then `lineage.py` from this directory using the shared scientific Python runtime; `lineage.py` reads the frozen parent gas arrays by relative path. `validate.py` checks the callable cycle, saved closures and source hashes. `results.json` holds pressure/cycle results; `regional.json` holds the optional topology; `lineage.json` holds exact cohort integrals. The [separate passive replay](../../warm-mechanical-cleanup/cold-interface-correction/README.md) provides corrected transient pressures. The earlier gas, buffer and source-routing packets remain unchanged.
