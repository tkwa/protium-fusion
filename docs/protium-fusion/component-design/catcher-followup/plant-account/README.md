# ✨ Catcher services in the complete plant account

The preferred warm-buffer interface can preserve approximately **20.7 mW/kg** in the existing fixed-rate ledger, provided its stated hardware and operating requirements are met. A direct cold local port gives **19.3 mW/kg** under the same ambient-cleanup allocation. These are conditional component-integration cases. The graphite transport, complete hot pumps and fuel-capsule source still require qualification.

[calculate.py](calculate.py) produces [results.json](results.json). The isolated [station_catcher.py](station_catcher.py) adds an explicit return of auxiliary dissipation to the hot loop to the archived station solver. No archived file or selected operating point is modified. The coherent delivered-B1/C13H2 source128 family, isotope flows, pulse rate, actual dispatch and 350 kg/m² graphite remain fixed.

## Reference interfaces

The following values are per receiver area:

| Item | Warm-buffer case | Direct cold-port comparison | Evidence/status |
|---|---:|---:|---|
| Shared gas-conduit hot-to-cold heat | 0.01411 W/m² | 5.66089 W/m² at the solved cold temperature | Warm-buffer conduit geometry versus conductance-sized direct cold mouth; hot radiation and solid conduction charged |
| Warm bearing/support heat bypass | 0.5 W/m² | 0.5 W/m² | Component requirement |
| Ambient warm-pump electrical input | 1 W/m² | 1 W/m² | Total motor/drive/standing-loss requirement, including molecular windage |
| Local warm booster input | 0.13 W/m² | 0 | Requirement informed by current pump-area comparisons; high-temperature apparatus remains to design |
| He active sorber area | 4.8667×10⁻⁵ m²/m² | Same | Retained conservative working area/capacity; alternating isolated banks |
| Refrigeration fraction of Carnot | 0.30 | 0.30 | Subjective central estimate; 80% interval 0.15–0.50, conditional on the specified thermal cycle |
| Further downstream processing/control electricity | 0.005 W/m² | 0.005 W/m² | Explicit reserve, including chemical separation work; apparatus still needs an allocation |
| New service hardware allocation | 300 kg/m² | 300 kg/m² | Acceptance allocation within the retained service mass |
| Credited warm motor-heat recovery | 0 | 0 | Conservative reference; 50% and 90% routing sensitivities are calculated separately |

The warm-buffer line is intercepted before the cold collector. Its very small admission mouth does not eliminate the retained 200 K view of the 5 K sorber, adsorption capacity, bank cooling, getter regeneration or gas compression. A fraction 10⁻⁴ of the admitted hot radiation is separately allowed to reach 5 K. The full ideal-gas H₂/He sensible-enthalpy drop from 1000 K to the radiator is also charged as a conservative hot-to-cold term, about 3.2 µW/m². Energetic first reflections remain inside the sphere until recapture; their heat is delayed rather than automatically sent down the cold gas line. The shield view, actual solved radiator temperature and sorber cycle are evaluated together. The active bank captures while the other warms electrically to the radiator temperature, remains isolated, and cools over a specified final 100 s before handoff. Constant areal heat capacity 1000 J/(active m² K) over 5 K to standby is an explicit conservative scenario, with actual Cp(T) still to replace it. Inactive radiation during cooling and finite local peak refrigeration are included through the [two-bank interface](../gas/cold-interface-followup/README.md). The full electrical sensible-warmup charge grants no reduction for passive warming radiation.

The proposed 1000 K buffer has 0.5 m³ of volume per receiver m². The actual H₂/He dispatch gives a full-incident-gas capacity response near 0.680–1.691 Pa after the passive cold-panel correction, with separately qualified source peaks up to 4.58 Pa. The [buffer calculation](../warm-mechanical-cleanup/BUFFER.md) includes finite source/plenum response, reverse molecular flux and precharge inventory. The retained cold panels give a 300 K reservoir near 0.000594 Pa mean and 0.001057 Pa peak, rather than a separately imposed pressure independent of their speed. The full-incident buffer source is a capacity comparison: a physical combined ledger partitions atoms between immediate local return, cavity stock and later recapture. It does not add globally recaptured atoms a second time. The global cleanup foreline stays isolated through a separate small booster. Its benchmark permits up to 10⁻³ Pa **partial pressure for each of H₂ and He**; the corresponding total can be 2×10⁻³ Pa. It is never connected directly to the several-pascal primary plenum.

The hot box's explicit mechanical screen is 3.15 kg/m². A current whole-pump package scale gives about 143.7 kg/m² for global cleanup, with another 35 kg/m² for a conservative particle baffle. The new rework factory requirement is 14.288 kg/m² at the inherited capacity coefficient. The tiny primary booster package comparison is 0.01165 kg/m². Together these named scales use about 196 kg/m² of the 300 kg/m² allocation. Liner, seals, valves, modified hot rotor/stator construction, extra stages, refrigeration, controls and installation must fit the remaining allocation or add mass. This comparison does not certify a bill of materials.

## Refurbishment and retained inventory

The simultaneous reference charges are:

- Entire 100 nm face renewed after every pulse, at 2 MJ/kg input: **0.052894 W/m²**.
- Every nominal E_d=25 eV NRT proxy event serviced at 10 eV electrical input: **0.189408 W/m²**.
- Entire 150 µm active layer reworked monthly while hot, at 2 MJ/kg: **0.228154 W/m²**.
- Returned H/He given a specified 100 eV useful free-energy upgrade per nucleus at 50% efficiency: **0.004938 W/m²**.

These distinct operations deliberately include overlapping maintenance options. A qualified service plan can replace redundant work. They do not assert that the full layer needs monthly replacement or that all NRT events survive. Broad material handling and the one-carbon-atom-per-proxy capacity equivalent together require **14.28798 kg/m²** at the inherited 10⁸ kg/(kg/s) manufacturing-capacity coefficient; this is new explicit machinery inside the service allocation. The much smaller existing generic routine factory remains charged separately.

The 100 eV per-nucleus selective-work comparison for the returned nitrogen stream is about 0.00214 W/m² at 50% efficiency, inside the downstream 0.005 W/m² reserve if assigned there. It supplies a work allowance, not a demonstrated separator. Captured H/N/O-bearing molecules require a chemical sorting path before helium sorption. They do not become fresh-carbon consumption merely by leaving a graphite lattice.

Monthly extraction needs rate control. A fully retained mean month of H/He has 3.293 Pa·m³/m² of gas at 1000 K, enough to raise a 0.5 m³/m² buffer by 6.587 Pa if released instantaneously. That cohort replaces its earlier prompt source. Local valves, sufficiently shared headers or longer release periods must keep the cold banks within their load and cycling capacities. Globally staggered service cannot cure an overloaded isolated local bank without a physical connection.

## Complete balance and specific power

| Case | Net mW/kg | Radiator temperature, K | Source peak-hardware multiplier |
|---|---:|---:|---:|
| Archived allocated reference | 20.7837 | 200.000 | 1 |
| Warm-buffer interface requirements | **20.6576** | **200.189** | 1 |
| Same, 50% more source peak hardware | 19.6066 | 200.189 | 1.5 |
| Same, twice source peak hardware | 18.6573 | 200.190 | 2 |
| Direct cold local port | **19.2729** | **202.226** | 1 |
| Direct cold port, 50% more source peak hardware | 18.2921 | 202.227 | 1.5 |

The source multipliers are hardware sensitivities. Doubling a pellet's radius doubles acceleration length at the same surface field and gradient and halves aggregate current for the same transported mass. It does not prove that source hardware doubles; the fraction that scales with length needs a launcher design. The separate orbital-support guides are not rescaled as fuel launchers.

The nominal 1 cm C13 capsule has a narrow solid-core thermal margin. A 2 cm comparison supports a larger gas-column budget but still needs source preparation, coating and slow-packet treatment. Consequently the first row of the new warm-buffer account is a **plant-interface requirement**, rather than an engineered operating point that establishes unchanged capsule survival. The [gas acceptance calculation](../gas/source-leakage-followup/) retains these dependencies explicitly.

The inherited ordinary auxiliary allowance is 10% of the gross process margin. Named existing and new loads consume that allowance first. New services are not debited again while they fit. The warm reference leaves about **1.00 W/m²** of auxiliary allowance and **378.606 kg/m²** of the previously unassigned service mass after allocating the new 300 kg/m². Additional loads raise the aggregate auxiliary fraction; additional mass beyond the retained allowance enters the support, spare and manufacturing feedback. Increasing selected repair work from 10 to 100 eV per proxy event lowers the warm case to **20.33 mW/kg**.

The radiator initially uses the full stipulated 90% of shell area. Any new net cold load therefore triggers a complete finite-exchanger solve at a slightly higher radiator temperature, with the same dispatch. Existing guide cryogenic duties are rescaled with that temperature. The source rate is not silently reduced. Net energy, all cold heat components, conversion equipment, radiator, support, reserve, spares and routine replacement close numerically for every reported complete row.

Internal warm dissipation is electrically charged once, added to the engine's hot input when a return fraction is selected, and subtracted from direct auxiliary cold rejection. The model retains the global nuclear energy balance. At 60% heat-to-electric efficiency, returning warm dissipation can reduce its net electrical penalty, while changing converter and heat-exchanger duties. The reference grants no such credit.

## A performance frontier, rather than one mandatory idle target

For the warm-buffer case and other reference parameters, keeping net output above **20, 10 or 1 mW/kg** allows total ambient-pump electrical input up to approximately **3.4, 24.9 or 44.3 W/m²**, respectively, with no warm motor-heat credit. The corresponding solved radiator temperatures are about 201.16, 214.46 and 224.62 K. A selected 1 W/m² input is therefore a useful high-performance target, rather than a universal viability threshold.

The small measured bearing-pump experiment scales to approximately 221 W/m² receiver at the chosen global aperture. Reaching the respective output thresholds would require roughly 65-fold, 8.9-fold or 5.0-fold lower areal input than that small-module comparison, or a changed architecture. Larger modules, reduced standing losses, colder compression stages with carefully limited heat paths, lower forepressure and useful hot dissipation are design variables. Neither extrapolating the measured small-module value unchanged nor assigning a future low-loss value by assertion constitutes a completed pump.

The machine-readable file also gives explicit 50%/90% warm motor-heat return, service-mass, source-hardware, repair-work, heat-leak and refrigerator sensitivities. A 0.52 W/m² local-booster allocation, representing a fourfold apparatus-input sensitivity for the narrower-pitch collection variant, still fits the retained auxiliary allowance. That row does not qualify the changed material or compressor geometry. No joint probability distribution or total-reactor credible interval is assigned. The refrigerator interval is a conditional engineering judgment based on primary industrial helium-refrigeration evidence; the other rows are deterministic performance requirements.

## Reproduction

Use the project's Python 3.10 research environment with NumPy, SciPy and pynucastro available, then run `python calculate.py` from this directory. Baseline reproduction and mass, electrical, nuclear and cold-heat closures are checked in the script. Input and code hashes are saved in the output. The [independent accounting review](INDEPENDENT-REVIEW.md) records the reconstructed balances and input hashes. The small full-bank-cycle update is reviewed as a delta after the main solver audit.
