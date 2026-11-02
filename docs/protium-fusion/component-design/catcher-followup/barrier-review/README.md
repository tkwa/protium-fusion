# ✨ Pressure and barrier review of the shallow graphite collector

**The physically better-supported configuration is a continuously evacuated collection layer feeding isolated downstream storage.** Holding an entire hot pulse at kilopascal pressure directly behind the irradiated cap makes the permissible connected leakage area roughly six orders of magnitude smaller than in the pumped comparison below. That is a design distinction, not a demonstration that the cap must fail. The unresolved quantity is its evolving gas transmission under irradiation, including solid hydrogen permeation and connected defects.

This independent review uses the geometry team's current variables:150 µm aligned active graphite,100 nm continuous caps over100 nm-wide vertical slits at20 µm pitch, a50–100 µm rear plenum, and separate gas conduits. The cap-mouth area is0.005 of receiver area. The full350 kg/m² graphite inventory remains. `input-snapshot.json` pins the exact delivered-B1/C13H2 source, cap damage proxies, all-charged-ion peak fluxes and geometry revision. No capture efficiency or material healing rate is assigned.

## Thermal gas flux and equivalent leakage area

For Maxwellian gas atT, mean speed is `vbar=sqrt(8kT/pi m)` and one-sided wall flux is `J=P/sqrt(2pi m kT)=n vbar/4`. A molecular aperture of effective open areaφ per receiver area gives conductance `Cfront=φ vbar/4` per receiver area. These follow kinetic theory; finite apertures require a transmission factor depending on geometry and wall interactions, as treated in [NASA/TM-2011-216970](https://ntrs.nasa.gov/citations/20110009985).

Hereφ includes geometric open area times molecular transmission. It is an **equivalent effusive area**, not necessarily a measured pore-area fraction. A long narrow pore passes less gas than an equal-area zero-thickness opening. For gas-permeating solid material, an equivalentφ can be defined at a particular pressure and temperature, but need not remain constant as conditions change. Molecular flow is appropriate for sufficiently small connected defects; large apertures at higher pressure require transitional/viscous conductance instead.

The exact mean returned H rate is1.26816×10¹⁴ nuclei/(m²s). If transported asH₂, its molecular rate is half that. For an illustrative allowable return to the cavity of0.001 of this H₂ inventory, a *maintained H₂ partial pressure* gives:

| H₂ pressure at1000 K, Pa | Thermal molecular flux, m⁻²s⁻¹ | Permissible equivalentφ |
|---:|---:|---:|
| 0.001 | 5.8686×10¹⁹ | 1.0805×10⁻⁹ |
| 1 | 5.8686×10²² | 1.0805×10⁻¹² |
| 1000 | 5.8686×10²⁵ | 1.0805×10⁻¹⁵ |
| 100000 | 5.8686×10²⁷ | 1.0805×10⁻¹⁷ |

At1000 Pa the last-but-one area corresponds to about344 unit-transmission1 nm-radius openings per receiver square metre. At100000 Pa it is3.44 such openings. A cap confined to0.005 of receiver area must have its area-averaged transmission below2.16×10⁻¹³ or2.16×10⁻¹⁵ respectively. These are acceptance curves for the chosen pressure and inventory budget, not predictions of the number of holes produced.

AtomicH andHe rows are in `results.json`. They are distinct gas scenarios, not extra simultaneous inventory. The reference comparison usesH₂+He; atomicH changes thermal speed and gas-particle count. Species conversion and recombination remain part of gas processing.

## Closed storage versus continuous extraction

If an entire pulse'sH₂+He is stored at1000 K, the pressure follows its actual inventory:

| Stage | Pressure in50 µm void, Pa | Pressure in100 µm void, Pa |
|---|---:|---:|
| A | 7993 | 3996 |
| B1 | 101.8 | 50.9 |
| B2 | 1725 | 862 |
| C13 | 2522 | 1261 |

These are full-pulse **void-volume sensitivities**. The final50 µm geometric plenum has25% posts, leaving37.5 µm local void-equivalent depth; A would reach10.66 kPa in that volume alone. Its slits add0.75 µm, while the longer gas-via network adds much more connected volume. If all model voids were isolated together and equilibrated at1000 K, their1.424 mm volume-equivalent depth would lower the A full-pulse pressure to~281 Pa. The actual valve location, conduit temperatures and equilibration determine which inventory/volume relation applies. Retained material lowers the gas inventory; overlapping gas from earlier pulses raises it. They are not the pressure during continuous removal.

For an isolated isothermal plenum of depthh with one front escape path, surviving gas is `exp[−φ vbar t/(4h)]`. In50 µm ofH₂ void, keeping front loss below0.001 requiresφ≤6.17×10⁻¹¹ for a1 s hold,6.17×10⁻¹³ for100 s, and8.02×10⁻¹⁵ for7698 s. Thus “open the extraction valve after the pulse” is a quantitative residence-time choice. An initially damaged cap may leak during that hold even if it is repaired before the next shot.

The current geometry estimates channel-plus-via conductanceβchannel=7.7393×10⁻⁶, where conductance isβvbar/4. The gas agent separately includes downstream pump speed. With its comparison valueSp=0.001 m³/(s·m² receiver), series combination gives:

| Species | Effective pumping speed, m³/(s·m²) | Equivalentβnet | φ permitting0.001 plenum return |
|---|---:|---:|---:|
| H₂ | 0.00086246 | 1.0644×10⁻⁶ | 1.0655×10⁻⁹ |
| He | 0.00081651 | 1.4201×10⁻⁶ | 1.4215×10⁻⁹ |

For linear molecular conductance into a negligible-pressure sink, eventual front fraction is exactly `φ/(φ+βnet)`. It is independent of pulse inventory and plenum volume; those determine pressure and emptying time. The mixedH₂+He mean pressure is0.001477 Pa. If all gas emerges promptly following the source arrival rates, the sums of individual species peak-pressure contributions are4.580,0.0580,3.270 and2.340 Pa forA/B1/B2/C13. These are upper sums for that immediate-release model, not limits on a later concentrated thermal-desorption burst.

The final geometry also places200 µm stopping pads over doglegged via entrances (0.00628 kg/m²) to block radial ion paths. Its conductance estimate retains a straight-channel approximation and omits some plenum, dogleg, bend and entrance resistance. Finite pump ultimate pressure and regeneration backflow also need to be included in the realized system. The earlierβ=3.35×10⁻⁶ slit-only discussion is preserved in the input as a sensitivity; its numbers are not silently substituted for the current finite-pump case.

The reference gas temperature here is1000 K. The latest scheduled material histories peak around1494 K in the150 µm skin; the1668 K conservative envelope must not be held for hours. At fixed stored inventory, molecular leak speed increases only as√T (about22% from1000 to1494 K). Activated solid permeation can change much more rapidly, so it requires the actual material-temperature history.

## Budget direct front release separately

Material that exits the front before reaching a collection channel has already entered the cavity. If that fraction isd and the subsequent plenum return fraction isb, total cavity admission is `d+(1−d)b`. For example, the geometry team's preliminary B1He100 nm cap-only floord=0.000662 would leaveb≤0.0003382 under a hypothetical total0.001 budget. With the finite pump above, theHe equivalentφ allowance would then be about4.80×10⁻¹⁰, before any additional direct losses. Parent routing calculations determine the fulld; the gas agent determines acceptable total admission from pellet and chamber limits. **The0.001 value in this review is a comparison point, not a universal plant requirement.**

## Radiation-created paths and repair timing

The all-charged-nucleus sampled peak fluxes are8.05585×10¹⁷,6.61324×10¹⁵,6.37869×10¹⁷ and5.75863×10¹⁷ m⁻²s⁻¹ forA/B1/B2/C13. The independent32→64 quadrature change is below7.1×10⁻⁵. These resolve the within-cell ballistic model; they are not a bound on unresolved physical substructure. The all-species flux interface avoids a previous draft's mistaken use of theH/He-only dictionary for an all-ion damage rate.

Definey as the number of through-connected gas paths produced per ion hitting a cap,Ah as each path's effective transmitting area, andτ as its mean open lifetime. For exponentially repaired defects, `dφ/dt = mouth_fraction × y × Ah × ion_flux(t) − φ/τ`. With a bounded source, the additionalφ is bounded by both the cumulative created area and the peak production rate timesτ. The initial defect/seam/permeation budget must be added separately.

For1 nm-radius **unit-transmission equivalent paths**, theH₂/He reference allowanceφ=1.0655×10⁻⁹ gives these conditional requirements forA:

- No within-shot healing: y≤8.05×10⁻¹⁰ per cap-hit ion.
- Peak-rate condition: yτ≤8.42×10⁻⁸ s. Ify were10⁻⁶, this would requireτ≤0.0842 s; ify were10⁻⁹,τ≤84.2 s. These are sensitivity examples, not assumed yields or repair capability.
- The annual-mean condition is much weaker and cannot replace the pulse condition.

A straight100 nm-long,1 nm-radius cylindrical pore has approximate diffuse molecular transmission8r/(3L)=0.0267 in the long-tube limit; its allowedy orτ is then37.5 times larger for the same actual radius. Atomic-scale chemistry, specularity and a crack's shape can change this, so the code reports the unit-transmission case and the geometric transmission scenario separately.

The cap's25 eV-threshold NRT displacement proxies per pulse are0.0295,0.000319,0.0316 and0.0197 dpa, totaling3.52 proxy dpa/year. These count collision damage opportunities. They do not count surviving vacancies, aligned holes through~300 graphitic layers, cracks or permeable seams. A connected1 nm-radius cylinder through100 nm contains roughly3×10⁴ carbon atoms; one displacement event cannot be equated with such an aperture. Conversely, collective cracking can connect existing damage without removing a cylinder atom by atom. An actual track/cascade/annealing model must establishy and its history. Neither inevitable puncture nor automatic perfect healing follows from the present proxy.

Pressure stress across a supported slit scales asP(width/thickness)². Here width and cap thickness are both100 nm, so even10⁵ Pa gives a stress scale around10⁵ Pa. The geometry's support and thermal-expansion calculation still applies, but this pressure load alone does not create a mechanical contradiction. Gas transmission is the sharper unresolved property.

## Hydrogen can permeate without an open hole

[Sun et al.2020](https://doi.org/10.1038/s41586-020-2070-x) measuredH transmission through monolayer graphene with approximately√P pressure dependence and1.0±0.1 eV activation energy over295–323 K. Its bilayer control was impermeable within the experiment's sensitivity at323 K. These results identify a chemical transport mechanism; they do not supply a1000 K,100 nm irradiated-laminate permeability.

A measured or calculated laminate lawJH(P,T) must therefore be charged separately from connected molecular apertures. At the reference meanH₂ pressure, the0.001 hydrogen-inventory budget corresponds to1.2682×10¹¹ H nuclei/(m² receiver·s). That is a direct acceptance scale for a whole-receiver equivalent permeability measurement. Conversion to per-cap-area flux divides by0.005 only when the caps are demonstrably the sole transmitting surface.

As an explicit pressure-law sensitivity, suppose pump throughput is proportional toP and solid leakage to√P. For one isolated gas cohort, definez as pump/leak throughput ratio at its initial pressure. Integrated front fraction is `2[z−ln(1+z)]/z²`. A ratio1000 gives0.001986 loss;2000 gives0.000996. Thus the constant-conductance fraction formula must not be used unchanged for a chemically permeating barrier. This mathematical result assigns no permeability coefficient to the candidate material.

## Useful design choices

1. **Pump during collection, then isolate compressed stock downstream.** This reduces pressure-time exposure at the irradiated cap. The high-pressure reservoir and its valves can be behind the backing and outside direct ion irradiation. Source arrival spans and fast mobile release make pumping during the pulse more useful than relying solely on a later valve opening. Delayed trapped release can continue into the same low-pressure collection network.
2. **Add a pumped guard space if one barrier is insufficient.** For linear conductancesC₁/C₂ and guard pumpingS, transmitted conductance isC₁C₂/(C₁+C₂+S). Two equal unpumped barriers merely halve conductance; their transmissions cannot automatically be multiplied. Pumping the intervening space produces the stronger suppression while returning intercepted gas to process stock.
3. **Separate sacrificial surface repair from the storage pressure boundary.** A graded or replaceable aligned surface may control implantation and channel routing, while protected downstream barriers control compressed storage. Extra front layers also change implanted-gas direct release, so their benefit must be evaluated with the routing calculation.
4. **Use hot-isolated cycling only for a specified retained population.** Holding strongly boundH and releasing it later into a closed extraction volume is possible in principle. MobileHe and readily releasedH cannot simply be assigned the same hold behavior. Cooling, valve timing and trapping must close the actual inventory balance.

Under the user's molecular-manufacturing premise, repair energy is not automatically a reason to reject the cap. Reconstructing100 nm only over the0.005 slit mouths after every event handles1 mg/m²/event and uses0.000132 W/m² at1 MJ/kg input. Reconstructing the entire100 nm face gives0.02645 W/m² at that same declared input. Neither is a required full-layer replacement or an atomization charge; targeted repair may handle less. The unknown is whether the morphology can be restored fast enough while meeting its gas-transmission budget, not an invented efficiency penalty.

The current calculations identify what must be demonstrated without assigning a failure probability. Continuous low-pressure collection with isolated storage has a plausible physical basis. Its quantitative performance still depends on the hot damaged laminate's gas transmission, direct-release routing, actual downstream conductance and repair timing. Present evidence does not justify a narrow80% credible interval for those material properties.

## Reproduction and checks

`python3 calculate.py` regenerates `results.json`; `python3 validate.py` checks Maxwellian flux identities, finite-hold limits, pressure/inventory closure, conductance partition, pumped barriers, the nonlinear pressure-law integral and source alignment. All89 checks pass. No hydrodynamic reruns or modifications to the geometry, gas or archived research packets were made.
