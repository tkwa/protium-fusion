# ✨ B1 source entropy, hydrogen ratio and N15-rate robustness

A controlled contact-entropy profile removes the low-N15-rate loss of ignition in the modeled H/N=2 pulse. The alternative H/N=1.5 source also burns strongly. The result supports comparing **source preparation and mixture together**, rather than requiring a recipe change from the original prescribed-flow sensitivity alone.

All six calculations use 128 radial cells, the full 27-species network, FastEOS, compression-cap viscosity, 0.001 m²/kg grey opacity, equal total prepared energy (3.19501678×10²³ J) and the original B1 physical duration. Each recipe has exactly the same initial profile for its carbon-only, low-N15 and high-N15 controls. The carbon-only label means the two median LUNA carbon-capture updates with the original N15 rate; the low/high cases additionally substitute the paired NACRE-II N15(p,alpha) evaluation limits. These limits have no probability level assigned here.

| Controlled source | Prompt nuclear/prepared gain | C12 per initial N15 | Peak temperature |
|---|---:|---:|---:|
| H/N=2, median carbon update | 8.082346 | 0.809695 | 542.89 MK |
| H/N=2, same source + low N15 | 8.259892 | 0.827415 | 575.42 MK |
| H/N=2, same source + high N15 | 7.901018 | 0.791599 | 534.09 MK |
| H/N=1.5, median carbon update | 8.013658 | 0.779284 | 598.37 MK |
| H/N=1.5, same source + low N15 | 8.176408 | 0.795053 | 646.80 MK |
| H/N=1.5, same source + high N15 | 7.847790 | 0.763210 | 551.31 MK |

All complete their requested physical durations. The 128-cell prescribed-flow H/N=2 low-N15 control elsewhere yielded Q=6.21285 and reached the selected station occupancy limit. That result remains valid for its own profile. The source adjustment changes the coupled compression/heating history enough to retain strong burn under the same low rate. Lower reaction rate increases integrated gain in these adjusted source examples; a one-parameter rate multiplier cannot predict their full hydrodynamic response.

`source-variants.json` names the exact unchanged records, profiles, trace files, durations, initial isotope abundances and hashes. Its trace audits use the stored EOS/energy states, with reaction-flow reconstruction disabled because the runtime network differs from the archived original rates. Input free-text scopes inherited from the contact construction describe the 64-cell preparation step; the actual result mesh and rate revision are the explicit top-level metadata. Whole-plant comparisons are supplied separately by the nuclear-rate integration packet, preserving the other three stages as named controls.

## Matched plant consequence

The independent four-row integration in `mixture-design/rate-uncertainty/b1source-integrated.json` holds A H3, B2 H2 and C13 H4 at their nominal128 median-carbon states. It changes only the named B1 source/recipe/rate endpoint, retaining the common 99% capture assumption. All rows pass 64 integration checks and are radiator-limited; main-flight occupancy lies between 0.720 and 0.747.

| B1 source control | 5000 kg/m² receiver | Separate 2000 kg/m² case |
|---|---:|---:|
| H2, carbon update | 6.626477 mW/kg | 12.483871 mW/kg |
| H2, low N15 | 6.640488 mW/kg | 12.510255 mW/kg |
| H1.5, carbon update | 6.620868 mW/kg | 12.473309 mW/kg |
| H1.5, low N15 | 6.633942 mW/kg | 12.497927 mW/kg |

H2 narrowly leads in both matched source comparisons. The present rate evidence therefore does not compel the additional cryogenic-nitrogen feed of H1.5. The small power difference leaves source simplicity, repeatability and later multidimensional results relevant to the choice. These rows are a controlled B1-only source/rate sensitivity, separate from the coherent four-source128 family and from qualification of the thinner receiver.

## Completed envelope with all four source profiles

`source128-envelope.json` instead fixes the other three stages to their **actual source128** states. This coherent source-profile comparison is the relevant follow-on to the four-row nominal-other-stage control above. The original REACLIB N15 curve is retained as the control; no NACRE-II adopted-source run is implied.

| B1 source | N15 control | NACRE-II low | NACRE-II high |
|---|---:|---:|---:|
| H2, 5000 kg/m² | 6.617345 | 6.631337 | 6.602431 |
| H1.5, 5000 kg/m² | 6.611747 | 6.624802 | 6.597908 |
| H2, separate 2000 kg/m² | 12.466631 | 12.492978 | 12.438548 |
| H1.5, separate 2000 kg/m² | 12.456089 | 12.480674 | 12.430031 |

Values are mW/kg under the same conditional plant assumptions. Every row remains radiator-limited; maximum main-flight occupancy is 0.76848 at 5000 kg/m² and 0.76900 at 2000 kg/m². H2 narrowly leads H1.5 for every rate curve. The two high-curve trajectories pass 28 trace checks, the original four pass 60 trace/profile checks, and the integrated six-row envelope passes 82 accounting/replay checks. This is a B1-local rate sensitivity with the other stages held fixed; the curves are not probabilistic credible-interval endpoints. No profile was retuned for high or low rates.

## The source adjustment

Both recipes use the same design rule: choose one contact-shock entropy near the nominal hot profile at 90% of its temperature, then map that entropy to each parcel's original hot density. This produces spatial temperature ratios roughly 0.87–0.92, rather than an artificial uniform multiplier. Radial kinetic energy is adjusted to preserve total prepared energy. The collision heat comes from controlled normal components of the incoming streams' stored kinetic energy; there is no new external flash-heater power debit.

For H/N=2 the contact conversion is 1.23943×10¹¹ J/kg, 4.3837% of the prepared pulse energy. For H/N=1.5 it is 1.13862×10¹¹ J/kg, 4.0272%. The latter gives contact density 4211.7 kg/m³, temperature 12.449 MK and pressure 2.9863×10¹⁴ Pa in FEOS. Each of a symmetric stream pair has normal speed 477.2 km/s; conditional half-angle mass quantiles are 8.97°, 12.18° and 17.31°. About 17.6% of the pulse mass rarefies after this contact state; the rest compresses. These angles use a stated local parcel-energy allocation and do not specify remote launch aiming.

The H2 profile was selected before the low-N15 sensitivity revealed its prescribed-flow weakness. The H1.5 profile applies the same contact-entropy rule. Neither profile is retuned between its two nuclear-rate calculations. They are constructive model-level source candidates, with three-dimensional contact, pressure work, wave propagation and packet homogenization still to be evolved.

## A physical H/N=1.5 feed

Per catalyst nitrogen atom, 0.375 hydrazine molecules plus 0.125 N2 molecules supply 1.5 hydrogen atoms. With N15 these are **77.27% hydrazine and 22.73% nitrogen by payload mass**. An NH3/N2 alternative would require 45.45% of its payload mass as N2. The hydrazine-rich option therefore reduces the separate cryogenic-nitrogen burden. Neither mixture implies a pure N2H3 material.

Two physical preparation candidates are a sealed heterogeneous solid composite near 50 K, or separately conditioned hydrazine and N2 packet populations mixed during contact. The second leaves the hydrazine warmer but adds a larger composition-homogenization task. The first cools the whole pellet and must manage heterogeneous thermal expansion, weak cryogenic inclusions, capsule pressure and surface-charge stress. The FEOS reference at 1000 kg/m³ and 250 K is an effective homogeneous contact surrogate; it does not establish the molecular starting phases or their actual Hugoniot.

Nitrogen's triple point is 63.14 K and its critical temperature about 126.19 K. It therefore needs a cryogenic or explicitly compressed storage path. Natural hydrazine melts near 274.63 K; its listed liquid density at 298 K is 1004 kg/m³. These phase properties exclude treating both constituents as an ordinary ambient homogeneous solid. [NIST nitrogen phase data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=4), [NASA property table](https://ntrs.nasa.gov/api/citations/20120003055/downloads/20120003055.pdf).

The selected 1000 kg/m³ bulk-density target is an engineering input. Illustrative constituent densities 1062.5 and 1000 kg/m³ give additive-volume density 1047.6 kg/m³, leaving a few percent for void space or dimensional adjustment. Measured confined nitrogen cryocrystals have molar volumes of order 27–29 cm³/mol, supporting the density scale, but their stated pressures and temperatures must not be relabeled as a measured zero-pressure N15 composite. [Konstantinov et al., 2005](https://nasplib.isofts.kiev.ua/bitstreams/6dfa3e52-8cf2-48d1-9002-389dd9cdc4cc/download).

Using the existing 1000 K cavity / smooth-aluminium absorptivity reference, a 1 cm-radius, 1000 kg/m³ packet travelling at 2.4 Mm/s for the full 1.1×10⁹ m receives about 170 kJ/kg. For a **chosen** 30 kJ/kg nitrogen warming allowance, the alternatives are Planck absorptivity below 0.00385, exposed flight below 1.94×10⁸ m, or packet radius about 5.7 cm with the same coating. These are conditional thermal design comparisons, not a universal pellet-size requirement. More warming may be acceptable with pressure/headspace and trajectory changes explicitly modeled. `recipe-screen.json` preserves the formulas and assumptions. It also distinguishes bare fuel fractions from capsule mass, which must be recovered or included in the pulse's impurity/energy ledger.

Hydrazine decomposition to N2 and H2 has a thermochemical scale about 1.14 MJ/kg of this mixture, roughly 10⁻⁵ of the modeled contact heat. Chemical energy cannot supply or substantially replace the required warm-contact energy. Its timing can still affect fragmentation and phase space during delivery, so small total chemical energy does not prove packet stability.

The next physical discriminator remains the contact-to-hot trajectory with actual density, chemistry and finite packet geometry. The reported gains establish tolerance of these specified hot entropy profiles and rate controls. They do not yet establish their complete cold-source realization or reactor lifetime.
