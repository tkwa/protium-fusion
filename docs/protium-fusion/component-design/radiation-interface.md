# ✨ Radiation at the fuel and receiver interfaces

Checkpoint, 2026-10-31. New component-design work; the earlier receiver packets remain frozen. Completed screens cover incoming infrared heating, receiver photon deposition, equilibrium vapor, return-gas columns, opacity and initial Coulomb energy scales.

## Cold packets see the cavity radiation field

The reference receiver has radius 1.1×10⁹ m and a 1000 K hot surface. An approximately closed, uniform-temperature cavity tends to blackbody spectral radiance even when its wall has low emissivity: reflection preserves radiation while emission and absorption bring it toward equilibrium. Reducing wall emissivity slows equilibration; it does not by itself establish a colder steady radiation bath. Actual cold fronts or infrared escape apertures require an energy balance. This is separate from the accepted shading of direct sunlight outside the reactor.

For a spherical grain of radius a and density ρ moving through this field, the absorbed power per mass is 3 ασT⁴/(ρa). The Planck- and angle-weighted absorptivity α refers to the **1000 K incident spectrum**. At v=2.4×10⁶ m/s, the full radial flight takes 458.3 s. For a=1 mm and ρ=1000 kg/m³, the black-absorber energy exposure is 77.97 MJ/kg. Its energy acceptance is therefore stringent even though the grain’s total mass is large compared with individual atoms.

[Chang’s primary polyethylene calorimetry](https://pmc.ncbi.nlm.nih.gov/articles/PMC6565414/) gives 269–302 kJ/kg to warm the two measured polyethylene samples from 20 to 300 K. The corresponding full-flight acceptance is α≤0.00345–0.00387. These two values are material cases, not a credible interval. The calculation treats whole-grain enthalpy as the acceptance test; surface temperature, differential thermal expansion, charge stability and loss of volatile components can fail sooner.

A reproducible optical reference uses [Rakić’s 1995 aluminum table](https://doi.org/10.1364/AO.34.004755), with its publicly displayed [complex-index data](https://opg.optica.org/ao/abstract.cfm?uri=ao-34-22-4755). Fresnel absorption is averaged over the isotropic incidence distribution 2cosθ sinθ dθ and then over the 1000 K Planck spectrum. Smooth, unoxidized, optically thick aluminum gives α=0.02182. This reference uses the source’s ordinary-temperature optical data; cold, stressed, charged, coated grains have not been qualified by that measurement.

| Grain radius | Required α for polyethylene 20→300 K | Al-reference absorbed dose | Time to 300 K, including self-emission | Maximum unshielded distance |
|---|---:|---:|---:|---:|
| 0.1 mm | 0.000345–0.000387 | 17.01 MJ/kg | 7.26–8.14 s | 17.4–19.5 million m |
| 1 mm | 0.00345–0.00387 | 1.701 MJ/kg | 72.6–81.4 s | 174–195 million m |
| 10 mm | 0.0345–0.0387 | 0.1701 MJ/kg | 726–814 s | 1.74–1.95 billion m |

The calculation therefore supports larger coated polymer grains or protected injection paths as concrete candidates. It does not qualify a full 458 s transit for a 1 mm grain with this optical reference. Cryogenic reductions in metallic absorption, multilayer reflectors, shorter free-flight paths and a genuinely cooler inner radiation field remain design variables. Their costs belong in the delivery apparatus.

A 100 nm aluminum coat adds approximately 0.081% to the mass of a 1 mm-radius, 1000 kg/m³ grain. In the tabulated optical band its normal-path bulk transmission screen is below 2×10⁻⁵ of the incident Planck energy, before reflection is credited. This supports using the optically thick reference as an initial energy estimate; oxide layers, finite-film interference, defects and core penetration still require a film calculation. The data range omits only a small, explicitly bounded Planck tail.

The coating’s temperature must be treated separately from bulk enthalpy. A thermally bonded coat can conduct into the colder grain. A thermally isolated exterior foil tends toward the 1000 K radiation temperature, above aluminum’s melting point. A supported liquid layer, a refractory exterior, or colder nested reflective layers is a physical design choice. Passing the bulk-dose test does not establish coating survival or uniform grain temperature.

The N:H=1 feed needs an actual stable material recipe and sealed-packet design before assigning a temperature limit. A numerical 30 kJ/kg heat budget is retained as a parametric example in the results: at 1 mm it would require α≤0.000385. It is not a measured enthalpy or phase-survival allowance for an unspecified nitrogen mixture. Bare volatile grains also require a vapor-pressure check.

## Evidence and reproduction

`packet_radiation.py` computes the angular and spectral integrals, enthalpy-limited flight times, radius scaling and coating mass. `packet-radiation.json` includes the transcribed primary optical values and source URLs. Run it with NumPy installed; `--output` can direct results to scratch. Six algebra, normalization and bound checks pass. The optical reference and calorimetry are real primary evidence. The cavity, grain dimensions, speed and allowable final temperature are engineering assumptions. This calculation adds no new hydro or fuel-yield result.

The receiver comparison below uses the saved corrected grey histories with conditional Planck spectral shapes. It preserves prompt and slow-particle timing and converts vapor sources into finite local-capture acceptance budgets.

## Metal and oxide collection surfaces

The added `receiver_radiation.py` screen uses the saved corrected 64-zone histories. Each escaped-grey-energy increment receives the earlier conditional Planck shape at its outer-cell temperature. The [NIST aluminum](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html), [tin](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z50.html) and oxygen mass coefficients then determine energy transfer versus depth. Mixture coefficients are mass weighted. Bulk μ_en represents local charged-particle energy retention; it does not resolve secondary-electron escape, fluorescent transport, or multiple scattering.

The sum ∫F_E μ_en(E)/c_p dE is a conservative no-conduction **surface sensible-heat equivalent** for the resolved 1–1000 keV spectrum. Temperature-dependent phase change occurs before any physically impossible extrapolated temperature. For aluminum, c_p=1177 J/kg/K follows the JANAF liquid table; alpha-alumina uses its 1000 K value, 1224 J/kg/K. Tin’s 250 J/kg/K remains an explicit screening parameter.

| Saved pulse | Aluminum surface equivalent | Alumina surface equivalent | Tin surface equivalent |
|---|---:|---:|---:|
| A | 28 K | 19 K | 376 K |
| B1 | 136 K | 96 K | 1200 K |
| B2, corrected64 | 40 K | 27 K | 523 K |
| C13 | 134 K | 90 K | 1626 K |

Depositing the unresolved sub-keV Planck tail entirely at the surface and convolving its original time bins with conduction adds only 0.97 K for C13 alumina and 1.65 K for B1 alumina, using k=10 W/m/K and ρ=3900 kg/m³. Halving/doubling the assigned color temperature is retained as a separate sensitivity. These results favor the lower-Z Al/Al₂O₃ interface over tin for flash tolerance. They are conditional spectral screens, not a frequency-dependent escaping-radiation calculation. The massive underlying receiver can also heat the layer; its slow history and thermal contact remain separate terms.

Pure molten aluminum presents a different problem even before a pulse. [JANAF Al(g)](https://janaf.nist.gov/tables/Al-005.html) gives p_sat=5.10×10⁻⁶ Pa at 1000 K. With unit evaporation coefficient the gross flux is 3.67×10⁻⁹ kg/m²/s, or 5.58×10¹⁰ kg/s across the reference sphere. This is circulation rather than irreversible consumption. A saturated 1000 K aluminum atmosphere would nevertheless have radial column 0.0182 kg/m² and would disrupt cold fuel delivery.

For a grain moving at v, a conservative bulk heating screen is Δu=3Σv²/(8ρa), where Σ is encountered gas mass per projected area and all incident kinetic energy is retained. Allocating 10 kJ/kg to gas heating gives Σ≤4.63×10⁻⁸ kg/m² for a 1 cm-radius grain, or 4.63×10⁻⁹ kg/m² for 1 mm. A local hot collector therefore needs an effective hot-area escape divided by effective global cold-capture area below 2.54×10⁻⁶ for the 1 cm case. Aluminum condensate exposed to the chamber must remain below about 747 K to keep its equilibrium vapor column below that same budget; the 1 mm value is 715 K. Those are fixed-criterion thresholds, not achieved engineering performance. The ideal bottom-center escape probability from a cold-walled narrow pocket scales as 1/[1+(depth/radius)²]; its required aspect ratio is about 627 for the 1 cm budget. Finite source size, specular reflection, reemission, gas collisions and wall temperature determine whether an actual pocket reaches this value.

The slow temperature pulse also matters. Reusing the separated backed-carbon histories gives aluminum excess evaporation of about 0.0172 kg/m² for B2 and 0.00971 kg/m² for C13. The maxima are **1387 K during slow arrival and 1713 K during a separate microsecond flash**; 1965 K was a sum of noncoincident bounds and is unsuitable as a 100 s evaporation temperature. The new diagnostic retains the separation. The resulting metal vapor is a strong reason to use a local enclosed recovery pocket or another surface rather than expose the entire hot metal layer to the delivery cavity.

Alpha-alumina avoids this pure-metal baseline. An independent seven-species equilibrium calculation includes Al, AlO, AlO₂, Al₂O, Al₂O₂, O and O₂, using their JANAF formation Gibbs energies and unit alpha-Al₂O₃ activity. The outgoing effusive atomic flux is constrained to O:Al=3:2. Total pressures are 1.19×10⁻¹⁹ Pa at 1000 K, 1.80×10⁻¹⁰ Pa at 1400 K and 2.00×10⁻⁶ Pa at 1700 K. These are thermal-equilibrium references for the oxide and its own vapor. Hydrogen/carbon reduction, radiation chemistry, preferential sputtering and defect production are additional source terms.

Integrating that oxide reference over the actual separated carbon histories, with a further 10 K interface allowance, gives 4.94×10⁻¹² kg/m² for B2, 1.41×10⁻¹² for C13 and 5.42×10⁻¹⁴ for A. Their saved-cadence mean is about 69 kg/s across the sphere. The largest isolated-pulse radial-column envelope, using the previous 3.26 geometric factor, is 1.61×10⁻¹¹ kg/m²: approximately 3.5 J/kg for the 1 cm grain criterion. Even this larger thermal source than graphite does not by itself force a global membrane. Actual O/O₂ capture needs a getter or pumping pathway; it cannot be inferred from aluminum condensation alone.

## Hydrogen, sputter products and dust

The saved cycle returns 1.859 postburn protons per produced helium nucleus, approximately 1.18×10⁶ kg/s at its reference rate. This excludes extra carrier hydrogen from a future richer mixture. If it all enters the central cavity at 1000 K, the unit-first-return-capture mean radial column is 7.43×10⁻⁸ kg/m² as H or 1.05×10⁻⁷ as H₂. A 1 cm grain with the 10 kJ/kg allocation therefore requires escaped fraction/effective global capture fraction ≤0.624 or 0.441 respectively. If effective global capture is 1%, at least 99.56% of the H₂ must be retained locally or routed to a gas line. Hydrogen does not condense on a 1000 K oxide. The required pumping or chemical storage is a concrete recovery interface, although these numbers do not establish a large plant-power penalty.

For dust, a 10.5 μm-radius, 3000 kg/m³ particle at 2.4 Mm/s carries enough relative kinetic energy to use the 10 kJ/kg allocation of one 1 cm-radius fuel grain. Requiring fewer than 0.1% of grains to receive such a hit gives a critical-size dust column below 4.62×10⁻¹¹ kg/m². A uniform 100 m/s dust source with unit first-return capture would reach this column at about 16 kg/s. For 1 mm fuel grains the corresponding dust radius and source are 1.05 μm and 1.6 kg/s. This is a useful containment specification, not an inferred sputter yield. It counts the fraction of individual grains harmed; it does not demand that an enormous macropulse encounter no dust at all.

Atomistic sputtering, ceramic spall and fracture can populate different parts of this spectrum. Atomized material obeys the gas-dose criterion, while larger fragments need the rare-interception criterion. The current data do not quantify either source. A locally enclosed, cold-walled recovery cassette can be evaluated against these finite targets, alongside the ion agent’s stopping/nuclear-loss and chemistry results. The thermal-vapor calculation supplies no general prediction of irradiation erosion.

## Useful uncertainty judgments

For a **specified thin Al reflector on a cold polymer grain**, an initial subjective central effective absorptivity is 0.022 with an 80% credible interval of roughly 0.007–0.06 for the 1000 K spectrum. This broad engineering judgment brackets the primary smooth-metal reference while allowing improved cold-metal conductivity and degraded/oxidized practical surfaces; it is not a fitted posterior or a qualified coating specification. With all other assumptions fixed it gives a 1 mm polymer maximum unshielded path near 180 million m, with an approximate 80% interval of 65–560 million m. Testing the manufactured surface’s infrared response would narrow this much more effectively than another unchanged hydrodynamic run. A 1 cm grain passes the reference energy test but still merits optical and thermal-gradient qualification.

The 10 kJ/kg gas/dust allocation and 0.1% damaged-grain fraction are acceptance choices. They are not measured reactor tolerances. The calculations report 1, 10 and 100 kJ/kg gas cases so a subsequent assembly/yield calculation can choose the allowable margin. No credible interval is assigned to an unspecified sputter or dust source; the finite source/column specifications identify what must be demonstrated.

## Mixture changes

With fully stripped catalyst (A,Z), H:catalyst number ratio x and fixed mass density and temperature, n_cat=ρ/[m_u(A+x)], n_e=ρ(Z+x)/[m_u(A+x)]. The initial capture-rate density scales as ρ²x/(A+x)², and nonrelativistic free-free emissivity as ρ²(Z+x)(Z²+x)/(A+x)². Thomson opacity scales with (Z+x)/(A+x). Consequently x=A maximizes initial capture density at fixed ρ, while x=Z^(3/2) minimizes the frozen-temperature free-free/capture ratio. These relationships motivate a mixture grid; pressure, internal energy, degeneracy, evolving composition and spectral transport determine actual gain. More unburned carrier hydrogen also increases the return-gas requirement above.

The current solver has composition-dependent electron density in its fully ionized electron/positron EOS, but a fixed grey transport opacity and an imposed LTE radiation energy aT⁴. The primary-opacity sensitivity screen below supplies physically motivated comparison scales for the new mixture runs. The previous grey escape output supplies a conditional receiver load; it does not establish the true emergent spectrum.

## Actual loss through openings

Gross evaporation, temporary surface retention and imperfect first-pass extraction remain circulating inventory until a terminal route is identified. The old [port geometry](../launcher-reference/ports-and-shielding.md) provides explicit opening scenarios, not a selected required muzzle size. For 1000 permanently open, radial 1 km-radius ports, the projected open fraction is 2.07×10⁻¹⁰; for 100 km-radius ports it is 2.07×10⁻⁶. A central isotropic outward pulse intersects that fraction on its first pass. More outlets at fixed total projected opening do not change it.

For diffuse well-mixed gas, if a fraction s of non-port wall encounters permanently captures the molecule into the return system, repeated flights give eventual chamber-exit probability f/[f+s(1−f)]. At s=0.01 the two example port sizes give approximately 2.07×10⁻⁸ and 2.07×10⁻⁴. A downstream catcher with recovery ζ reduces true station loss by another factor (1−ζ); functioning closed receiver shutters reduce the open fraction during their closed intervals. Chamber-facing sources supplied by bent service passages need no equivalent permanent full-speed hole. Actual directed gas near a port needs its own geometric transport calculation.

If s=0, a nonreacting gas eventually exits through any permanent openings; that limiting result describes the absence of recovery, not an unavoidable reactor property. The local gas-return system must supply s>0 or route material before it reaches the central cavity. No gross oxide or metal vaporization in this study is automatically charged to fresh carbon or catalyst consumption.

## Radiation transport sensitivity for mixture optimization

`opacity_interface.py` reconstructs density, temperature and evolving composition from the four saved corrected64 traces, then evaluates 48 representative initial, maximum-density, maximum-temperature, radiation-release and terminal states. The electron/positron source EOS supplies actual degeneracy and pair populations. The [Poutanen 2017 Compton-opacity fit](https://doi.org/10.3847/1538-4357/835/2/119) provides a primary comparison in its 2–300 keV fit range. It predicts roughly 0.0059–0.0211 m²/kg in the sampled hot states, compared with the existing constant 0.001. The same paper distinguishes Rosseland diffusion opacity from free-streaming spectral response. The code evaluates electron and positron contributions separately as a diagnostic extension; the sampled positron population reaches about 7.8% of net electron number in C13.

A classical, unsuppressed, g_ff=1 free-free calculation, using the explicitly stated [frequency-dependent absorption expression](https://doi.org/10.1051/0004-6361/202451835), supplies an additional scale. In the hot burn it is usually around 10⁻⁵–10⁻⁴ m²/kg, below scattering. It rises strongly in the cold prescribed initial skin. Its density/temperature dependence is useful; treating that classical result as a validated dense-plasma opacity would overstate it. [Relativistic Gaunt-factor calculations](https://arxiv.org/abs/1502.06788) and [dense-mixture average-atom calculations](https://doi.org/10.1103/PhysRevE.102.043209) describe corrections missing from this screen.

A useful bounded hydro comparison therefore retains 0.001 and adds **0.01 and 0.03 m²/kg** for promising mixture cases; 0.1 is a further opaque sensitivity if needed. These constants test sensitivity to the physically motivated hot-state scale. They do not replace a composition-, density- and temperature-dependent opacity calculation. Increased opacity can retain more energy while changing where pressure and radiation are deposited, so the old lower opacity does not establish the sign of every gain or source-preparation change.

The saved traces span approximately 2.2–1647 MK and 92–2.1×10⁸ kg/m³. Their initial outer zones are only about 0.19–0.26 keV. Fully ionized ideal-ion EOS assumptions and missing Coulomb/ionization terms merit a separate warm-dense comparison there. The diagnostic records degeneracy, ion-coupling and plasma-frequency scales. During much of the sampled burn and outer escape, the classical effective thermalization depth across one outer cell is large; late B1 approaches order unity. This supports an opaque interior while leaving the spectral photosphere and its boundary layer unresolved. A fixed LTE aT⁴ radiation energy and fixed grey flux limiter do not demonstrate that the outgoing spectrum equals an outer-cell Planck function.

## Initial nonideal pressure and energy scale

`coulomb_initial.py` reconstructs the initial density, temperature, mass and kinetic-energy profiles, then applies the rigid-background ion-sphere scale U_C/V=−0.9[e²/(4πε₀a_e)]Σn_i Z_i^(5/3), with P_C=U_C/(3V). The strong-coupling term is a diagnostic approximation within [ion-liquid EOS theory](https://doi.org/10.1103/PhysRevE.62.8554). At moderate coupling its attraction is stronger than the complete liquid correction; this is a useful scale rather than a rigorous universal error bound. Electron exchange/polarization, bound states and quantum ions remain additional terms.

| Baseline stage | Ion-sphere energy / total prepared energy | Ion-sphere energy / initial internal energy | Mass-weighted magnitude of pressure correction / EOS pressure | Maximum local pressure ratio |
|---|---:|---:|---:|---:|
| A, H1 | −0.91% | −13.0% | 8.0% | 21.6% |
| C13, H1 | −0.78% | −11.2% | 6.9% | 18.6% |
| B2, H2 | −0.70% | −10.0% | 6.2% | 16.6% |
| B1, H1 | −0.18% | −2.1% | 1.35% | 3.4% |

The denominator of the pressure ratios includes the source EOS's Fermi electrons, ideal ions and LTE radiation. The whole-pulse prepared-energy correction is below 1% in these cases because most initial energy is convergent kinetic energy. Local pressure corrections can still shift compression or ignition near a threshold. H2/H4 comparisons retain the same density/temperature profiles and recalculate electron/ion density and internal energy; the JSON reports both unchanged-flow prepared energy and the original fixed-prepared-energy denominator. A thermodynamically consistent nonideal sensitivity for final contenders is more informative than interpreting Γ≈8 alone as a catastrophic energy error. Eight source-energy reconstruction checks pass.

The actual saved electron table extends to10⁵ K, so the prescribed≈2 MK outer skin lies inside its interpolation domain. This says nothing by itself about complete ionization or Coulomb accuracy. Table domain and constitutive physical approximation are separate questions.

## Methane offers a different cold-feed tradeoff

CH₂/CH₄ combinations can supply H:C ratios between2 and4 without free condensed H₂. Methane still requires a specified sealed-core pressure/expansion path. [Younglove's primary calorimetry](https://pmc.ncbi.nlm.nih.gov/articles/PMC6742817/) measured dense liquid near439 kg/m³ with C_V≈34 J/mol/K. Along an almost fixed-density measured path,102.4→117.4 K raises pressure4.84→28.25 MPa while accepting only about 32 kJ/kg. A100 nm optical skin on a1 cm sphere cannot supply those containment stresses. Expansion space, a structural capsule or protected flight is therefore part of this particular recipe. Saturated-liquid C_σ is a different thermodynamic path and should not be integrated as sealed-core C_V.

Using the earlier ordinary-Al optical reference, a1 cm,439 kg/m³ sphere in the 1000 K cavity absorbs roughly 0.39 MJ/kg during the full458 s flight. The above32 kJ/kg interval is used in about 38 s, corresponding to approximately90 million m. This is an illustrative near-isochoric interval, not a total methane survival allowance: colder solid preparation, fusion, expansion and other recipes change the usable heat budget. The result supports a concrete cold-feed design comparison rather than an unavoidable20 K hydrogen penalty for every richer hot mixture.

## Follow-on mixture and gas results

[Warm-return gas](return-gas.md) adds He product and the A-H2/B1-H2 source case, with 1/10/100 kJ/kg grain-heating allocations. [Mixture photon deposition](mixture-radiation.md) compares both outer-temperature and flux-equivalent spectral proxies for C13 H2/H4 and the completed 0.01-opacity runs. The softer alternative materially changes oxide flash heating while leaving graphite a live baseline. Neither a universal separating film nor a mandatory radius increase follows from the current evidence.
