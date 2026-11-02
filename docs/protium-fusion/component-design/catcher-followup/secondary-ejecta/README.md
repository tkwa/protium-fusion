# ✨ Carbon ejecta, local recapture and collector chemistry

The closed sphere preserves material inventory. Its large radius nevertheless gives emitted carbon days in transit, during which it contributes to the column crossed by incoming fuel. **The present physical-sputtering extrapolation looks manageable; thermal radiation-enhanced sublimation and volatile chemical erosion are consequential unresolved surface processes.** A protected local collector can change this conclusion by returning carbon to a nearby solid surface before it enters the main cavity. Collector capacity and chemical regeneration work are small at the relevant mass fluxes.

This packet uses the delivered-B1/C13H2 128-zone family and the exact 350 kg/m² dispatch. It provides response coefficients, two explicitly unvalidated erosion diagnostics, a spherical ballistic transport calculation, and conditional collector accounting. No measured reactor yield or permanent sticking coefficient is assigned. Existing archives are unchanged.

## Primary evidence and its applicable scope

| Source, accessed in full | Relevant evidence | Application here |
|---|---|---|
| [Eckstein et al., *Sputtering Data*, IPP9/82 (1993)](https://pure.mpg.de/rest/items/item_2131245/component/file_3587307/content), pp8,42,46 and accompanying tables | Normal-incidence graphite measurements and binary-collision calculations. The tabulated C/He/O experiments mostly cover tens of eV to 10 keV. At 10 keV the reported C experimental yield is about 0.507 carbon atoms/ion, versus a calculated 0.210. Low-keV oxygen erosion includes substantial chemistry. | Its physical-yield fits can motivate a MeV screen; they do not validate the present spectra. Nitrogen fit parameters below are an explicit C/O interpolation. Surface orientation and roughness matter. |
| [Roth, Vietzke & Haasz, *Erosion of graphite due to particle impact* (1991)](https://www-pub.iaea.org/MTCD/Publications/PDF/Pub23_web.pdf), pp63–67 | Distinguishes physical sputtering, chemical erosion and radiation-enhanced sublimation (RES). Gives an empirical RES temperature/flux expression for low-keV irradiation, with flux data roughly 10¹³–10¹⁷ ions cm⁻² s⁻¹. RES is associated with mobile irradiation defects and predominantly atomic carbon. | The actual pulse temperatures enter the temperature range where RES deserves attention. Its incident-energy dependence at MeV energies remains unresolved. |
| [Nygren et al., *Radiation Enhanced Sublimation of Graphite in PISCES Experiments* (1989)](https://www.osti.gov/servlets/purl/5244031) and [*Low Energy Data…* (1990)](https://www.osti.gov/servlets/purl/7184809) | Helium-plasma experiments at 30–300 eV and 900–2000°C. The discussion distinguishes thermal RES carbon from energetic physically sputtered carbon, and treats redeposition explicitly. | Supports a thermal-emission transport case for RES. Their plasma fluxes, energies and temperatures do not set the present yield. |
| [Bossion et al., A&A692,A249 (2024)](https://doi.org/10.1051/0004-6361/202452362), TableA.1 | ReaxFF molecular dynamics on ten amorphous-carbon structures, 740 incident particles per temperature. Thermal-C chemisorption is 0.247±0.018 at1000K, 0.295±0.016 at1250K, and 0.299±0.020 at1500K; quoted errors are the authors’ sampling errors. CO has negligible sticking in their warm cases. | Supports nonzero short-time carbon capture. The 0.45ns observation and mean bulk binding-energy estimate do not establish permanent retention on an irradiated collector. |
| [Lehtinen et al., PRL91,017202 (2003)](https://doi.org/10.1103/PhysRevLett.91.017202) | Spin-polarized DFT gives an isolated C adatom on clean graphene adsorption energy1.40eV and migration barrier0.47eV. | Shows why a pristine basal surface and an edge-rich collector need separate retention models. |
| [NIST-JANAF SiC](https://janaf.nist.gov/tables/C-100.html), [TiC](https://janaf.nist.gov/tables/C-107.html) | Formation Gibbs energies at1000K are −63.517 and −172.653kJ/mol respectively; both remain negative at1500K. | Carbide capture and regeneration have favorable inventory/work scales. These equilibrium data supply no sticking rate or passivation lifetime. |

The earlier [nitrogen chemistry packet](../nitrogen-chemistry/README.md) supplies additional low-energy N/H chemical-erosion evidence. None of these sources establishes current CH₄, HCN, CO or dust fractions. Displacement counts and NRT damage energy are not converted into sputter yields.

## Exact incident flux and reusable yield coefficients

Let `Y_i` be the carbon atoms emitted into the chamber per incident nucleus of group i, including any local capture/routing before cavity entry. Source atomic counts are projected from the saved endpoint with finite-flight weak decay. Oxygen here means oxygen at the wall, including radioactive oxygen; its eventual decay chemistry is separate.

| Incident group | Nuclei/s | Carbon source at Y=1, kg/s | Gross carbon removal at Y=1, mg/m²/year |
|---|---:|---:|---:|
| H | 1.92828×10³³ | 3.84238×10⁷ | 79.746 |
| He | 4.15089×10³² | 8.27127×10⁶ | 17.166 |
| C | 4.99709×10³² | 9.95744×10⁶ | 20.666 |
| N | 8.10638×10³² | 1.61531×10⁷ | 33.525 |
| O | 2.36997×10³² | 4.72252×10⁶ | 9.801 |

The whole sphere has area1.52053×10¹⁹m². A common C/N/O yield of one therefore gives **63.992mgC/m²/year**, equivalent to31.996nm/year of graphite at the stipulated density2000kg/m³. Redeposition can return this mass; this is gross removal, not net material consumption. Reaction processing and external carbon supply should charge only actual permanent export, unrecovered material or deliberate replacement.

For each stage, the code integrates the individual isotope yields over its velocity/energy distribution. Median N energies range3.50–7.74MeV, and the C13-stage carbon median is26.98MeV. The physical-yield diagnostic extends the published fit far beyond most of its experimental calibration. The alternate fit uses the published He/C experimental-fit coefficients with the same O calculation and another explicit N interpolation. It is a sensitivity, not a confidence interval.

| Source diagnostic | Carbon kg/s | Gross thickness-equivalent µm/year | Atomic-C steady column at3eV, kg/m² |
|---|---:|---:|---:|
| Calculated physical-fit extrapolation | 2.37162×10⁵ | 0.0002461 | 1.23836×10⁻⁸ |
| Alternate He/C-fit extrapolation | 3.57478×10⁵ | 0.0003710 | 1.86660×10⁻⁸ |

This low erosion rate does not bound irradiation damage, release of previously trapped material, or chemical products.

## Thermal erosion is the consequential source uncertainty

The explicit RES diagnostic is

`Y_RES = 54 A^1.18 exp[−0.78eV/(kT)] (Γ/10^16)^−0.1`,

with Γ in ions cm⁻²s⁻¹. Each stage uses its warmest exported150µm-skin temperature waveform and actual sampled total charged-ion flux. Peak waveform temperatures are1352K,1190K,1393K and1494K for A/B1/B2/C13. The formula is also evaluated with contributions below1200K removed, purely to show the importance of the low-temperature extension. Neither evaluation is a prediction or a rigorous bound. Most B1 pulses have cooler backgrounds; B1 flux is entirely below this expression’s quoted flux range. MeV recoil-depth distributions, defect recombination, mixed chemistry and irradiation history require a different calibration.

The all-temperature and ≥1200K-only diagnostics produce72.34 and51.57million kgC/s respectively. With thermal1000–1500K emission and permanent capture on the next global wall, their steady columns are respectively(1.287–1.576)×10⁻⁵ and(0.917–1.124)×10⁻⁵kg/m². Finite-dispatch peaks span(2.054–2.384)×10⁻⁵ and(1.676–1.924)×10⁻⁵kg/m². Gross removal is only0.0751 or0.0535µm/year, and the largest individual diagnostic pulse removes about1.12nm-equivalent.

These figures establish that a high-temperature release coefficient deserves treatment, even when structural erosion is modest. Determining its actual value needs a surface model incorporating the near-surface recoil source, migration, recombination, chemical reactions and strong trapping; total displacement energy alone cannot settle it. A collection geometry can reduce the sensitivity to this uncertain coefficient.

## Atomic carbon: transient column and recapture

The transport calculation assumes uniform full-sphere emission, a cosine angular distribution, a broad flux-Maxwell speed distribution, and straight collisionless trajectories until the next wall encounter. The illustrative physical-ejecta mean energies are1,3 and10eV. Thermal RES cases use mean emitted energy2kT at1000 and1500K. These speed distributions are modeling choices.

For mono-speed emission of areal mass σ, `x=vt/R` gives radial column `σ h(x)`, where

`h(x) = [(1+x²)ln(1/|1−x|) − (2x−x²)/2] / x²`, for `0<x<2`.

The perfectly monochromatic, simultaneous case has an integrable focusing singularity. Integrating over the broad speed distribution gives a finite single-pulse maximum3.250σ. Its time integral reproduces `4R/v`; averaging over the speed law gives steady column `4R q/v̄`, where q is emitted kg/m²/s. Mean first-flight residence is12.78days for thermal1000K C or3.06days for the3eV case. The finite schedule starts empty, uses the actual dispatch, and treats the few-minute material arrival as an impulse on this days-long scale.

For common C/N/O yield Y and unit next-wall capture:

| Emission | Steady column/Y, kg/m² | Finite370-day peak/Y, kg/m² |
|---|---:|---:|
| Thermal1000K | 6.7170×10⁻⁶ | 9.0627×10⁻⁶ |
| Mean1eV | 2.7886×10⁻⁶ | 5.3160×10⁻⁶ |
| Mean3eV | 1.6100×10⁻⁶ | 4.3328×10⁻⁶ |

The existing fast99.9%-speed,2cm C13 lower-enthalpy example allows3.754×10⁻⁷kg/m² above ideal H/He. Allocating its entire remaining steady margin to carbon corresponds to Y=0.0559 for thermal1000K or0.233 for3eV. Finite-pulse maxima give smaller allocations. Other sources share that margin, and an actual pellet-family redesign can choose a different margin; these are comparison points, not universal yield requirements.

Permanent capture probability s matters. If unretained atoms reemit thermally at1000K, the steady result becomes

`Σ = Σ_first + [(1−s)/s] Σ_thermal_single_flight`.

For the calculated physical-fit source with3eV initial ejecta, this gives1.24×10⁻⁸,6.40×10⁻⁸,1.67×10⁻⁷ and4.77×10⁻⁷kg/m² for s=1,0.5,0.25,0.1. Short-time chemisorption is only one part of s. With an illustrative10¹³s⁻¹ prefactor and the clean-graphene DFT energies, isolated-adatom residence is1.14µs at1000K and5.06ns at1500K; equal migration/desorption prefactors and a0.142nm hop give RMS exploration distances31.3 and5.18nm. These are conditional kinetic estimates, not a required collector texture.

## Local collection and small chemical inventories

The geometry team’s radial honeycomb places protected collector faces near the emitting floor. Its escape fraction must include both repeated encounters and exposed-tip emission. Rough, edge-rich carbon can offer incorporation sites. Refreshable Si/Ti surfaces provide a second route through carbide formation. Protected collectors can deliberately differ from the basal-aligned gas barrier.

For Y_CNO=1, one-to-one carbide capture requires149.77mgSi or255.26mgTi per receiver m² annually. At stipulated metal densities2330/4500kg/m³ this is64.28/56.72nm of metal per unit collector area; forty times the collector area reduces it to1.607/1.418nm/year. A10nm Ti layer over40 times the receiver area weighs0.0018kg/m² and holds0.4513gC/m² stoichiometrically, about7.05years at Y=1 before competingN/O chemistry. Accessibility through a passivating carbide and refresh cadence determine useful service life. The useful optional specification is a protected, warm, renewable reactive collector face, with both initial capture and permanent incorporation established separately.

Reversing carbide formation to metal plus graphite at1000K has reversible Gibbs work0.6583eV/C forSiC or1.7894eV/C forTiC. With the user’s50% minimum chemical-manufacturing efficiency, these imply21.47/58.35µW per receiver m² at Y=1. A stipulated positive processing requirement2eV/C gives65.2µW/m² at50% efficiency;100eV/C gives0.003261W/m². The larger RES diagnostic scales this to about0.00765W/m². These numbers leave room for targeted capture, separation and refresh; they do not price a complete service machine or assume atomizing all recycled receiver material. Competing N/O/H chemistry, passivation and collection kinetics remain the substantive component issues.

The sampled peak carbon source divided over40 times the receiver area is0.00647/0.0000551/0.00833/0.00836 atoms perinternalnm²/s for Y_CNO=1 in A/B1/B2/C13. The RES diagnostic gives0.0183/0.0000979/0.0261/0.0416. These are area-normalized handling scales. First-collision bands near the cell floor can have higher local flux. The largest value averages oneC per24s pernm², leaving a low mean atom-handling rate even for that adverse source scenario.

## Volatile carbon and dust need their own balances

Stable CO, CO₂, CH₄ and nitrile scenarios use whole-molecule mass and a warm-pump capture g, rather than unit atomic-carbon deposition. At1000K and g=0.01, the carbon-emission yield into the cavity needed to add10% of the ideal H/He column is4.66×10⁻⁶C/H forCH₄,1.64×10⁻⁵C/O forCO,8.32×10⁻⁶C/O forCO₂,5.06×10⁻⁶C/N forHCN or3.79×10⁻⁶C/N forC₂N₂. Each includes liberation and cavity routing. These steady stable-molecule sensitivities are not predicted chemical fractions; radioactive O decay can change a molecule before recapture. Carbon-rich volatile products should join the gas extraction/separation ledger and return their atoms to processing.

Dust production is independently parameterized. For a2cm,80%-filled C13 core at the fast4.151Mm/s quantile, one carbon grain of radius1,10 or100µm deposits a kinetic-energy equivalent of2.66,2662 or2.66million J/kg-core if fully coupled. A uniform column giving1% probability of at least one such impact is respectively6.70×10⁻¹⁴,6.69×10⁻¹¹ or6.63×10⁻⁸kg/m². These are geometric/energy comparisons, not fracture predictions or mandatory failure rates. At1ppm of incoming mass turned into dust, unit first recapture and speed10m/s give steady column1.19×10⁻⁹kg/m², but the global steady state would take years to establish. Slow grain trajectories also depend on forces and spatial concentration. Full-wall average columns cannot bound a local service plume. Enclosed repair and local fragment collection address this without requiring zero erosion.

## Reproduction and checks

Run `calculate.py`, then `sticking_screen.py` and `getter_screen.py` with NumPy available. `results.json` records exact source/trace/temperature/schedule hashes, isotope quadrature, model formulas and45 checks. The32→64 incident-quadrature comparison changes physical-fit emission by less than10⁻⁵; the ballistic integral agrees with4 to better than10⁻⁶ relative;256→512 angular/speed quadrature differs by less than2×10⁻⁴ in the normalized impulse response. Halving dispatch time bins1800→900s changes peaks by less than1%. Source uncertainty and surface chemistry dominate these numerical errors. `evidence.json` and `manifest.json` identify downloaded primary material and all new artifacts. No reactor-yield credible interval is manufactured from this sparse evidence.
