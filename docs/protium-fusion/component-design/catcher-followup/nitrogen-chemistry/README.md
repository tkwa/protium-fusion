# ✨ Nitrogen chemistry in the hot catcher

The catcher needs an explicit nitrogen-bearing volatile stream alongside retained C–N sites. Published implantation and hot-graphite experiments support both populations. They do not determine their proportions in this reactor. Assuming that all returned nitrogen remains chemically fixed is particularly consequential because a small cavity-directed branch can affect the incoming fuel: the current gas comparison reaches a 10% increase in total residual-gas column when **7.2506 ppm of nitrogen return enters the cavity as N₂**, at global N₂ capture probability 0.01.

This calls for species-inclusive collection and recovery, with the volatile fraction and front routing modeled separately. It does not establish a need to discard the graphite receiver or replace captured matter. Nitrogen retained in the laminate, collected behind it, or recaptured from the chamber remains available to the material cycle.

## Exact source and concentration scales

[calculate.py](calculate.py) reads the frozen [inputs.json](inputs.json) and produces [results.json](results.json). The source is the delivered-B1/C13H2 family and representative 350 kg/m² cadence. Atomic masses are used for process flows. N14, N15, O14 and O15 at the wall all contribute to eventual nitrogen; arriving N13 becomes carbon.

| Stage | Eventual N dose, nuclei/cm² per pulse | N mass, µg/m² per pulse | Fraction arriving as O14/O15 | Arriving H / eventual N nuclei |
| --- | ---: | ---: | ---: | ---: |
| A | 2.63487×10¹⁵ | 647.700 | 30.392% | 2.197 |
| B1 | 5.33119×10¹² | 1.32791 | negligible | 5.941 |
| B2 | 2.45235×10¹⁵ | 570.235 | 20.794% | 0.509 |
| C13 | 2.78726×10¹⁵ | 670.602 | 22.001% | 0.655 |

The mean nitrogen return is **24.4528 million kg/s**, spread over 1.52053×10¹⁹ m². The annual dose is **2.10565×10¹⁷ nuclei/cm²**, or **50.750 mg/m²**. The exact nominal CSDA stopping bins give the largest no-removal concentration in the 2–5 µm bin: about 36–53 N atoms per million carbon atoms after one large pulse, and **0.00421 N/C after a year without extraction**. These are bin averages from ballistic stopping, before diffusion, rework and chemical release. Straggling, aligned channels and unresolved Bragg-profile structure can change local concentrations.

Important arriving species have CSDA median depths of roughly 2–7 µm and 99th-percentile depths up to about 17 µm. The packet preserves individual-species quantiles rather than combining their percentiles. Laboratory nanometre-thick films with several atomic percent nitrogen are therefore chemically different even when their areal implantation dose resembles one shot.

**23.355% of eventual nitrogen initially stops as radioactive oxygen.** Oxygen bonding, volatile oxygen-bearing products and subsequent nuclear transmutation must be included when assigning initial chemical populations. A nitrogen-only implantation experiment cannot specify these sites. The source audit preserves the atom counts; it does not calculate postimplant recoil or molecular branching. Chemical products containing these nuclei still belong in the same recovery balance.

## Primary evidence and its applicability

The full primary files and their hashes are indexed in [sources.json](sources.json). PDF methods and numerical tables were checked against page renders. Raw papers remain in the research output directory.

| Primary source | Measured or calculated result | Transfer to this catcher |
| --- | --- | --- |
| [Li et al., HTR2014-41336](https://nucleus.iaea.org/sites/htgr-kb/HTR2014/Paper%20list/Track4/HTR2014-41336.pdf), full conference paper | IG-110 implanted with 100 keV N⁺ at 5×10¹⁶–5×10¹⁷ cm⁻² had a desorption peak assigned to N₂ at 953.75–960.65 K during a 10 K/min ramp. | Relevant baseline temperature and graphite; substantially greater single-dose loading and shallower implantation. The m/z=28 channel also receives CO. Peak evidence supports a release channel, while the integrated amount needs the normalization audit below. |
| [Kumar, Ganguly and Papakonstantinou, 2012](https://arxiv.org/abs/1205.3155), full paper, DOI 10.1088/0953-8984/24/23/235503 | Vertically aligned graphene flakes received 2 keV N⁺, 10¹⁵ cm⁻². XPS N concentration fell from 9.59 to 5.56 atomic percent after heating to 1073 K; graphitic N became dominant. | Comparable areal dose, very different depth, local concentration and edge density. XPS establishes redistribution of bonding populations and a surface concentration decrease; it does not measure an absolute N₂ release fraction. |
| [Kondo et al., 2012](https://arxiv.org/abs/1206.6593), full paper, DOI 10.1103/PhysRevB.86.035436 | HOPG exposed to 200/500 eV N ions retained pyridinic and graphitic N after 900±50 K for 300 s. XPS concentration changed 3.1→2.7 atomic percent. | Supports a stable population near operating temperature. The low STM dose is not the XPS dose: the authors did not measure the XPS implantation current. Gas identity was inferred, not measured. |
| [Schlüter, Hopf and Jacob, 2008](https://doi.org/10.1088/1367-2630/10/5/053037), full paper | Atomic H enhanced erosion of a-C:H under 30–900 eV N₂⁺ bombardment. H/ion flux was about 380–400; investigated temperatures were 110–350 K. | Establishes coupled H–N–C chemistry, with ion damage and hydrogen passivation. Reactor H/N is roughly 0.5–6, temperatures are higher, and deposition is deeper. Its erosion yield cannot be applied to our pulse. |
| [Bystrov et al., 2013](https://www.differ.nl/bibcite/reference/7800), full paper, DOI 10.1063/1.4822166 | Fine-grain graphite lost carbon under ≤35 eV nitrogen plasma. Yield depended on temperature, energy and flux; the high-flux empirical temperature fit peaked near 943 K. | Confirms hot graphite can participate in nitrogen-driven chemical erosion. Flux was roughly 4×10²⁰–4×10²³ m⁻²s⁻¹, versus reactor all-ion sampled peaks below 10¹⁸. Mass loss was measured; CN/C₂N₂/HCN product signals were below useful sensitivity. No product fractions transfer. |
| [Deng and Souda, 2002](https://doi.org/10.1016/S0040-6090(01)01785-0), publisher abstract/sections accessed | Carbon-nitride films made by nitrogen implantation into graphite emitted N₂ and C₂N₂ on heating, using gas analysis and secondary-ion measurements. | Supports multiple volatile products. Full numerical inventory and kinetics were unavailable; no reactor retention fraction is taken from this source. |

A further [2026 quantum-chemistry/kMC study](https://doi.org/10.1021/acs.jpcc.6c03128) treats nitrogen migration, pairing and N₂ evolution from graphitic sites at dilute concentrations. The accessible primary table and conclusion concern ideal graphitic nitrogen and emphasize defect/edge extensions. Its direct full-text retrieval returned HTTP 403, so no kinetic prediction in this packet depends on it. It is a useful follow-up model for the strongly retained population, rather than a substitute for the damaged, H/O-containing material.

Li's integrated-signal normalization does not close as an implanted-N fraction from the stated geometry. Ten 2×2×1.5 mm cuboids have at most 2 cm² of external area. Even assigning implantation to every face, the lowest dose contains at most 0.001861 mL STP-equivalent N₂; Table 4 lists 0.01303 mL. The factor-seven difference could involve batch/area normalization, background or other nitrogen inventory. The table alone cannot resolve it. We use the peak location and preserve the numerical discrepancy, rather than treating the listed volume as a measured retained fraction.

## Time-dependent model and what remains unidentified

A usable model separates at least three populations: labile implanted/edge sites, substitutional or strongly bound C–N sites, and already formed mobile molecules. O14/O15-origin sites require their own initial branching and decay history. H exposure can change both bonding and volatile products. Neither a single bond energy nor a TDS peak supplies all these rates.

For specified first-order populations, the supplied `release_kernel(times, temperatures, frequency, barrier)` implements

`F_j(t) = 1 − exp[−∫ ν_j exp(−E_j / kT(t)) dt]`.

It can be convolved with each implanted source cohort and multiplied by its population weight. A separate helper implements homogeneous recombination `dn/dt = −2 k(T)n²`; the initial concentration explicitly affects that release. The helpers have analytic checks, but **no nitrogen rates or population weights have been fitted**. Their output becomes a reactor prediction only after supplying a material-specific network and transport to an actual free surface or collection channel.

As a conditional scale, a first-order population with ν=10¹³ s⁻¹ releases 10% in 100 s at 1000 K if its barrier is 3.170 eV; the corresponding barrier at 1500 K is 4.755 eV. These are acceptance diagnostics under a stated kinetic law. Actual temperature histories must be integrated: the latest laminate/post material peak is about 1496 K, while the earlier 1668 K value was a conservative envelope. Prompt photon bounds last microseconds. Holding either maximum for hours would invent a release load.

The main unidentified quantities are the initial bonding distribution after the collision cascade, defect-assisted N motion and pairing, H/O-dependent molecular products, and transport from their actual 2–17 µm stopping region to the nearest channel. A 17.5 cm bulk-graphite diffusion assumption is unnecessary. No defensible 80% credible interval for a reactor N₂ release fraction emerges from these experiments; the scenario fractions below are deliberately unweighted.

## Gas acceptance and constructive recovery

Let `r` be the fraction of returned nitrogen liberated as N₂ during the accounting interval, and `b` the eventual fraction of that liberated stream entering the cavity. At the frozen gas model's `g_N2=0.01`, the 10%-total-column comparison is

`r × b ≤ 7.25056×10⁻⁶`.

| N liberated as N₂, r | Maximum cavity fraction of that N₂, b | N₂-bearing N handled, kg/s |
| --- | ---: | ---: |
| 0.01% | 7.2506% | 2,445 |
| 0.1% | 0.72506% | 24,453 |
| 1% | 0.072506% | 244,528 |
| 10% | 0.0072506% | 2,445,279 |
| 100% | 0.00072506% | 24,452,786 |

A 1% liberated population with 0.1% front escape would increase the baseline column by 13.79%; it warrants treatment. Conversely, efficient sealed collection can accommodate a large liberated fraction. Changing global N₂ capture or the permissible pellet column changes the comparison directly. Time dependence changes the transient gas load; storing and later emitting the same amount does not remove its steady-state throughput.

The constructive modification is to extend the existing channels, isolated manifold and warm cleanup to nitrogen-bearing species. Mobile nitrogen can be removed during or after the pulse; strongly bound nitrogen can be transferred through enclosed selective rework. N₂, CN-bearing compounds and any H-bearing products need distinct transport/separation coefficients. Their isotope and C/H inventories must return to the source-preparation ledger. Molecular formation is not automatically fresh-carbon consumption. Keeping nuclear daughters enclosed during decay also prevents their chemical history from being silently replaced by a pure-N input.

The atom-scale energy coefficient is small: **one eV of positive process free energy per returned N atom, supplied at 50% efficiency, costs 2.138×10⁻⁵ W/m²**. Actual reagents, temperature, separation and dissipation must define that work. Activation barriers are not a charge paid independently for every atom, and recycling the surrounding graphite does not require atomization. Whole-layer handling and apparatus heat leaks can still dominate, as the geometry packet already shows.

A one-day mean nitrogen hold-up is 2.113×10¹² kg globally but only 0.139 mg/m². This is retained processing inventory. A mean residence-time allowance cannot replace capacity for an individual large pulse, whose dose is roughly 0.6–0.7 mg/m².

The next useful computational work is a branching calculation for N/H/O collision-damaged sites, followed by a kinetic network on the actual temperature and dose histories, coupled to side-channel collection. It should retain a graphitic-N population and explicit volatile products. It can discriminate a need for faster selective extraction or stronger cavity cleanup before a large experiment. Material-specific rate and branching measurements remain necessary to assign a narrow quantitative release fraction.

## Reproduction and boundaries

Run `python3 calculate.py` in this directory. Fourteen checks cover independent atom/mass closure against the cold material ledger, depth-bin conservation, analytical kinetic helpers and the gas threshold. Parent source hashes and raw-paper hashes are recorded. This packet adds no new hydrodynamic simulation, assumed sputtering yield, release fraction, carbon loss, or refurbishment cadence. It preserves the archived source and catcher studies.
