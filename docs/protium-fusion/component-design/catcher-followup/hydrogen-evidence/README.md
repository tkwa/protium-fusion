# ✨ Hydrogen recovery from the catcher: evidence and conditional kinetics

This packet concerns the delivered-B1/C13H2 source128 family, a 1000 K receiver with excursions approaching 1668 K, and hydrogen implanted within roughly the first 0.1–100 µm. The current 350 kg/m² graphite backing is 17.5 cm thick. **The relevant recovery distance is the distance to a connected collection surface, not the backing thickness.** Implantation, movement within the material, desorption, collection, and return to the fuel plant are separate material transfers. None is automatically a loss from the closed sphere.

The evidence supports an engineered shallow collection layer and selective treatment of trapped hydrogen. It does **not yet select a numerical release fraction or an irradiated-layer diffusion tensor**. A passive bare surface, a capped collection channel, and targeted extraction have different gas loads even when their nuclear stopping is identical. This distinction can change the gas-handling design substantially; the chemical work of handling the hydrogen itself is small.

## Primary evidence and applicability

Full papers read are retained in `/home/tkwa/src/protium-fusion/output/catcher-followup/hydrogen-evidence/`; `sources.json` records access and hashes. This table distinguishes a full-paper measurement from a model and an accessible abstract.

| Source | Evidence actually available | Relevance to this receiver |
|---|---|---|
| [Atsumi et al.2013](https://doi.org/10.1016/j.jnucmat.2013.03.041), requested paper | Publisher abstract, introduction and section excerpts; **full text was not obtained**. Gas-charged graphites/CFCs show several D₂ TDS peaks; high peaks around1300–1400 and1600–1700 K are assigned effective activation energies3.48 and6.93 eV. CD₄ follows the higher peak. | Demonstrates distinct populations. The material was gas charged without implantation damage. Peak locations give neither the populations in our receiver nor absolute isothermal release rates. The stated~20% CD₄/D₂ ion-current ratio is not a calibrated molecular yield. No retention fractions are inferred here. |
| [Muto & Morita1992](https://doi.org/10.3327/jnst.29.980), full experiment | Saturated keV H implantation in IG-110U; isothermal annealing723–873 K. ERD finds an initial decrease followed by a slow tail. A local free-H recombination/retrapping model fits the profiles. Methods say3 keV H₂⁺, abstract5 keV: the discrepancy is retained. | Supplies the audited concentration-dependent kernel below. Implantation dose3×10²²/m² and nanometre depths differ strongly from this reactor. Applying the fitted law above1000 K is an explicit sensitivity calculation. |
| [Hasebe, Sonobe & Morita1992](https://doi.org/10.3327/jnst.29.859), full experiment | Previously implanted H/D released by1.5 MeV He; room temperature to523 K. Release depends on initial trap occupancy; effective H detrapping/recombination cross-section2.9×10⁻²² m². | Shows radiation-assisted release and retrapping. The present B1 He dose gives dimensionless2σΦ=1.2189×10⁻⁴. That does not establish complete cleaning, especially at low occupancy. Our energies, depths and temperatures differ. |
| [Tsuchiya & Morita1995](https://doi.org/10.1016/0022-3115(94)00593-1), abstract | Above1000 K, implanted-H profiles develop depth gradients; reported activation energies3.1±0.6 eV for diffusion and3.8±0.6 eV for atomic ejection. | Supports a change from the lower-temperature local-recombination regime. Verified prefactors and full profiles are unavailable here; these barriers are not converted into absolute residence times. |
| [Reichart et al.2006](https://doi.org/10.1016/j.nimb.2006.04.014), full experiment | 2.25 MeV protons in HOPG produce a narrow retained-H peak at46 µm. No resolved c-axis broadening. Microspots have18–36% fewer detected H atoms than implanted; lateral migration and local heating are discussed. | Directly relevant implantation depth. Missing local inventory is not measured escape. Broad-implant dose uncertainty was~50%; this is not a precise prompt reflection measurement or a1000 K transport coefficient. |
| [Deng et al.2018](https://doi.org/10.11889/j.0253-3219.2018.hjs.41.050502), full experiment | Three centimetre-scale graphites gas charged at623 K and desorbed up to973 K; gas release extends over2–8 days, with substantial grade dependence and residual inventory. | Rules out universal instantaneous release from arbitrary graphite at~1000 K. Long specimen/grain paths and gas charging differ from a shallow channel skin. Abstract/body interchange the smaller two grade-specific fractions, so those labels are not used. |
| [Warrier et al.2006 author manuscript](https://www.mv.helsinki.fi/home/knordlun/pub/NuclFus2006.pdf), full multiscale simulation | Atomic transport through grains, surface paths and pores; porosity/jump lengths change effective diffusion by orders of magnitude. Trace-H simulations omit molecular formation. A surface-channel coefficient used is1.2×10⁻⁴exp(−0.9 eV/kT)m²/s. | A conditional fast-channel coefficient, not the diffusivity of all implanted H. Trap escape, access to the channel, channel connectivity and molecular chemistry must be supplied separately. |
| [Herrero & Ramírez2010](https://doi.org/10.1088/0022-3727/43/25/255402), full tight-binding MD | Ideal interlayer atomic-H diffusion fit D∥=2.0×10⁻⁷exp(−0.38 eV/kT)m²/s; intact-sheet crossing has a calculated barrier~5 eV. Simulations cover300–2000 K. | Supports a physically possible lateral collection route in aligned material. It does not quantify the effective anisotropy after implantation, seams, grain boundaries or cracking. |
| [An et al.2019](https://doi.org/10.1002/smll.201901722), full primary DFT/MD | Graphite neutralizes a proton, and H binds to carbon; no facile H-transfer pathway was found under the investigated conditions. | An independent reason not to promote one ideal-crystal mobility calculation into an engineered coefficient. Trapping, molecular state and model conditions matter. |
| [Sun et al.2020](https://doi.org/10.1038/s41586-020-2070-x), full experiment and theory | Monolayer graphene transmits H₂ at295–323 K through proposed ripple-assisted dissociation/atom flipping, with measured activation energy1.0±0.1 eV. Bilayer controls showed no detectable transmission during3 days at323 K. | A He-tight monolayer does not establish a hydrogen barrier. Neither the monolayer nor the bilayer experiment qualifies an irradiated1000–1668 K cap. No extrapolated hot leak rate is credited. |
| [Franzen & Vietzke1994](https://doi.org/10.1116/1.578829), abstract and institutional metadata | Atomic H release becomes detectable around900 K; molecular and atomic thermal-release peaks differ. | Gas handling must allow atomic H as well as H₂ and hydrocarbons. Absolute species fractions remain unspecified. |

The Atsumi full-text access gap remains explicit. Publisher text-mining access returned metadata only, OpenAlex listed no open copy, and the author/institutional searches did not provide one. The model uses complete independently retrieved implantation studies rather than pretending to have extracted Atsumi’s curves.

## Exact hydrogen load and mass convention

`input-snapshot.json` pins the source spectrum by SHA-256. Counts use nuclei. Hydro mass uses A×mᵤ; the separate neutral-protium mass uses1.00782503223mᵤ. These are two mass conventions for the same atoms.

| Stage | H nuclei per m² per pulse | H baryon mass, µg/m² | Period, s | Mean H:C number ratio if spread uniformly over10 µm |
|---|---:|---:|---:|---:|
| A | 5.78923×10¹⁹ | 96.1324 | 810809 | 5.768×10⁻⁵ |
| B1 | 3.16746×10¹⁷ | 0.525969 | 7698.20 | 3.156×10⁻⁷ |
| B2 | 1.24919×10¹⁹ | 20.7432 | 988256 | 1.245×10⁻⁵ |
| C13 | 1.82663×10¹⁹ | 30.3319 | 11207481 | 1.820×10⁻⁵ |

These averages are concentration scales, not the Bragg-peak concentration or trap occupancy. The average return is6.64551 mg/m²/year in the baryon convention,6.69751 mg/m²/year as neutral H; whole-sphere rates are3.20198×10⁶ and3.22704×10⁶ kg/s respectively. The companion stopping calculation now gives the following proton projected mean-depth proxies (quantiles over the energy distribution):

| Stage | Median depth, µm |99.9% depth, µm |
|---|---:|---:|
| A | 2.783 | 32.814 |
| B1 | 2.008 | 9.696 |
| B2 | 9.044 | 29.335 |
| C13 | 11.095 | 87.380 |

`stopping-snapshot.json` pins these values. They are evaluated stopping/depth estimates, not measured retention profiles or a longitudinal straggling distribution. Direct integration of the updated NIST proton stopping table differs from its printed range column by~7% near0.3 MeV and~1.2% at1 MeV; the companion packet preserves both. C13 therefore extends beyond the earlier40 µm reference, while most H stops within a few to~11 µm.

## A time-dependent model that can be used now

`calculate.py` provides conditional release kernels with prescribed temperature histories. The populations are disjoint initial cohorts: prompt release, mobile diffusion, a user-specified first-order trap, the audited local recombination/retrapping kernel, or held material. Their fractions are inputs. No default reactor mixture is invented. Released H is divided between cavity and local collection, then atomic H/H₂/CH₄, conserving H nuclei exactly and recording accompanying methane carbon.

For mobile material starting a distanceL from an absorbing connected plane, define Θ(t)=∫D[T(t)]dt. The one-dimensional half-space first-passage result is F(t)=erfc[L/(2√Θ)]. It assumes a spatially uniform scalarD, no retrapping and an absorbing boundary. A distribution of implantation depths can be integrated with this kernel. The half-space has an unbounded interior; a finite reflecting layer has different long-time behavior.

With Warrier’s **conditional surface-channel** coefficient,100 µm paths have50%/90% passage times3.15/90.6 s at1000 K,0.0968/2.79 s at1500 K. A10 µm path is100 times faster. This says that already-mobile H with access to a connected channel need not wait for bulk17.5 cm diffusion. It does not make grain-trapped H mobile.

The Muto kernel uses initial occupancyq=n₀/C₀ and surviving fractionN:

`Gq(N)=(1/q)(1/N−1)+2 ln N+q(1−N)=∫8.2 exp[−1 eV/(kT)]dt`.

At1000 K, the extrapolated half-release time is31.5 minutes forq=.95 and154 days forq=.001. At1500 K those values become39.5 seconds and3.23 days. The large difference is the calculated effect of retrapping/concentration, not an uncertainty interval for this reactor. The measured fit range is723–873 K; the high-temperature primary evidence above prevents treating the extrapolation as an operating prediction. Repeated pulses change occupancy, so independently superposing fixed-q kernels would be incorrect.

For an eventual coupled calculation, evolve mobile atomic/molecular concentrations and every trap occupancy explicitly: `dn_i/dt=kcap_i m(C_i−n_i)−kdet_i n_i`, with matching opposite terms in the mobile equation, diffusion, recombination and surface fluxes. The pore gas then has its own conductance and collector/cavity boundary conditions. TDS peak energies alone do not specify C_i, capture rates, prefactors, molecular branching or these boundaries. Physical uncertainty in those quantities is larger than numerical integration error.

The present code integrates piecewise-linearT with256-panel Simpson quadrature per interval and inverts the implicit survival relation by120 bisections. `validation.json` contains84 successful numerical/ledger checks, including1024-panel comparison through a1000→1700 K ramp. This validates the implementation and conservation, not the unknown reactor coefficients. Temperature histories must be material temperatures along the implanted path; a radiative effective temperature alone is insufficient.

## Conditional aligned-laminate screen

A50–100 µm aligned skin with collection channels at20–100 µm pitch offers lateral paths around10–50 µm. At1000 K, the ideal-crystal Herrero coefficient is2.43×10⁻⁹ m²/s. The same first-passage model then reaches90% over50 µm in~32.6 s. Only~1.03×10⁻¹¹ m²/s is needed for90% passage over50 µm during one B1 period. These comparisons justify investigating the geometry; they do not establish its collection fraction.

For a diffusion tensor, the competition between front escape and side collection depends on `x sqrt(Dparallel/Dnormal)/pitch`. The geometry packet evaluates this using actual implantation depths. No evidence-supported hot, damaged-layer numerical anisotropy was found. Pristine barriers, thermal-conductivity anisotropy and room-temperature barrier tests cannot supply it. Grain edges, vacancies, connected pores, cap seams and cracks can provide vertical paths. Conversely, long trapping delays need not change the eventual side/front probability if the same tensor and boundaries remain applicable; they change hold-up and timing.

A continuous1–5 µm cap over buried channels can allow energetic protons to cross into a closed plenum and stop at its floor while retaining thermalized gas. That is a two-part stopping/sealing geometry, whose joints and damage evolution need calculation. The full backing remains the final catcher. Channel opening, front leakage and trap-release efficiency should be varied separately, rather than replacing all three with one “recovery efficiency.”

## Targeted extraction and refurbishment

At the exact nucleus throughput, one eV of positive chemical upgrading free energy per H costs2.03182×10⁻⁵ W/m² ideally, or4.06364×10⁻⁵ W/m² input at50% thermodynamic efficiency. For example,10 eV/H at that efficiency is4.06×10⁻⁴ W/m². Activation barriers are rate parameters; they are not automatically irreversible work per recovered H. Chemical endpoints, coproducts, recoverable heat and separation determine net work. This honors the user’s molecular-manufacturing premise without assigning zero energy to an unspecified process.

By contrast, broad-layer handling after every pulse processes much more graphite:

| Treated depth each event | Mass per event, kg/m² | Annual handled mass, kg/m² | Power at1 MJ/kg input, W/m² |
|---|---:|---:|---:|
| 1 µm | 0.002 | 8.346 | 0.2645 |
| 10 µm | 0.020 | 83.46 | 2.645 |
| 100 µm | 0.200 | 834.6 | 26.45 |

These are replaceable processing scenarios, not required replacement or atomization. Selective access to traps, controlled gas recombination/extraction, damage-site repair, and less frequent surface restoration have different throughputs. Carbon carried as CHₓ can be recovered and returned; retained H is process inventory. If η is the efficiency of positive upgrading, use ΔG/η plus explicitly specified handling and separation loads. Do not charge the cohesive energy of all recycled graphite on every shot.

The highest-value next discrimination is whether the proposed irradiated skin routes most released gas into its collection channels while maintaining sufficient lateral mobility. Geometry and transport calculations can establish the required tensor, cap leakage and service rate now. Existing measurements do not determine the realized hot irradiated tensor or population fractions; those remain material properties to characterize. There is no basis here for a narrow80% credible interval on them. The reported scenario ranges are transparent conditions, not probability estimates.

## Reproduction

Run with ordinary Python3:

```sh
python3 calculate.py --out /tmp/hydrogen-kernels.json
python3 validate.py
python3 calculate.py --config example-mobile-cohort.json --out /tmp/hydrogen-history.json
```

The example is a normalized already-mobile cohort, not the whole incident H population. No article, previous snapshot, stopping result or plant ledger is modified by this packet.
