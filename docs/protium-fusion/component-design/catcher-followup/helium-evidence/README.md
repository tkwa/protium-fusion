# Helium release and local extraction evidence ✨

The present evidence supports designing for a mobile helium population that can return rapidly to the exposed graphite surface. It does not establish the fraction remaining trapped in the proposed irradiated receiver. A 17.5 cm solid plate does not automatically deliver implanted helium to its rear plenum. A credible next design is a characterized near-surface diffusion/channel layer with separately measured front leakage, feeding isolated collection and compression equipment.

This packet preserves the exact delivered-B1/C13H2 source128 family at radius1.1e9 m and receiver areal mass350 kg/m². [The spectrum interface](../stopping/surface-spectra.json) gives B1 helium2.10153e17 atoms/m² =1.39587 µg/m² per pulse, a7698.20 s mean period, and number-median incident energy0.6517 MeV. The other stages contribute negligible helium by comparison. Those are modeled incident quantities; stopping, reflection and subsequent release are separate operations.

## Primary evidence and access

|Source|Conditions and observation|What it constrains|
|---|---|---|
|[Li et al2019](https://doi.org/10.1080/00223131.2019.1624654), full article body read|ETU-10/15 densities1750/1840 kg/m³, porosities15/11%;200 keV He,1e19–1e21 m⁻²,implant below373 K. Overnight400 K bake, then0.5 K/s TDS to1200 K. Release starts above400 K, peaks near550 K; signal negligible above800 K. Vacancy clusters persist above1200 K.|Releasing and surviving defect populations differ. No calibrated absolute retained fraction was obtained. The proposed hot mixed-ion history differs.|
|[Jung1992](https://doi.org/10.1016/S0022-3115(09)80070-6), abstract only|Pyrolytic graphite/SiC,5–28 MeV alpha irradiation,typically1e17 m⁻²,up to1100 K. Highest graphite release for implantation parallel to basal planes. Reported mobile-reference D=5e−5exp[−1.17 eV/kT]m²/s over750–1050 K; SiC1.1e−6exp[−1.14 eV/kT] over800–1050 K.|The coefficient supplies a conditional mobile-path calculation. A slower retained component is reported but unquantified. No measured diffusivity tensor or bulk retained fraction is assigned.|
|[Alimov et al1995](https://doi.org/10.1063/1.360664), abstract/metadata only|Graphite orientation, intergranular pathways and blistering affect reemission.|Motivates orientation-specific testing. Quantitative retention results are deferred until full text is obtained.|

Li's methods, results and discussion were read through the indexed publisher full article; Table1 was read in the indexed PDF. The direct publisher PDF returned403 and the Delft repository timed out. Figure6's image/calibrated integral remains unavailable. Jung's full paper was not obtainable through the public publisher/institutional searches; OpenAlex reports no open copy. Alimov's public MPG record contains metadata without an attached manuscript. These access limits are explicit in [evidence.json](evidence.json); no retention fraction has been reconstructed from an abstract.

## Time-dependent mobile-population calculation

For a helium atom initially at depthx in a homogeneous half-space with an absorbing exposed surface, define A(t)=∫D[T(t)]dt. The image solution to the diffusion equation gives released fraction

`F(x,t)=erfc[x/(2 sqrt(A(t)))]`.

The median first-passage time at constantD is1.0991 x²/D. Its long tail matters: the mean first-passage time in the unbounded half-space diverges. A finite trapping/collection structure removes that idealized long tail. `helium_release.py` implements this kernel and does not select the mobile fraction. If a fractionf is mobile and the remainder remains trapped over the interval, total release isfF; detrapping adds a time convolution. Temperature-dependent trapping, radiation detrapping, cracks and pore transport require separate measured inputs.

At1000 K the graphite reference coefficient is6.34497e−11 m²/s:

|Depth|Median surface-arrival delay|Time for95% of this mobile population to emerge|Fraction emerged after7698 s|
|---:|---:|---:|---:|
|0.1 µm|0.173 ms|0.0200 s|99.9919%|
|1 µm|0.0173 s|2.00 s|99.9193%|
|10 µm|1.73 s|200 s|99.1927%|
|100 µm|173 s|20041 s|91.9411%|

These are conditional mathematical results, not measured receiver escape fractions. For100 µm paths, Arrhenius continuation gives median1.88 s at1500 K and0.754 s at1668 K. Both temperatures exceed the coefficient's measured domain. The JSON also reports a sensitivity holdingD fixed at1050 K; neither extrapolation is a calibrated high-temperature material model. A measured receiver temperature history can be inserted throughA(t) without changing the kernel.

Keeping95% of mobile helium at1 µm depth for one mean B1 period would instead require effectiveD≤1.69e−17 m²/s,3.75million times below the reference. At10 µm the corresponding bound is1.69e−15 m²/s. These are storage-design requirements, not evidence that all graphite fails them.

A single thermally activated trap with attempt frequency1e12–1e13 s⁻¹ would need binding/escape barrier3.41–3.61 eV for95% survival over7698 s at1000 K. Surviving a100 s hold at1668 K requires5.06–5.39 eV under the same first-order model. Neither the trap spectrum nor the mobile fraction has been measured for this catcher. The1.3 eV Frenkel-recombination barrier discussed by Li is not a measured helium-detrapping barrier and is not used here.

### Exact incident-spectrum convolution

[helium-spectrum.json](helium-spectrum.json) applies the same kernels to the completed [stopping interface](../stopping/surface-stopping.json). B1 He projected-depth proxies have1st/50th/99.9th number percentiles0.266/1.942/5.264 µm. For the mobile reference at1000 K, the ensemble50/90/99% post-implantation delay is0.0653/2.118/213.8 s. Convolving actual arrival times gives cumulative50/90/99% release200/367/977 s after ignition. By the next mean B1 ignition,99.8302% of this incident mobile population has emerged;1.10e−5 of the incident population has not yet arrived. MultiplyingD by1e−3 or1e−6 changes that release fraction to94.6705% or13.1235%, respectively. Those multipliers are material sensitivities, not credible intervals.

Each packet uses the evaluated mean projected range at its energy, not a stochastic implantation profile. Both erfc(x/c) and the strip's positive exponential series are convex for positive depth. Thus, within these idealized boundaries, unresolved depth straggling increases average front release relative to substituting the mean depth. The current result is an optimistic front-retention screen; prompt backscatter is also separate. The stopping agent independently confirmed this mean-depth interpretation.

## Where the released atoms go

For isotropic diffusion in a slab with perfectly absorbing front and rear faces, eventual rear-exit probability isx/L. A10 µm implant inL=0.175 m has only5.71e−5 rear probability; most exits are through the nearby front. This boundary-value result establishes why bulk through-wall diffusion cannot be credited as local pumping. Gas enters whichever connected surface/pore is nearest in the actual transport tensor; release into a sealed extraction channel is recoverable inventory, while release to the reactor cavity becomes a chamber-pumping load.

The mechanical agent proposes a50–100 µm aligned graphite skin, basal planes parallel to the exposed front, with side-connected collection channels at20–100 µm pitch beneath a continuous cap. For uniform lateral implantation, perfect side absorption and a semi-infinite depth, the front-exit probability is

`Pfront=(8/pi²) sum(n odd) exp[-n*pi*x*sqrt(eta)/p]/n²`, where `eta=Dparallel/Dnormal`.

A collecting rear boundary reduces this probability. Finite channel sticking, return from the plenum, damaged seams and cracks change the boundary conditions. Atx=1 µm andp=20 µm, front probability is77.4%,16.9%,0.564% for eta1,100,1000 respectively. Targets5%,1%,0.1% require eta315,783,1818. A tenfold shallower path requires100times larger anisotropy. These percentages are selectable evaluation thresholds, not acceptance requirements imposed on the reactor.

At the mobile in-plane referenceD, mean side-exit time in the strip ignoring front loss isp²/(12D):0.525 s at20 µm pitch and13.1 s at100 µm. There is no retrieved measurement establishing eta1000 after the proposed repeated mixed-ion exposure. The in-plane interpretation of Jung's coefficient is an inference from the stated implantation orientation, not a verified tensor decomposition. Therefore the channel calculation supplies a property specification to test, not demonstrated collection efficiency.

With the complete B1 mean-depth proxy distribution,20 µm pitch requires eta186/982/11020 to hold eventual front leakage below5/1/0.1%;100 µm pitch requires25times those values. At eta1000 and20 µm pitch the optimistic front fraction is0.982%. The shallow tail therefore matters beyond the median implant depth.

The strip calculation places absorbing sides all the way to the exposed front. Buried channels under a continuous cap need a separate cap-domain solution. Before reaching the channel depthh, an atom implanted atx<h has front-exit probability at least1−x/h for homogeneous normal diffusion, even if the entire bottom of the cap were perfectly collecting. Weighting the actual B1 depth proxies gives lower bounds0.00176%,0.0662%,4.963% for cap depths10 nm,100 nm,1 µm. These lower bounds are not additive to the strip loss without solving the connected domain. The current mechanical candidate uses100 nm cap depth,100 nm channel-mouth width,20 µm pitch and150 µm collection-layer depth; this is a nanofabricated design proposal. Finite side capture, cap damage and gas backflow remain explicit boundary conditions.

There is primary support for the basal-barrier mechanism itself. [Bunch et al2008](https://sites.bu.edu/bunch/files/2015/06/Impermeable-Graphene-Membrane.pdf) bounded thermal-He transmission through micron-scale graphene devices near room temperature below1e−11 per incidence; observed device leakage was consistent with glass/interface paths. [Leenaerts et al2008](https://arxiv.org/abs/0810.4056) calculated GGA monolayer penetration barriers of11.69 eV for pristine graphene,6.12 eV for Stone–Wales defects,5.75/3.35 eV for two divacancy geometries, and1.04/0.44/0.10 eV for tetra/hexa/decavacancies. Both full papers were retrieved; the barrier table was visually checked. These results make directional collection physically plausible while identifying connected larger defects as a specific danger. They do not measure the transport tensor, crack connectivity or lifetime of the proposed hot irradiated multilayer skin. A pristine-sheet barrier cannot be used as its effective diffusivity.

## Helium collection after release

[Day et al1999](https://www.osti.gov/etdeweb/servlets/purl/20047879) tested approximately1 mm activated-charcoal panels near5 K. Initial helium sticking fell from0.64 at6 K to0.17 at7.5 K. Figure2 gives roughly13 m³/(s m²) at low pure-He loading, then strong loading-dependent decline. At6.5 K the reported helium equilibrium pressure is around1e−6 Pa;85 K regeneration removes the helium capacity. The full paper and Figures1–2 were inspected. A loaded panel at that equilibrium pressure cannot be treated as a zero-pressure sink for a cavity targeting much lower pressures.

[Day et al2004](https://scipub.euro-fusion.org/wp-content/uploads/2014/11/EFDC040504.pdf) documents dedicated coolant/heater circuits, temperature sensing, conductance corrections and transfer pumps during regeneration. Its large hydrogen/tritium capacities are not transferred to helium. [Hauer/Day's FZKA6745 report](https://publikationen.bibliothek.kit.edu/270052792/3814277) supplies4.21–40 K He isotherms chiefly at much higher pressure; its non-Langmuir behavior precludes an unsupported extrapolation to the reactor's extremely dilute cavity.

A side-channel plenum can operate at higher pressure, with shielded cold collection, staged compression and valves isolating regeneration. A continuously exposed final sorption stage needs its own low-coverage equilibrium and conductance model. The proposed geometry does not yet choose the pumping train.

Own thermodynamic/throughput calculations, using the incident B1 helium inventory:

- One B1 pulse supplies7.925e−4 Pa m³ of helium at the273.15 K throughput convention per receiver m². With cryopanel area1e−4 of receiver area, that is7.93 Pa m³ per panel m² per shot. An illustrative1000 Pa m³/m² loading is reached after126shots. This illustrates inventory capacity only; it does not establish low-pressure performance.
- Ideal constant-pressure cooling of200 K helium to5 K removes1.013 MJ/kg; its minimum reversible work with a200 K heat sink is2.821 MJ/kg. Adsorption heat, compression, regeneration and finite machine efficiency remain additional loads.
- A5 K panel with a black200 K view receives90.7 W/m² of panel area. Its reversible refrigeration power alone is3538 W/m² of panel area; at area fraction1e−4 this is0.354 W per receiver m². A direct1000 K view instead adds56.7 kW/m² at the cold panel. Multistage baffles and cooling are therefore part of the collector, not optional free surroundings.

## Decisive next material tests

1. Measure absolute retained and front/side-released He inventories during the actual energy distribution and1000 K irradiation, followed by the modeled temperature excursion. Include the C/N/H irradiation that creates and modifies traps; match cumulative damage as well as He dose.
2. Measure in-plane and cross-plane transport and channel collection on the proposed capped lamellar skin, including seams, evolving porosity, prompt backscatter and desorption. An isolated pristine-crystal diffusivity is insufficient for that boundary condition.
3. Record front leakage versus cycle count, temperature and regeneration while demonstrating an inventory-closing extraction train. Refurbishment must return the retained He and implanted catalyst as ordinary process stock.

A numerical80% credible interval for total hot-wall retention or eta would be fabricated from the present evidence. The useful conditional results are the calculated transport/property requirements above. The proposed material and its defect history need to be specified before a meaningful statistical estimate can be assigned.

## Reproduction

Run `python helium_release.py` and `python helium_spectrum.py` with numpy/scipy. [helium-release.json](helium-release.json) records input/code hashes and57 checks; the spectrum convolution adds20 checks. These include the analytic median, exact strip boundary, independent series and monotonicity. Raw retrieved files remain under `/home/tkwa/src/protium-fusion/output/catcher-followup/helium-evidence/`; no earlier model packet was changed.
