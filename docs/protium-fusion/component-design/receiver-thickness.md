# ✨ Receiver mass is a design variable

The existing 5,000 kg/m² graphite receiver is a conservative allocation, rather than an established minimum. A 2,000 kg/m² receiver remains a serious candidate after the present photon, thermal-storage and local-load screens. A 1,000 kg/m² receiver is also worth carrying as a sensitivity; its photon containment needs the tighter transport calculation. Neither alternative changes the fixed-denominator mixture comparison in `plant-scan-003`.

The reproducible thermal and mechanical calculation is [receiver_thickness.py](receiver_thickness.py), with [results](receiver-thickness.json) and an [independent validation packet](receiver-thickness-validation.json). Photon transport is handled separately in [photon-containment/README.md](photon-containment/README.md). All cases here keep the same radius, pulse profiles and cadence: A-H3, B1-H2, B2-H2.5 and C13-H4, at radius 1.1 billion metres and 4.4734×10³² helium nuclei per second.

## Heat storage and heat extraction

For the selected mixture, one complete A/B1/B2/C13 event delivers approximately 31.58/0.1795/49.65/59.18 MJ per square metre if all non-neutrino energy is collected. These totals include the prepared pulse energy as well as nuclear heat. They do not create a second electrical credit: the station ledger already debits preparation. The mass-conserving velocity interpolation used for thermal delivery differs from the exact nodal kinetic-energy ledger by about 0.1%; both conventions are recorded without rescaling the mass profile.

At 1,000 K graphite has a heat capacity close to 1,800 J/kg/K: the NIST-hosted thermodynamic table gives 21.610 J/mol/K. Historical screens used a constant 1,000 J/kg/K. Both are retained as controlled comparisons. Through-thickness conductivity is grade- and assembly-dependent; 25, 50 and 100 W/m/K are explicit scenarios. High-temperature graphite measurements establish substantial grade differences, so a room-temperature value cannot qualify an irradiated assembled receiver. [NIST thermodynamic table](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=921062), [high-temperature grade measurements](https://doi.org/10.1016/j.jnucmat.2020.152176).

The following comparison uses density 2,000 kg/m³, heat capacity 1,800 J/kg/K and conductivity 50 W/m/K. It assumes an intact receiver, distributed extraction and the existing one-day guard between large pulses.

| Graphite areal mass | 1,000 kg/m² | 2,000 kg/m² | 5,000 kg/m² |
| --- | ---: | ---: | ---: |
| Thickness | 0.5 m | 1 m | 2.5 m |
| Through-thickness diffusion time L²/α | 5 h | 20 h | 125 h |
| C13 energy / whole-slab heat capacity | 32.9 K | 16.4 K | 6.58 K |
| Conservative sum of all-stage bulk temperature swings | 78.1 K | 39.1 K | 15.6 K |
| Mean front-to-back gradient at full retained load | 1.26 K | 2.51 K | 6.28 K |
| C13 peak heat pickup with back fixed at 1,000 K | 6.08 kW/m² | 1.52 kW/m² | 0.243 kW/m² |
| Surface-temperature screen with finite receiver storage | 2,041 K | 2,002 K | 1,984 K |

The full-retention mean heat flux in this thermal model is 125.5 W/m². Consequently, holding the back face exactly at 1,000 K can demand much larger peak heat-exchanger capacity than the mean-rated converter package. The displayed back-face values are for surface deposition; volumetric gamma heating changes their timing. They are an illustration of the sink requirement, rather than a containment or exchanger qualification.

The alternative uses the receiver itself as storage and extracts the mean heat rate. For each periodic pulse family, the range of stored energy about its mean-removal trajectory is at most one complete pulse energy. Therefore the sum of these energies divided by the receiver heat capacity bounds the bulk-temperature swing without assuming a free heat bath. The calculation shifts each stage's bulk-energy contribution so its minimum is zero, placing the combined volume-average bulk above a selected 1,000 K floor. The actual pickup temperature differs from that average and must be matched to the variable-temperature engine cycle. A dispatch schedule can reduce that conservative all-phase excursion. Variable electricity production is another option, with its converter peak rating and source-store interaction accounted for explicitly.

The surface-temperature entries add a finite nonuniform conduction response to this bulk storage. They place all heat at the surface and preserve the separation between the microsecond prompt flash and the much later debris arrival. They do not add the prompt and debris maxima as though simultaneous. The radioactive heat calculation uses the actual endpoint parent inventory and a simple delivery-time envelope. Coupled positron/photon transport and local hot spots remain separate qualifications; those surface numbers are conditional thermal screens. The bulk energy-storage bound depends on total energy and remains useful when temporal delivery is spread differently.

Deep graphite does little to the early skin peak: the thermal disturbance has penetrated only a small fraction of even the 0.5 m slab at that time. Conductivity, photon deposition depth and a compliant renewable surface matter more. At conductivity 25 W/m/K the same conservative surface screens exceed 2,390 K at every thickness. Extra depth is therefore an inefficient remedy for poor surface heat transport. Thermal expansion, cracking and irradiation damage need a surface/joint design; routine displaced material remains inside the closed collection system for reuse.

## Photon containment is a separate requirement

The current spectrum contains substantial 511 keV annihilation radiation and 2.3126 MeV oxygen-14 daughter de-excitation radiation. An old uncollided-attenuation example at 20 MeV does not set the necessary thickness for that spectrum. The companion photon packet distinguishes uncollided removal, energy deposition and multiple scattering through the wall and cavity.

Its broad analytic local-cascade bound also puts all annihilation-plus-positron kinetic energy below 5.143 MeV and the saved grey prompt energy at 20 MeV. With a conservative 1.10 Doppler factor and the entire bounded direct hot-photon bypass additionally charged as lost, its closed-shell, inner-source and localized-secondary assumptions bound photon loss by about 4.457%, 0.8390% and 0.00859% of the complete heat boundary at 1,000/2,000/5,000 kg/m². The 2,000 kg/m² case leaves 0.16095% of that boundary for all other losses if total capture is to reach 99%. The tighter Compton transport sensitivity is encouraging for 1,000 kg/m², but omits other interactions and charged-secondary transport. The plant must use a justified total capture fraction rather than inherit 99% automatically. Radioactive processing stock outside the absorber requires its own enclosing geometry; moving a decaying source behind the shield changes the calculation.

## Local supports and capital material

Thinner slabs need closer local support. A simply supported 100 m square of 0.5 m graphite, using a selected 10 GPa elastic modulus, has a 9.7-second first bending-mode period. Under the selected pulse momentum, the absolute sum of all linear mode displacements has an impulsive bound of about 25 mm; the finite pulse duration generally reduces that response. The same slab spanning 1 km has a displacement bound around 2.5 m. This supports a subdivided receiver with local edge/subframe support. The edge support, joints and their connection to the guide deck still need an actual layout inside the service allocation.

These are local plate calculations. Global shell support still uses the guide/stream architecture, and predictable common motion belongs in its trajectory model. The plate screen does not qualify thermal-shock strength or supply an unsupported global shell.

Graphite capital stock at the three thicknesses is 1.52×10²², 3.04×10²² and 7.60×10²² kg. Reducing the receiver from 5,000 to 2,000 kg/m² directly removes 4.56×10²² kg of construction carbon before any dependent plant reallocations. That is a substantial resource and specific-power opportunity. The operational fresh-carbon-to-hydrogen ceiling remains a separate account; it does not constrain initial construction stock. The maintenance calculation must continue charging actual processing throughput and useful work while retaining recovered material. No blanket atomization charge is introduced.

## Checks and next decision

The analytic slab kernels agree with an independently assembled finite-volume model to 0.0062% at 1,024 cells, with convergence from coarser meshes. The back-face impulse kernel integrates to one; halving the actual pulse time step changes the slow surface peaks by at most 0.0413%. The finite late-tail screen retains the advected internal-energy t⁻² contribution separately from the kinetic t⁻⁴ contribution.

The next integrated receiver case should specify the collected photon fraction, receiver property curve, allowed bulk-temperature range, surface recovery architecture and support pitch together. The 2,000 kg/m² branch is sufficiently promising to carry forward as a conditional design. A thinner receiver should receive a new plant calculation only after those interfaces are stated; no higher W/kg result is assigned by this memo alone.

## Packaged 2,000 kg/m² candidate

[receiver-candidate-2000.json](receiver-candidate-2000.json) records the separate conditional candidate and links its thermal, photon and support inputs. The service allowance remains 1,000 kg/m², with the source and electrical store retained. The component comparison keeps the archived mixture and rate; it does not substitute a refined mesh, different mixture, or a newly optimized cadence. A 99% capture value remains a performance target subject to the photon packet's qualifications.

The lower mass changes the balance of self-gravity and solar tides. [receiver-support-family.json](receiver-support-family.json) evaluates the effect using the station's updated mass ledger at the same rate. At 2,000 kg/m², the conditional installed mass is 6.4785×10²² kg and a/n²=5.122, compared with 9.638 for the original mass. The present three-family stream construction requires a/n²>4, or total mass above 5.0594×10²² kg at this radius and orbit. All examined 30/100/300-day attitude-control cases remain non-growing in the nine-mode model with the retained 10⁸ N/m extensional mesh stiffness. The largest local stream-to-receiver mass ratio falls to 0.669 of its old value, within the range examined by the local guide-feedback screen. The fixed guide package is retained.

The 1,000 kg/m² comparison falls below that particular three-family threshold. Uniform membrane prestress of at least 0.221 MN/m would restore nonnegative family densities algebraically, adding 0.000402 Pa inward pressure. A positive design margin and a new geometric-stiffness/control check would be required; this is a support-layout change rather than a universal mass floor. The 2,000 kg/m² candidate needs no such change.
