# Titanium catalyst recycling audit ✨

The subsequent [complete allocation with explicit recycling](NETWORK-PLANT.md) fits these processors into the shared magnet and counts their kinetic, potential, radiation, extraction and converter costs. Its conditional 140 T reference gives 1.474 mW/kg after replacement. This page retains the standalone isotope/stock audit.

The principal titanium side products have a quantitatively testable return network. A sixteen-processor calculation using TENDL-2025 closes its isotope and atomic-mass energy ledgers while explicitly allowing a small waste stream and fresh titanium feed. **This replaces much of the earlier unspecified nuclear-restoration problem with named reactions. It has not yet replaced the 50 MeV side-cost assumption in the complete plant.** Auxiliary confinement, extraction, isotope separation and neutron transport still need to be integrated into that allocation.

At 1.6 MeV ion temperature, the network consumes 0.02561 neutrons in catalyst-restoring captures per primary Ti49(p,n) event and retains 1.26958 neutrons for light-fuel processing. Auxiliary proton reactions produce the additional neutrons. About 1.98% of the source's useful neutron absorptions go to catalyst restoration in the ideal capture model.

## Return paths

The main path remains Ti49(p,n)V49, followed by cold external electron capture back to Ti49. Neutrons captured on hydrogen make deuterium; the closed D cycle supplies helium. Cold neutral-atom decay rates are applied only after extraction and electron restoration.

Several side paths can be restored with low-energy neutrons:

- Ti49(p,d)Ti48 and Ti49(p,np)Ti48: Ti48(n,gamma)Ti49. The first produces deuterium directly. The second returns the emitted neutron to the catalyst and is an energy-circulation loss rather than an extra hydrogen-burning cycle.
- Ti49(p,2n)V48: cold decay to Ti48, then one neutron capture back to Ti49. One neutron remains available for fuel production.
- Ti49(p,alpha)Sc46: cold beta-minus decay to Ti46, followed by three neutron captures through Ti47 and Ti48 to Ti49.
- Ti49(p,gamma)V50: V50(n,gamma)V51, then V51(p,alpha)Ti48 and Ti48(n,gamma)Ti49. Competing V51(p,n)Cr51 events return through cold electron capture and make useful neutrons while the vanadium waits for its alpha-emission exit.

The vanadium processor itself has side products. Cr50 captures a neutron to Cr51, which decays back to V51. Cr52 can return through (p,alpha)V49 and cold decay to Ti49. Its proton-capture branch makes Mn53. **Mn53 is sent to a proton-processing channel:** Mn53(p,alpha)Cr50 feeds the Cr50 return path, and competing Mn53(p,n)Fe53 events return through Fe53's short cold decay. This avoids relying on Mn53's 3.7-million-year electron-capture lifetime for this branch. A similar treatment uses Ni59 as a processor when it appears farther along the network.

The calculation follows further chromium, manganese, iron, cobalt and nickel branches. It includes Mn54's small beta-minus branch to Fe54, whose subsequent neutron capture and cold Fe55 decay reach Mn55. The dominant Mn54 electron-capture branch reaches Cr54. NUBASE's branching entries are percentages; its `0.93e-4` beta-minus value corresponds to a probability of 9.3e-7.

These paths require separate isotope handling and suitably chosen ion channels. Their existence and balance do not establish the rates, losses or mass of the extraction hardware.

## Data and interpretation

The main input is the [TENDL-2025 proton library](https://tendl.ic.ac.uk/tendl_2025/tendl2025.html). Sixteen standard ENDF files, their source URLs and SHA256 hashes are archived in `beam-data/tendl2025-proton-recycling.zip`. TENDL supplies evaluated model data; these tables are not sixteen new measurements. The parser uses ENDF reaction definitions checked against the [ORNL SCALE documentation](https://scale-manual.ornl.gov/6.3.4/XSLibAppA.html).

The standard files give explicit reaction channels below 30 MeV laboratory energy and an MT5 remainder above that boundary. An earlier attempt to combine alternate-format MT5 cross sections with their product yields failed conservation checks. For example, its thermally weighted Ti49 heavy-product sum was 1.145 per MT5 event, and the Mn53 value was 1.267. Those results are retained only as a diagnostic in the research history; they were not used to repair or rescale the physical result. The standard files also leave small sum-rule discrepancies in the high-energy residual yields. The default calculation therefore discards the whole MT5 event, credits no useful products or nuclear heat from it, and reports a separate electrical reserve of 1000 MeV per such event. This is an explicit waste-and-makeup convention. It provides an accounting allowance for poorly resolved events; it does not assume hardware can identify every collision by its incoming energy.

The main titanium neutron-production rates agree reasonably across the inputs at 1.6 MeV:

| Evaluation | Ti49 total pn reactivity, m3/s |
|---|---:|
| PADF-2, 2026 | 3.74718e-22 |
| TENDL-2023 | 3.66204e-22 |
| TENDL-2025 | 3.66305e-22 |

The side branches have larger differences. Ti49 proton capture relative to pn increases from 0.001046 in TENDL-2023 to 0.001544 in TENDL-2025. V51 alpha emission relative to pn decreases from 0.009378 to 0.005505. These changes materially affect the required auxiliary processing. Cross-version differences are sensitivity evidence, not a calibrated error bar.

For Cr52, the PADF ground-state Mn52 yield alone is 1.01500e-23 m3/s. Including its metastable state gives 6.23564e-23, close to the TENDL-2023 total of 6.26418e-23. Comparing a ground-state yield to a total would falsely suggest a sixfold disagreement. Mn52m's beta decay and internal transition both ultimately return Cr52.

## Steady flows and energy ledger

The solver fixes one primary Ti49(p,n) event and solves the stationary heavy-isotope balance for the other processors. Their proton reaction rates use a common specified temperature; their density, confinement volume and equipment mass have not yet been allocated.

| Processor | pn events per primary Ti49 pn, at 1.6 MeV |
|---|---:|
| Ti49 | 1 |
| V51 | 0.245069 |
| Cr52 | 0.0155342 |
| Mn55 | 0.00827185 |
| Mn53 | 0.00459181 |
| Cr54 | 0.00128778 |
| Fe56 | 0.000422591 |
| Fe57 | 0.000243186 |
| Fe58 | 0.0000346881 |
| Remaining selected processors | below 0.000003 each |

The auxiliary total is 0.27546 pn events per primary event. A larger reaction count is not a free improvement: it burns additional hydrogen and requires additional ion-processing volume. Some processors have substantially lower reactivity than Ti49.

The downstream light-particle ledger includes neutron capture on hydrogen, net 2D-to-He4 burning, tritium completion through DT plus neutron recapture, and He3 neutron capture. It supplies no external D, T, He3 or manufactured neutrons. At the reference temperature:

| Quantity per primary pn | Value |
|---|---:|
| Net hydrogen atoms consumed, including tiny waste input | 2.555858 |
| Helium-4 atoms produced | 0.638964 |
| Catalyst-restoring neutron captures | 0.025614 |
| Net source neutrons available for light fuels | 1.269576 |
| High-energy events discarded | 1.60707e-6 |
| Fresh Ti49 atoms required | 1.60773e-6 |
| Neutrino-energy allowance | 0.903471 MeV |
| Nuclear heat after excluding positive makeup-transmutation credit | 16.17665 MeV |
| Additional 1000 MeV/discard-event electrical reserve | 0.001607 MeV |

The neutrino allowance assigns all ground-state atomic weak-decay Q values to escaping neutrinos. It adds 0.005892 MeV per primary event as a conservative excitation allowance for known weakly decaying isomers, assuming their most energetic listed state is populated whenever applicable. Actual decay branching and gamma/charged-particle energy spectra still require a fuller energy-partition calculation. IT-only isomers return excitation electromagnetically. Cold inventory estimates presently use ground-state lifetimes.

The 16.17665 MeV number is **nuclear heat**, not electricity or surplus after ion preparation. It excludes the tiny positive energy contribution of converting makeup titanium into discarded heavy products. The main-cycle energy remains hydrogen-to-helium conversion. The model's source heating, rotation, potential handling, end barriers, D equipment and conversion efficiencies must be recomputed for the whole network before comparing with 1 mW/kg.

At 1e24 primary pn events/s, the conservative high-energy purge requires 4.12 kg/year of enriched Ti49 makeup. At the earlier 140 T plant's primary rate, it would require 15.56 kg/year. Trace copper-63 is the largest explicitly identified boundary residue, at about 6.6e-10 atoms per primary event; high-energy disposal dominates the makeup requirement. The spent material and replacement feed remain in the ledger. The network therefore makes no claim of mathematically zero catalyst consumption.

At 1.0 MeV, auxiliary processing is much larger and the copper boundary reaches about 5.8e-4 per primary event. At 2.2 MeV, the high-energy discard rises to 9.6e-5. The preferred temperature cannot be chosen from nuclear heat alone.

## Cold neutron-capture inventory

`capture_bank_screen.py` estimates isotope stock from N=R_capture/(phi_2200 sigma_2200). Here phi_2200 is the 2200 m/s equivalent thermal-neutron flux, using a 1/v-style screening convention. It is not automatically equal to the scalar flux of a Maxwellian spectrum.

Most cross sections come from the [NIST NCNR absorption table](https://www.ncnr.nist.gov/resources/n-lengths/list.html). Assigning total absorption to the desired capture makes the resulting stock a lower bound; competing reactions require transport and depletion calculations. For V50, the screen uses the measured radiative-capture cross section of 21(+4,-2) barn from [Michaelsen, Lieb and Robinson (1991)](https://doi.org/10.1007/BF01295765). The TENDL-2023 pointwise activation file gives 40.70 barn, while the older NIST total-absorption entry is 60(40) barn. These are not interchangeable precision measurements. A raw resonance ENDF file's MF3 background alone is also insufficient to obtain the thermal capture value.

Using the earlier 140 T plant's rate and water mass purely as a scale gives phi_2200=1.167e14 m-2 s-1 and about 1.905e11 kg of listed capture stock. That is 0.0140% of the earlier plant mass. Ti46, Ti48 and Ti47 dominate. Their inferred mean individual capture times range from about 0.35 to 4.6 million years; the smaller Ti50 inventory has a roughly 15-million-year mean capture time. Initial loading of stable isotope stocks can establish their required capture rates immediately. A long residence time does not by itself require waiting that long to operate.

The stock screen covers all but 1.21e-8 capture steps per primary pn. Its homogeneous neutron balance omits spatial flux depressions, resonance self-shielding, fast-neutron damage, parasitic absorptions, chemical carriers and irradiation of already-produced return isotopes. The auxiliary ion channels have not been fitted into the earlier plant geometry. This is an inventory scale, not a transport result or a revised plant mass.

For orientation, processing the listed isotope inventory with a 1000-year mean turnover and an assumed 100 MJ/kg electrical cost would require about 6.04 kg/s and 604 MW at that scale. The corresponding simple competing-rate estimate gives a Ti49 recapture probability of 0.00081 before extraction. These are sensitivity parameters, not an isotope-separation process design; carrier-fluid processing and equipment remain open.

## Reproduction and next discriminator

Run `catalyst_network.py`, `capture_bank_screen.py`, `compare_recycling_data.py`, and `validate_recycling.py` with the research dependencies plus `endf-parserpy==0.18.0`. The archived data need no network access. `recycling-comparison-input.json.gz` contains the selected earlier MF3 tables with source hashes and URLs. The raw-yield diagnostic is separate from the default disposal cases.

The 172 checks compare the small MF3 reader against the IAEA parser, Gaussian quadrature against independent incomplete-gamma integrals, stationary isotope balances, fresh-feed versus waste counts, and an independent atomic-mass energy ledger. Across the five temperature cases, baryon residuals are below 1e-15 per primary event and energy-ledger residuals below 1.1e-11 MeV. This verifies the stated calculation, not the extraction or neutron-capture engineering.

Next, allocate the auxiliary ion populations and their extraction/heating/rotation costs into the shared solenoid, replacing the generic 50 MeV parameter. Then model finite extraction and electron clearing, followed by multigroup neutron/gamma transport and conductor damage. No article chronology has been selected, and this checkpoint alone does not update the earlier conditional 1.37 mW/kg result.
