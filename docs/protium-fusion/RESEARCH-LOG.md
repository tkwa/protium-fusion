# Protium fusion research log ✨

Status: exploratory calculations, 28 October 2026. No reactor has been validated and no new future-history date has been selected.

## Scope and decisions

Thomas asked the present agent to take this article out of the automatic writer's queue and write it personally. The manual reservation is live. He requested high scientific accuracy, an investigation of bremsstrahlung mitigation, economic implications, a system smaller than approximately one Jupiter mass, and preferably at least 1 mW/kg. The investigation uses net useful power divided by complete installed system mass, including operating inventories and heat rejection. The ten-million-token allowance is a ceiling.

The natural hydrogen-burning stellar mass limit does not establish an engineered reactor's minimum mass. Conversely, optical thickness or a positive local reaction rate does not establish net useful power. Both mistakes must be avoided.

## Existing canon

- The Venus and starlifting briefs place the early 2040s before practical protium fusion. No subsequent breakthrough date is fixed.
- Eta Carinae's 10480–11360 CE programme uses the target core's gravity to support a very large hot-CNO apparatus. Its rejection of freestanding burners concerned that programme's throughput; it is not a universal impossibility theorem.
- Mature black-hole conversion is a competing technology in the far-future economy. Compare delivered energy, infrastructure, accessible feedstock and recoverable helium, not fuel abundance alone.

## Calculations currently present

1. `nuclear_screen.py`: steady closed CN-I/hot-CN-I catalyst residence times, reaction power and elementary EOS/radiation diagnostics. At high temperature the two hot-cycle beta waits total approximately 278 seconds. The resulting specific power agrees approximately with the published beta-limited value, 5.6e11 times catalyst mass fraction W/kg. The code excludes leakage and degeneracy corrections and is not a confinement model.
2. `beam_screen.py`: measured narrow 14N(p,gamma) resonances versus neutral-nitrogen stopping. Even granting a whole CNO cycle for every capture, stopping exceeds fusion credit by roughly 0.5–3 million at the examined peaks. This addresses these cold-target resonances only. The quoted natural-width reference frame needs verification; a 7% frame difference cannot change the order-of-magnitude conclusion.
3. `opaque_screen.py`: classical free-free plus electron-scattering Rosseland opacity, a dimensional diffusion size, pressure and deliberately optimistic confinement mass floors. Opacity can be obtained far below Jupiter mass. Required pressure, stability, boundaries and driver work remain unresolved. Its ideal-EOS region must be enforced.
4. `staged_capture_screen.py`: two prescribed-state finite-fuel pulses, separated by cooling and beta decay. Stage A is 14N+p -> 15O, followed by cold decay to 15N. Stage B is 15N+3p -> 14O+alpha, followed by cold decay to 14N. This relocates the long weak-decay waits outside the hot plasma. It does not eliminate their time or catalyst inventory. No implosion, expansion, driver or separation model is included yet.

## Important emerging distinction

At 300–500 MK, a dense staged pulse can have a capture time much shorter than a beta-decay time. Preliminary fixed-state thermal-energy ratios and sound-crossing masses therefore justify further investigation. They do not establish ignition, target gain or plant gain. The energy deposited during burning substantially changes temperature, pressure and reaction branching. Follow the complete reacting material, including partially processed catalyst; do not credit its energy while silently discarding its regeneration cost.

The first alpha-leak rate check found 13N(alpha,p)16O / 13N(p,gamma)14O rate-coefficient ratios approximately 3.8e-5 at 300 MK, 7.3e-4 at 500 MK and 8.0e-3 at 800 MK, before abundance weighting. 14O(alpha,p)17F also grows rapidly. A larger network and self-heating calculation are required before selecting a temperature.

## Other live branches

### Expansion-model result and revision

`pulse_model.py` adds a 27-isotope, 194-rate network to a homogeneous homologously expanding sphere. It uses the scalar virial acceleration and PdV work, includes self-gravity, and starts from an already compressed hot state. Its radiation escape prescription is approximate diffusion with a fixed opacity. It is not a spatial radiation-hydrodynamics or implosion calculation.

The initial combined stage B fails its intended chemical output at rho=1e8 kg/m3, initial T=300 MK and radius=1 km. Self-heating reaches about 1.55 GK. Only about 8% of catalyst ends in the desired 14O state, while most ends as 16O. Nuclear energy alone therefore overstated useful catalytic performance. This is a substantive rejection of that particular two-pulse parameter choice.

A revised three-pulse loop separates 15N+p -> 12C+alpha (B1) from 12C+2p -> 14O (B2), removing helium between them. Preliminary results preserve approximately 99% or more of the desired output in each of these individual idealized pulses. Their energy release divided by initial internal energy is approximately 13–22, depending on stage and parameters. These figures exclude driver losses and complete recovery/reprocessing. Minor unwanted nuclei require an explicit recovery ledger. They are candidates for further modeling, not a validated plant.

`fast_eos.py` tabulates the pynucastro electron/positron EOS and evaluates ions and radiation analytically. A first 20-point comparison gives maximum relative errors of about 3.3e-6 in pressure, 2.6e-6 in internal energy, 4.2e-5 in heat capacity and 1.2e-4 in adiabatic index. Integration energy residuals must be inspected separately. The rate evaluator reproduces the library at sampled temperatures to relative error below 4e-12 and exactly conserves nucleon number in every included reaction.

- Nonneutral accelerator cycles: 13C(p,n)13N followed by beta decay restores the carbon isotope and converts a proton to a neutron. Subsequent neutron capture and a closed deuterium-burning network can in principle complete hydrogen-to-helium conversion. Verify all cross sections, ion scattering, space charge, electron contamination, catalyst-destroying channels, neutron capture, energy recovery and plant mass. Neutral-target stopping results cannot be reused as a universal bound on this architecture.
- Important 13C side channels include (p,d)12C, (p,alpha)10B, (p,pn)12C and inelastic excitation. A neutron emitted by destroying a catalyst nucleus is not automatically a catalytic conversion of hydrogen.
- A 2023 centrifugal-electrostatic confinement calculation provides a primary-source example of nonneutral fusion theory, for DT. Its radial model and assumed axial confinement do not validate a CNO machine.
- Steady magnetic/self-pinched opaque burning, low-temperature cold CNO, pure pp, and small-black-hole-assisted confinement remain incomplete. The last branch requires an explicit seed-manufacture assumption and separation of fusion from accretion power.

## Antecedent attribution

A 2020 deepfuturetech essay reproduced by Atomic Rockets discusses both staged CNO burning and nonneutral (p,n) cycles. These are useful idea antecedents, not physical validation. Its claims about arbitrary thick-shell stress reduction and confinement must be independently checked. Record any adopted idea's antecedent and derive its feasibility from primary nuclear/plasma sources.

## Next discriminating checks

1. Enlarge the staged network, include leakage/reverse reactions, and test screening sensitivity.
2. Couple nuclear energy release to a consistent EOS and expansion. Check optical depth and diffusion throughout the pulse.
3. Build a closed repeated-cycle isotope and energy ledger, including incomplete burn, cooling, separation, reheating and catalyst loss.
4. Derive driver and capture requirements and a complete plant-mass budget rather than extrapolating local specific power.
5. Obtain quantitative 13C beam-channel data and compare ion-only scattering/thermalization with nuclear yield.
6. Evaluate economic competitiveness and only then select any future-history chronology.

## Spatial expansion, catalyst processing and new branches

The spatial model and its verification are documented in `RADIAL-MODEL.md`. The main nitrogen pulse's bulk conversion converges near 92% at 300 m radius and rho=1e8 kg/m3. It starts compressed; successful expansion is not evidence of an achievable implosion.

Splitting B1 onto a lower-density trajectory (rho=1e6 kg/m3, 200 MK, radius 1500 m) lowers its peak to about 443 MK. The uniform model diverts approximately 1.2e-4 of catalyst into oxygen and much less into heavier nuclei. The 64-cell calculation gives 99.17% desired carbon output; its rare branches require tighter absolute tolerances.

O16 recovery has a possible proton-capture route through F17 and Ne18. Cold decays lead to O17 and O18; subsequent proton/alpha reactions can return N14 or C12. Dense O16+2p pulses at rho=1e10 kg/m3, 500 MK and radius 30 m convert about 49% to Ne18 in the uniform model. The 100 m case hit the 2 GK guard before expansion completed and must not be treated as a completed pulse. Stage-specific temperature/density choices are important: a dense O18+p pulse produces substantial unwanted heavy nuclei, whereas a dilute O18+2p pulse favors C12+2alpha.

`two_proton_screen.py` examines omitted O15(2p,gamma)Ne17 and Ne18(2p,gamma)Mg20. Applying the next hotter published upper tabulated rate throughout the whole pulse gives leakage estimates below 3e-7 for the main nitrogen case and about 0.0015 for the dense cleanup case. These are screens within the cited rate envelope, with conservative constant density/proton abundance, not a complete expanded network. The 2018 O15 study narrows the older uncertainty interval at low temperatures.

A simpler auxiliary hypothesis is to manufacture replacement C12 from a fraction of helium ash, rather than regenerate every leaked isotope. A pure-helium uniform expansion at rho=1e10 kg/m3 and 1 GK gives C12 mass fractions 0.00191 and 0.00654 for 30 m and 100 m initial radii. The larger case requires about 0.41 GeV of initial fuel energy per produced carbon atom. This is unfavorable as a standalone power source, but a catalyst loss fraction of 1e-4 would charge only about 0.041 MeV of initial fuel energy per main cycle, before driver/recovery efficiencies. Recycled unburned helium and manufacturing infrastructure must be counted. The 30 m case has much smaller pulse inventory but higher energy cost per carbon atom. No spatial manufacturing or compression calculation has yet verified these screens.

The next decisive question is how to compress the fuel with acceptable input energy and entropy. A prescribed compressed initial state does not answer it. A spatial implosion calculation should include initial inflow energy or external boundary work, nuclear reactions during compression, and nonuniform stagnation. Direct inductive recovery of expanding plasma is a published concept for other fuels, but its efficiency and geometry cannot simply be assumed for this plant.

## Compression checkpoint

A controlled convergent-flow calculation now includes inward kinetic energy in the input ledger. Its initial velocity/density/temperature profiles are constructed by reversing a nonreacting adiabatic expansion. The 512-cell nitrogen-stage result burns about 83% of catalyst and releases nuclear energy 12.63 times the initial thermal-plus-kinetic energy, with a relative total-energy residual below 1e-6. This is progress beyond starting directly at maximum compression. It remains a prescribed inflow, whose preparation, three-dimensional stability and driver efficiency are not demonstrated. See `RADIAL-MODEL.md` and `compression-summary.json`.

Screening diagnostics from the installed Chugunov–DeWitt implementation give forward-capture enhancement factors about 1.11 for N14+p at 300 MK, 1.02 for dilute N15+p at 200 MK, 1.09 for C12+p at 300 MK and 1.46 for the dense O16+p cleanup at 500 MK. These are initial-state values. The library returns the natural logarithm of enhancement; `screening_screen.py` explicitly exponentiates it. No screened integration is claimed.


## Closed-cycle and collector checkpoint, 29 October

The full-network 128-cell A compression calculation completed successfully: conversion 0.811296, matching the reduced result to a few parts per million. It took about 978 seconds with the earlier abundance-clipped rate evaluator. Removing clipping of tiny negative Newton iterates in the polynomial mass-action rates substantially improved stiff solver progress; physical positive-abundance rates are unchanged, and accepted-state negativity remains at numerical tolerance. The 64-cell full result gives conversion 0.786403 and energy residual about 8.5e-6 of input.

B1 compression gives desired C12 output 0.760645 at 64 cells and 0.801618 at 128 cells. B2 gives O14 output 0.880011 and 0.908248. These differ from earlier uniform or initially compressed expansion models. B1 oxygen leakage remains around 1e-4 per processed catalyst. Further convergence and rare-branch tolerance checks are warranted.

A low-density C13 cleanup trial (200 MK, rho=1e6 kg/m3, radius 1500 m) barely burns and is retained as a failed choice. The higher-density 300 MK, rho=1e8 kg/m3, radius 100 m convergent-flow case gives N14 output about 0.83185 and nuclear/input-energy ratio about 16.0. Its energy residual is about 5.1e-5 of input.

The spatial helium factory at rho=1e10 kg/m3 and 1 GK gives carbon mass fraction 0.003496 for radius 30 m and 0.011567 for radius 100 m at 64 cells. The 100 m, 128-cell result gives 0.013691 and a higher peak temperature; neither yield nor local peak is fully converged. The factory is an auxiliary energy consumer. Its complete-cycle contribution is small only when catalyst losses are small.

`cycle_ledger.py` now solves a stationary four-isotope processing balance, including cold decay, unburned returns, heavier ash and replacement carbon from helium. All default stages use specified convergent-flow results. The older uniform-auxiliary hybrid is preserved as `cycle-ledger-hybrid.json`. The current ideal-collection case consumes about 4.002 ordinary-hydrogen atoms per net helium atom, gives available nuclear heat about 24.68 MeV and requires prepared fuel energy about 1.78 MeV. Baryon and atomic-mass-defect energy residuals are at floating-point precision. With assumed 40% thermal recovery and 30% driver coupling it leaves about 4.66 MeV per helium before other loads.

Cold N13 produces C13, which may share the carbon mixture rather than need an isotope-pure cleanup stream. A 64-cell B2 calculation with initial carbon isotope fraction 0.0411 C13 nearly matches its own steady return, with residual about -7.3e-5 C13 atoms per carbon-stage attempt. A 5% mixture is also saved. This simplification is not selected: the closed ledger retains separate C13 processing and exact stationarity.

`collector_screen.py` computes ballistic arrival and a constant-property surface conduction response. At radius 2e10 m, the sum of separate peak rises for illustrative carbon remains around 525 K or less among current cases. Prompt radiation dominates B1 and is less numerically settled than bulk burn. Photon energy is treated as surface heating for this temperature screen; gamma transport and wall survival are not calculated.

`plant_budget.py` combines the ledger with explicit component allowances, active support losses, pulse storage and gross-rated conversion/cooling. Its base allocation is about 6.78e25 kg (0.036 Jupiter masses) and 6.78e22 W net, exactly 1 mW/kg by the chosen operating rate. It assumes 30% driver efficiency, 40% thermal conversion, 99% energy capture, 99.99% catalyst collection per stage, the smaller helium factory and a 10% auxiliary reserve. Its 1e25 kg driver and area-based processing allowances are not engineered designs. Some worsened efficiency/collection cases fail. See `PLANT-BUDGET.md`.


## 29 October: nonneutral carbon-13 branch

Archived the current PADF-2 C13 evaluation (EVAL-MAR26, distributed April 23) and parsed all 34 MF10 production channels. The data distinguish useful N13 production from inclusive neutron yield. Independently evaluated D and C12 production differ even where they should match; the discrepancy is recorded, not hidden in a conserving network. Older ACSELAM estimates substantially understate competing deuteron production relative to the newer evaluations.

`nonneutral_screen.py` closes the reduced neutron-catalysis/deuterium cycle, including p,d restoration by neutron capture. It charges consumed proton kinetic energy, Coulomb recoil and minimum inelastic excitation, and separates charged energy from neutron/photon heat. At 4.5 MeV with 0.1 MeV RMS width and ln Lambda=25, recoil costs 22.64 MeV per useful event; 90% acceleration / 90% direct conversion / 40% heat conversion loses 0.94 MeV, whereas 95% / 95% / 40% leaves 1.85 MeV before other plant loads.

The same Coulomb collisions produce about 65 rad² of accumulated small-angle variance per useful event. An explicit angular-diffusion escape/reinjection model charges 8.72 MeV even for a 90-degree acceptance cone at 95% / 95% recovery, eliminating that surplus. This does not apply to all-angle confinement or prove an absolute recovery limit. It identifies the missing phase-space treatment.

A necessary cold Brillouin/field-energy screen can exceed 1 mW/kg of allocated magnetic structure at 50–100 T under the preliminary positive energy balance. It omits angular recycling and the rest of the plant, so it is not a successful complete-system example. Prompt N13 extraction avoids storing inactive catalyst in the expensive field volume. A conducting-pipe cold-fluid calculation shows why an infinite homogeneous two-stream argument cannot by itself reject every narrow beam; the stated narrow pipe also fails unrestricted pitch-angle acceptance.

`validate_nonneutral.py` checks isotope closure for arbitrary p,d substitution, exact energy accounting (maximum residual about 1.8e-12 MeV), impact-parameter integration, the small-width limit, and the ideal cold-pipe dispersion threshold. See `NONNEUTRAL-BEAMS.md` for equations, numerical results and unresolved physical assumptions.


## 2026-10-29 — Thermal ions and a finite radial equilibrium

Added the thermal C13 neutron-catalysis screen and two-species rotating Poisson–Boltzmann calculation, with independent mathematical checks and a four-panel sensitivity figure. Ion–ion thermalization carries no false beam-stopping debit. Product extraction instead removes reaction-selected kinetic energy, and the trap adds rotation, potential-work circulation and a specified full-tail wall loss. Vacuum ion radiation is a separate sensitivity because conducting boundaries alter emission.

The 100 T structure-only result is poor under 40% thermal recovery and 50 MeV side-product processing, but improves substantially when those assumptions change. The documented range reaches 1.2 mW/kg only in an extreme 80% heat-conversion/zero-regeneration-cost case; no full plant is established. Actual catalyst recycling, stronger structural materials, magnet current/stress limits, end plugs, neutron capture, deuterium-plant loads and heat-engine mass are the next discriminators. A fixed 1 MJ/kg strength or 40% heat efficiency must not be turned into a physical impossibility claim.

Validation: incomplete-gamma vs Gaussian nuclear integrals 1.1e-9 relative; independent Gauss-law check 2.4e-10; tighter-equilibrium overlap 2.4e-10; Monte Carlo reaction-energy bias, dilute Gaussian and field-scaling checks passed. Reduced-cycle nuclear energy closure <2e-12 MeV. Physics assumptions and inclusive-channel inconsistencies remain visible.


## 2026-10-29 — Broad catalyst survey and heavy-ion confinement

Parsed 61 PADF-2 proton-target files and used NUBASE2020 to identify 44 initial stable-target/ground-state weak-return candidates. This replaces the earlier implicit focus on carbon with an explicit comparison. All raw evaluations are archived in the original ZIP; selected parsed curves are compressed reproducibly. The only existing-parser change is Latin-1 comment decoding, needed for non-ASCII author names.

Ti49, V51, Cr53 and Mn55 have favorable lower thresholds and smaller destructive-channel burdens than carbon in the tested range. A generic main-cycle thermal ledger and charge/mass-dependent Poisson–Boltzmann solver now compare them. The generic side-cost treatment must not overwrite the dedicated carbon p,d regeneration model. Heavy residual production is inclusive; the side-cost sum remains a specified diagnostic, not a closed factory.

At 100 T, 1 MJ/kg strength, factor 5 field mass, 95% heating/direct recovery, 99.9% potential recovery each direction and 50 MeV side cost, Ti49 reaches0.163mW/kg field-only at1.4MeV/40% heat conversion and 0.375 mW/kg at 1.7 MeV / 60%. Other tested catalysts perform worse. At 1 MeV the Ti main-cycle surplus 5.651 MeV is reduced to 4.222 MeV by the modeled trap costs; potential work circulated is 81.2 MeV per reaction. Complete magnet, end/extraction, cooling, neutron capture, D-reactor and isotope-factory mass/load are still absent.

Cold decay inventory is small for V49 despite its 330-day half-life. Mn53's 3.7-million-year half-life makes startup a distinct problem even though its steady inventory can fit within the reference field mass. Electron capture occurs after extraction and electron restoration; neutral decay rates are never applied to fully stripped trapped ions.

Validation checks 352 quadrature moments against incomplete-gamma integrals to 8.22e-9 relative, reproduces the existing carbon equilibrium to 6.22e-15 absolute, and passes independent Gauss-law, solver convergence, charge/mass separation, field/strength scaling and main-cycle energy closure. The four-panel figure was visually inspected and its footer spacing repaired. Mathematical agreement does not validate the unresolved side network or actual reactor hardware.

Read the TUNL Li7 XLS records; 1976SE05 is a pn dataset with a 0.578 barn peak near 2.259 MeV. This corrects an earlier mistaken interpretation of the link context. The separate Taskaev alpha data are not yet combined into a lithium model. Primary Ti49 measurement metadata are identified, but publisher/EXFOR access failures prevent a fresh numerical-table check through those routes. The investigation continues with titanium restoration and a complete magnet/cooling budget, alongside the opaque-pulse driver gap. No future date or final article has been selected.


## 2026-10-29 — Shared solenoid, whole-mass allocation and economic sensitivity

Refined 504 Ti49 equilibria at 1.7 MeV and 60% heat recovery. The strongest retained profile improves the earlier field-only result from 0.375 to 0.569 mW/kg at 100 T. Tight numerical convergence confirms the overlap and net-energy values. A distributed conductor/hoop winding replaces the constant field-mass multiplier in a separate complete allocation. It counts field energy outside the bore, axial compression and buckling screens, finite on-axis field reduction, warm water, walls under electric and water pressure, cold decay inventory, conversion of circulated energy, the D plant, end allowances, refrigeration, radiators, replacement factories and spare inventory.

The base 120 T, 1.6 MeV, 16 GPa hoop / 8 GPa compression case has 1.410 mW/kg bus output and 1.080 mW/kg after replacement manufacture at 100 MJ/kg and 0.01% per year. Total allocated mass is 7.049e14 kg. It assigns 50 MeV per inclusive side-restoration event and 0.2 MeV per useful pn event to extraction. Lower recovery or higher process costs can erase its target margin.

The finite-end check exposed an additional requirement: a 40 MV proton barrier presents 920 MeV to V23+. An explicit branch now charges that coherent circulation, its recovery loss and converter mass. At 120 T and 60% heat recovery, 99.9% transfer each way reduces output to 0.590 mW/kg; 99.99% gives 0.987 mW/kg. At 140 T, 70% heat recovery and 99.9% transfer, the allocation reaches 1.367 mW/kg at 1.360e15 kg. A selective path with lower circulation needs its own physical design before the base extraction-only allocation can be preferred.

The analytic infinite-current-density winding has energy per length proportional to exp[B²/(mu0*sigma)]-1. Together with B^4 reaction power, it gives a structural optimum near 179 T for 16 GPa, matching the broad numerical turnover. Independent analytic/quadrature, virial, finite-sheet field, source-convergence and ledger checks pass. The final validation has 72 checks. Ion–ion dipole bremsstrahlung is about 0.000166 MeV/pn at 1.6 MeV in the classical screen; electron production and clearing remain open.

Primary reading includes the July 2026 CNT paper reporting 12.5 GPa breaking strength in thin ribbons. This updates the earlier 8 GPa specimen without turning either into a reactor working allowance. The thin-ribbon and rolled-fibre moduli are different. REBCO data inspected through 45 T do not establish 100–180 T conductor operation. Magnetoelastic coupling and neutron damage remain distinct constraints. The attempted local Aachen PDF download returned HTML; web full text was used.

The base reference embodied-energy payback is about 2247 years at bus output; the explicit 140 T example gives about 1865 years. These are model results under the specified manufacturing energy. The next decisive work is a conserving Ti side network and continuous extraction/electron management, with transport and lifetime checks. The opaque-pulse driver gap remains active. No final article or future date has been selected.


## 2026-10-29 — Titanium recycling and independent reaction data

Retrieved standard TENDL-2025 proton evaluations for sixteen Ti/V/Cr/Mn/Fe/Co/Ni processors and archived the original files with hashes. A small MF3 reader agrees with the IAEA parser. Alternate2023 MT-file residual yields failed sum-rule checks when combined with their MT5 cross section; those trial results were rejected. The 2025 standard files expose explicit channels below 30 MeV. The high-energy MT5 tail is conservatively removed from useful production and assigned material makeup plus a 1000 MeV/event electrical reserve.

The new stationary network accounts for p,d; p,np; p,2n; p,alpha; p,gamma; additional exclusive channels; cold decays; and isotope-selected neutron captures. Mn53 and Ni59 are proton-processed, avoiding their long electron-capture storage paths. Mn54's small beta-minus branch is included. At 1.6 MeV the auxiliary pn total is 0.27546 per primary Ti pn, catalyst capture demand 0.025614, useful source neutrons 1.269576, and fresh Ti49 feed 1.60773e-6 atoms. Makeup is 4.12 kg/year per 1e24 primary events/s. Main nuclear energy comes from hydrogen; positive trace makeup-transmutation heat is excluded. The known weak-isomer margin is 0.005892 MeV/primary, including Mn52m.

The ground-state/total-state distinction explains an apparent sixfold Cr52 pn discrepancy: PADF's Mn52m yield dominates its ground state. Ti49 pn rates agree within about 2.3% between PADF 2026 and TENDL-2025; titanium capture and vanadium alpha side ratios differ more across TENDL versions. These differences require sensitivity treatment, not averaged certainty.

The cold-bank screen uses 2200 m/s equivalent flux and NIST absorption values as capture upper bounds, with the measured V50 radiative-capture value 21(+4,-2) barn. It gives 1.905e11 kg of listed stock at the earlier 140 T rate/water scale. This is an inventory lower bound, not neutron transport or a revised complete plant. Stable precursor pools can be preloaded. Likewise, a long-lived radioactive decay seed would produce its steady activity immediately if manufactured in advance; its production energy and stock must be counted rather than assuming a mandatory wait of several half-lives.

Validation: 172 independent-reader, incomplete-gamma integral, material-balance, makeup/waste and atomic-mass energy checks pass. The five temperature cases close baryon accounting below 1e-15 per primary event and energy below 1.1e-11 MeV. The shared-solenoid model still uses its original 50 MeV side parameter; next allocate the auxiliary ion channels and replace that parameter with this network's complete costs. No final article/date or new 1 mW/kg claim was selected.


## Explicit isotope processors in the complete plant

`NETWORK-PLANT.md` replaces the generic restoration cost with the conserving TENDL-2025 network. MF6 spectra determine neutron/photon energy; all fifteen nonzero-flow ion processors receive their own equilibrium and field volume. The finite-family allocation uses 2161 central configurations with eight wall choices. Auxiliary channels take 41.0% of active volume. Radial/potential/end energy circulation, converter capacity, stock processing, wall displacement, cooling and replacement are counted.

The 140 T, 1.6 MeV, 70%-heat, 99.9%-end-transfer reference gives 2.160 TW at 1.466e15 kg, or 1.474 mW/kg after replacement. The extra electrical margin to 1 mW/kg is 1.11714 MeV per primary Ti event. Lowering coherent end transfer to 99.8% gives 0.942 mW/kg. Ion heating at 90% gives 0.680, and 8/4 GPa hoop/compression gives 0.369. Finite extraction, electron clearing, transport and materials remain unvalidated. This is no completion date or demonstrated reactor.

At the stated 100 MJ/kg manufacturing energy, construction equals roughly 1770 years of bus output or 2150 years of net output. Hydrogen consumption is 563 tonnes/year and delivered efficiency is 0.135% of fuel rest mass. These are this conditional reference's economics, separate from the mature-fusion canon benchmark.


## Electron access and static clearing-field response

Saved a relativistic invariant-orbit model and an independent Lorentz integrator for the selected cylindrical source. A wall electron requires at least287.3MeV to reach the titanium99%-reaction radius. A30MeV electron is limited to the outer1.499mm by the necessary all-direction access bound. All2560 sampled conservative orbits return to the wall. The independently integrated trajectories agree in penetration below2.3e-9 relative. Perturbative radiation is small per tested pass but does not establish recapture of arbitrarily grazing particles.

Linearized Poisson–Boltzmann response gives median-radius axial-wall-bias transmission0.001269 and reaction-weighted0.007273. A simple1kV wall gradient cannot be treated as1kV across the core. Constant-density modified-Bessel solutions and finite nonlinear equilibrium perturbations independently verify the calculation. The72wall and35bias checks pass; results, input hashes and reproduction times are saved. The figure was visually inspected.

This checkpoint leaves the conditional1.474mW/kg plant unchanged. Bulk births, strong-magnetization collisions, cyclotron radiation, positron extraction, atomic preparation and opposite-charge instabilities remain to be modeled. NRL and primary strongly magnetized plasma work identify where ordinary unmagnetized collision coefficients need qualification. No final reactor, article or breakthrough date has been selected.


## Trace-lepton collision/radiation benchmark

Built a local unmagnetized Landau test-particle energy closure for the main titanium source, retaining relative ion/E-cross-B flow and free-space cyclotron radiation. Five reaction radii, prescribed Coulomb logarithms 5, 15 and 30 and flow on/off give 30 cases. At the median radius and logarithm 15, mean energy is 3.468 keV and radiation 16.86 MeV/s per lepton. The local tangential drag implies an electron outward force drift 12.14 mm/s under this approximation; positrons have the opposite force-drift sign.

A self-adjoint velocity-correlation solve retains energy changes during direction decorrelation. It gives median homogeneous axial diffusivity 7.343e9 m2/s, versus 8.545e9 from the fixed-speed average. A homogeneous absorbing diffusion interval would imply 4.65 ms mean residence, a benchmark rather than the actual channel's escape prediction. An independent Ornstein–Uhlenbeck case verifies the diffusivity solver. All86 mathematical checks pass; figure visually inspected; source/output hashes and run times saved.

The source has rms electron gyroradius about 0.001–0.002 of its local ion static screening length. Thus the unmagnetized coefficients require replacement or bounds. The next check is a gyro-averaged velocity simulation with finite electron escape, followed by applicable strong-field encounter/response physics. Bulk birth spectra, positron removal, atomic processes, stability and global energy/angular-momentum accounting remain open. The 1.474 mW/kg plant allocation, article and future date are unchanged.


## Gyro-averaged electron escape and source requirements

Replaced the isotropic closure with a stochastic velocity/pitch/axial-position benchmark, retaining unmagnetized Landau encounter coefficients and fixed median-radius source parameters. An Euler pilot caused angular numerical heating and a Milstein pilot produced negative-energy proposals; both were rejected. The retained positive radial-energy/energy-conserving angular scheme passes independent Lorentz, Maxwellian and radiation limits.

At20.235km,4096cold electrons have mean residence4.719±.081ms and radiation77.806±1.354keV; halving the step gives4.665±.082ms and76.583±1.344keV. All sampled electrons escape, with no negative-energy proposals. Mean occupation energy3.469keV agrees with the old energy closure. Free-space diffusivity is higher: two controls give10.42±.35 and9.27±.32 billionm²/s. Half/reference/double-length escape gives a quadratic diffusion term10.34±.99billionm²/s and a linear boundary correction about31% of the reference mean time. Thus the earlier4.65ms diffusion agreement involved partial cancellation, rather than validating the isotropic transport closure.

Mean diagnostic outward displacement60.79µm stores51.63keV of electrostatic potential energy. Laboratory collision work134.86keV closes against83.23keV in the drift frame plus that potential transfer. Actual charge return and tens-of-MeV end circulation remain to be designed before assigning net electrical cost. A separate parameterized table translates an assumed complete cost into allowable fresh hydrogen release; it does not assign actual pressure, desorption or birth rates.

All72checks pass. Seven production/control jobs completed; per-trajectory data, source/output hashes and timings are archived. The figure was visually inspected. Ady2015 provides primary wall-coating evidence under4keV synchrotron irradiation, not reactor nuclear-radiation yields. The titanium neutral-collision energy must include laboratory rotation, about.78MeV/u at the median radius; low-energy charge-exchange tables cannot be transferred unchanged. The next discriminator is applicable strong-field encounter physics and then physical birth/collection, with positron removal separate. No reactor allocation, article or breakthrough date has been changed.


## Magnetic collision scales and an independent cold-response limit

Mapped source and encounter scales at five titanium reaction radii. At the median radius, the electron rms gyroradius is1.182µm and local ion static screening is565.47µm; proton/titanium thermal gyroradii are1.333/.422mm. The titanium self-coupling parameter is5.75e-7. Close encounters occupy a large geometric logarithmic interval, with a separate distant interval that spans many electron turns. Spatial gradients and ion turning limit a universal homogeneous approximation.

Derived the linear velocity response of a cold electron to a prescribed straight-ion Coulomb pulse using its zero-frequency parallel force and cyclotron-frequency transverse Fourier components. Modified-Bessel functions produce exponential suppression of distant transverse heating while retaining an axial impulse. Direct Lorentz ODEs, Fourier quadrature and the zero-field NRL limit agree. A drifting-Maxwellian cutoff screen gives median Ti initial heating1558.88MeV/s versus2114.28withoutB for the selected weak-inner/short-outer range; this is an initial cold-electron response, not the warm-population cost.

All67checks pass and the figure was visually inspected. Scripts, outputs and input hashes are archived. Read the2023Jose–Baalrud accepted primary paper's method and parameter ranges; its attractive examples atΓ=1and heavy-ion/electron-bath ordering cannot supply this source's tensor directly. The next calculation must retain finite electron gyroradius, consistent drift/diffusion and matching across encounter regimes. The plant allocation and article chronology remain unchanged.


## Finite-energy distant collisions and continuous helical response

Derived the weak distant guiding-center diffusion tensor and its consistent Ito drift, including first-order ion recoil/polarization. Independent second-order encounters with a smooth screened potential verify covariance and mean force below 3.2e-8 relative and the per-encounter energy identity below 2.0e-10 normalized. Direct Maxwellian quadrature, NRL recovery, cold-force-pulse normalization, positive diffusion and local balance bring the suite to 71 passing checks. A low-speed error-function cancellation discovered by derivative testing was fixed with a nine-term series.

At the median Ti radius with prescribed log interval 4, axial heating changes from +375.03 MeV/s for a stationary electron to -22.63 MeV/s at 1 keV. Diffusion alone incorrectly gives +15.36 MeV/s at that energy. A prescribed Gaussian parallel temperature 2.312 keV averages +11.34 MeV/s and 10.50 mm/s outward drift; these are probes of the operator and do not replace the earlier stationary/escape result.

Built a finite-gyroradius spectral prototype with straight Maxwellian ions, static Yukawa screening, a smooth wave-number filter and explicit cyclotron harmonics. It retains derivatives of the resonance and Bessel amplitude plus bath polarization. All 28 initial checks pass, including strong/weak-field limits, independent state finite differences, cold-perpendicular boundary and quadrature refinement. An independent full-orbit encounter comparison and source-specific short/long matching remain the next discriminator. No new escape job, birth rate, plant power, article or breakthrough date was selected. Figure visually inspected; code, results and hashes retained.


## Independent full-Lorentz finite-gyroradius verification and source pilot

Verified the Yukawa/Gaussian real-space potential, force and gradient by inverse Fourier integration. First/second-order full Lorentz encounters reproduce the helical spectral covariance and mean force below 5.68e-6 relative in two refined ensembles, including ion recoil and both field signs. Work identities close below 4e-8 normalized. Refining impact/gyroangle quadrature gives monotonic convergence. All 39 checks pass.

Applied continuous short/long matching to a cold state and a representative warm state at the median Ti radius, with the same prescribed total log 15 as the old benchmark. Warm perpendicular/parallel/total energy rates at filter factor1 change from +3.474/-29.451/-25.977 MeV/s to +3.568/-29.328/-25.760 on quadrature refinement. Perpendicular heating therefore remains about 2.7% sensitive to resolution. Initial-grid filter factors2/4 give totals -25.583/-25.301. No source convergence, stationary distribution or reactor load is inferred. The source data and cold parallel-limit suite passes 24 checks.

All seven source pilots and the independent orbit jobs are terminal. Their code/data, run times and hashes are retained; figure visually inspected. Next resolve the narrow angular resonance using coordinates aligned with the relative guiding-center velocity, refine Bessel oscillations, and then establish filter convergence before updating transport. The plant allocation, article and future chronology remain unchanged.


## Resolved source quadrature and streamed helical evaluation

Separated radial and angular numerical errors using a resonance-aligned integral and a magnetic-axis moment reduction. The two methods agree below1e-7relative at filter factor1. Polar refinement64to192already converges, while radial refinement corrects the previous warm pilots. At factor2the old perpendicular rate was about52%high; its much smaller total error came from cancellation. Retain the old pilots as historical numerical diagnostics.

Resolved filter factors1/2/4/8give perpendicular heating3.512814/3.591389/3.618456/3.627076MeV/s and total-25.777738/-25.776864/-25.776068/-25.775725MeV/s. Further radial refinement at factors2and4and all-three-grid refinement at8agree to numerical precision. The4to8matching change is.238%perpendicular and.00133%total. These are single-state values, not a stationary distribution or physical screening qualification.

A stable Bessel recurrence uses direct special functions in the evanescent region, a contour bound for negligible tails and two-order streaming. It reduces memory and run time while matching the direct table. The28prior mathematical checks and128new independent-coordinate, recurrence, source-resolution and positive-covariance checks all pass. The source/validation/plot jobs are terminal; the figure was visually inspected; code/data/hashes and timings are saved.

The next physical ordering check must account for the saved xi=.5rigid rotation: magnetic and Coriolis terms cancel in Ti's co-rotating thermal equation. Electric/centrifugal curvature, gradients, collective screening and electron-frame corrections remain. Lepton births and finite removal determine whether transport refinement matters economically. No new reactor power, abundance, residence time or date was selected.


## Actual rotating-frame trajectories and distant-range sensitivity

Derived the thermal Lorentz equation in the saved rigidly rotating source. Atxi=.5the Ti magnetic and Coriolis terms cancel exactly; protons retain7.188e9rad/s. The nonlinear electric/centrifugal force and its gradient remain. Independent laboratory and rotating orbits, theH-wP_phi invariant and the equilibrium density gradient pass58checks, with maximum normalized residual1.02e-10. At median radius, maximum sampled Ti relative-path curvature is.589%of a screening length for the warm probe and6.26%for cold. This does not bound slow encounters or collision-rate corrections.

Spatial gradients remain a separate limitation: the Ti logarithmic-density change across one static screening length is.101atmedianradius,2.84at99%and10.5at99.9%reactionradius. A consistent long-range-removal test raiseskappaby1/2/10while loweringtotal logbyln(f). Warm axial cooling changes29.37→26.02→18.24MeV/s; perpendicular heating stays3.591–3.593. Thus numerical convergence does not remove screening uncertainty.

An analytic exponential-integral formula verifies the spectral logarithm and the1.16e-5finite-filter correction to the complementary log. All48screening checks pass, including radial refinement and independent distant-operator comparison. All jobs are terminal; scripts/results/hashes are retained. Next prioritize actual lepton supply and finite collection so further transport work is tied to the power margin. Plant output, article and date remain unchanged.


## Prompt internal pairs and finite-residence weak-decay ledger

Folded the evaluated photon spectra through ordinary BrIcc pair coefficients, retaining all six pure multipolarity scenarios and explicit out-of-table intervals. In 1.1–8 MeV, the gamma-only cascade scenarios give 0.000253–0.000906 pairs/primary. Unrecovered 40 MeV positron work would cost 0.0101–0.0362 MeV/primary, 0.91–3.24% of the remaining electrical margin. E0, actual multipolarities and stripped-ion cascade changes remain unresolved. The MT5 photon contribution is separately retained and small.

The exclusive photon energy reproduces 3.29596427 MeV/primary with 2.15229285 photons/primary. Coefficient-grid halving changes integrated pair scenarios by at most 9.53e-6 relative. Direct-residual neutral-rate weak-decay proxies at an assigned 1 ms residence are 6.36e-9 ground-state and 1.49e-8 when every nuclide is placed in its fastest known weak state. These include EC as if it supplied a lepton and leave charge-state changes explicit.

All 111 numerical/data checks pass; all jobs are terminal; figure visually inspected. Code, source hashes, coefficient grids, detailed ledgers and validation are retained in PROMPT-PAIR-SUPPLY.md and its manifest. Next quantify gas/injection/external-photon supply and finite positron collection without excessive ion loss or RF heating. Plant allocation, article and future date remain open.


## Transverse clearing, collective crossings and heat/torque accounting

Computed ion transverse spectra in the saved nonlinear rotating equilibrium. Exact magnetic splitting, harmonic circular modes, an independently forced oscillator, adaptive nonlinear ODEs, Parseval and the velocity spectral moment verify the integrator and circular-response normalization. Production ensembles contain 16,384 particles per species; a matched 4096-particle half-step control is second order. The initial sparse proton tail failed the 12% seed-comparison criterion at 21.1%. Refinement reduces that difference to 4.54%; the largest refined broad-band difference is 5.82%. Narrow-bin precision is not asserted.

Derived the rigid two-species centroid force from density overlap and the grounded cylinder's image dipole. Independent displaced-density integration verifies the mutual-force coefficient below 7.3e-7 relative, improving fourfold on displacement refinement. Equal-q/m controls recover common translation and internal-force cancellation. Four real closure poles include two crossings of the lepton drift sweep at 5.858 and 20.548 mm. These identify coupled-response requirements, not kinetic instability or a reactor loss.

The clearing field can extract ordered laboratory energy while heating the rotating distribution. Maintaining rotation and crediting useful replacement heat must share one energy ledger. At the representative drive frequency, ideal restoration and recovery are 2.031 and 1.031 times rotating heat. Gross converter ratings and RF losses remain separate. All 100 retained mathematical/broad-band checks pass; figure visually inspected; all jobs terminal; code/data/hashes saved. No revised plant power, selected collector, article or date is claimed.


## Construction energy, replacement and alternative supplies

Added a quantitative energy-capital comparison. The nonneutral reference's 100 MJ/kg construction repays in 2150 years; charging its stored field at an assumed 95% adds 19.47 years. Replacement already consumes 17.70% of bus output. Ideal continuously reinvested construction grows at only 0.0465%/year before field charging. Discounted energy payback and transportation-energy scales expose the importance of long service life and local manufacture.

Found an explicit-accounting gap in the opaque allocation: its auxiliary reserve did not supply a separate replacement ledger. Reclosed mass, support, manufacture, spares, conversion and radiators. Replacement of 0.01%/year at 100 MJ/kg requires 31.75% more pulse throughput to preserve 1 mW/kg. The revised conditional mass is 0.035765 Jupiter and aggregate pulse interval5.460days; radiators occupy29.64%of collector area. The0.1%-replacement scenario exceeds the available spherical radiator surface and is flagged. Its numerical mass closure is not geometric qualification.

Compared both DD fuel cycles with neutron recapture distinguished, ideal isotope-separation entropy, local versus bulk fuel lifting, solar areal-mass scenarios and published thin-disk accretion energy yields. At historical Jovian abundance, the reference power needs401kg/s of hydrogen feed for the assigned30%-efficient DD comparison; the protium reference accesses22500times more energy per initial hydrogen reservoir. Solar and easier fuel access can dominate the economics long before hydrogen fuel exhaustion. No arbitrary helium revenue or market price was assigned.

All46checks pass, including independent mass closure, conservation and chemical-potential integration; figure visually inspected; all jobs terminal; code/data/input hashes retained. Primary Galileo, IAU and Thorne texts inspected. Article architecture, future date and final feasibility judgment remain to be synthesized.
