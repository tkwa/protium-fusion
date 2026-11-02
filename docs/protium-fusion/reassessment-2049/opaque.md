# Opaque protium route: reassessment for a superintelligent civilization ✨

29 October 2026. This memo supersedes the earlier selection of a 2049 gain of 2.6 and an unsuccessful opaque campaign. That selection was invented history. The incomplete cold-fluid calculations do not establish that future engineers would suffer that loss.

**Recommendation:** write the opaque route as a successful late-2040s development of controlled catalytic pulses, followed by an adaptive, much smaller collection architecture. Repeated high-gain shots and completion of the catalytic sequence by 2049 are plausible extrapolations in this canon. A first experimental power station around 2049–2052 is also physically credible if the article explicitly supplies efficient launchers, smaller recurring pulses, regenerative capture, and a shared catalyst-manufacturing industry. The old eleven-Earth-mass thermal allocation is neither a minimum reactor mass nor the design that competent engineers must build. The strongest reasons for limited deployment are capital allocation, recovery/maintenance engineering and competition with solar and deuterium, rather than an established impossibility of creating the compressed state.

The dates, successful integrated hardware and any eventual measured output remain selected future history. The calculations below bound what that history must accomplish; they are not predictions of a future experiment.

## 1. Controlled inflow is a design problem with an identifiable route

Use composition-controlled cold parcels, launched on a radial velocity/arrival schedule, discharged after acceleration, and assembled through a graded density profile with a designed weak-shock or ramp sequence. Neighboring parcels should be nearly co-moving when their voids close. The design should be optimized forward with a nonideal EOS; it need not reproduce the exact time reversal of the old numerical expansion.

The earlier symmetric-jet merger screen illustrated the cost of intersecting streams at substantial angles. It did not establish an inevitable entropy floor. In a smooth radial flow, the difference between neighboring velocities scales approximately as `delta_v ~ v ell/r`. At `v = 2.4e6 m/s`, `r = 20 km` and parcel spacing `ell = 0.1 m`, that is about 12 m/s. Its kinetic energy scale is tens of joules per kilogram, compared with gigajoules per kilogram for ionization and terajoules per kilogram for the whole inflow. Errors, porosity, nonuniform packing and coherent distortions can be much larger; the important point is that full-speed head-on collision is avoidable.

A prospective parcel design can use frozen nitrogen/hydrogen mixtures for nitrogen stages and hydrocarbon compositions for carbon stages, with thin conductive, load-bearing surfaces whose mass enters the reaction network. A 1 GV/m surface field gives 4.43 MPa of Maxwell pressure. For a separate spherical shell carrying that pressure, a 1 GPa allowable stress and shell/core densities of 2000/1000 kg/m³ give 1.33% shell/core mass; at 10 GPa the ratio is 0.133%. These are material requirements, not demonstrated macroscopic shell properties. Atomically precise manufacture makes careful composite optimization more credible than merely assuming an uncoated cryogenic solid can carry the stress.

A strong first shock is not the only way to obtain a useful adiabat. Shaped shock sequences already change the compression/stability trade in ordinary ICF. MacPhee et al. experimentally reduced perturbation growth with an adiabat-shaped drive; their particular DT capsule performance is not transferred to this system. It establishes that entropy generation and instability growth can be design variables. [MacPhee et al., 2015](https://doi.org/10.1063/1.4928909).

The cold EOS is also tractable scientific work, rather than an unbounded unknown. First-principles nitrogen EOS calculations already cover molecular, warm-dense and plasma regimes, including ionization-related compression structure; experiments constrain shock-compressed nitrogen to 800 GPa. Mixtures, surface coatings, loading history and extreme compression still need their own models and validation. [Driver and Militzer, 2016](https://arxiv.org/abs/1601.05787); [Nitrogen shock measurements, 2022](https://doi.org/10.1103/PhysRevLett.129.015701). The parallel transport reassessment finds that replacing the assumed prepared internal energy by additional cold-flow kinetic energy requires roughly 7% more kinetic energy, or 3.7% more speed. That is an energetic scale, not proof that the optimum entropy history has been found.

**Physical risk that remains:** low-order asymmetry, premature shocks and reaction-driven expansion could prevent enough mass from reaching useful density, even after optimization. Opacity alone cannot resolve this. However, the current guard stops and viscosity dependence are evidence that our surrogate has not decided that risk; they are not evidence that a superintelligent research programme fails it.

## 2. Discharge and high driver efficiency are credible engineering choices

For the illustrative small nitrogen pulse, retaining all accelerating charge gives a collective Coulomb barrier near 175 km. A uniform-sphere example requires more than 99.87% compensation to limit the residual charge-energy penalty at 30 m to 1% of incoming kinetic energy. The proper architecture therefore discharges or compensates parcels after acceleration, while their trajectories are still correctable. Neutralized drift compression has laboratory precedents, although those ion beams are very different from a complete macroparticle pulse. [Seidl et al., 2007](https://nonneutral.pppl.gov/pdfpapers2007/SeidlHIF06.pdf).

The compensating electrons for that example have only about 0.185 kg of rest mass. The isolated grains' stored electric-field energy is approximately 13.3 kJ/kg, about four billionths of their several-TJ/kg preparation energy. Neither the charge supply nor its unavoidable discharge energy forces poor efficiency. Macroparticle Larmor radiation is likewise negligible at the specified acceleration. Current transport, field emission, leakage, switched-field energy recovery, wall interception, metrology and failed-pellet rejection determine the real efficiency.

A mature 50–70% electrical-to-prepared-pulse efficiency is a reasonable **future engineering assumption**, and is not ruled out by these energy scales. An illustrative accounting makes the assumption concrete. If 93% of prepared energy is bulk motion and 7% is heat, the combined efficiency is

`eta_prepared = 1 / (0.93/eta_motion + 0.07/eta_heat)`.

Motion-channel efficiency 0.6 with heat-channel efficiency 0.2 gives 52.6%; 0.7/0.3 gives 64.0%; 0.8/0.4 gives 74.8%, before any source loads not included in those channels. The channels must include power conditioning, applied-field coupling, interception and deliberate heating; plant refrigeration and processing still require their own accounting. This is a component requirement, not a measured projection. Reversible electrical acceleration of bound parcels has no intrinsic laser-ablation-like efficiency penalty.

Use segmented induction/electrostatic acceleration, local pulse storage, charge-controlled steering, precise optical timing and recovery of unused field energy as the named mechanisms. Heavy-ion accelerator studies establish related engineering concepts; they do not establish these macroparticle efficiencies. The old 30% figure remains a conservative scenario, not a physical bound or necessarily the median of an optimized future design.

## 3. The complete cycle has considerably more margin than the failed-shot story implied

I recomputed the unchanged ledger with `collection_loss=1e-4` and the existing 30 m helium factory. Before collector losses, per net helium atom it gives:

| Quantity | MeV per net helium |
| --- | ---: |
| Available nuclear heat H | 24.71663 |
| Prepared input U | 2.11879 |
| Final bulk kinetic energy, including recovered input | 19.63614 |
| Final hot internal energy | 1.44016 |
| Delayed decay heat | 5.75912 |
| Helium-factory preparation included in U | 0.36847 |

The saved plant's `cycle_heat` field of 24.448 MeV already incorporates its capture-loss convention. Applying another capture fraction to that field would double-count some losses. Use the pre-collection values above for new comparisons.

For complete-cycle collection fraction `f`, ordinary thermal conversion `eta_h`, and electrical-to-prepared efficiency `eta_d`,

`W = f eta_h (H + U) - U/eta_d`.

Input energy is recovered once. If a sensitivity changes preparation to `lambda U` while preserving reaction yields, its break-even multiplier is `f eta_h H / [U(1/eta_d - f eta_h)]`. That sensitivity is an accounting envelope, not a prediction of how a changed implosion burns.

| Driver | Thermal | Collection | Net MeV/He at reference yields | Break-even preparation multiplier |
| ---: | ---: | ---: | ---: | ---: |
| 30% | 40% | 99% | 3.564 | 1.573 |
| 40% | 45% | 99% | 6.658 | 2.530 |
| 50% | 45% | 99% | 7.718 | 3.343 |
| 60% | 50% | 99% | 9.752 | 4.928 |

Finite heat-engine temperature approaches and actual collector temperature history must support the selected thermal efficiency. These rows should not substitute a Carnot ceiling for an engineered cycle.

Direct recovery can enlarge the margin further, but apply it only to the bulk kinetic component. At 50% driver efficiency, 70% conversion of the computed kinetic component, 45% conversion of hot internal plus delayed-decay heat, and 99% collection, the net is 12.578 MeV/He. If additional preparation energy becomes heat, the same purely accounting comparison permits a multiplier of 4.819. It does not assume direct conversion of every decay photon.

The four catalytic stages and helium-to-carbon makeup need linked material balances. Weak decays taking minutes do not impose a stellar-size machine: products expand, cool and decay in a separate processing phase. Recovery and isotope sorting can occur before the next pulse. Catalyst loss at the 1e-4 level is a real economic requirement, but not a thermodynamic lower bound. Nested catch surfaces, offset ports, regenerative coatings and recovery from secondary gas are appropriate design responses. Correlated leaks and transmutation of implanted catalyst must be included.

## 4. Smaller pulses permit a much smaller thermal station

A full absorber at radius 1e9 m and areal mass 1e4 kg/m² weighs **1.257e23 kg**, 2.6% of Venus. That is 400 times less absorber than the old 2e10 m sphere. The favorable 30 m nitrogen calculation contains 1.13e13 kg and gives immediate capture gain 8.182; its opacity remains ample. This is enough to justify investigating a smaller recurring-pulse station, rather than treating the old full-scale allocation as mandatory.

At fixed temperature and density, reducing the pulse's lengths by ten reduces mass/energy by 1000 and hydrodynamic time by ten. Moving the collector inward by twenty then gives about 0.4 times the old fluence for an otherwise identical outgoing spectrum. A simple conduction similarity gives temperature-rise factors of about 1.79 for the ballistic-arrival component, 1.26 for a prompt component whose duration follows the hydrodynamic time, and 0.4 for a component tied to unchanged beta-decay time. These are scaling diagnostics, not predictions for new burn trajectories. They show that smaller catcher mass does not automatically imply an impossible thermal transient.

We have not calculated the smaller B1, B2 and C13 pulses. Their conversion cannot be copied from A: the shorter residence time and sequential reactions matter. More recycling and a higher driver efficiency are available responses. This is the nearest inexpensive model discriminator for the small-station branch, followed by an actual full-cycle ledger and catcher temperature history.

**The helium factory is the hidden scale mismatch.** Its existing 30 m radius uses density 1e10 kg/m³, and hence mass 1.13e15 kg, one hundred times the small nitrogen pulse. Its delivery needs a distinct helium parcel or neutralized-ion architecture; a carbon coating would itself alter the factory's reaction network. Ordinary dense nitrogen/hydrocarbon parcels do not automatically solve helium preparation.

Even the existing factory is not an immediate thermal contradiction at the smaller catcher: rescaling its saved thermal screen from 2e10 to 1e9 m gives approximately 567 K of ballistic temperature rise plus 153 K of prompt rise in the constant-property carbon example. Nonlinear response and volumetric photon absorption still matter. Its maximum electrical pulse falls from 1.36e29 J at 30% driver efficiency to 6.79e28 J at 60%; at the old installed storage assumption of 0.1 MJ/kg, that alone requires 6.79e23 kg. Thus the smaller shell does not make the source or storage free.

For a station of order 1e24 kg at 1 mW/kg, an unchanged large factory would be needed only on an approximately decadal cadence, depending on the improved cycle margin. Duplicating a mostly idle factory and its storage at every station is poor architecture. A shared regional helium-to-carbon works, with its energy and apportioned equipment counted, or a smaller factory with lower yield and more frequent batches, is a credible adaptation. The factory currently costs only 0.368 MeV of preparation per net helium atom; severalfold extra factory work need not destroy a cycle with efficient drivers. Transported replacement carbon is tiny compared with the circulating catalyst inventory.

The old 1e25 kg driver allowance also cannot be imported unchanged as a minimum. Peak source power, field hardware, storage discharge rate and parcel manufacturing must be designed for the chosen pulse schedule. Conversely, simply multiplying every piece of equipment by the pulse-mass ratio would be unjustified.

## 5. A magnetic catcher is credible, but a second design branch

A magnetic piston can intercept the expanding conductive plasma, transfer its work into segmented pickup coils, and send slowed products to thermal/chemical processing. Segmentation limits local voltage and allows recovery into distributed storage. A primary 1992 analysis studied this mechanism with segmented coils and obtained high analytical efficiency for a cylindrical expansion; its greater-than-80% figure is not evidence for a spherical protium catcher. Plasma direct converters have also operated experimentally at much smaller power in different geometries. [Mima et al., 1992](https://wx1.ans.org/pubs/journals/fst/a_30054); [Barr and Moir, 1983](https://wx1.ans.org/pubs/journals/fst/a_20820).

The scale estimates are favorable enough to retain the branch. A 0.6 T field occupying a spherical equivalent volume of radius 1e8 m stores about 6e29 J. A structural energy/specific-stress estimate with geometry factor 3 and allowable specific stress 3 MJ/kg is about 6e23 kg. Those are explicit hypotheses, not a coil design or mass lower bound. They leave room below a Venus mass for a catcher, but they do not allocate conductors, thermal protection, storage, source, processing or support automatically.

Four constraints decide whether that redesign works:

1. **Recombination and ionization.** An expanding pulse cannot be presumed fully ionized at the catcher. Single-ionizing nitrogen costs order 1e8 J/kg; complete ionization costs order 1e10 J/kg, much less than its several-TJ/kg input and tens-of-TJ/kg outgoing energy. An ionization layer can therefore be energetically affordable. Its optical thickness, photon delivery, recombination kinetics and neutral escape require a transport design. At 1 T and 1e7 m/s, singly charged N14 has a gyroradius of order 1.5 m, far below a 1e8 m catcher: particle magnetization is easy once ionized, while collective plasma interception is the real issue.
2. **Topology and leakage.** One dipole is not a closed spherical catcher. Use overlapping mirror/cusp sectors, offset secondary capture paths and thermal backstops for escaping material. Finite direct recovery of 60–75% of intercepted kinetic energy is a sensible future design target; exact efficiencies require nonlinear MHD and circuit validation. Plasma interchange, reconnection, field penetration and cusp losses cannot be erased by a field-energy comparison.
3. **Momentum and pulsed loads.** A 5e29 J outflow near 1e7 m/s has a scalar impulse of order 1e23 kg m/s. Spherical cancellation of vector momentum does not remove local loads. A 1e24 kg assembly receives a radial impulse scale of 0.1 m/s; the fields, coils and supports must store and return that motion. Asymmetric escaping material also causes recoil and requires recovery or stationkeeping.
4. **Delayed heat.** About 5.76 MeV/He in the reference ledger is delayed decay heat. Slowing the products does not eliminate it. Shielded distributed decay reservoirs, thermal conversion and radiators are needed. An illustrative 1e29 J heat load requires 1e23 kg at a usable 1 MJ/kg thermal-storage swing; local photon deposition and mixing can be more restrictive than that bulk capacity. This term is one reason not to equate magnetic field structure with total station mass.

A sub-Venus magnetic plant is therefore a credible advanced-engineering branch, rather than an established design. For the first history, a smaller thermal catcher with progressively added inductive recovery is simpler than claiming a perfect magnetic sphere solved every collection problem simultaneously.

## 6. Classification of the major uncertainties

| Item | Primary class | What would settle it, and what the bound actually says |
| --- | --- | --- |
| Nuclear mass defect, neutrino loss, energy conservation | Fundamental bound | Fixed ledger; useful heat is below total mass defect. No design removes neutrino escape. |
| Opaque radiation retention | Fundamental transport condition, already favorable in the modeled regime | Absorption/transport depth and radiation-hydrodynamic evolution. Bremsstrahlung is retained energy when truly trapped, not automatically a lost luminosity. |
| Minimum useful burn residence time | Fundamental reaction kinetics plus design | Full-network evolution along optimized trajectories. It constrains pulse scale; it does not mandate Jupiter mass. |
| Cold EOS, ionization and phase transitions | Computable material physics with experimental constraints | First-principles free-energy EOS, Hugoniot/ramp tests, mixture validation. Current ideal-EOS failure is not a reactor result. |
| Packing, shock sequence, asymmetry and premature burn | Design/optimization task with experimental validation | Multidimensional forward simulation and increasingly integrated shots. Co-moving parcels avoid an assumed head-on-merger loss. |
| Grain strength, conductivity, retention and field emission | Experimentally limited combined material properties | Select an actual composite and qualify it under charging, acceleration and radiation; count its mass and composition. |
| Charge neutralization and trajectory metrology | Design/control task | Complete charge/current circuit and finite transport; very little electron mass or compulsory energy is involved. |
| Electrical-to-prepared efficiency | Design/optimization task | Measured component and integrated energy account. 30% is not a fundamental ceiling; 50–70% is a marked future assumption. |
| Catalyst collection, sorting and replenishment | Process design plus experimental surface yields | Nested capture and full isotope balance; chemistry, implant transmutation and correlated leakage matter. |
| Magnetic-catcher efficiency and topology | Design task with experimentally constrained nonlinear plasma behavior | Coupled MHD/kinetic/circuit model and staged capture experiments; field volume alone is insufficient. |
| Absorber lifetime and thermal cycling | Experimentally limited material/maintenance properties | Resolved transient deposition, erosion, recoverable coatings and replacement inventory. |
| A complete station's specific power | System design and operating schedule | Include all source/catcher/storage/radiator/processing mass and recalculate at chosen cadence. 1 mW/kg was a target, not a ceiling. |
| Completion date and number of stations | Construction schedule and economics | Source material, freight, receiving works, opportunity cost and desired market. Missing calculations do not by themselves set a completion date. |

## 7. Recommended actual history and its performance boundary

The strong canon changes the expected research pace. By the 2030s the civilization has atomically precise manufacture, a mostly complete solar swarm, and industrial magnetic guidance, braking and separation. Venus reclamation supplies nearly 4.9e24 kg in 2039–46. These capabilities support substantial research infrastructure and a sub-Venus first station; they do not automatically provide a validated implosion or catcher. The established absence of practical protium power in the early 2040s can remain intact.

A coherent selected history is:

- **2043–46:** Moravec's programme replaces opposing hot-jet concepts with composition-controlled, co-moving parcels and forward-optimized compression. Material/EOS measurements and small hydrodynamic tests determine the first-shock and ionization sequence. Neutralizers operate after the acceleration/steering stages. Launcher efficiency rises into the 50–60% range through field-energy recovery and reduced interception.
- **2047–48:** The open range demonstrates repeatable opaque nitrogen and carbon burns. The headline is repeatability, measured input energy and recovered isotopes. It is reasonable to choose gains of order several to ten for the appropriately sized capture stages, with stage-specific values set by the small-stage checks; the earlier 2.6 is not retained as a predetermined failure. Their microscopic opacity result survives the real assembly process because the optimized flow reaches the required column density.
- **2048–49:** Linked catalytic campaigns include cold decay, isotope return and carbon makeup. The programme chooses smaller frequent main pulses, builds out the 1e9 m thermal collector in modules, and supplies makeup catalyst through a dedicated regional works rather than duplicating its largest rare pulse everywhere. Inductive recovery is introduced on selected sectors and retained where it reduces the measured plant burden.
- **Around 2049–52:** A first specialist power station commissions, with a whole-system mass of order 1e24 kg rather than the inherited 6.8e25 kg allocation, subject to the actual source/storage/capture design. A defensible chosen operating boundary is complete-cycle nuclear/prepared gain above roughly 4–5, 50–60% driver efficiency, about 40% thermal conversion, and high catalyst recovery. For example, at 60% driver, 40% thermal and 99% energy collection, break-even is Q > 3.209 before auxiliaries; Q = 5 gives 0.709 delivered-energy units of margin, and Q = 7 gives 1.501. A reported net output at or above 1 mW/kg must then come from the complete mass, cadence, maintenance and auxiliary ledger. These are narrative performance requirements, not outputs of a completed new plant calculation.

My preference is a successful first station near the end of that window, with an already successful reactor development by 2049. A 2049 first station is not excluded by energy or material supply if the story explicitly allocates part of Venus and fast receiving/transport works. Extending commissioning a few years should be attributed to that chosen construction programme and integrated lifetime qualification, not to our inability to finish a cold-fluid simulation.

The strongest credible stall is economic: a technically working pulse system needs costly collection, storage and catalyst recovery while solar and deuterium remain attractive. At 120 AU, a 30%-efficient solar collector with complete areal mass 10 kg/m² gives about 2.84 mW/kg; a 1 kg/m² version gives 28.4 mW/kg, before additional storage/transmission. A first opaque station is thus a research and strategic dark-location capability, not automatically the cheapest local electricity. Its later market is long-lived installations far from stars, deuterium-depleted feedstocks, or a preference for access to the full hydrogen inventory.

Nor must every mature plant remain at exactly 1 mW/kg. Increasing pulse cadence can amortize fixed mass, while raising heat-transport, source, recovery and maintenance requirements. The old millennial construction-energy payback at that particular specific power is not a universal economic theorem for optimized opaque reactors.

## Scope and reproduction

This reassessment ran no new hydrodynamics and modified no shared article or prior results. Numerical energy-component values were obtained by importing the existing `cycle_ledger.ledger(collection_loss=1e-4, factory_file="compression-He-64-r30-burning.json")`; all remaining figures are the stated algebra or saved thermal/driver results. The key upstream records are `PLANT-BUDGET.md`, `ECONOMICS.md`, `collector-screen.json`, and `phase-two/driver/`. Further small-stage checks can refine the preferred scale without changing the distinction between an incomplete calculation and an actual future failure.
