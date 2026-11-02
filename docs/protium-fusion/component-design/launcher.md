# Launcher and fuel assembly: developed component checkpoint ✨

The preferred starting architecture is **1,000 organizational source complexes containing spatially distributed accelerator columns**, with condensed fuel grains passing through charge-controlled acceleration, final measurement/correction, neutralization, and a protected exit into the chamber. This provides a concrete alternative to treating each complex as a compact point gun. The number of independently positioned outlets remains a design variable: 1,000, 10,000 and 100,000 are being compared at fixed pulse inventory and prepared energy. A lower-yield simpler assembly may be preferable to the most uniform one.

This checkpoint includes nonideal merger entropy, a finite space-charge envelope, a neutralizer flux/energy design, and completed full-network tests of two selected-mixture source profiles. The ammonia profile retains 99.65% of reference gain after a small controlled-preheat adjustment; the methane profile retains 96.75%. A complete three-dimensional cold-to-burning trajectory, qualified capsules and the accelerator circuit remain to be demonstrated.

## Fixed requirements and chosen architecture

Physical requirements are conservation of energy/momentum and charge, attainable material states, a useful catalytic yield, collection consistent with 0.007 kg fresh C per kg fresh H, and whole-station net output above 1 mW/kg. The current 5.4 mW/kg account has room to trade yield for simpler hardware. Retaining 95% of the reference yield is an initial comparison criterion, not a required physical threshold.

The following are choices to test, not laws: the 1,000 organizational complexes; condensed spherical grains; 1 GV/m positive grain surface field; 1–10 MV/m effective accelerating gradient; roughly 100 m column beam radius; reflective grains and shielded source tubes; and spatially distributed outlets. The earlier 10 m spot/4 microsecond example supplies no mandatory tolerance.

The concrete sequence is:

1. **Make stage-specific condensed fuel.** CH₂ and CH carbon-stage cores allow the required C/H atom counts, with isotope-specific carbon. A conductive carbon skin can be included in the intended carbon inventory by adjusting the core's hydrogen loading. Nitrogen/hydrogen chemistry requires its own recipe, below.
2. **Charge and accelerate in distributed columns.** Positive surface charging avoids treating freely emitted electrons as the acceleration charge. Time-varying staged cells provide repeated work; the integrated accelerating voltage is not placed across one terminal. Alternating-gradient optics confine the charged beam. Injector current, longitudinal bunching, beam loading and lens stability require a circuit/envelope design.
3. **Measure and correct while charge remains.** Position/velocity pickup and optical measurements feed electrostatic timing and transverse correction. The control law uses the station's predictable orbital and support motion. A known common displacement is not random focus error.
4. **Neutralize after the final charged correction.** A co-moving electron supply discharges each grain and compensates collective beam charge before dense convergence. Surface-field energy, image-current energy and electron recoil belong in the thermal/phase-space account. Any later charged steering requires an explicit recharge step; neutral grains otherwise coast.
5. **Coast through a specified radiation/gas environment.** Reflective skins, cold shielded exit paths and packet size provide thermal margin. Slow tails can use separately timed seed delivery or longer protected storage/routing, rather than silently sharing the main flight time.
6. **Merge at controlled density.** Neighbouring channels have specified positions and velocities. The phase transition, grain breakup, shock formation and subsequent compression must produce the desired fluid profile. Merger heating is part of prepared fuel energy; adding it to an already heated initial state would double-count energy.

## Grain size, acceleration and charge

The calculation uses `q/m = 3 epsilon₀ E_surface/(rho a)` and `L = v²/[2(q/m)E_acc]`. For a 1,000 kg/m³ core and the C13 RMS speed of 2.339 million m/s:

| Grain radius | q/m at 1 GV/m | Effective gradient | Acceleration length |
|---|---:|---:|---:|
| 0.1 mm | 0.266 C/kg | 1 MV/m | 10.3 million m |
| 1 mm | 0.0266 C/kg | 1 MV/m | 103 million m |
| 1 mm | 0.0266 C/kg | 10 MV/m | 10.3 million m |
| 1 cm | 0.00266 C/kg | 10 MV/m | 103 million m |

The 1 cm/10 MV/m and 1 mm/1 MV/m alternatives have the same acceleration length and time. Larger grains improve radiative tolerance and reduce charge/current, at the cost of higher local accelerating gradient or a longer column. A centimetre grain's integrated accelerating voltage is about 10¹⁵ V; a repeated-cell machine is essential.

At 1 GV/m the outward electrostatic pressure is **4.43 MPa**, and the initial stored surface-field energy is **13.3 kJ/kg** for this density. A core capable of carrying that tension can support a thin conductor. If the skin carries it by itself, a spherical membrane with allowable tensile stress 1 GPa requires `t/a≈0.00221`, giving roughly **1.33% carbon skin mass/core mass** for density ratio two. At a genuinely qualified 10 GPa allowable, the ratio would fall tenfold. The material's irradiated, defective, temperature-dependent allowable must be supplied; ideal graphene strength is not automatically an assembled-capsule specification.

Laboratory electrostatic dust accelerators demonstrate charged conductive and coated grains, measurement and selection, but at very different throughput and velocity. Colorado's facility reports micron/submicron particles approaching 100 km/s, generally of order one selected particle per second or less. This establishes mechanisms to build on, without certifying the required many-grain source. [IMPACT facility](https://impact.colorado.edu/dustaccelerator), [NASA facility description](https://assets.science.nasa.gov/content/dam/science/cds/researchers/research-facilities-for-planetary-science/2024/Facilities_Quad_Chart_Colorado_Dust_Accelerator.pdf).

Free molecular/cluster beams are alternatives, but internal thermal motion becomes beam emittance. A 100 K nitrogen-mass particle has about 244 m/s one-axis thermal speed and spreads roughly 115 km over the nominal flight if no further focusing/collisions intervene. A rigid grain carries this atomic motion internally until it releases its material. Plasma jets permit easier electromagnetic control, while requiring their thermal expansion and early merger entropy to be paid explicitly. The current preference for condensed grains follows this phase-space distinction.

## Finite distributed outlets

At **100,000 independently positioned outlets**, each organizational complex owns about 100 columns. Their typical spacing on the receiver is **12,300 km**, spread over a complex's approximately 69,500 km-radius angular territory. This is an actual large footprint, not extra packets emitted from nearly the same direction.

For 1 cm grains, a 700 s large-shot supply window and 100,000 outlets, the mean charged current is about **430 A per outlet**. A uniform 100 m-radius exiting beam has a self-field of about **33 kV/m**. The corresponding field at a 100 km/s injector is about 0.77 MV/m at the same current and radius. Distributed lenses and injector shaping are required, but these fields provide a quantitative envelope target rather than an unspecified collective-charge cure.

A 100 m-radius, 103 million m-long column has about 6.5×10¹⁰ m² of cylindrical side area. At 10 kg/m², 100,000 such skins total approximately **6.5×10¹⁶ kg**, below the existing 4.22×10²¹ kg source/conditioning/storage allocation. This is a geometric skin floor. Conductors, insulators, power feeds, lenses, thermal shields, mounts, controls, repair and radiator duty must still be assigned. Peak power remains of order 10¹⁸ W per outlet; dividing the source does not remove that demand.

A 200 m-diameter exit can fit between the support team's kilometre-pitch rails. Wide injectors can sit behind the receiver plane. Any source count comparison must charge repeated end hardware and the actual distributed column supports, in addition to the common electrical capacity.

## Nonideal merger calculation

Primary FPEOS tables are archived under `launcher-data/`, with their original headers and hashes. CH and CH₂ include molecular/atomic interactions and ionization in the sampled warm-dense states. The code interpolates pressure and internal energy inside the table's convex hull and solves the Rankine–Hugoniot relations. It does not extrapolate to absent cold states or claim a thermodynamically consistent free-energy surface from independent interpolation. Nitrogen tables are also retained for subsequent mixed N/H work. [FPEOS database](https://militzer.berkeley.edu/FPEOS/), [hydrocarbon calculations](https://doi.org/10.1103/PhysRevE.96.013204), [nitrogen shock measurements](https://doi.org/10.1103/PhysRevLett.129.015701).

For an illustrative equal-area sector carrying parallel incoming velocity, the transverse energy available on radial merging is

`q_transverse = v_rms²/N × (1−2/(3N))`.

This diagnoses a geometry's energy scale. Actual oblique interfaces, the fraction of material initially shocked, occupancy and engineered fanning can change its deposition. The table below gives the equivalent normal-shock state if this energy is deposited in condensed material; it is not a complete spherical merger.

| Independent directions | B2 transverse energy | CH₂ shock T | C13 transverse energy | CH shock T |
|---|---:|---:|---:|---:|
| 1,000 | 6.69 GJ/kg | 846,000 K | 5.47 GJ/kg | 872,000 K |
| 10,000 | 670 MJ/kg | 120,000 K | 547 MJ/kg | 129,000 K |
| 100,000 | 67.0 MJ/kg | 12,000 K | 54.7 MJ/kg | 13,200 K |

The associated compression ratios are approximately 4.2, 3.5 and 2.6 for CH₂, and 4.3, 3.6 and 2.8 for isotope-13 CH. At one million directions the diagnostic shock falls below the available table's low-temperature coverage; no extrapolated result is printed.

These results explain why effective angular coverage is consequential. They do not establish that 1,000 directions fail or that 100,000 are optimal. Heating can supply intended preheat, and simpler hardware can tolerate some lost burn yield. Its location and density determine the entropy delivered to the hot phase. The next physical calculation must carry these contact states through compression rather than equating a small fraction of bulk KE with a small entropy penalty.

## First actual yield sensitivity

`launcher_yield_probe.py` increases the initial hot-profile temperature uniformly and reduces inward kinetic energy by exactly the corresponding EOS internal-energy increase. Total prepared energy is fixed. Corrected compression-only artificial viscosity, original EOS, radiation and the reduced nitrogen capture/reverse network are retained. All four runs complete four reference sound times.

| Initial temperature multiplier | Prompt gain | N14 capture fraction | Capture/reference |
|---|---:|---:|---:|
| 1.0 | 8.19634 | 0.539367 | 1.000 |
| 1.1 | 7.96035 | 0.523838 | 0.9712 |
| 1.2 | 7.73584 | 0.509064 | 0.9438 |
| 1.5 | 7.09637 | 0.466983 | 0.8658 |

There is no abrupt cliff over this range. The results support evaluating source complexity against whole-cycle output, rather than automatically buying enough outlets to retain 95% yield. A temperature multiplier is not yet a tolerance on a launch time, angle or grain size. The full-network and stage-specific controls, plus a physical contact-to-hot mapping, are needed before using it that way.

## Radiation, nitrogen chemistry and scheduling choices

The radiation-interface agent finds a smooth unoxidized-Al optical model gives a 1,000 K Planck absorptivity near 0.022. With that reference, a 1 mm grain traversing the full cavity can reach the 300 K polyethylene enthalpy allowance in 73–81 s; a 1 cm grain passes the same integrated-energy screen for the nominal flight. Slow material and internal gradients still need a thermal history. This favors larger grains, better spectral reflectors or a shorter exposed path. The optical and thermal calculations belong to the companion radiation packet; these are coupled design inputs.

For N:H=1, two candidate recipes remain:

- A sealed cold N₂/NH₃ mixture, with equal molecule counts giving equal total N and H atoms. Its volatile core requires a real pressure barrier and thermal path, even at low temperature.
- An isotope-selected ammonium-azide core, [NH₄]N₃, which has exactly the required atom ratio in a solid. Primary measurements document sublimation and pressure/temperature transformations. Decomposition, fragmentation and chemical energy release must enter the launch/contact entropy account. A convenient stoichiometry does not prove survival. No synthesis procedure is needed for this assessment. [Sublimation calorimetry](https://doi.org/10.1016/0040-6031(85)85381-8), [high-pressure transformations](https://pmc.ncbi.nlm.nih.gov/articles/PMC7560398/).

These nitrogen recipes are provisional because the new mixture optimization may favor a different H:catalyst ratio. Their carbon or metal skins also change isotope feed purity. The purification study is testing carbon admixtures in the nitrogen stages; coating mass must be included there and in material recovery.

A promising slow-tail change is a separate set of faster seed packets that collide near the focus to establish the central hot material. Their kinetic energy would supply the already-budgeted central internal energy, rather than being added to it. This could remove the hours-long low-velocity launch tail. Central seed mass/energy, arrival order and its pressure effect on the main implosion require a dedicated test; the pulse/schedule agent is extracting those inputs.

## Uncertainty and next discriminators

The most likely changes concern grain thermal protection/material choice, effective angular coverage, and the cold-to-hot merger. Electrical source efficiency and apparatus mass are meaningful optimization variables, not physical constants. There is insufficient evidence yet for a credible interval on the minimum outlet count: the actual contact-to-burn acceptance function is still being built. The 1,000–100,000 direction comparison is a design sweep.

For the condensed-grain architecture, a **working grain-radius estimate of 1 cm with a subjective 80% range of roughly 1 mm–5 cm** is useful only under the stated 1 GV/m surface field, 1–10 MV/m effective gradients, cold/reflective exit options and substantial manufacturing capability. It expresses which region is worth developing, not a measured survival interval. Nitrogen chemistry and changed optimized fuel ratios may move it.

Priority continuation: (1) couple nonideal contact states to compression and compare permissive lower-yield designs; (2) specify source electrical cells and a beam-envelope/neutralization sequence; (3) integrate seed/tail delivery with the recurring schedule; (4) derive timing and differential pointing tolerances from actual yield and return changes. Whole-plant mass/output and coating contamination must be checked alongside each improvement.

`launcher_screen.py` and `launcher-screen.json` contain the geometry, grain/beam alternatives, primary-data hashes and 20 current algebra/root checks. These are mathematical consistency checks, not experimental validation of the apparatus. `launcher-yield/summary.json` retains the completed hot-profile sensitivity cases and trace identities. All older source, pulse, article and station packets remain unchanged.


## Second checkpoint: actual neutralization and circuit sequence

The selected sequence for further design uses an electron reservoir linked to the charging station, conductive return paths and image-charge electrodes around each accelerator. Charge removed from fuel remains in the local circuit until the final neutralizer returns it. A source made of positively charged grains plus an otherwise uncharged apparatus would build an unwanted station-scale electric field. The beam, surrounding electrodes and return current must be modeled together.

For the illustrative 1 cm, 1,000 kg/m³ grain, the charge is 11.13 microcoulombs, or 6.94×10¹³ missing electrons, and its isolated potential is 10 MV. A co-moving electron supply avoids adding the full beam-relative electron kinetic energy. A 1 eV, 10¹⁵ m⁻³ electron population has an unenhanced isotropic geometric current of 26.8 A/m². Dividing the surface charge by that current gives **0.330 milliseconds**, corresponding to a **793 m neutralizer** at 2.4 million m/s. At 10¹⁴ m⁻³ the same calculation gives 7.93 km. These are explicit flux cases, not predicted sheath solutions.

At the 100,000-column mean current, a 100 m-radius neutralizer must supply 429 A per column; the advective electron flux through that cross-section is about 12 million A at the chosen density. Each grain consumes the electrons in about 0.0694 m³ of plasma, a sphere of radius 0.255 m, compared with a mean grain spacing of 12.5 m in this mean-duty example. Gross particle availability therefore gives no immediate veto. The geometric interception of an equally dense stationary hydrogen plasma over 793 m would add at most about 286 J/kg in directed impact energy; a co-moving supply reduces the relative energy. This estimate omits collective fields and sheath-mediated scattering.

The Debye length is 0.235 mm, much smaller than the grain. A small-grain orbital-motion-limited formula does not justify the actual collection rate. The initial 10 MV grain potential, nonlinear sheath, primary-electron penetration, secondary emission and final floating potential require a kinetic calculation. Published ion-beam neutralization studies demonstrate capture followed by slower cold-electron accumulation, and sensitivity to source position; positively emitting dust can also depart from ordinary orbital-motion-limited theory. [Lan and Kaganovich, 2019](https://doi.org/10.1063/1.5128521), [Delzanno and Tang, 2014](https://doi.org/10.1103/PhysRevLett.113.035002), [Kaganovich et al., 2010](https://nonneutral.pppl.gov/pdfpapers2010/Kaganovich.PoP.2010.pdf). Their geometries do not validate this macroscopic-grain neutralizer.

The available electric field energy is **13.3 kJ/kg**. It must appear once as deposited heat, escaping electrons or radiation. It is small compared with launch energy but appreciable beside a cold capsule's thermal allowance. The largest incoming electrons can reach MeV energies in the grain potential; using only the 1 eV supply temperature to estimate material heating would miss this energy. A final surface potential of magnitude 5 V would leave only 5×10⁻⁷ of the initial grain charge. At that explicitly assumed residual, the uniform-sphere charge-energy diagnostic at a 2 km assembly radius is only 1.87×10⁻¹¹ of bulk kinetic energy. The required result is low net charge, not a particular choice of floating potential.

Individual grains receive their last measurement and steering before the neutralizer, with the mean deflection of the neutralizer included in the control law. An exit pickup measures the remaining charge and ensemble velocity distribution. A correction after neutralization requires another charged stage or a specified mechanical/optical force. The neutralizer is not itself a promise of arbitrary post-launch steering.

Beam current also constrains the injector separately from merger geometry. For a uniform charged beam, `Eedge=I/(2π epsilon₀ v b)`. The 429 A mean example has roughly 33 kV/m self-field at 2.4 million m/s and radius 100 m, rising to 0.77 MV/m at 100 km/s. At 10 km/s it needs about 770 m radius to stay at 1 MV/m, or multiple parallel injection branches that combine after acceleration. Current peaks require their own margin. Thus a few thousand effective angular regions may still use many more physical low-current source columns. The number of accelerator columns and the number of independently converging angular regions should be separate optimization variables.

## Second checkpoint: contact entropy connected to the hot profile

`launcher_feos.py` adds an explicit, thermodynamically consistent comparison using the published Frankfurt FEOS 16.7 package. It combines Thomas–Fermi electrons, the Cowan ionic model and a bonding correction, plus equilibrium radiation here. The public source archive is preserved unchanged. We supply illustrative polymer reference density, 300 K temperature and 2 GPa bulk modulus; no liquid–vapor Maxwell construction or fitted soft-sphere model is used. This model has broad coverage but misses discrete shell structure and detailed low-temperature chemistry. [Faik, Tauschwitz and Iosilevskiy, 2018](https://doi.org/10.1016/j.cpc.2018.01.008), [public GPL source](https://doi.org/10.17632/6vjsv6v48p.1).

The calculation uses a clearly limited geometry: equal mass per solid angle in parallel-sector streams. Their transverse kinetic energy is `q=K[2/N−4/(3N²)]`, where **K is specific bulk kinetic energy in J/kg** and N counts effective angular regions. A normal shock thermalizes that energy at the chosen condensed-material reference density. The resulting FEOS entropy is then conserved while each parcel is compressed to the corresponding density in the prescribed hot starting profile. This constructs an optimistic entropy comparison with no later dissipation. It does not produce the required spatial density, momentum profile, grain breakup or three-dimensional symmetry. Changing where the merger occurs changes the result.

| Original reference stage | Effective angular regions giving the original mean entropy | Shock energy | Post-shock density | Post-shock temperature |
|---|---:|---:|---:|---:|
| B2, C12+2H | 1,941 | 3.45 GJ/kg | 3,761 kg/m³ | 0.440 MK |
| C13, C13+H | 1,676 | 3.27 GJ/kg | 4,604 kg/m³ | 0.521 MK |

The implied carbon charge is about 3.2–3.3, not six, and hydrogen is about 0.83–0.85 ionized in this model. The earlier fully stripped contact assumption was therefore qualitatively unsuitable. Across N=100–100,000, FEOS shock temperatures differ from the primary FPEOS interpolation by about −11% to +11%, and compression by about −7% to +1%. This agreement is useful for an initial entropy screen, not an error bar on the subsequent isentrope. The six sampled first-law/entropy derivative tests agree to better than 10⁻⁷ relative error; they check the implementation's thermodynamic consistency, not the EOS's physical accuracy.

At N=1,000, the mapped mean initial temperature is **1.70× reference for B2 and 1.50× for C13**. At N=2,000 the ratios are 0.98 and 0.88; at N=10,000 they are 0.30 and 0.29. A lower-entropy stream needs additional deliberately placed preheat to reproduce the reference. The controlled heat would come from the already budgeted prepared energy. A higher-entropy stream can instead be accepted with a different compression and burn yield. The completed A sensitivity calculation above establishes that such a trade can be gradual, although its yield numbers cannot be transferred directly to B2 or C13.

This result supports testing a few thousand effective directions before selecting a much larger number. It does not establish a minimum number of launcher outlets. A useful next simulation is a small converging sector with several resolved neighboring streams, condensed or nonideal contact EOS, and prescribed total energy; measure the entropy distribution and coherent momentum reaching the hot starting density. Its output can feed the existing full-network burn calculation. The two independent steps are then tested rather than presumed identical.

## Cold-source burden when the hydrogen ratio changes

The parallel mixture study varies H/catalyst; the recipes above must vary with it. For nitrogen, mixtures of sealed NH3 and N2 can realize H/N up to three without free hydrogen: the fraction of nitrogen nuclei in ammonia is H/N divided by three. This can remove the need for an energetic ammonium-azide recipe in the currently promising H/N≈2–3 region. Both phases still need a mechanically qualified container and a thermal history.

For carbon, a simple accounting recipe uses CH2 plus catalyst-rich material below H/C=2 and CH2 plus condensed H2 above it. This is a scalable illustration, not the best possible chemistry. CH4-rich and CH2/CH4 mixtures can cover H/C=2–4 without free 20 K hydrogen and deserve comparison if the hot mixture optimum lies there. Molecular and isotopic homogenization must be achieved by the contact/compression stage; separate inclusions do not automatically give a homogeneous nuclear mixture.

The additive-volume screen takes a chosen 1,000 kg/m³ carbon matrix, a chosen 850 kg/m³ nitrogen matrix and liquid hydrogen at 70.8 kg/m³ near 20.3 K. The last value is a physical property; the matrix values are explicit screening inputs. [NASA cryogenic-storage demonstration](https://www.energy.gov/sites/default/files/2021-10/lh2-storage-handling-demonstrations.pdf).

| Illustrative recipe | Excess H2 mass fraction | H2 share of volume | Overall density | IR dose for the same 1 cm Al-coated example |
|---|---:|---:|---:|---:|
| C12:H=1:2 | 0 | 0 | 1,000 kg/m³ | 170 kJ/kg |
| C12:H=1:3 | 6.67% | 50.2% | 533 kg/m³ | 319 kJ/kg |
| C12:H=1:4 | 12.5% | 66.9% | 379 kg/m³ | 449 kJ/kg |
| N14:H=1:2–3 | 0 | 0 | 850 kg/m³ | 200 kJ/kg |
| N14:H=1:4 | 5.56% | 41.4% | 527 kg/m³ | 323 kJ/kg |

The table uses the earlier room-temperature optical reference, full 1.1 billion m flight, a 1,000 K radiating cavity and 2.4 million m/s speed. It is absorbed energy, not a predicted fuel temperature. Dense H2 heated above its low-pressure phase range requires pressure containment or permitted expansion, so a universal cold-grain temperature cannot close the design. Larger grains, lower-absorptivity skins, protected flight and alternative molecular recipes are competing choices. The H/C=16 boundary would have density about 132 kg/m³ and 93% H2 volume, but current hot-grid results do not favor selecting that recipe.

Lower density does not simply make acceleration harder. At fixed grain radius it raises q/m. Taking surface field proportional to the square root of density holds stored field energy per fuel mass constant; the acceleration length then scales as the square root of density at fixed accelerating gradient. On the other hand, a skin that alone carries a fixed electrostatic pressure becomes more expensive per kilogram of low-density core. Thermal survival, mechanical containment, optical skin impurities and usable source density should be optimized together.

The second-checkpoint artifacts are `launcher_neutralizer.py` / `launcher-neutralizer.json` (156 algebra and range checks) and `launcher_feos.py` / `launcher-feos.json` (50 shock/root, positivity, local-energy closure and thermodynamic checks). Input/archive/code hashes are recorded. The upstream FEOS source emits legacy compiler warnings, including mismatched array allocation/deletion; these bounded state calculations do not constitute numerical-library qualification. Existing hot trajectories and the concurrently optimized fuel-ratio family remain separate evidence.

## Accelerator cells and a finite space-charge envelope

A long, continuously filled beam cannot gain 10¹⁵ volts of electrostatic energy while the apparatus merely resets an unspecified local voltage. The concrete alternative is **phase-bunched acceleration in independently timed gaps**, with drift shielding during the decelerating part of each cycle. The familiar linac mechanism permits repeated work without one terminal at the integrated voltage. At 2.4 million m/s, illustrative 100 m cells imply kilohertz timing; the exact cell spacing, synchronous phase, bunch length and beam-loading circuit must be solved together. The 1–10 MV/m numbers denote an effective gradient after the duty/phase factor. An assumed 50% longitudinal bucket filling doubles instantaneous current. [CERN introduction to RF linear accelerators](https://cds.cern.ch/record/865936/files/p53.pdf).

`launcher_envelope.py` now integrates a specific transverse model, rather than assigning ideal focusing. It solves the cold, uniformly filled elliptical-beam envelope in alternating electrostatic quadrupoles, including linear space-charge defocusing. For half-period focusing/defocusing gradients ±G, the equations are `a″=−kappa a+2K/(a+b)` and `b″=+kappa b+2K/(a+b)`, with `kappa=(q/m)G/v²` and `K=(q/m)I/(2π epsilon₀ v³)`. The script finds a periodic solution and perturbs all four envelope initial coordinates to check its local Floquet stability. This is the zero-emittance, constant-speed slice of a KV-type model; it does not establish halo control, longitudinal phase stability, grain-size/charge tolerance, stability of the actual particle distribution or the full accelerating transition. [Qin et al., generalized KV theory](https://arxiv.org/abs/1610.00126), [Lund and Bukh, envelope stability](https://doi.org/10.1103/PhysRevSTAB.7.024801).

At G=10,000 V/m² and `kappa L²=4`, the illustrative 1 cm source gives:

| Speed | Full quadrupole period | Maximum envelope radius, 429 A mean | Maximum radius at 8× current |
|---|---:|---:|---:|
| 10 km/s | 3.88 km | 1.08 km | 3.06 km |
| 100 km/s | 38.8 km | 343 m | 969 m |
| 500 km/s | 194 km | 153 m | 433 m |
| 2,400 km/s | 931 km | 69.9 m | 198 m |

The 8× case is a deliberately explicit example: four times the mean pulse mass current and 50% bucket filling. It has not been fitted to the complete source waveform. The final bare beam needs about 396 m diameter before aperture clearance; a roughly 450–500 m outlet is a candidate for the current guide geometry. Its quadrupole field at the beam edge is about 2 MV/m. Upstream injection requires a much wider region or many small branches, which can sit behind the receiver rather than passing through the same outlet hole. Quadrupole electrodes would be segmented mechanically and electrically; the model's long periodic focusing region is not a claim that one solid rod spans its full length.

The computed periodicity residual is below 10⁻⁸ in normalized coordinates, and the envelope perturbation eigenvalues lie on the unit circle within the saved numerical tolerance. Eighteen checks pass. These results establish a consistent transverse-envelope example at the proposed current and aperture scale. They provide the next circuit/particle-transport study with an actual starting geometry. The missing acceleration ramp, grain charge distribution, RF bucket acceptance and halo are substantial design tasks; the calculation does not convert the old source mass allocation into a manufactured accelerator.

## Updated preparation priorities

The most consequential unresolved step is now the three-dimensional contact/merger, followed by the coupled cold-capsule thermal and mechanical design. The new scalar entropy result makes a huge mandatory outlet count less credible, while the source envelope shows why physical transport columns can still outnumber the angular regions. A central seed retiming concept also has a concrete energy budget: the inner 0.1% of each large reference pulse contains about 1.13×10¹⁰ kg, with 3.55–4.35×10²¹ J of prescribed internal energy. Opposing packets at 0.79–0.88 million m/s can supply that heat from kinetic energy and arrive in 21–23 minutes. These values come from the companion schedule/EOS calculation. The seed's density, pressure and timing still must match the incoming main pulse; no extra energy is required merely by replacing the slow central flight.

The methane option illustrates why a chemical formula alone is insufficient. Primary fixed-density calorimetry near 439 kg/m³ shows that warming methane from about 102.4 to 117.4 K raises pressure from 4.84 to 28.25 MPa while absorbing only about 32 kJ/kg. A 100 nm optical film cannot contain that pressure on a 1 cm capsule. A structural shell, expansion volume, colder/protected flight or larger grain must be included. At a genuinely qualified 10 GPa allowable carbon shell stress, the thin-shell pressure formula would need a few percent shell mass for that pressure scale; that carbon changes H/C and belongs in the fuel recipe. At 1 GPa the mass is much larger. This is an available design trade, not proof that methane is excluded. [Younglove, 1974 methane calorimetry](https://pmc.ncbi.nlm.nih.gov/articles/PMC6742817/).

No narrow credible interval on minimum outlet count or acceptable arrival-time error is justified yet: we have a new thermodynamic screening scale and an envelope solution, but no calibrated three-dimensional yield response. The earlier 1 cm grain-size working estimate remains conditional. Quantities that can change the architecture with an appreciable chance include a higher-than-expected merger entropy, capsule heating/fragmentation, RF phase/charge spread causing unacceptable halo, and a neutralizer heat/charge tail. Each now has a concrete next calculation: sector merger; sealed-core thermal/mechanical trajectory; longitudinal/transverse particle tracking through the accelerating lattice; and nonlinear electron/grain neutralization. The user's optimized H/catalyst ratios should feed all four, rather than freezing an obsolete molecular recipe.


## Shell-dependent entropy and an actual mapped-profile burn

The pulse-mean entropy calculation is useful but hides the radial speed distribution. The refined screen assigns each cell its consistent share of nodal kinetic energy, `K_i=(v_inner²+v_outer²)/4`, and uses `q_i=K_i[2/N−4/(3N²)]`. Cell-weighted kinetic energy exactly reproduces the original total. The inner 0.1% mass is reserved for the separate seed concept; its original temperature is retained.

At fixed N=2,000, the C13 mapped initial temperatures have mass-weighted 10th/50th/90th ratios **0.434/0.809/1.355** to the reference. The B2 ratios are **0.465/0.899/1.537**. Uniform angular coverage therefore does not make uniform entropy. Matching the old reference entropy instead would require a time-dependent effective angular count: roughly 611/1,502/2,956 over those C13 mass percentiles, and 705/1,741/3,429 for B2. These use the existing hot-profile speeds as a delivery proxy; the actual cold velocity must also supply compressive work.

The next discriminator has now been run: **full-network, corrected-viscosity C13, 64 radial cells, N=2,000 mapped temperature profile, fixed total prepared energy, original mass/density profile, opacity 0.001 m²/kg and four original sound times**. FastEOS evaluates both original and new internal energy; a uniform 1.009247 velocity multiplier restores total prepared energy. The exact isotope abundances, mapped profile and time convention are archived.

| Quantity | Original corrected C13 | Mapped entropy profile |
|---|---:|---:|
| Nuclear/initial-energy gain | 15.538979 | 15.529991 |
| Desired N14 per initial C13 | 0.8166214 | 0.8177223 |
| Physical simulated duration | 70.0036 microseconds | 70.0036 microseconds |

Both calculations complete. Gain changes by **−0.0578%**, and desired-product yield by **+0.1348%**; the mapped run's energy residual is −3.66×10⁻⁶ of input. Seven direct comparison/termination checks pass. The small gain/yield difference should not be interpreted as a precisely optimized improvement at this resolution. The substantive result is that this particular broad entropy profile retains the reference burn performance in the radial model. Reproducing one uniform target adiabat is therefore not established as necessary.

The result still starts with the prescribed hot density and radial momentum shape. It does not show that finite streams produce that density, preserve symmetry, or pass through the source/neutralizer. `launcher_mapped_yield.py` and `launcher-mapped-yield/` retain the calculation and inputs; the raw trace path and hash are in its summary. A three-dimensional sector assembly remains the most valuable preparation calculation.

## Molecular recipes for the current mixture shortlist

The concurrent process optimization's provisional H/catalyst choices admit specific molecular feed candidates: **NH3 for N14:H=1:3; N2H4 hydrazine for N15:H=1:2; C4H10 butane for C12:H=1:2.5; and CH4 methane for C13:H=1:4**. These are source candidates, not selected qualified capsules. The optimized ratios remain subject to mesh and process checks.

Hydrazine is especially relevant because it gives H/N=2 at a density near 1,000 kg/m³ without requiring liquid hydrogen or nitrogen. It freezes near 275 K. NASA's actual freeze/thaw experience warns that phase-volume changes can damage confined hardware; a one-use capsule's initial fill volume and phase trajectory must be designed explicitly. Its liquid formation enthalpy is about +50.74 kJ/mol, corresponding to about 1.58 MJ/kg released to N2+2H2 at the same reference conditions. That is small beside several GJ/kg of intentional merger heat, but premature decomposition can fragment the grain and spoil its trajectory. Surface compatibility, irradiation/electron capture and acceleration are therefore real qualification tests. [NASA freeze–thaw technical bulletin](https://www.nasa.gov/general/risk-of-hydrazine-use-following-freeze-thaw-exposure/), [Argonne Active Thermochemical Tables](https://atct.anl.gov/Thermochemical%20Data/version%201.222/reaction/?c1=-1&c2=-1&c3=1&species_number%5B%5D=2028&species_number%5B%5D=2363&species_number%5B%5D=24), [NIST hydrazine thermochemistry](https://janaf.nist.gov/tables/H-090.html). No chemical preparation procedure is needed for this comparison.

For H/C between two and four, alkane/blend choice avoids treating the CH2+H2 illustration as a compulsory penalty. Butane gives exactly 2.5 before coatings. A carbon structural skin on methane shifts the whole-grain H/C below four; the actual isotope-specific skin belongs in the nuclear mixture optimization. Nitrogen capsules need an explicit impurity allowance if their skins contain carbon or another element. Thermal protection, capsule structure and the full fuel recipe must consequently be chosen together.

## Selected-mixture source profiles and a constructive preheat adjustment

The same calculation has now been repeated for the provisional optimized C13:H=1:4 and N14:H=1:3 mixtures, using their actual initial abundances, hot density profiles, prepared energies and physical durations. The new methane EOS uses the C13+4H inertial mass of 17 u. Its chosen cold density is 439×17/16=466.4 kg/m³, an equal-molecular-volume isotope scaling of the measured methane density; its reference temperature is 102.4 K. The original FPEOS methane table has a denser cold reference, so the comparison resets its initial volume to this chosen density while retaining its molecular ground-state energy. The omitted cold cohesive/pressure-energy adjustment is small on the GJ/kg shock scale but is recorded. At the mean N=2,000 shock, FEOS predicts a temperature 5.2% below the direct methane FPEOS interpolation and a compression 8.2% lower.

The ammonia mapping uses a chosen 850 kg/m³, 150 K reference with 2 GPa bulk modulus. There is no direct NH3 FPEOS table in the downloaded archive. A separate additive-volume comparison uses pure nitrogen and hydrogen first-principles tables at common pressure and temperature. At the mapped mean contact state, FEOS pressure is 12.8% above that mixture approximation and thermal energy 5.5% above it. This supports a bounded exploratory sensitivity; it does not qualify ammonia's cold chemistry or replace a direct mixed-material EOS.

For N=2,000, the C13:H4 mapped temperature ratios have mass 10th/50th/90th percentiles 0.490/0.969/1.677. The full 27-isotope burn completes at the unchanged prepared energy and 70.0036 microsecond duration. Its gain is **22.88617**, compared with **23.65567** for the optimized reference, a 3.253% decrease. Final N14 per original catalyst rises from 0.222027 to 0.252914. The hydrogen-rich network proceeds beyond N14, so this final product is not a scalar substitute for total capture or whole-cycle performance. The complete isotope endpoints must pass through the material ledger.

The initial ammonia mapping exposed a useful failure of the **chosen profile**, rather than a physical temperature limit: it reached the existing 2.8 GK numerical/EOS guard after 33.604 microseconds of a requested 68.395 microseconds. The hottest cell enclosed 0.0103% of total mass, and 0.278% of mass exceeded 2 GK. Its partial gain is preserved separately and is excluded from the selected operating family.

The constructive next case adds controlled initial preheat wherever the mapped temperature lies below the original reference. This affects 70.56% of mass, with all energy debited from coherent radial motion. It adds 74.15 GJ/kg compared with the cooler uncontrolled mapping; relative to the original optimized reference, final internal energy rises by only **4.47 GJ/kg**, about **0.15% of total prepared energy**. Initial radial speeds fall to 0.9992075 of reference. Such an adjustment corresponds to varying merger angles or placing additional heating in the underheated shells. It is an initial source prescription, with no clipping of simulated hot temperatures and no raised numerical guard.

The revised ammonia case completes the full physical duration. Gain is **10.25252**, compared with **10.28882** for the optimized reference, a **0.3528% decrease**; desired O15 changes by the same relative amount. Its energy residual is 4.29×10⁻⁷ of input. This demonstrates one inexpensive adjustment that avoids the mapped profile's excessive core compression within the model's existing domain. Actual angular-control hardware and three-dimensional assembly remain to be demonstrated.

The selected family is now supplied as `selected-source-stage-map.json`, with the complete A:H3 controlled-preheat and C13:H4 mapped cases plus unchanged B1:H2 and B2:H2.5 reference cases. It embeds exact starting profiles and abundances, records physical duration and raw-trace hashes, and gives portable worktree-relative result paths. `selected-source-paths.json` is the flat loader-compatible map. Thirty-eight checks independently reproduce duration, composition, baryon conservation and the energy ledger. Compact profiles/results are archived; raw trace paths are explicit. The guarded A run remains a documented incomplete numerical continuation, not a failed reactor design.

## Capsule design still couples size, charge and thermal transport

A larger methane capsule can trade source length for radiative tolerance. At the 466.4 kg/m³ reference density, full-flight absorption under the existing Al optical proxy is about 365 kJ/kg for a 1 cm radius, 73 kJ/kg for 5 cm, and 45.6 kJ/kg for 8 cm. At 1 GV/m surface field its stored electrical energy is another 28.5 kJ/kg if fully retained on neutralization. Reducing surface field to 0.5 GV/m lowers that term to 7.1 kJ/kg, but doubles the required accelerator length at fixed grain radius. For example, an 8 cm capsule at 0.5 GV/m and 10 MV/m effective acceleration needs roughly 0.81 billion m of accelerating path at 2.4 million m/s. This is a concrete alternative cost, not a free thermal improvement.

A partly frozen capsule with expansion space can use sensible and latent heat before reaching the dense-liquid pressure path. Its required enthalpy allowance, headspace and local temperature gradients must be demonstrated with the actual isotope and skin. Bulk absorbed energy alone is insufficient: the conductive path from the skin to the interior can bottleneck heat uptake, and melting must preserve containment and charge control. Under the current optical proxy, the surface flux is about 1,240 W/m². A quasi-steady internally warming sphere has surface-to-center temperature difference approximately `q_surface a/(2 k_eff)`; an 8 cm capsule needs effective radial conductivity near 5 W/(m K) to hold that difference to about 10 K. This is an engineering requirement, not a measured conductivity of methane. Internal conductive structure or a different thermal path must be included in both mixture and mass accounts if needed.

The finite source-envelope calculation was independently checked: Gauss's law gives its stated elliptical-beam force and dimensionless normalization. For the peak-current example, the quadrupole potential scale at a 198 m beam radius is about **196 MV**, with a 2 MV/m edge field and half-period focusing region of about 466 km. These are substantial electrical layout requirements even though the local field is below the chosen effective accelerating gradient. The independent check supports the mathematical transport example; it does not establish a manufactured, phase-stable accelerator.
