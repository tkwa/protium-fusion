# Engineering hypotheses to test ✨

These are dimensional leads for further calculations. None is a selected plant design or an established performance claim.

## Pulse capture and recovery

The computed main pulses release several 1e29 J. Their expanding gas converts most initial internal energy and nuclear heat into bulk kinetic energy before radiation escapes. A remote magnetic catcher is therefore worth comparing with a thermal absorber. Published inductive-conversion studies for other fuels motivate this branch; they do not establish its efficiency or stability here.

An equivalent spherical volume of radius 1e8 m with RMS magnetic field 1 T stores about 1.67e30 J. This is an energy-volume comparison, not a realizable spherical field topology. Coil geometry, magnetic escape paths and plasma interchange modes need explicit treatment. A structural virial estimate scales as field energy divided by allowable specific stress, multiplied by a geometry factor. At 1 MJ/kg specific stress and an illustrative factor of five, the field structure alone would be about 8e24 kg, below Jupiter mass. An actual mass budget must also count conductors, shielding, pulse storage, drivers, support, heat rejection and processing.

At such masses, reaching 1 mW/kg might require only one enormous cycle every few months or years. Cold beta-decay residence times would then be a small part of the turnaround time. Pulse absorption is much faster than the average power suggests: gas moving near 1e7 m/s reaches a 1e8 m catcher in about ten seconds. Field/current transients, conductor voltage, local mechanical response, plasma stability and irreversible losses must be calculated against that arrival time. An average generator rating cannot substitute for pulse handling.

For a hypothetical cycle with deposited nuclear energy Q and fuel preparation energy U, thermal conversion efficiency eta_h and wall-plug-to-fuel efficiency eta_d give ideal net work `eta_h*(Q+U)-U/eta_d`, before processing and other auxiliaries. A direct converter would need separate efficiencies for bulk kinetic energy, residual heat and decay photons. Do not assign all components a published efficiency from a different geometry.

## Whole-system specific power

Mass that scales with average output should enter self-consistently. If fixed pulse equipment has mass M0, output-dependent equipment has mass k times net power P, and the target specific power is s, then `M=M0/(1-s*k)` when `s*k<1`. If `s*k>=1`, that assumed balance of plant cannot meet the target regardless of fixed-equipment mass. Radiator mass follows waste-power ratio, emissivity, temperature and areal density; conversion and driver masses follow their gross-power ratings, not net output alone.

The repeated-cycle ledger is now in `cycle_ledger.py`. `plant_budget.py` gives a conditional requirements allocation, described in `PLANT-BUDGET.md`; it does not derive the driver or processing facility from an engineered design.

## Preparing the convergent flow

The compression experiment specifies an initially inward-moving dense plasma. Plausible physical implementations to compare include convergent liners, synchronized neutralized plasma streams, and staged compression with a late heating shock. Their entropy, energy coupling, fuel mixing and asymmetry could materially worsen the ideal result. External beams cannot be assumed to heat the center of an already opaque 300 m fuel sphere uniformly. A future driver must create the state through a physically described compression/heating history.

## Catalyst replacement

Replacement carbon could be manufactured from helium ash. The initial uniform pure-helium screens cost about 408 MeV of fuel preparation energy per carbon atom for the 100 m case and about 1401 MeV for the 30 m case. Specified convergent-flow calculations at 64 cells improve these to approximately 231 and 763 MeV. The 128-cell 100 m calculation increases carbon yield further; its peak temperature and yield are not yet converged. The cycle ledger explicitly recycles unburned helium and charges manufacture. Trace loss rates, real separation, inflow preparation and infrastructure remain to be checked.

## Enlarging the collector to limit transient heat flux

A thermal collector at radius 1e10 m deserves an explicit transient calculation. A 5e29 J pulse spread over 4*pi*R^2 deposits roughly 4e8 J/m2. A broad velocity distribution near 1e7 m/s spreads arrival over approximately 1e3 s, suggesting flux of order 4e5 W/m2. Actual arrival flux must be derived from the saved velocity/energy profiles. Cold-decay gamma energy arrives sooner, over the beta-decay times rather than the gas flight time, and needs a separate deposition calculation.

At an illustrative areal mass 1e4 kg/m2, such a collector has mass 1.26e25 kg. A per-kilogram heat capacity of 1 kJ/K would give a mean pulse temperature rise of about 40 K if the energy were distributed through its whole thickness. It is initially deposited near the surface, so that average is insufficient: solve surface conduction for ions and volumetric gamma deposition using attenuation data. A high-conductivity refractory absorber might handle the calculated flux, but this has not yet been checked.

A freely suspended rigid spherical shell at this scale would buckle under its own gravity. At fixed areal mass Sigma, self-gravitational pressure is 2*pi*G*Sigma^2, about 0.042 Pa for Sigma=1e4 kg/m2. Support must be explicit. Candidates include magnetic forces with suitable coil geometry or actively supported structures using circulating mass streams. For total supported mass M at radius R, the stream kinetic-energy scale is G*M^2/(4*R), and the idealized support mass at stream speed v is G*M^2/(2*R*v^2). Losses per circuit become a persistent parasitic load. Existing Infinite Wiki guidance already treats active orbital rings as speculative momentum-conserving engineering requiring containment, bearings, power and redundancy (Route 66 brief).

Ion implantation, collector erosion, catalyst extraction, contamination and recovery fraction matter as well as heat. Increasing area reduces surface fluence but does not by itself guarantee material recovery or economical construction.

## Thermal collector checkpoint

`collector_screen.py` now reconstructs a piecewise velocity distribution from the spatial pulse results. All examined outgoing profiles have increasing radial velocity, so they do not immediately form ballistic shell crossings. It normalizes the reconstructed kinetic energy to the hydro result and convolves the arriving heat flux with the constant-property semi-infinite conduction kernel. It treats remaining hot internal energy as following the same arrival distribution, separately computes prompt radiation from the saved escape history, and treats cold beta-decay energy as surface heat. Photon transport, erosion, temperature-dependent properties and nonlinear material response remain absent.

At radius 2e10 m, the illustrative carbon constants (k=100 W/m/K, rho=2000 kg/m3, cp=1000 J/kg/K) give a sum of separate peak temperature rises of about 245–286 K for A, 525 K for the 128-cell B1 pulse, 440 K for B2, and 246 K for the larger carbon factory. This sum is an upper superposition screen, not the actual temperature history. The prompt B1 contribution dominates its result and remains sensitive to radiative boundary resolution. The smaller 30 m carbon factory produces a much gentler pulse, at higher preparation cost per carbon atom. At radius 1e10 m the same screen gives rises above 2000 K for some stages.

The collector's inner surface could exchange thermal radiation across a closed cavity. Its outer surface needs insulation and designated heat rejection panels: freely radiating the entire absorber at 1000 K would overwhelm the cycle's average energy supply. The base allocation includes 100 kg/m2 for insulation but does not design it.

## Active support and the surrounding gravitational field

For a thin self-gravitating sphere, W=-G*M^2/(2R). Circulating support streams can supply 2K=-W, with K=G*M^2/(4R), while guides transfer their outward reaction to the stationary absorber. This is an ideal virial and momentum balance. It requires realizable ring geometry, local magnetic bearings, momentum control, collision avoidance and redundancy. A rigid unstressed shell does not provide the same support.

Do not automatically charge K/(material specific strength) as a second global containment mass: gravity supplies the opposing force in this ideal active structure. Conversely, the ideal virial balance does not prove that the allocated guide mass works. Support losses are explicitly charged per stream circuit in the budget.

The current allocation assumes an isolated machine. A solar-mass star's tidal field would dominate this very diffuse object near 1 AU. Its simple Hill scale is several AU; placing it tens to hundreds of AU away reduces stellar tides, but a real orbital arrangement still needs analysis.

## Driver delivery is the principal unresolved requirement

A power budget for average electricity is insufficient. The fuel driver must deliver the prepared pulse energy with the appropriate timing, entropy and angular symmetry. Its peak power, energy storage, focusing and distributed mass must be counted. Locating an Earth mass of conventional apparatus within a radius of a few thousand kilometres would create a separate self-gravity problem; the large fixed driver allocation is assumed distributed over the outer structure.

Possible delivery approaches include solid pellets that remain mechanically bound during flight, or neutralized heavy-ion beamlets with large focusing apertures. A freely streaming thermal gas launched across 2e10 m cannot be assumed to focus into a few hundred metres: its transverse emittance would spread it. Bound pellets avoid that particular thermal spreading, but acceleration, vaporization, mixing and collision entropy become critical. Neutralized beams need actual phase-space transport and collective-stability calculations. Neither route has yet been shown to create the specified convergent-flow initial condition at 30% electrical-to-fuel efficiency.

The source state can also be extended backwards to a larger, colder converging sphere, rather than postulating uniform remote heating of the opaque hot core. Its ionization, shock formation and entropy must be modeled. An exact velocity reversal of a numerical expansion is a diagnostic, not a physical driver prescription.

## Collection and processing

The base allocation requires recovery of 99.99% of circulating catalyst after each stage. At collector radius 2e10 m the incident catalyst mass per unit area is only a few milligrams/m2 per main pulse. Ions stop near the surface, so a recoverable coating, reactive getter or circulating liquid layer could be preferable to processing the whole absorber. The materials and recovery process remain unspecified. Re-emitted neutral gas would cross the chamber much more slowly than the original MeV ions; storage of such gas and its interaction with later pulses must be counted.

Thermodynamic isotope-separation work is small compared with MeV nuclear energies, but heating and moving a vastly larger carrier mass can be costly. The 10% reserve in `plant_budget.py` is an allowance to test against a process design, not evidence that processing fits within it.

## Alternative cold-envelope catcher

An engineered sub-Jupiter hydrogen/helium envelope could provide gravitational containment of pulse debris and thermal storage while remaining too cold for continuous stellar burning. A dimensional example, M around 1e27 kg and R around 1e9–3e9 m, has binding energy much greater than one 1e29–1e30 J pulse. This might replace a huge active empty collector. It needs a hydrostatic envelope/opacity calculation, transport and power extraction, a way to drive the central pulses, and catalyst capture before dilution into the envelope. It is an untested branch, not a selected architecture. Natural bulk gas would have a different construction-energy cost from manufactured collector hardware.
