# ✨ Warm local ion pumps for residual helium

A warm, electrically active cleanup layer has a concrete gas-transfer candidate: local electron-impact ionization, followed by controlled-energy implantation into the routed graphite laminate and collection in an isolated plenum. The ionization and acceleration energy are modest. **Electron confinement, protected admission and the magnetic field's interaction with outgoing plasma remain the decisive open requirements.** A static broad superconducting field changes the ballistic receiver; its switching or localized return geometry must be included before claiming a working cleanup system.

The calculation identifies useful engineering targets: confine an externally supplied electron population for roughly0.1–1 s or recover its loss energy, provide local magnetic/electrostatic fields with low leakage into the chamber, and form sufficiently strong helium traps or transfer the implanted atoms promptly. Reproducing an ordinary historical sputter pump at this scale would be excessively massive. A future integrated module would need a different mass and conductance design.

## Required cleanup rate

The exact delivered-B1/C13H2 plant emits4.150894×10³² He atoms/s before local routing. For a residual fraction f=.001, a1000 K cavity, and geometric wall-removal probability g=.01:

`Γ = 3 g v̄/(4R)`, `nHe = f QHe/(Γ V)`, `PHe = nHe kT`.

The result is nHe=4.74778×10⁹ m⁻³, PHe=6.55502×10⁻¹¹ Pa, and mean residence6.37699×10⁷ s. The required net capture is2.72990×10¹⁰ atoms/(m² s), equivalent to4.37 nA/m² and5.74984 m³/(m² s) pumping speed. Whole-sphere flow is2758.88 kg/s; its size reflects the1.52053×10¹⁹ m² receiver area.

These reproduce the gas agent's ledger. The pressure is a stationary mean; pulsed release and finite transport need the separate gas model. Temporary implantation followed by reemission into the chamber does not supply the required Γ.

## What existing pumps establish

The [NASA1977 experimental report](https://ntrs.nasa.gov/api/citations/19770011041/downloads/19770011041.pdf) studied a1.5 kg flight pump with0.11 T field. At10⁻⁶ torr it measured helium speed0.6 L/s and~0.35 pumped atoms per discharge ion; helium speed sometimes fell by a factor two after nominally similar conditioning. Loading, surface treatment and low-pressure discharge behavior mattered. It does not measure1000 K performance. A dimensional replication of0.6 L/s per1.5 kg would require**14375 kg of pump per receiver m²** for the present g target, versus2196 kg/m² for the entire current plant ledger.

The [CERN JUAS lecture](https://indico.cern.ch/event/471931/contributions/1149644/attachments/1226014/1794773/JUAS_Vacuum_Technology.pdf), pp.163–176, describes3–7 kV Penning cells,~0.1 T fields and noble-gas implantation/burial. Its nominal He/N₂ speed ratios are0.1 for diode and0.3 for triode pumps. A speed ratio is not a wall sticking probability or a hot-retention efficiency. It also describes release of previously implanted gas and pressure-dependent performance.

The retrieved PDFs/text and inspected evidence pages are hashed in[sources.json](sources.json). The NASA table and tungsten-trap page were visually checked. The citations initially failed in the web viewer; direct primary downloads succeeded.

## One concrete candidate path

A grounded, warm conductive entrance and labyrinth shield each local pump from the chamber. Gas enters recessed Penning-like cells; their hollow anodes, two collection cathodes and local magnetic circuit remain behind the shield. Externally injected electrons replenish the ionizing population at pressures where a self-sustained discharge cannot be assumed. A floating injector prepares the intended electron energy relative to the cloud potential; energetic cathode secondaries need their own loss/recovery route.

Helium ions accelerate into electrically negative collector surfaces inside the shield. Those surfaces are thin refractory cathode cassettes or locally exchangeable patches. A loaded patch is isolated from the chamber, then transferred or switched into a regeneration enclosure feeding a closed compression train. This supplies an actual destination for recovered gas. The return electron from neutralizing He⁺ and the emitter/collector currents close through the local electrical supply. The scheme does not require charging the entire receiver or introducing an electron cloud into the source-packet flight volume.

The cathode may use implantation-created or stabilized vacancies in a refractory coating. It must qualify both initial trapping and eventual transfer. Pump entrances can be recessed from direct outgoing MeV trajectories; their conductance, pulse protection, contamination and thermal loads must be counted. If protection requires closing them, live duty enters the capture requirement explicitly.

The later routed-graphite cathode calculation below is the preferred bounded transfer candidate because it shares the primary catcher's material qualification. The refractory-trap calculation remains an alternative, with its own capture and regeneration requirements.

## Ionizing column, fields and electrical circuit

The bounded model is a planar layer with independently thermalized incoming/outgoing normal velocities. For a flux-weighted1000 K Maxwellian, x=vnormal/√(kT/mHe) has density x exp(−x²/2). With β=ne⟨σve⟩L/√(kT/mHe):

`Pion = 1 − [∫₀∞ x exp(−x²/2−β/x) dx]²`,

`g = a C D η Pion`.

Here a is physical active area fraction, C entrance conductance/admission, D live duty, and η net transfer probability conditional on ionization. The two-crossing geometry is a column discriminator, not a solved field/trajectory model of the proposed cells. A measured residence-time distribution can replace it.

Using σ=3.5×10⁻²¹ m² near100 eV gives⟨σve⟩=2.07584×10⁻¹⁴ m³/s. The cross-section scale follows the primary report and[NIST helium ionization work](https://www.nist.gov/publications/cross-sections-singly-differential-and-total-ionization-helium-electron-impact). A monochromatic100 eV distribution is an imposed preparation condition; a real electron distribution needs its own averaged coefficient.

For a=.1, L=.1 m, η=.5 and D=1:

|Entrance admission C|Required Pion|ne, m⁻³|Chosen B, T|Ideal magnet mass diagnostic, kg/m²|Electron replenishment at1 keV/0.1 s, W/m²|
|---|---:|---:|---:|---:|---:|
|1|0.20|6.72×10¹⁶|0.200|23.9|1.08|
|0.5|0.40|1.64×10¹⁷|0.259|40.2|2.62|
|0.25|0.80|5.96×10¹⁷|0.495|146.5|9.55|

These are conditional design requirements. The quoted B satisfies the cold single-species necessary bound `ne≤ε₀B²/(2me)` at≤50% of that limit; finite warm-cloud stability and losses remain uncalculated. At C=1, a1 mm cloud radius gives100 eV gyroradius0.169 mm and axis-to-edge space-charge potential304 V. This is compatible in scale with kilovolt electrodes, while requiring an actual three-dimensional field and injection solution. No m=1 plasma response from the earlier ion reactor is imported here.

The ideal magnet diagnostic uses field energy U=B²aL/(2μ₀) and the ideal load-line relation `Mmag≥2U ρmag/(BH)max`, with **assumed**(BH)max100 kJ/m³ andρmag7500 kg/m³. At C=1, U=159 J/m². Return flux, fringing, mechanical support and cooling add mass; there is no claim this energy product is available at1000 K. Magnets could occupy a thermally shielded cooler part of the receiver while the gas-facing parts remain warm. Individually winding every1 mm cell with warm resistive wire gives~61 MW/m² in the illustrative winding model and is excluded for that geometry. Shared magnetic circuits, cooled permanent magnets or superconducting windings need their own engineering ledger.

The net ionization threshold costs only2.15×10⁻⁷ W/m² atη=.5. Accelerating every formed He ion through10 kV costs8.75×10⁻⁵ W/m²; the ion target current is8.75 nA/m². In a prescribed1 cm,10 kV gap, the ballistic transit is28.8 ns. These small numbers do not describe the complete pump power.

The larger electron inventory is6.72×10¹⁴/m² at C=1. Its100 eV kinetic energy is0.0108 J/m², and its internal self-field contributes another0.00818 J/m² before external-gap energy. At the gas density here, an electron's mean ionizing collision time is~10⁴ s. Electron loss, not useful ionization, can dominate replenishment:

`Prefill = D Ne Eeffective (1−recovery)/τloss`.

Eeffective includes the actual injector/collector circuit work and is scanned at100 eV,1 keV and10 keV. Merely assigning100 eV to each replacement electron does not close a kilovolt circuit. At1 keV and no loss-energy recovery, keeping this term below1 W/m² requiresτloss≥0.108 s at C=1 or≥0.262 s at C=.5. The calculation does not establish those lifetimes.

Shielding is central. A uniform spherical wall voltage creates no interior extraction field. An uncompensated volume electron population sized for this ionization rate would require a huge space-charge potential; a quasineutral volume plasma would instead alter charged-grain collection and focusing. Local electrodes with a grounded outer shield avoid relying on either arrangement. Mesh leaks, fringe fields and electron escape must be bounded for the source-packet charge tolerance.

## Geometry discriminator: broad shallow cells

A narrow, deep tube does not supply the planar layer's path distribution. [capillary.json](capillary.json) checks diffuse molecular access using the cylindrical Knudsen coefficient `DK=2r v̄/3`. For uniform ionization frequency ν, a reflecting rear boundary and a partially admitting entrance:

`Sint=√(DK ν) tanh[L√(ν/DK)]`, `Seff=(1/Sint+4/(C v̄))⁻¹`.

The net rate is `g=4 a η D Seff/v̄`. This homogenized diffusion model is a bounded geometry check; the1 cm case has aspect ratio10 and needs a later molecular-trajectory treatment.

The electron cloud's cold rigid rotation also contributes particle energy. Its slow branch obeys `ω(ωc−ω)=ne e²/(2ε₀me)`. Holding the illustrative edge rotation below40 eV prevents silently assigning a100 eV ionization distribution to a much faster rotating cloud.

For1 mm radius,10 cm depth and10% area, neutral depletion drives ne to6.49×10¹⁷ m⁻³ even with open admission. The rotation-compatible field reaches1.59 T; its ideal magnet diagnostic is1503 kg/m². With50% admission, it rises to4696 kg/m². These results reject the sparse, deep geometry as the favorable interpretation of the earlier planar example.

A **broad-area,1 cm-deep candidate**, with50% entrance admission,η=.5 and full live duty, instead needs ne=5.924×10¹⁶ m⁻³. At0.2 T, the self potential is268 V and edge rotation24.3 eV. Its electron inventory is5.924×10¹⁴/m², replenishment0.949 W/m² at1 keV/0.1 s, and ideal gap-magnet diagnostic23.9 kg/m². It suggests integrating the ion cells with shallow warm catcher channels. Opaque warm baffles, source interception and gas access must share a real layout.

A25% live-duty version needs a larger instantaneous population: the chosen field reaches0.741 T and ideal magnetic diagnostic328 kg/m². Duty cannot be omitted from the hardware sizing. Field-return volume can also dominate the active gap: placing magnets far behind the17.5 cm backing does not preserve the1 cm field-volume assumption. Magnet placement, thermal insulation and flux return are therefore the next mechanical closure, rather than a free24 kg/m² component allocation.

## Actual duty, placement and hydrogen

[placement.json](placement.json) unions the exact4,228 scheduled events over370 days. Closing from1 s before focus until10 s after the99.9%-mass main outgoing arrival leaves **79.7006% live time**. Incoming launches are allowed during cleanup, apart from their local source apertures. This is a main-pulse duty assumption; it does not dispose of the slower material. With this duty, the broad shallow model needs ne=7.5662×10¹⁶ m⁻³, B=0.203823 T, electron inventory7.5662×10¹⁴/m² and342 V self potential. Replenishment at1 keV per electron,0.1 s loss time and no recovery is0.966 W/m² averaged over the schedule.

The warm1 cm cells would be on the chamber-facing side of the150 µm collection skin. A physical grid/baffle assembly needs its own primary interception, heating and gas-return map; the old flat-skin response cannot simply be reused for surfaces placed in front of it. The factor C=.5 describes molecular admission, not50% opaque projected area. Retracting or shielding electrodes during the main pulse is a proposed hardware action. It does not turn off a magnet.

At the same B, an ideal permanent-magnet diagnostic based on0.2/0.3/0.5 m effective field-return volume is496/744/1240 kg/m², before actual yokes, fringing, support and thermal margins. Cooler magnets behind the17.5 cm wall require such a placement account. No1000 K permanent-magnet rating is assumed.

Using the primary [NIST H2 BEB parameters](https://physics.nist.gov/cgi-bin/Ionization/table.pl?ionization=H2) and [atomic-H parameters](https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?differential=Y&id=HI&initial=&ision=I), the100 eV ionization cross sections are9.6898×10⁻²¹ and5.4598×10⁻²¹ m² respectively. In the same cloud and with net transferη=.5, the access/ionization model gives gH2=.018224 and gH=.007968. The corresponding η needed for equal g=.01 is.2744 or.6275. For equal escaped hydrogen mass, the faster atomic-H thermal speed makes its steady mass column0.887 times the gas packet's H2/g=.01 reference; H2 gives0.549 times. These are conditional rate comparisons. Dissociation, reflection, cathode transfer and backstreaming remain separate.

The cleanup g acts on gas already in the cavity. Capturing freshly desorbed gas before it escapes changes the upstream residual fraction f; the same atom cannot also earn a second cleanup credit.

## Controlled-energy routed graphite cathode

[cathode.py](cathode.py) imports the current NIST stopping model and the existing finite-volume cap/channel solver without altering either archive. A300 kV negative collector, behind a grounded entrance, gives these mean-range results:

|Incident ion|Energy per stopping fragment|Mean projected depth|Conditional mobile side collection at Dparallel/Dnormal=1000|
|---|---:|---:|---:|
|H+|300 keV|2.533 µm|99.9997%|
|H2+|150 keV/H|1.205 µm|99.7830%|
|He+|300 keV|1.082 µm|99.5990%|

The domain has a100 nm cap,20 µm channel pitch and150 µm depth, with an absorbing side channel and cavity face and reflecting rear. The99% mean-depth thresholds are107 kV for H+,214 kV for H2+ and231 kV for He+. At anisotropy100 rather than1000,300 kV He collection falls to85.0%; the geometry/material requirement matters. Molecular fragments are assigned equal velocity after entry. Entry charge equilibration, backscatter, straggling, damaged-barrier permeability and trapped populations are not supplied by this mean-depth calculation. A broad depth distribution can increase front leakage relative to evaluating its mean alone.

Using the gas packet's routed residual loads,2415.36 kg/s He and3326.18 kg/s H, one300 kV pass costs0.00429 W/m² for He+ plus H2+, or0.00743 W/m² if H is atomic. Two attempts cost0.00858/0.01486 W/m². These are source-work requirements; recovered impact heat could enter the thermal ledger only once. They exclude escaped secondary electrons and electron-cloud replenishment. Each secondary electron that escapes a negative300 kV cathode adds300 keV to the circuit, so a recessed suppression electrode and a measured effective escaping yield are necessary.

The positive-ion path is electrostatic. A300 kV,0.3 m gap has1 MV/m mean field and158 ns cold-start He transit. The plane-diode Child–Langmuir scale is49.9 A/m² for He, while the residual global current is3.83 nA/m² before retries. Space-charge throughput is therefore generous, even with sparse cathodes. High-voltage insulation, local gas, field emission, focusing and secondary return require actual geometry. The accelerated He gyroradius at0.204 T is0.774 m; the electron-cloud field is not itself an ion guide through a centimetre cell.

Ionizer electrons remain in their own source/collector circuit behind grounded shielding; positive ions are extracted into the negative cathode. The laminate then supplies the same side-channel gas destination as the main catcher. This removes a mandatory new tungsten storage requirement. It does not remove the net optical transfer factor: to obtainηnet=.5 with the calculated99.599% mobile He collection requires at least50.20% ion delivery, before reflection, trapped populations or other losses. Temporary implantation remains distinct from transfer into the isolated gas circuit.

## A65 K coil-array requirement bound

[coil_array.py](coil_array.py) supplies an explicit return-field geometry. Alternating rails carry±I, spaced p/2, closing into finite racetrack circuits. For the ideal periodic straight sections at distance z from their plane:

`|B(x,z)| = μ0 I/p / sqrt[sin²(2πx/p)+sinh²(2πz/p)]`.

The minimum front field fixes `I=Bmin p cosh(2πd/p)/μ0`. Flux returns between adjacent rails; it occupies both sides of the current plane. For a uniform-current round cable of radius r, the stored energy per receiver area is

`U/A = μ0 I²/(2πp) [ln(p/(πr))+1/4]`.

The last term is internal field energy for that specified current distribution. A second agent independently checked the coefficient and field formula. Actual superconducting current distribution, finite ends, forces and protection remain separate. I is an aggregate ampere-turn requirement; multiple series turns can reduce supply current without eliminating cable volume or field energy.

Cable mass uses a required complete-cable operating Je of30 or100 A/mm² at65 K and the calculated local field. These are acceptance scenarios, not universally measured REBCO ratings. [Senatore et al.](https://arxiv.org/abs/1512.01930) measured substantial dependence on temperature, field, orientation and manufacturer. A full cable includes stabilizer, substrate, insulation and current margin.

The model places the entire coil, cold buffer and local graphite halo behind the17.5 cm backing plus5 cm clearance. Their radius changes the field distance and therefore the current. The measured [65–70 K polyethylene heat capacity](https://nvlpubs.nist.gov/nistpubs/jres/77A/jresv77An4p395_A1b.pdf), with a5% capacity reserve, supplies2257.9 J/kg pulse storage; effective density900 kg/m³ is assumed. The largest whole-pulse photon reserve is12.069 MJ/m² for B2, obtained by assigning all non-neutrino delayed heat to≤20 MeV photons with10% Doppler allowance, plus saved prompt radiation, behind the350 kg/m² main wall. The actual cycle-average cold-equipment heat is7.494 W/m². The halo attenuation uses the minimum [NIST carbon mass energy-absorption coefficient](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html) through20 MeV,0.001198 m²/kg, with local-secondary and joint-streaming qualifications.

For a100 m pitch, Je30 A/mm² and2000 kg/m² radial halo, the coupled solution requires16.41 MA-turn per rail,0.417 m cable radius,7.87 T peak cable field, and a coil center2.448 m behind the hot face. Cable including simple end returns is88.2 kg/m², local graphite433.1 kg/m² and PE buffer68.0 kg/m²: **589.4 kg/m² subtotal**, before cryostats, supports, thermal spreaders, electrical feeds or protection. This is26.8% of the entire existing plant mass. A20 kA winding would need821 turns per rail in this reference.

Only the halo's actual intercepted solid angle contributes to the cold dose; no full cold plane is exposed to7.494 W/m². In that case mean photon load at65 K is0.0954 W/m². At15% of Carnot COP for65→200 K, photon refrigeration costs1.32 W/m². A200 K shield with effective emissivity.01 adds1.91 W/m² electrical refrigeration; an assumed effective1000→200 K emissivity10⁻⁴ gives0.870 W/m² hot bypass. Pulse thermal contact, polymer irradiation, joints and current-sharing margin remain unqualified. This conservative annular whole-pulse-buffer layout becomes much heavier at10–50 m pitch. Shared shielding, remote storage and time-resolved refrigeration could alter it; the calculation does not exclude those other layouts.

## Field interaction and switching: the broad array is not a drop-in receiver

At0.203823 T the magnetic pressure is16.53 kPa, while sampled source ram-pressure peaks range0.000492–0.115 Pa. Quasineutrality does not remove magnetic stresses. The charged component of the outgoing plume is strongly deflected by the broad field; switching off only the electron cloud leaves this effect intact. The100 m array's field decays on a15.9 m scale into the cavity. Its transverse field integral on the ray above a rail is8.43 T m, compared with0.144 T m rigidity for a1 MeV proton.

[switching.json](switching.json) therefore preserves this as a **switching or field-localization requirement**, not a recommended continuously energized design. The100 m case stores2.470 MJ/m². Cycling once per scheduled pulse circulates653 W/m² of combined charge/discharge energy. A loss of0.1% of stored energy per complete cycle costs0.327 W/m²;1% costs3.27 W/m². If that loss occurs at65 K, refrigerator work adds13.85 times the dissipated cold heat. Recovery-electronics loss, conductor hysteresis/AC loss and cold loss location require separate qualification.

A100 s ramp needs at least24.7 kW/m² converter rating and changes the main-mass live fraction to77.06%. Keeping the original cloud gives g=.009668; restoring g=.01 requires ne=7.8497×10¹⁶ m⁻³ and Bmin=.21066 T in the same rotation screen, followed by a revised coil/thermal sizing. Thus the base-field ramp numbers are requirement estimates, not a closed g=.01 operating point. At1 kW/kg they add24.7 kg/m² converter mass; an independent100 kJ/kg magnetic-energy buffer adds another24.7 kg/m². The existing driver store is not credited without a joint schedule. Field cancellation also matters: a chosen0.01 rad kick limit for the1 MeV proton reference requires the residual broad field to fall below1.71×10⁻⁴ of its on amplitude along that ray. This is an illustrative transport criterion, not a universal tolerance. Eddy currents in conducting backing, shields and supports must follow their actual segmented geometry.

Waiting until the archived `full_resolved` timing benchmark leaves only0.012% live time over the finite schedule, essentially the startup interval. This is a finite-grid velocity benchmark; interpolation toward the zero-velocity center can contain still later, smaller parcels. B1's26.1 ks benchmark already overlaps its7.70 ks cadence. Reopening after99.9% main mass leaves approximately0.125% of B1 He nuclei and0.0829% of its H nuclei in later packets. Those particles remain in the closed station; they are not automatically lost. Their magnetically redirected collection, gas release and heat path nevertheless need their own calculation. A field-guided slow-tail collector could be constructive, while retaining the ballistic reference for the main pulse only.

A compact local magnetic return path is the other live alternative. For the illustrative0.01 rad/1 MeV proton criterion, the transverse field integral must stay below1.44×10⁻³ T m—about7.2 mm at0.2 T. Lower-energy ions impose tighter integrals. A normal field aligned with the incident trajectory and explicit collection optics can use a different criterion. Cooler magnets or pole structures must then be close enough to deliver the centimetre-scale field, with their insulation, primary exposure, gamma heat and temperature-dependent magnetic properties accounted. Placing an unshielded cold plane at the front, or silently using a1000 K superconductor, does not supply that geometry. No completed local return-circuit design is asserted here.

## Retention and regeneration at1000 K

There is a material mechanism worth testing. [Nguyen-Manh and Dudarev](https://arxiv.org/abs/1408.0630) calculate4.63 eV binding and~4.70 eV dissociation for one He atom in a tungsten vacancy, comparing with a4.6 eV experimental value. Further He added to an already occupied vacancy has substantially weaker binding. Their calculation supports a specific trap mechanism, not the fraction of ions captured by a hot coating.

For a first-order trap, `S=exp[−νt exp(−E/kT)]`. At1000 K andν=10¹²–10¹³ s⁻¹,95% survival for1000 s requires E=3.23–3.43 eV; for one day it requires3.62–3.81 eV. Under that isolated-barrier model,4.6–4.7 eV traps retain He strongly at1000 K, while heating to1800 K yields95% release in2.27–4.32 s. The3.1 eV comparison has a420 s mean residence at1000 K forν=10¹³ s⁻¹. These are conditional rate calculations, not a measured hot-cathode release curve.

Empty hot vacancies can migrate or anneal; a reservoir of stable, vacant traps is not automatically available. Low-energy He can escape before meeting one. Implanting at10 keV permits up to~835 eV transfer to a W atom in a single collision, so creation of local damage/traps is kinematically possible, with capture yield and defect evolution still unresolved. Surface reflection, near-surface escape, repeated occupancy, bubbles, sputtering and co-implanted hydrogen need measurement or kinetic/material simulation.

A100 nm W coating over10% of the wall carries0.000193 kg/m². This collector area is an independent focusing/geometry assumption; coating the full active face multiplies its mass and whole-coating processing figures by ten. If its *effective usable* single-He capacity were10⁻⁴ per W atom, one day's residual flow occupies3.73% of that capacity. This is a stated capacity requirement. Heating just that coating through800 K at150 J/(kg K) once per day costs0.000268 W/m² before support and heat recovery. Whole-coating processing at2–100 MJ/kg would instead cost0.00447–0.223 W/m² per day; at1000 s intervals it rises to0.386–19.3 W/m². Selective extraction need not process every atom in the coating. No blanket atomization cost is imposed.

Regeneration must occur after isolating the patch from the chamber; otherwise it simply reverses pumping. A trap barrier is not itself the electrical work of selective atom handling. Capacity, residence and restoration work remain separate inputs.

## What this checkpoint resolves

The residual throughput and acceleration energy do not prohibit warm cleanup. A controlled-energy routed graphite cathode provides a concrete inventory destination under the same laminate qualification as the primary catcher. The unresolved coupled requirements are electron confinement, protected molecular access, field localization or recovered-energy switching, slow-tail collection and net transfer. The broad static magnet model cannot be inserted unchanged into the ballistic plant. Its mass and cooling requirements are material compared with2196 kg/m² plant mass and45.64 W/m² export; the radiator is already fully allocated, so additional heat needs explicit reallocation or derating.

A plausible next acceptance test would demonstrate measured C, Pion, η and electron loss power together at1000 K, through a protected entrance, while closing the helium inventory over regeneration. This is more informative than demonstrating ionization alone. Current evidence cannot support a numerical80% credible interval for the complete future pump efficiency; the parameter grid is an explicit engineering sensitivity family.

Run `screen.py`, `capillary.py`, `placement.py`, `cathode.py`, `coil_array.py` and `switching.py` with the archived scipy runtime. The gas agent independently checked the capture normalization, Maxwellian integral, capillary admittance and periodic-array field/energy formulas. The new files preserve their exact source hashes. No plant, article or earlier catcher record is modified.
