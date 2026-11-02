# Opaque-pulse engineering reassessment through 2049 ✨

29 October 2026. Independent bounded reassessment of the opaque catalytic route. Existing pulse trajectories and source ledgers were read; the only new calculations were algebraic energy, compression and recovery checks. No new hydrodynamic simulation was run. The previous selected Q=2.6 history is not evidence for this assessment.

**Assessment:** a superintelligent civilization with the established 2040s industrial capabilities could plausibly develop useful catalytic pulses, including an A-stage capture gain approaching the favorable model's order of magnitude. No thermodynamic or radiation-loss prohibition has emerged. Preparing the converging flow is a serious coupled engineering problem, and exact Q=8.18 is not established. The strongest justified history is therefore allowed to contain technical success, extensive optimization and limited deployment. Persistent low gain should be attributed to a demonstrated surviving mechanism if selected, rather than to our unfinished cold-fluid solver.

## What the favorable result establishes

The saved 64-cell small A pulse contains 1.131e13 kg, receives 3.493e25 J of kinetic plus internal energy, and gives prompt capture Q=8.182 over 68.4 microseconds. It converts 53.84% of the initial N14, peaks at 1.149 GK, and reaches sampled peak density 1.879e8 kg/m³. It starts from a prescribed favorable hot inflow and uses a reduced three-species network. The full material's minimum sampled outer radius is 71.4 m; its 30 m parameter is a reference scale. These are a constructive spherical burning example, subject to the documented numerical qualifications, rather than a driver demonstration or a theoretical optimum. [Saved calculation and qualifications](../phase-two/driver/COLD-COMPRESSION.md); [independent transport extraction](../phase-two/pulse-transport/small-pulse-check.json).

The nuclear energy budget has substantial room for capture gain. At the same delivered specific energy, complete conversion would give prompt Q about 15.2. Q around 8 therefore does not require complete burning. It does require substantially more successful compression than an arbitrary head-on collision. Target shape, radial mass flux, velocity schedule, entropy profile, ignition timing, hot fraction and pulse scale remain optimization variables. Localized ignition or a propagating burn might improve the energy investment further, but no such improvement is credited here.

## Thermodynamics of preparing the flow

The earlier cold probe did not solve a physically faithful inverse preparation problem. It backmapped positions at constant velocity, assigned a uniform cold temperature, and retained the favorable hot flow's kinetic energy. That preserves a ballistic schedule while omitting the pressure work needed to arrive at the same hot state.

For A, the favorable total energy is **3.089 TJ/kg**, of which the initial convergent kinetic energy is about **2.872 TJ/kg**. The cold probe therefore starts with approximately 7% less total energy. Restoring the favorable total energy as incoming kinetic energy would increase speed by roughly **3.7%**, before solving for a proper density and velocity history. This is an accounting correction, not a claim that that speed adjustment produces the favorable trajectory.

An independent read of the existing electron/ion/radiation EOS at 1e8 kg/m³ gives:

| State | Specific internal energy | Pressure |
|---|---:|---:|
| 100,001 K, fully ionized degenerate comparison | 1.110 TJ/kg | 7.22e19 Pa |
| 300 MK, original A reference state | 3.089 TJ/kg | 1.93e20 Pa |

The low-temperature entry is a dense-plasma energy scale, not an EOS for cold pellets. The available incoming energy is sufficient to pay for compression and a hot state; the route must avoid turning too much energy into entropy before sufficient density is reached. Starting with genuinely cold material does not impose a minimum entropy above the desired hot state. Work can produce compression, and controlled dissipation can then establish the chosen entropy. Ionization, Coulomb terms and phase changes alter that path and must be treated consistently.

The prior contact audit puts even an isolated-atom upper estimate of complete ionization below 0.5% of whole-pulse input, while finding that it can exceed the ideal contact internal energy severalfold. This means the original cold adiabat is unreliable. It does not establish a large unavoidable loss of total driver energy. A physical chemical/free-energy EOS, verified shock and ramp solutions, and convergence studies would distinguish the favorable path from excessive early heating. [Microphysics audit](../phase-two/pulse-transport/README.md).

Laboratory evidence supports controlling compression entropy as an engineering variable. Smith et al. measured ramp compression of diamond to 5 TPa and discussed the need to shape pressure loading to prevent shock formation. Their regime is far below the pulse's final pressure and is not a demonstration of this reactor, but it directly refutes treating all rapid compression as one fully dissipative shock. [Smith et al., Nature 2014](https://duffy.princeton.edu/sites/g/files/toruqf616/files/smith_nature_14.pdf).

## What an engineered assembly could actually do

The useful candidate is a **programmed converging material flow**, with a controlled transition from separate cold projectiles or streams to an ionized fluid. Three distinctions matter.

First, the approximately 2400 km/s bulk speed is not automatically thermalized when neighboring grains merge. In a smooth radial field, neighbors have similar directions and velocities. Dissipative heating depends on their relative velocity. For streams meeting at half-angle theta, the available merger heat is `v² sin²(theta)/2` per unit mass. Segmentation, correlated launch positions and velocities, and radial timing can reduce that quantity. Grain number does not force the merging angle to equal the large angular separation of a few launcher stations; a station can supply a structured family of trajectories.

A scale example illustrates the distinction without specifying a complete apparatus. For a local velocity gradient around 2,000/s, neighboring material separated by 2 cm differs in velocity by about 40 m/s, releasing about 200 J/kg if two equal-mass parcels mix inelastically. That is vastly below the TJ/kg bulk energy. A poorly phased collection of intersecting streams could instead dissipate a substantial fraction of the bulk energy. The continuum limit and the bad-merger limit are physically different initial conditions.

Second, “engineered preheat” should name its energy-delivery mechanism. Uniformly illuminating all ordinary grains after the cloud becomes optically thick does not work automatically; heating them early enough to remain visible creates a vapor-expansion problem. That rejects one simple external-heating recipe. It does not reject programmed weak shocks, staged compression waves, or other distributed mechanical heating after the material becomes contiguous. Specific heating of order 2 GJ/kg corresponds to a dissipative velocity scale around 63 km/s relative to the common flow, only a few percent of bulk speed; two equal streams would differ by twice that velocity. The appropriate value must come from a physical EOS and an optimized trajectory; arbitrary uniform heating is not interchangeable with entropy control.

Third, the fluid must form at a planned density and time. A cold grain's geometric swept-mass column is of order `4 rho_grain a/3`; millimetre-to-centimetre grains have a macroscopic stopping and fragmentation problem distinct from the stopping of individual ions. Relative collisions can vaporize and ionize material while preserving most of its common inward momentum. Ion beams require their own stopping, charge-neutralization, emittance and beam–plasma calculations. An assumption that particles penetrate a dense target without heating it would be inappropriate. Conversely, energy deposition at controlled interfaces is useful when it creates the intended compression wave.

Neutralized beam compression and focusing are real control mechanisms. NDCX-II produced roughly millimetre-radius, nanosecond ion pulses using a velocity ramp, plasma neutralization and final focusing. Its particle number and geometry are immensely smaller than this application, so this is precedent for the mechanism, not for its throughput. [Seidl et al., 2016](https://w3.pppl.gov/~ikaganov/pdf_files/Short_Pulse_Compressed.pdf).

A credible successful history can therefore specify optimized radial phase space, many locally smooth streams, charge compensation before final convergence, and a measured controlled ionization/merger stage. Source calibration and launch-schedule optimization are plausible routine work for the established civilization once the material and plasma models are qualified. Qualification of the material transition, collective instabilities and irradiation-dependent behavior remains substantive physical development. Superintelligence should accelerate that development; it should not be invoked as a missing force or an uncounted energy source.

## Radiation, stopping and instability

The smaller trajectory is strongly opaque where it burns. Its energy-release-weighted median Compton opacity is 0.00720 m²/kg. The 1st-percentile effective cell thermalization depth is 17,007. A sensitivity reducing both true absorption and scattering by 100 still leaves that depth at 170 and the 99th-percentile absorption-time/local-evolution-time ratio at 2.95e−4. The existing qualified free-free and Compton treatments require refinement, but these margins do not suggest bremsstrahlung immediately escapes from the burning bulk. The core supplies its own opacity independently of any external collection shell. [Independent extraction](../phase-two/pulse-transport/small-pulse-check.json).

Radiative precursors and shocks can still alter the incoming flow. Balancing diffusion with advection gives a characteristic precursor length `c/(3 kappa rho v)` where diffusion applies. This is a local scale, not a justification for ignoring a shock. At the hot dense interior it is small compared with the target. At dilute interfaces, the real opacity, absorption and geometry must determine whether premature heating is important. The bulk opacity argument should not be promoted into a claim that every driver interface is benign.

Conduction is small on the saved resolved hot profiles, while the source calculation does not resolve microscopic shock widths. Thin-interface transport, radiation-mediated entropy production and mixing belong in a refined local model. These are reasons to compute and optimize the actual transition, not reasons to assume the entire pulse loses its energy.

Spherical convergence amplifies errors and makes low-mode asymmetry, angular momentum, shock corrugation and mixing real risks. An inward flow whose layers reach the focus over a long interval can stagnate and launch an outward shock before most of the fuel compresses; temporal focusing is essential. A smooth filled flow also differs from a thin heavy ablator surrounding a light hotspot, so instability limits from one ICF geometry cannot simply be copied to this one.

The experimental literature supports both risk and control. Yates et al. measured strong nonuniformity when merging plasma jets shocked, and improved symmetry after balancing the jets. MacPhee et al. showed that a shaped four-shock drive could suppress measured instability growth while retaining a low fuel adiabat. Neither establishes stability at this pulse's conditions; together they support treating merger geometry and shock sequencing as design variables. [Yates et al., 2020](https://arxiv.org/abs/2002.03006); [MacPhee et al., 2015](https://doi.org/10.1063/1.4928909).

The cold-probe record supplies no physical instability ceiling. Its initial artificial-viscosity stress exceeds physical pressure by factors of thousands, its gains change strongly with resolution and viscosity, one intermediate-resolution solve fails, and higher-resolution burns stop at a temperature guard. Good total-energy conservation does not cure uncontrolled entropy generation. Neither the low coarse-grid gain nor a guard-stopped high gain should decide the history. [Cold-probe record](../phase-two/driver/COLD-COMPRESSION.md).

## Q near 8, electricity and the remaining stages

For prompt capture gain Q, collection fraction f, electrical-to-delivered efficiency eta_d and conversion efficiency eta_h, the immediate comparison is

`E_net/E_delivered = f eta_h (1 + Q) - 1/eta_d - auxiliary loads/E_delivered`.

The incoming energy is recovered once. If the source loses energy outside the target, that loss is already represented by eta_d. If preparation requires additional delivered energy, the denominator and outgoing energy must both be updated.

At Q=8.182 and full collection, the immediate-only comparison gives:

| Driver efficiency | Conversion efficiency | Break-even Q before auxiliaries | Net electricity / delivered input |
|---:|---:|---:|---:|
| 30% | 40% | 7.333 | 0.339 |
| 38% | 44% | 4.981 | 1.408 |
| 50% | 50% | 3.000 | 2.591 |

The first row is a narrow margin: about 10.2% of source electricity, or room for roughly 10.4% reduction in prompt nuclear yield at unchanged input. The 38%/44% pair appears in a published heavy-ion fusion point design; it was a design calculation, not a demonstrated performance transferable to macroscopic grain acceleration. The 50%/50% row is a future sensitivity, not a selected qualified specification. Efficiencies are therefore consequential engineering choices, not immutable physical constants. [HIF point-design paper, indexed primary text](https://www.osti.gov/servlets/purl/815523).

Delayed energy materially changes the comparison. The existing O15 ledger assigns about **1.754 MeV of available decay heat per capture**, after its 1 MeV neutrino allowance. Applied to the small A conversion, that is about **1.967 delivered-input units** of additional heat. If collected, total available nuclear Q becomes approximately **10.15**, and the first row's net increases from 0.339 to **1.126 delivered-input units**, before auxiliaries. This is an existing decay-ledger estimate, not a new beta-spectrum calculation. Delayed positrons and annihilation or nuclear photons cannot all be assigned a high bulk-kinetic conversion efficiency. [Cycle ledger](../cycle-ledger.json).

Collection must be integrated over those components. At the immediate-only Q=8.182, the minimum collection fractions for the three rows are approximately 90.8%, 65.1% and 43.6%. Small diagnostic panels cannot establish net range electricity. A thermally absorbing enclosure, moving or replenished absorbing material, or a suitable charged-product recovery architecture could have different mass and conversion. The old 10,000 kg/m² absorber and multi-Earth-mass allocation are design assumptions rather than a lower bound. Hard decay-photon absorption, impulse, heat distribution, catalyst return and maintenance must remain in any lighter replacement.

A captures a proton on N14; its O15 decay produces N15. B1 and carbon stages must return the catalyst and produce helium before this is a closed ordinary-hydrogen cycle. The existing ideal-recovery cycle has about 24.68 MeV of available nuclear heat and 1.777 MeV of prepared input per net helium. Replacing only A's original conversion with 53.84%, while holding all other stage properties fixed, increases prepared input by approximately 0.301 MeV to **2.078 MeV/He**. This is a useful accounting sensitivity, not a newly solved small-scale cycle. All scaled B1/B2/C13 yields still require their own calculations; shorter hydrodynamic residence times cannot be ignored.

Stage geometries need not share a radius, density, temperature or pulse cadence. The helium-to-carbon factory is a makeup route whose load depends strongly on catalyst loss; existing carbon and nitrogen inventory can supply early experiments. A completed plant must account for eventual makeup, separation and processing, including cold beta-decay heat. Thus successful A gain can reasonably lead to a useful catalytic sequence and then an integrated plant, while each milestone has its own measured account.

## Implication for a history through 2049

The established industrial setting can afford parallel material tests, local merger experiments, independent radiation-hydrodynamic codes, pulse-shape searches and repeated large shots. It should have substantially stronger evidence than this bounded investigation. The canon's absence of practical protium fusion in the early 2040s remains compatible with successful late-2040s capture stages and a limited integrated installation. [Civilizational capability audit](capabilities.md).

It is reasonable to select a successful order-8 capture-stage result after the history names a controlled assembly method and the measured output/input convention. It is equally necessary to distinguish that selection from a prediction that our particular initial profile will survive exactly. A claim that gain remains permanently near 2–3 requires additional physical or institutional reasons; our inability to finish a solver supplies neither. A claim of cheap mass deployment likewise requires a favorable complete plant and alternatives comparison, not just target gain.

The engineering sequence that best fits the present evidence is: qualify cold-material EOS and controlled contact; establish reproducible converging flows and stage gains; complete the catalytic sequence and delayed-energy account; integrate high material return with adequate energy conversion; then decide whether the installation's specific output, replacement burden and service value justify replication. Those activities can overlap extensively. A limited successful installation with weak commercial demand is a coherent outcome. A smaller open range remains useful for development, but it need not be the maximum achievement of the entire civilization in 2049.

No new hard obstruction to the opaque route was found in this reassessment. The central unresolved question is whether a robust, affordable driver/collector combination can realize enough of a favorable pulse family. That is a physical engineering question with several plausible control mechanisms, rather than a date inferable from one unsuccessful numerical probe.

## Calculation and access notes

New arithmetic used the existing FastEOS table and saved A/cycle JSONs. The dense-state energy check used rho=1e8 kg/m³, Ye=8/15 and Yi=2/15 at 100,001 K and 300 MK. The late-heat estimate scales the existing N14 stage's cold available heat by its conversion. The A-only cycle sensitivity is `0.4801715923 × (1/0.5384000787 − 1.2311411282) = 0.3006900422 MeV/He`. These are bounded checks; no new fluid solution or complete plant optimization is claimed.

Primary papers were read for ramp compression, shaped-drive instability control, merging jets and neutralized ion compression. The HIF point-design table was available in indexed primary text; direct OSTI retrieval returned 502, and the institutional landing page presented a robot check. Neither access problem was bypassed. Its efficiency pair is used solely as an explicitly labeled design assumption. All other research files remain unchanged.

## Subsequent bounded small-stage check

The follow-up [small-stage packet](small-stages/README.md) completes full-network A, B1, C13 and a 3 m helium factory. B2 reaches gain 10.52 and conversion 0.8193 before an identified low-temperature network/outer-EOS failure; the stored endpoint is close to capture freeze-out but remains incomplete. Two explicitly conditional cycle ledgers quantify the 30 m versus 3 m factory tradeoff. These results supersede this memo's earlier statement that no small B1/B2/C13 diagnostics had been attempted. They preserve the distinction between numerical closure, a plausible engineered route and the selected future history.
