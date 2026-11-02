# Collection apparatus: first component checkpoint ✨

The strongest current baseline is **material interception coupled to the existing 60% thermal plant**, with a chemically separable, recoverable catcher surface. The source calculation does not require removal of 99.6% of kinetic energy before impact. Staggered annual or multi-month extraction is now the default process comparison. An electromagnetic braker remains an optional response to a demonstrated loss or contamination problem, and an optional route to higher electrical efficiency.

This checkpoint calculates stopping-depth distributions and nuclear-loss acceptance tests for the actual four-stage, finite-flight isotope mixture. It has not yet closed isotope extraction, all N/O nuclear channels, or the local vapor return circuit. The proposed hardware is therefore a conditional design with explicit acceptance tests.

## Important update: batch extraction removes a premature throughput requirement

The pulse interval is **not an extraction deadline**. Fuel may remain implanted while a larger circulating inventory supplies subsequent shots. At the reference cadence, the mean gross fuel loading is only **9.3×10⁻⁵ kg/m²/year**. Processing a 10 µm graphite skin annually requires about9.64×10⁹kg/s across the whole receiver and0.00127W/m² at the conditional2MJ/kg processing cost. A 50 µm alumina skin annually costs0.0124W/m² on that same assumption. A synchronized full year of gross fuel stock is1.41×10¹⁵kg; a uniformly staggered collection age distribution holds about half that amount on average. This is a small fraction of the existing station mass.

Consequently, the expensive whole-skin-per-pulse case is a sensitivity, **not a requirement of direct capture**. Annual or multi-month graphite recovery remains a serious baseline, with the chemical and isotope plant processing accumulated material in staggered cassettes. Lifetime, isotope leakage, secondary reactions in the accumulated fuel and actual release chemistry remain acceptance tests. The higher throughput of a new coating should be justified by better separation or vapor behavior before replacing graphite.

## What the stopping calculation changes

`collection.py` reuses the archived carbon stopping model and extends it with measured IAEA electronic-stopping rows for aluminium, tin and hydrogen. It retains individual isotope mass numbers and the 512 velocity packets per pulse. It integrates electronic plus screened elastic stopping in the continuous-slowing approximation; elastic nuclear stopping does not itself destroy a nucleus. New aluminium coverage contains 135 carbon, 108 nitrogen and 154 oxygen measurements, covering the main source energies. Solid-state stopping is used to estimate molten-metal areal ranges; temperature, charge equilibration, scattering and straggling are further qualification requirements. [IAEA stopping database](https://www-nds.iaea.org/stopping/).

CNO number-weighted median / 99.9th-percentile path lengths are:

| Stage | Graphite, 2,000 kg/m³ | Aluminium, 2,300 kg/m³ |
|---|---:|---:|
| C13 | 3.95 / 21.88 µm | 5.10 / 24.10 µm |
| A | 2.99 / 9.08 µm | 3.89 / 10.86 µm |
| B1 | 2.47 / 4.86 µm | 3.16 / 6.23 µm |
| B2 | 3.51 / 17.70 µm | 4.55 / 19.64 µm |

These are CSDA path lengths rather than simulated projected implantation profiles or qualified retention probabilities. The fastest sampled proton stops in about 46 µm of aluminium in the analytic proton model; the fastest sampled CNO stops within 25 µm. A 100 µm aluminium layer (0.23kg/m²;3.50×10¹⁸kg over the receiver) supplies a useful starting thickness without requiring the whole layer to be replaced each pulse. Detailed range straggling and reflection determine its final thickness. The tin measurements cover much less of the relevant energy range, so its screening result has lower evidential weight.

## Nuclear events and useful-isotope conversion along the slowing path

For target isotope density per areal mass `n_t/ρ`, the relevant probability is

`P_j(E₀)=1−exp[−∫₀ᴱ⁰ σ_j(E) (n_t/ρ) dE/S(E)]`.

Each reaction then needs an isotope/product map. A capture or transfer into a usable CNO isotope can change processing and cycle routing without permanently losing catalyst. A neutron or proton emitted from the target also needs transport; no secondary cascade is included here.

A useful acceptance diagnostic is to replace the unknown **nonelastic-event** cross section by a constant solely to determine the required scale. Across the actual cycle, a 1 barn constant gives CNO reaction probabilities of32.3ppm in graphite,21.1ppm in aluminium and13.9ppm in tin. For aluminium a 3 barn test gives about63ppm, compared with373 ppm total cold-stage nonreturn. These are neither measured destruction rates nor assertions that the relevant cross section is everywhere bounded by those values. They show that a small stopping column can make even appreciable nuclear cross sections tolerable. The historical373 ppm return parameter should not be assigned to every processing inefficiency. In a closed system, unreleased material and retained tails increase inventory and residence time. They require makeup only if a terminal sink is explicitly introduced.

The graphite structural-carbon allowance is tighter than the catalyst allowance: approximately0.0004704 carbon atoms per completed helium are available from the0.0014 structural C:H reserve. All-projectile stopping exposure corresponds to0.0003147 target events per helium per constant 1 barn, or a1.495barn stopping-weighted destructive average if the entire material allowance were assigned to this mechanism. Most of that diagnostic exposure comes from proton and alpha trajectories, whose actual reaction channels and cross sections differ greatly from heavy-ion fusion.

For one genuine energy-dependent nuclear calculation, the measured C12+C13 fusion fit of Notani and colleagues has `Rb=7.4fm`, `Vb=5.7MeV`, `ℏω=2.1MeV` over4.9–8.0MeV centre-of-mass energy. Applied to the actual carbon-projectile subset, it gives about**9.28×10⁻⁷ fusion events per helium**, around0.2% of the full structural-carbon allowance. Below the fit range the script uses an explicitly matched Gamow tail; above it, the Wong continuation. This is an illustrative carbon-isotope model, not a complete nuclear bound. The low-energy continuation is not a measured fit and does not reproduce every resonance. It is not applied to nitrogen or oxygen. [Notani et al., PRC85,014607](https://doi.org/10.1103/PhysRevC.85.014607).

Relevant CNO-on-carbon fusion measurements exist, but their coverage is uneven relative to this source. Kovar and colleagues measured C12+N14 only from15.1MeV centre-of-mass energy, while much of our mixture slows below that interval. Filling this gap with a constant barn-scale cross section would erase the Coulomb-barrier dependence and would not justify mandatory electromagnetic braking. [Kovar et al., PRC20,1305](https://doi.org/10.1103/PhysRevC.20.1305). Published C12+Al27 residual-production data beginning around39MeV projectile energy likewise lie above the maximum carbon energies here; those data identify channels but cannot alone certify the receiver. [Sharma et al., PRC91,024608](https://www-nds.iaea.org/nrdc/india/ws2017/article/pr_c_91_024608_2015.pdf).

## Concrete collection route and present decision

The material route comprises locally supported receiving cassettes thermally attached to the existing graphite absorber; a reusable collection layer; a closed process manifold; cold traps and gas return; and chemical/isotope processing with retained product and tails. Its primary heat continues to enter the1,000 K receiver and60% engine. It transmits the archived momentum impulse to the shell. It requires no megavolt ion retarder and works with both neutral and charged incoming atoms.

Graphite is the simplest interception material, but implanted carbon-13 becomes diluted in target carbon. Skimming or vaporizing the whole stopping layer after every pulse can dominate processing throughput. This is a process issue, separate from the small carbon-fusion estimate. The frequent B1 pulses have a shallower range and should not automatically receive the same renewal schedule as the rare C13 pulse.

An abundant aluminium layer avoids target C/N isotope dilution. Incoming carbon and nitrogen can in principle be chemically tagged as carbide/nitride and later released into separate process streams; the material-loop study is evaluating actual reaction conditions and regeneration throughput. A 100 µm molten layer is a stopping-thickness candidate, not a requirement to pump or chemically regenerate0.23kg/m² after each pulse.

**Bare molten aluminium at1,000 K has a serious newly identified interface problem.** The independent contamination screen obtains an equilibrium vapor pressure near5.1×10⁻⁶Pa from primary thermochemical data. In a gigametre chamber, even that pressure can impose a large column on incoming cold packets. Global material conservation does not keep that column low. The molten route needs local vapor capture, shutters or a compatible lower-activity surface; its vapor topology is not yet closed.

A refractory alumina collection skin is consequently a priority alternative: it avoids carbon and nitrogen matrix dilution and can have much lower vapor release at1,000 K. Most radioactive oxygen-15 in the fuel becomes nitrogen-15, allowing chemical separation from stable target oxygen. The added measured alumina stopping calculation gives CNO median/99.9th-percentile depths of2.69/12.62µm for C13,2.16/5.29µm for A,1.77/3.23µm for B1 and2.47/10.30µm for B2 at3,900kg/m³. Its constant 1 barn CNO-event sensitivity is25.3ppm, and the3barn case about76ppm. Nuclear-channel and extraction closure remain separate. A 50 µm coating is a range-margin candidate. The independent radiation-interface calculation finds modest integrated equilibrium oxide evaporation using the separated prompt and slow temperature histories, rather than incorrectly applying their summed maxima for a full pulse; see[radiation-interface.md](radiation-interface.md). Tin remains a secondary comparison because of its larger barrier and possible carbon phase separation, but availability, prompt absorption, wetting and sparse stopping coverage must be counted.

The retention criterion is the complete material balance: recoverable nuclei in the coating, condensed vapor, trapped gas, secondary surfaces, process product and retained tails, minus permanently transformed or externally lost atoms. A rebound into a closed receiver is not automatically permanent loss; it can still be unacceptable contamination for the next launch. Conversely, an implanted atom is not yet a usable returned feedstock. Both boundaries must be demonstrated.

The5,000kg/m² backing absorber is approximately61,700times the largest graphite CSDA range among the retained source species (0.0811kg/m², a trace alpha component); even a factor-two stopping uncertainty leaves a huge margin. Passage through a thin collection coating therefore chiefly changes the location of retained atoms. It does not imply escape from the station. Fast neutrons and photons require separate transport, and externally open launch ports or actual permeation/ejection paths are the relevant physical escape mechanisms. None is replaced here by an assumed generic loss fraction. Conversion into nuclei outside the four-stage fuel set creates a separate ash or recycling stream rather than disappearance of matter. Deliberate disposal would be an operating choice, not an intrinsic consequence of a finite extraction yield.

## Charge and electron circuit

A conductive material catcher accepts the quasineutral plume directly. Incoming electrons, conductor electrons and local return paths neutralize deposited ions; secondary electrons and a changing surface sheath are handled within the cassette. The fully stripped peak positive-current scale is below0.6A/m² in the archived source. No global megavolt field is needed to hold the entire pulse charge apart. An adequate conducting path is a local surface-current requirement rather than a central conductor carrying the sum of every cassette's current.

For scale, assuming conductivity3×10⁶S/m and100µm thickness, a metal coating has sheet resistance0.0033Ω. Even a deliberately pessimistic0.6A/m² net current imbalance across a1m patch gives millivolt-scale drops. Actual sheath potential, secondary yields and collection asymmetry depend on electron temperature and surface state. The conductive-path example does not predict these quantities or add electron/ionization energy a second time to the source heat.

An electrostatic direct converter would instead need electron separation, stage-dependent `E/q` acceptance and local high-voltage loads, as developed in the cold-sink direct-recovery packet. It would also require charge preparation for neutral material. An inductive or MHD braker preserves quasineutrality but needs a conducting, coupled plasma and an external electrical load: static magnetic fields alone do no work. At the favorable fully stripped C13 peak, the prior calculation gives roughly1km ion inertial length and650m nitrogen gyroradius at a0.49mT ram-pressure field. These are apparatus-scale constraints, not proof of an efficient compact braker. Material capture presently avoids that entire preparation requirement.

## Priority and uncertainty

The following judgments guide the next work; they are subjective80% ranges conditioned on the saved source and advanced industrial fabrication, not measurements or calibrated forecasting intervals:

- Probability that implanted-fuel extraction forces a substantial change from a simple graphite wall: central35%, range15–65%, after recognizing annual batch extraction. The accepted criterion is full cycle return with tolerable inventory, regeneration energy and chemical yield. The earlier per-pulse throughput objection alone does not meet this criterion.
- Probability that an uncovered molten-Al surface needs a major vapor-control redesign: central90%, range75–99%. The criterion is preserving incoming packet entropy, not merely retaining all atoms somewhere in the chamber; the vapor-pressure calculation makes this a high priority.
- Probability that nuclear losses alone force near-complete electromagnetic braking after noncarbon catcher options are examined: central10%, range2–30%. This remains uncertain until N/O channels, tails and secondaries are integrated. The stopping columns and existing carbon subset favor testing simpler capture first.

No credible interval is attached to the exact geometry or arithmetic. A measured thickness/release efficiency is not invented. The next decisive work is the joint catcher-surface, chemical release and vapor-capture model, followed by energy-dependent N/O nuclear closure. Those tasks can change the design substantially; a large PIC calculation is not yet the most useful discriminator.

## Reproduction

Run `collection.py --root /home/tkwa/.codex/worktrees/protium-fusion/website` with NumPy and SciPy. `collection-stopping-data.json` preserves selected primary measurement metadata; `collection.json` records source and code hashes, all stage/species ranges, acceptance tests and positivity checks. The existing stopping files and article are unchanged. The extension moves the low-energy analytic matching point to0.5MeV/u for tin so that its high-excitation-energy Bethe expression remains positive; this explicitly modeled completion is another reason to give the tin comparison less weight.
