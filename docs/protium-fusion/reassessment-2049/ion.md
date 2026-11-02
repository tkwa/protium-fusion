# Ion reactor: reassessing the 2049 outcome ✨

This memo supports the article's choice of a future history. The calculations establish restrictions on particular clearing fields and a conditional energy envelope. They establish neither a complete selective collector nor a commissioning date. In the established ASDI/APM civilization, **a limited net-positive installation by 2049 is a defensible choice**, provided its history includes a successful nonlinear capture and recovery programme. Roughly 0.2 mW/kg with 1.75–2 MeV of additional processing per primary reaction is a more forgiving choice than approximately 0.4 mW/kg with 1.2 MeV. The [plant memo](plant.md) owns the final mass, sky-load, auxiliary and electricity calculation.

The important distinction is between an unsolved design and a physical obstruction. A reactor that failed to remove most of its trapped positrons would fail energetically; our incomplete calculation does not establish that future engineers could not build a different extractor. Conversely, a literature demonstration of autoresonance does not qualify this reactor. The selected history must supply a specific successful preparation, capture, transport and energy-recovery arrangement, whose realized performance is an invented future result recorded as such in the technical documentation.

## Source point and scope

The new screen uses the actual source selected at 60% heat conversion in `network-plant.json`. Changing its radiator from 30 to 50 K does not change that source selection. Its main Ti processor has B = 137.1491 T, wall radius 22.6163 mm and 1.318593×10^11 primary reactions per active metre per second. One MeV per primary is **0.0211262 W per active metre**. The old 70% source had a 20.7316 mm radius and a different rate; its detailed susceptibility coefficients cannot be transplanted unchanged.

[ion-cara-screen.py](ion-cara-screen.py) is a small standard-library script; [ion-cara-screen.json](ion-cara-screen.json) contains its inputs, cases, limitations and SHA-256 hashes. Run `PYTHONDONTWRITEBYTECODE=1 python3 docs/protium-fusion/reassessment-2049/ion-cara-screen.py` from the dedicated worktree. Checks verify retarded-time kinematics, equality of the two optimized variable costs and the stipulated geometric acceleration allowance. They verify the specified screen, not nonlinear particle capture. No large orbit kernels were rerun.

All work for this reassessment was written inside `/home/tkwa/.codex/worktrees/protium-fusion/website`, with temporary preparation files under `/tmp`. The main checkout and article were not edited by this agent.

## What short cells actually improve

Dividing the original approximately 20 km channels into shorter independently controlled cells gives distributed pumping, feed points, diagnostics, wave shaping and finite extraction regions. It reduces transport lengths and allows staged maintenance. The complete installed active length remains enormous: shortening the cells does not reduce a field's absorbed power per active metre.

These are physical ion cells with local collector/barrier regions, rather than RF feed subdivisions of an uninterrupted 20 km plasma. A positron accelerated in a long uninterrupted column would still need a timely through-path to a collector before radiation consumes its axial energy headroom. More ends add high-voltage apparatus, switching capacitance, leakage and support. End-field/fringe occupancy, transfer gaps, RF access, electrode support and new recovery circuits must fit the added mass and power reserves; the active-bore factor below does not count them. For a Maxwellian proton tail, retaining the same fractional end-loss rate when length changes requires approximately `Delta V = (kT/e) ln(L_old/L_new)`. At 1.6 MeV, reducing 20 km to 1 km adds about 4.8 MV to the barrier requirement, before the actual end distribution and fields are recalculated. Indiscriminately opening positive end gates also releases the ions. A successful selective scheme needs its own phase-space transport, while gate-field recovery and ion losses remain in the circuit.

The earlier [axisymmetric wall-wave calculation](../phase-two/ion-dynamics/README.md) gives a useful counterexample. In the old source, a 1 kV, k = 40/m, 0.1c maintained wave heats ions at approximately 51.5 kW/m. A one-microsecond pulse deposits 0.0515 J/m, corresponding to 1.496 MeV per primary at one pulse per second. It can be a small perturbation to the stored plasma while consuming a material share of the operating margin. The old extra 1 mW/kg margin was 0.0384 W/m, not several watts per metre.

Some absorption can replace planned ion heating. If H is absorbed ion heat and h is the amount that actually displaces existing preparation, the incremental electrical cost is `H/eta_RF - h/eta_heater - eta_heat*(H-h)`. The species, distribution and timing must make that substitution useful. A combined heater and isotope separator can therefore be better than insisting on negligible ion coupling. Canonical torque and its recovery must be accounted for separately whenever the field transfers angular momentum. The [transverse clearing analysis](../TRANSVERSE-CLEARING.md) supplies the relevant heat/torque identities and collective-response warning.

## A high-frequency accelerator: useful kinematics, inadequate simple capture

A circularly polarized electromagnetic wave travelling along the field can supply axial momentum while resonantly increasing lepton energy. This avoids the particular problem of heating a positron purely transversely to tens of MeV and waiting for it to reach a distant magnetic nozzle while it radiates. Cyclotron autoresonance acceleration, or CARA, is a real accelerator mechanism [1–3]. It has been demonstrated for a prepared energetic electron beam [4]. Those results do not supply broad capture of the reactor's initially trapped positrons.

For an ideal uniform-B plane wave and a particle initially at rest, use `u = p/(m_e c)`, `a = eE/(m_e c)` and `Omega = eB/m_e`. The invariant is `s = gamma-u_z = 1`; on resonance,

```
gamma = 1 + u_perp^2/2
u_z = gamma - 1
Delta eta = u_perp/a
t = (u_perp + u_perp^3/6)/a
z = c*u_perp^3/(6*a)
t - z/c = Delta eta
```

The RF pulse's retarded duration `Delta eta` is much shorter than the particle's laboratory transit. The RF power calculation must use the pulse duration, and the physical interaction region must accommodate the transit distance.

At 137.149 T the resting-particle resonance is 3.839 THz, wavelength 78.09 micrometres. A 60 MeV particle has gamma = 118.417, `u_perp = 15.324`, `u_z = 117.417`, pitch tangent 0.1305 and gyroradius 0.190 mm. The ideal acceleration screen is:

| Electric amplitude | Acceleration length | Particle transit | Retarded pulse | First-order synchrotron energy |
|---:|---:|---:|---:|---:|
| 0.1 MV/m | 3,065 m | 10.48 microseconds | 261.2 ns | 2.722 MeV |
| 1 MV/m | 306.5 m | 1.048 microseconds | 26.12 ns | 0.2722 MeV |
| 10 MV/m | 30.65 m | 0.1048 microseconds | 2.612 ns | 0.02722 MeV |

Radiation was integrated on the unperturbed orbit using `tau_perp = 3*pi*epsilon0*m_e^3*c^3/(e^4*B^2)` and `E_rad = (m_e*c^2)/(2*tau_perp*a)*(u_perp^3/3 + u_perp^5/10)`. The nonrelativistic perpendicular cooling time is 137.1 microseconds. This first-order work does not maintain resonance against radiation reaction. Dispersion, the radial source field, rotating equilibrium, finite ends, taper and chirp remain to be solved. In particular, the strong source field cannot be replaced with a uniform empty solenoid merely because the gyrofrequency is high.

There is a decisive limitation to the simple monochromatic version. For a general initial particle, `s = gamma-u_z` is conserved by the ideal travelling wave and resonance requires `omega*s = Omega`. A detuning `delta = Omega-omega*s` limits driven transverse momentum to `2*a/abs(delta)`. Reaching the target at 1 MV/m therefore requires fractional detuning below **3.174×10^-6**, approximately a 12.2 MHz half-width at 3.839 THz. A warm population spans much more than this interval.

For illustrative 100 eV, 1 keV and 10 keV axial Maxwellians, generously allowing two separately placed resonance bands gives maximum addressed fractions of 3.621×10^-4, 1.145×10^-4 and 3.621×10^-5. Applying these fractions to the conditional ledger below gives respectively **12.17, 21.47 and 38.01 MeV per primary**, before other handling. Real fixed-frequency repeated pulses need not sample fresh particles independently, so even this reuse of the stochastic-removal ledger can be optimistic. These temperatures are examples, not calculated reactor lepton distributions.

Increasing the amplitude broadens this band in proportion to E, but the empty pulse energy also rises in proportion to E. A frequency sweep alone does not repair the invariant: particles already locked to one ideal s are detuned by changing the resonant frequency. The attractive 10% acceptance case below requires a physically different capture stage. It cannot be called a result of this simple CARA calculation.

## A live nondegenerate capture architecture

A concrete research direction is a finite-cell bank with paired, independently phased counterpropagating circular modes in its preparation regions, followed by a predominantly travelling accelerating mode. The preparation fields form standing/beat structures; a spatial magnetic or electrostatic taper and frequency programme vary their resonances. This combination breaks the single travelling wave's exact `gamma-u_z` invariant. It offers nonlinear phase-space islands that can sweep through cohorts of the trapped distribution. The prepared cohorts are handed to a separately phased travelling accelerator and final recovery channel. Both polarizations, charge signs, finite field matching and induced ion response must be treated in the Maxwell/kinetic design.

This is an identifiable physical mechanism rather than an assumed unscreened electrostatic force. Spatial autoresonance theory [3] and experiments on autoresonantly driven trapped plasmas [5] provide starting points. Neither establishes its acceptance, RF fluence or stability in the Ti source. The programme must demonstrate that its genuinely nondegenerate drive captures a sufficient population and retains it as energy rises. Exporting particles to an external cold preparation trap would also help, but the transport into that trap is itself part of the unresolved extraction problem.

The phase space is processed in cohorts; Hamiltonian compression of an entire warm distribution into one arbitrarily cold packet is not assumed. Entropy leaves through radiation, rejected-particle handling and final collectors, all with assigned energy and apparatus. Empty-field loss, unsuccessful kicks and preparation pulses must be included in measured fluence. The permanent uncaptured population needs an explicit purge route. Merely citing repeated 10% attempts silently leaves an inaccessible subset accumulating forever.

The successful-history choice can be that this combined preparation-and-acceleration system was invented and empirically qualified. Its achievement is broad-population locking with acceptable total fluence, rather than discovery that an isolated monochromatic plane wave has broad capture. This is a substantial engineering invention. The present calculation leaves its existence and performance open; no conservation or phase-space argument here excludes it.

## The conditional RF envelope and its hardware

At 1 MV/m the uniform-bore circular mode has 0.014228 J/m field energy, 4.265 MW peak wave power and 0.111415 J per empty pulse. A free Gaussian beam with waist equal to the 22.6 mm bore radius has a Rayleigh length of only 20.6 m. A 306 m interaction therefore requires a guide or renewed injection. Corrugated hybrid-mode guides are a physical candidate [6], with loss and mode purity requiring qualification at 3.84 THz in the irradiated reactor geometry.

The screen assigns 1 km source cells, effective power attenuation alpha = 0.01/m (0.04343 dB/m), 50% RF generation efficiency and two pulses per complete attempt. Single-end transmission through 1 km would be only 4.54×10^-5; distributed feeds must replace local wall loss and beam loading. The laboratory guide results in [6] are at substantially lower frequency and do not demonstrate this attenuation requirement. Peak source, coupler, switch and cooling mass is additional to mean-power arithmetic.

For bore area S, `W_pulse = epsilon0*S*E^2*c*Delta eta`. Empty RF electricity per active Ti metre per complete attempt is

```
E_empty = G * n_directions * W_pulse * (alpha + 1/L) / eta_RF
```

The full network occupies 1.3588695 times the Ti-only active bore volume per primary. This factor G charges equal-field, equal-duty clearing across all fifteen processors; Ti alone is 73.59% of that volume. Stage-specific timing may improve the result only after the appropriate stage sources and extraction schedules are supplied. Exiting empty pulses are discarded in this envelope, with no heat or electrical recovery credit. Distributed beam loading is charged separately.

For constant birth b per primary, attempt frequency f and independent removal fraction p, the average susceptible residence is `(1/p-1/2)/f`. With retained-lepton radiation P drawing energy from the ions,

```
C_added = A*f + B/f + b*W_removed/eta_RF
A = E_empty/(R_primary * 1 MeV in joules)
B = b*P*(1/eta_heater-eta_heat)*(1/p-1/2)
f_opt = sqrt(B/A)
```

The diagnostic uses b = 0.000906 prompt pairs per primary, P = 100 MeV/s per retained positron, heater efficiency 0.95 and heat conversion 0.60. The birth rate is one recorded prompt-pair scenario, not a total-source upper bound. P is a sensitivity input, not a measured steady lepton loss. The total supplied RF work W_removed = 120 MeV per removed positron allows for the 60 MeV target, radiation and unsuccessful excitation; fitting the real kicks into that allowance is an explicit requirement.

At the **hypothetical** p = 0.10, the optimum is 1.112 complete attempts/s, mean susceptible residence 8.547 s and local ideal-pulse duty 5.807×10^-8. Empty-field electricity costs 0.35049 MeV/primary, retained-lepton heating replacement another 0.35049, and particle RF work 0.21744: total **0.91842 MeV/primary**. The average pulse success removes about 1.075×10^8 positrons per active metre per attempt, and supplies 0.002066 J/m of particle RF work. Beam loading therefore must be replenished along the path; comparing only the empty pulse energy with the whole particle beam would be inconsistent.

At the same attenuation, p = 0.05 costs 1.222 MeV, and p = 0.01 costs 2.486 MeV. At p = 0.10 but ten times the attenuation, the cost is 2.342 MeV. These sensitivities explain why broad capture and actual guide construction can decide net operation. The low duty factor alone does not demonstrate a cheap machine.

## Electrons, heavy products and the remaining circuit

Electron endpoints are more favorable electrostatically than positron endpoints: the positive ion-end electrodes attract them. The ions nevertheless screen weak applied axial wall voltages. A credible successful installation uses finite locally controlled collector regions and shaped guiding fields, with the co-born positive gas products routed back through recovery and preparation. It cannot replace this with an assertion that a low wall pressure guarantees an electron-free core.

For scale, a fresh-H2 wall throughput of 10^-10 Pa m^3 s^-1 m^-2 gives a complete-stripping ceiling of approximately 0.052 electrons per primary for the new Ti profile at 300 K. At an assumed complete electrical cost of 0.1 MeV/electron that is approximately 0.0052 MeV/primary; 100 times the gas throughput raises the figure to about 0.52 MeV. These are requirements/sensitivities, not irradiated-material measurements. The earlier gas screen found roughly 30% first-flight ionization, with much larger accumulation after repeated flights when sticking is weak. Local pumping improves transport without erasing that first encounter. Energetic secondaries, other gas chemistry and particle-induced desorption need qualification. Co-born positive and negative electrostatic energies initially cancel; the birth potential is a gross circuit scale, not an automatic net electrical loss.

The heavy-product route can combine resonant isotope separation with useful ion heating, steering daughter cohorts into shielded low-density recovery channels. Bare q/m differences help select a drive but do not establish selectivity in a rotating, interacting plasma; collective frequency shifts and overlap must be calculated [7]. The Ti inventory divided by its primary rate is about 210 s for the new source. A chosen 0.1 s daughter residence would imply an inventory fraction around 4.8×10^-4 before secondary reactions; its acceptability requires the corresponding nuclear loss calculation. The earlier 1 ms target was a diagnostic choice, not a universal physical requirement.

Coherent high-voltage circulation and random thermal energy require separate recovery accounts. The baseline already pays for its large daughter charge circulation and converter throughput. A new RF scheme should neither charge that same useful work twice nor omit new circulation, torque restoration, gate capacitance or converter mass. High beam-energy recovery in an ordered accelerator does not establish the same recovery of a warm ion distribution. Adiabatic direct-conversion theory identifies real geometric and dynamic limitations [8]; changing a quoted thermal efficiency is an additional machine design, not a free percentage adjustment.

A rare full-plasma purge may remove a tiny permanently inaccessible tail without becoming the normal operating cycle. Its actual new-profile inventory, recovery, downtime, converter peak power and complete tail clearing need closure. The script retains a clearly labelled old-70%-profile reset screen only as a scale diagnostic. It does not transfer that result into the selected 60% plant.

A selected 1.75 MeV total added-handling result leaves 0.8316 MeV beyond the illustrative 0.9184 MeV positron envelope. That remainder must cover every additional preparation field beyond the specified pulses, residual electrons, product-handling increments, rare purges, recovery imperfections and omitted radiation channels. It is not independently validated spare capacity. The strongest article treatment is a measured complete-system processing burden, with subsystem calorimetry and charge/energy closure as named commissioning achievements. The larger 2 MeV choice has more handling headroom but lower net output. The final plant ledger should hold previously allocated auxiliaries fixed and include the full sky load, as the plant review recommends.

## Defensible selected history and reasons for a later date

The established industrial chronology supports large parallel theory, manufacturing and empirical campaigns after the 2035 swarm. It supports a plausible sequence of choices:

- During the late 2030s and early 2040s, source cells, high-field supports, irradiated liners and isotope loops were qualified separately. Cold V stocks were prepared early enough for the approximately eleven-month neutral decay to be included in later complete-cycle campaigns.
- In the middle 2040s, finite-cell experiments mapped the actual lepton distribution, rejected simple narrow-band clearing and developed the nondegenerate preparation stage. The decisive progress was reproducible broad capture with recorded failed-kick energy, followed by guided electromagnetic acceleration and local recovery.
- In 2047–2048, distributed feeds, daughter sorting, useful heater substitution, electron collectors and rare purges were integrated. The complete Ti/V, neutron-capture and deuterium-burning return sequence was operated with charge and energy accounting.
- A 2049 commissioning result can then be selected at roughly 0.2 mW/kg, with a measured total processing burden near 1.75–2 MeV/primary and the additional plant apparatus charged explicitly. Limited deployment is compatible with that technical success and with its poor mass productivity relative to the established energy system.

These dates and achieved handling performance are proposed worldbuilding choices. They are not inferred from the ideal pulse or from present accelerator demonstrations. The article should describe the chosen achievements directly; this memo preserves the real unresolved physics and the alternative outcome.

A later date remains equally physically credible if broad capture requires too much empty-field fluence, guide/coupler losses exceed the allowance, isotope separation disrupts the source, irradiated surfaces exceed the gas budget, or actual high-field structures cannot sustain their assigned working conditions. An installation might then operate useful source and nuclear-cycle experiments while consuming net electricity. Those are concrete coupled engineering reasons to stall. Absence of a complete design in this investigation is not itself one.

The fundamental restrictions are conservation of energy, charge and angular momentum; phase-space and entropy accounting; magnetic stress; finite field bandwidth for a specified geometry; and the real weak-decay clock. Nonlinear capture, waveform design, recovery routing, stable distributions and module optimization are theory and simulation tasks with experimentally testable outputs. Surface yields, insulation, guide loss, irradiated conductor/support performance and replacement lifetimes remain empirical. ASDI and APM make a strong successful branch credible without abolishing any of those requirements.

## Primary sources and the limits of their use

1. K. H. Tsui, *Electron-beam acceleration by cyclotron-autoresonance interaction*, Physical Review E 51, 649 (1995), [DOI](https://doi.org/10.1103/PhysRevE.51.649). Accelerator theory; no reactor-population acceptance result.
2. Y. I. Salamin et al., *Feasibility of electron cyclotron autoresonance acceleration by a short terahertz pulse* (2015), [author preprint](https://arxiv.org/abs/1504.00990). Numerical vacuum-pulse proposal; no qualifying experiment for the present source.
3. L. Friedland, *Spatial autoresonance cyclotron accelerator*, Physics of Plasmas 1 (1994), [author PDF](https://lazarfriedland.huji.ac.il/sites/default/files/lazar_friedland/files/76.pdf). Spatially varied resonance is a concrete theoretical starting point for a genuinely different capture stage.
4. M. A. LaPointe et al., Physical Review Letters 76, 2718 (1996), [DOI](https://doi.org/10.1103/PhysRevLett.76.2718). High RF-to-beam conversion for an injected approximately 90–100 keV, multi-megawatt electron beam; this does not demonstrate dilute thermal-positron capture or complete wall-plug efficiency.
5. G. B. Andresen et al. (ALPHA), *Autoresonant Excitation of Antiproton Plasmas*, Physical Review Letters 106, 025002 (2011), [DOI](https://doi.org/10.1103/PhysRevLett.106.025002). Trapped-plasma experiment with distribution-dependent response; it establishes neither MeV-source selectivity nor the present energy budget.
6. E. A. Nanni et al., *Low-loss Transmission Lines for High-power Terahertz Radiation* (2012), [DOI](https://doi.org/10.1007/s10762-012-9870-5), [MIT manuscript](https://library.psfc.mit.edu/catalog/reports/2010/11ja/11ja053/11ja053_full.pdf). The reported 19 mm corrugated-guide measurements near 250 GHz are a physical precedent, not measured 3.84 THz reactor attenuation.
7. J. M. Dawson et al., Physical Review Letters 37, 1547 (1976), [DOI](https://doi.org/10.1103/PhysRevLett.37.1547), and M. Affolter et al., Physics Letters A (2014), [DOI](https://doi.org/10.1016/j.physleta.2014.06.015). Ion-cyclotron separation and interacting-plasma frequency shifts motivate selective handling while limiting vacuum-frequency arguments.
8. J.-M. Rax, E. J. Kolmes and N. J. Fisch, *Efficiency and Physical Limitations of Adiabatic Direct Energy Conversion in Axisymmetric Fields*, PRX Energy 4, 013007 (2025), [DOI](https://doi.org/10.1103/PRXEnergy.4.013007), [author preprint](https://arxiv.org/html/2410.15485v1). Geometric and dynamical conversion analysis, not an achieved recovery percentage for the proposed plant.

The independent CARA audit is maintained by the gas-review agent. Its two important corrections—thermal resonance width and full-network guide loading—are incorporated above and in the saved script.
