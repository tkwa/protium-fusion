# First-contact compression and a small-pulse diagnostic ✨

The bounded compression calculation does **not determine a physically delivered final gain**. A spherical fluid started at 1000 kg/m³ maximum density can reach densities above 10⁸ kg/m³ in this surrogate, but burning results depend strongly on grid resolution and artificial viscosity. The higher-resolution burning integrations terminate at the 2.8 GK temperature guard. The cold equation of state and granular-to-fluid transition remain unvalidated.

## What was initialized

`cold_driver_probe.py` takes the saved N14 + p incoming profile, coarsens its mass grid to 32, 64 or 128 cells, and rescales velocities to preserve the original 3.24887 × 10²⁸ J of kinetic energy. It moves each node backward at constant velocity until the densest cell reaches the chosen contact density. The principal reference is 1000 kg/m³; 2000 kg/m³ is a compressed-material scenario. This constructs a radial mass/velocity requirement. It does not establish a pellet material, source or contact process.

The initial temperature is set uniformly to either 1000 K or 100,000 K. The code then treats the material as a fully ionized fluid. Discrete contacts, chemical binding, partial ionization, Coulomb corrections and multiphase material behavior are absent. The fully ionized model's electron degeneracy energy at low temperature is not the physical internal energy of a frozen pellet. A physical cold EOS and shock Hugoniot must replace this surrogate before its entropy or compression is interpreted as a deliverable state.

The electron/positron table was extended to 100 K using the same `pynucastro` source EOS, with 24 additional temperature nodes per decade. All original hot nodes were preserved. The source's documented degeneracy-root upper bracket was increased from 10⁷ to 10¹⁰; no physical EOS equation was altered. The successful table build took 164 seconds using two workers. The independent interpolation checks pass 6/6, with pressure and energy errors below 7.7 × 10⁻⁷ at the checked source points. This validates interpolation of that surrogate, not its cold-material physics.

Each run uses the existing spherical Lagrangian solver, constant opacity 0.001 m²/kg, and the reduced N14 + p capture/reverse network when burning is enabled. Nuclear screening, conduction, multidimensional mixing and weak decay are absent. The outer boundary expands into vacuum. A 1000 kg/m³ maximum contact density corresponds to an outer radius of about 22.9 km for the full 1.13 × 10¹⁶ kg pulse; a uniform sphere of that density is not assumed.

## Results and termination status

The density maxima below are sampled from at most 160 saved times. They can miss brief extrema. Nuclear/input ratios in guard-stopped runs are **values at the stop**, not final yields, gain predictions or rigorous lower bounds on a corrected physical model.

| Cells | Initial T, K | Maximum initial density, kg/m³ | Nuclear evolution | Termination | Peak T, GK | Sampled maximum density, kg/m³ | Nuclear/input at stop |
| ---: | ---: | ---: | --- | --- | ---: | ---: | ---: |
| 32 | 1000 | 2000 | Off | Requested final time | 0.377 | 4.54 × 10⁶ | 0 |
| 32 | 100,000 | 2000 | Off | Requested final time | 0.373 | 4.40 × 10⁶ | 0 |
| 32 | 1000 | 2000 | On | Requested final time | 0.420 | 4.28 × 10⁶ | 0.141 |
| 32 | 1000 | 1000 | On | Requested final time | 0.346 | 2.59 × 10⁶ | 0.0274 |
| 128 | 1000 | 1000 | On | **2.8 GK guard** | 2.800 | 4.43 × 10⁸ | **2.81, incomplete** |
| 128 | 1000 | 1000 | On; half viscosity coefficients | **2.8 GK guard** | 2.800 | 4.78 × 10⁸ | **6.84, incomplete** |
| 128 | 100,000 | 1000 | Off | Requested final time | 1.076 | 1.82 × 10⁸ | 0 |

The default artificial-viscosity coefficients are linear 0.5 and quadratic 2. The reduced-viscosity run uses 0.25 and 1. At initialization the mass-weighted ratio of artificial-viscous stress to physical pressure is about 7556 for the 32-cell, 1000 K, 1000 kg/m³ run; it remains about 1396 at 128 cells. Halving the viscosity coefficients halves that ratio and changes the burning evolution materially. Artificial viscosity is intended to generate shock entropy, but its large initial value in the smooth convergent flow makes a coarse-grid suppression of compression a numerical possibility. The 32-cell low yields cannot be used as evidence of physical failure.

The 128-cell burning runs reach their guards after 4.922 and 4.968 ms. The mass fractions above 10⁸ kg/m³ at sampled times reach only 0.0348% and 2.49%, respectively. A very dense central cell is insufficient evidence that the whole pulse attains the original favorable stagnation profile. The nonreacting 128-cell run shows that nuclear energy is not required for a small central region to reach the specified density within this surrogate. Its initial 100,000 K is a numerical sensitivity input, not a demonstrated heating design.

Energy residuals relative to initial input range from about 2 × 10⁻⁹ to 9.2 × 10⁻⁶. Good energy accounting does not establish the correct entropy generation or converged hydrodynamics. Individual successful runs took approximately 3–22 seconds; no large parameter survey was launched.

## Numerical failures retained

Two runs have structured failure records rather than a physical outcome:

- The 64-cell, 1000 K, 1000 kg/m³ burning run fails with a singular factorization, following negative-temperature Newton-trial warnings. It has no accepted final yield.
- The 128-cell, 1000 K, 1000 kg/m³ nonreacting attempt crosses the 100 K EOS table floor during a solver trial; an outer trial temperature reaches about 2.04 K. No final accepted state or trace was saved. This is not a measured physical cooling trajectory.

These records are included in `cold-driver-summary.json`. The second record was saved from observed terminal output because it preceded structured exception capture. A separate 128-cell nonreacting run at 100,000 K completed, as reported above; it does not repair or erase the failed lower-temperature attempt.

## What can be learned without a physical experiment

The next useful computation is a shock-convergent fluid calculation with a partially ionized, nonideal free-energy EOS and material-consistent initial state. The solver needs a verified convergent-compression/shock test and a viscosity or Riemann-solver convergence study. Only then should the full network be coupled to assess premature ignition and expansion. The current reduced network and guard-stopped trajectories cannot determine whether merger-generated heat ignites a useful fraction too early.

The discrete-to-continuum transition is a separate computation: a local pellet-assembly or jet-merger calculation needs solid/vapor contact, ionization, radiation exchange, collision geometry and angular/radial nonuniformities. It must produce an entropy distribution that the larger fluid calculation can accept. A three-dimensional calculation could then test perturbation growth and the mass fraction that reaches useful burn conditions. Even an excellent simulation would depend on material EOS, charging, opacity and damage data; the present output is not a substitute for those inputs.

This checkpoint therefore leaves both substantial derating and improved compression open. It does not attach a probability or fixed loss factor to either. It specifically rejects interpreting a low-resolution cold run as a settled reactor gain.

## A 30 m prescribed-hot diagnostic

At the parent's request, `scaled_hot_probe.py` reduces every initial radius in the existing 64-cell favorable N14 + p inflow by ten. Density, temperature and speed are unchanged. Mass and input energy decrease by 1000, and the hydrodynamic time decreases by ten. This is a useful scale sensitivity of the old prescribed-hot condition; it does not close the delivery gap described above.

The completed run gives:

| Quantity | Value |
| --- | ---: |
| Reference compressed radius | 30 m |
| Initial outer radius | 170.86 m |
| Mass | 1.13097 × 10¹³ kg |
| Initial delivered energy | 3.49316 × 10²⁵ J |
| Deposited capture energy | 2.85798 × 10²⁶ J |
| Capture energy / initial energy | 8.18164 |
| N14 conversion | 53.8400% |
| Peak temperature | 1.14905 GK |
| End time | 68.395 microseconds |
| Energy residual / initial energy | 2.095 × 10⁻⁷ |
| Termination | Requested final time |
| Wall time | 5.27 seconds |

This is one 64-cell diagnostic with the reduced three-species capture/reverse network and constant opacity. It has no resolution sweep or full-network cross-check at this new scale. The parallel transport screen finds a large optical-depth margin after the tenfold spatial reduction; the shorter hydrodynamic time nevertheless reduces the converted fraction. The result permits a small-pulse exploratory scenario to remain under consideration. It does not predict its real gain or a late-2040s deployment date.

This stage converts N14 + p to O15. Subsequent O15 decay and the other catalytic stages must be accounted separately before claiming a closed hydrogen-to-helium cycle. The short-run capture-energy ledger excludes the later decay heat.

## Energy convention and the proposed fictional range

For comparison with a plant, define `Q = deposited nuclear energy / delivered pulse energy`. Recover the original delivered energy once. With thermal efficiency `eta_h`, electrical-to-pulse efficiency `eta_d`, and complete collection,

`E_net/E_delivered = eta_h (1 + Q) - 1/eta_d`.

At `eta_h = 0.4` and `eta_d = 0.3`, the break-even condition is `Q > 7.3333` before auxiliary loads. At 99% collection it becomes `Q > 7.4175`. Thus a chosen fictional observed `Q = 2.6` loses electricity even in a favorable hypothetical complete-collection comparison. The fraction of energy intercepted by a diagnostic open range is much smaller.

`prototype_range.py` and `prototype-range.json` record an expressly **fictional scenario input**, Q = 2.6, for a 1.13 × 10¹³ kg single-stage shot at a launch radius of 10⁹ m and heliocentric distance 120 AU. They do not turn that number into a computed physical result. Its immediate nuclear plus delivered energy is 1.258 × 10²⁶ J, giving an isotropic fluence of 10.01 MJ/m² at the launcher radius. Delayed O15 decay is separate. The actual spectrum, directionality and arrival distribution determine shield performance.

A 10 MJ/kg assumed sacrificial-material enthalpy would give a simple energy/enthalpy quotient of approximately 1 kg/m² per shot. That quotient is not an ablation-depth prediction: surface heating, ion stopping, radiation transport, shock impulse and ejected shield material can dominate damage. Shielded nodes and shuttered diagnostics are an explicit engineering requirement. A projected collector area of 10¹² m² intercepts only about 8 × 10⁻⁸ of an isotropic shot at this radius. The proposed installation is a compression test range, not a complete collector or a closed catalyst-recovery plant.

Direct sunlight at 120 AU is approximately 0.0945 W/m². A black isothermal sphere has an equilibrium temperature near 25.4 K. At the RMS inflow speed, the 10⁹ m flight takes 417 seconds. An illustrative 10 MV/m, 10⁸ m long accelerator with a 1 GV/m pellet surface field allows a maximum radius of about 9.25 mm at that speed and 1000 kg/m³ material density; direct solar absorption during flight is then only about 3.2 J/kg. The accelerator may extend outward from the range's launch radius. Nearby warm equipment, prior-shot radiation and residual gas still require a design. The hot closed-cavity survival requirement does not apply to an isolated open-range flight illuminated only by the Sun.

A charge-neutralization or grain-discharge system is also essential: uncompensated accelerating charge would give an illustrative uniform-sphere Coulomb barrier near 175 km, well outside the 30 m core. Keeping that model’s residual electrostatic penalty below 1% of incoming kinetic energy requires more than 99.87% charge compensation. The delivery report and `prototype-range.json` retain the calculation; compensating electron mass is small, while charge transport and phase-space control remain unmodeled.

The 1 GV/m field corresponds to 4.43 MPa of surface Maxwell pressure. If a separate conducting shell carries that load, a 1 GPa allowable membrane stress and 2000 kg/m³ shell density require shell mass about 1.33% of a 1000 kg/m³ core. A thin conducting coating can be much lighter only if the core and adhesion carry the stress. Shell mass, contaminants and hydrogen retention cannot be omitted from a physical pellet specification.

Solar differential acceleration across the range is about 4.59 × 10⁻¹¹ m/s², corresponding to 4 micrometres of relative displacement over that characteristic flight. Self-gravity is similarly negligible for pulse compression: a uniform 30 m sphere's binding energy is only about 15 J/kg. Orbit-following nodes rotate their coordinate frame enough to require a predictable aiming correction of roughly 63 m over the flight. These arithmetic checks reveal no gravitational need for a massive support shell. They also do not validate source pointing: an uncorrected 30 m footprint demands an angular error below 3 × 10⁻⁸ rad and timing error below 12.5 microseconds.

A test range can therefore be written without implying a full multi-Earth-mass enclosing collector. Its shot delivery, shielding, source power, material consumption and incomplete cycle remain explicit costs. The particular observed gain and merger losses in the article are worldbuilding decisions, not outputs of these simulations.

## Reproducibility

`cold-driver-summary.json` contains the seven saved traces and two numerical-failure records. `cold_driver_probe.py`, `diagnose_cold_driver.py`, `validate_cold_eos.py`, and the saved initial JSONs reproduce the assumptions and diagnostics. `scaled_hot_probe.py` reproduces the separate 30 m prescribed-hot run. `prototype_range.py` uses only elementary arithmetic and checks its energy and acceleration accounting.

Run the numerical scripts with the existing `/tmp/protium-python` dependencies, `PYTHONPATH=/tmp/protium-python`, `PYTHONDONTWRITEBYTECODE=1`, `OPENBLAS_NUM_THREADS=1`, `NUMBA_NUM_THREADS=2`, `NUMBA_CACHE_DIR=/tmp/protium-numba`, and `MPLCONFIGDIR=/tmp/protium-mpl`. The source hydrodynamic/EOS files remain unchanged. `MANIFEST.json` records packet and source hashes.

The underlying incoming-state construction is described in [RADIAL-MODEL.md](../../RADIAL-MODEL.md), and the transport qualifications in the [parallel transport report](../pulse-transport/README.md). The delivery report links the relevant accelerator and plasma-liner primary literature. The solar mass parameter and astronomical-unit convention follow the [IAU nominal-constant note](https://arxiv.org/abs/1510.07674).
