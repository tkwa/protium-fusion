# Driver delivery and first-contact compression ✨

> **Scope:** This study uses an earlier pulse family, often at larger mass than the current design. Its cold-material EOS and unfinished or failed runs remain useful physical and numerical diagnostics, but its operating points, gains and any illustrative future scenarios are not current reactor predictions. Compare against the [current mixture design](../../mixture-design/README.md) and [corrected 30 m regression](../../b2-completion/regression.md). Saved input hashes refer to the original research snapshot.

This checkpoint replaces the unspecified delivery system with explicit requirements for one candidate: a distributed accelerator launches cold, composition-controlled pellets from the collector, on a velocity and timing schedule that assembles the incoming radial mass profile. It has not established a working accelerator, pellet material, merger process, or complete driver. The subsequent cold-contact compression calculation is deliberately a surrogate test; its cold equation of state is not a measured material model.

## What the saved pulse actually demands

The present main pulses contain about 1.1–1.4 × 10¹⁶ kg. Their mass-weighted kinetic energy corresponds to speeds of 2.27–2.59 million m/s. At a launch radius of 2 × 10¹⁰ m, the median flight is 2.3–2.6 hours. For the nitrogen-14 stage, the middle 90% of the mass is launched between about 1.6 and 4.7 hours before the prescribed incoming state. Slower material near the center starts earlier. Launching different radial layers at different times avoids requiring one simultaneous several-megatonne-per-nanosecond discharge from the collector.

| Stage | Input mass, kg | RMS kinetic speed, km/s | Peak source electric power at an assigned 30% conversion, W |
| --- | ---: | ---: | ---: |
| N14 + p | 1.131 × 10¹⁶ | 2397 | 2.46 × 10²⁵ |
| N15 + p | 1.414 × 10¹⁶ | 2273 | 2.58 × 10²⁵ |
| C12 + 2p | 1.131 × 10¹⁶ | 2589 | 3.09 × 10²⁵ |
| C13 + p | 4.189 × 10¹⁴ | 2339 | 8.32 × 10²³ |
| Helium-to-carbon factory | 1.131 × 10¹⁵ | 6304 | 4.42 × 10²⁵ |

These powers come from finite-volume launch-time bins. They include delivered kinetic energy and the assigned efficiency, but not a physical accelerator's additional losses. A 512-cell nitrogen profile changes the peak by 0.75% and the median flight by 0.14% relative to 128 cells. The source-conditioning equipment therefore needs a peak-power rating substantially above the plant's approximately 2 × 10²³ W average driver consumption. The previously allocated 10²⁵ kg driver facility might accommodate some of this hardware, but its mass must be designed against these ratings rather than counted a second time without explanation.

The final compression concentrates stored projectile energy; it does not require remote electrical power to follow the millisecond burn instantaneously. The existing electrical-storage allocation still needs a charging and discharge schedule consistent with these hours-long source pulses and overlapping earlier launches.

## A coupled accelerator–thermal-survival requirement

For a conducting spherical pellet with radius `a`, material density `rho`, and surface electric field `Es`,

`q/m = 3 epsilon0 Es / (rho a)`.

A constant applied field `Ea` over length `L` can supply at most `(q/m) Ea L` of specific kinetic energy. The surface Maxwell pressure is `epsilon0 Es²/2`. These are elementary electrostatic requirements, not specifications that a cryogenic pellet has demonstrated. An induction accelerator can accumulate the integrated accelerating voltage in stages; the calculation does not require every stage to sit at the summed voltage relative to ground.

A pellet flying through an approximately isotropic cavity radiation field at temperature `Tw` absorbs, while remaining much colder than that field,

`e_abs = 3 epsilon_P sigma Tw⁴ t_flight / (rho a)`.

Here `epsilon_P` is its Planck-mean absorption efficiency for that radiation. It is not automatically equal to its optical reflectivity or to the absorber's emissivity. The calculation permits an explicit 1 MJ/kg of absorbed energy before the assumed cold state is lost. It does not assign that value to actual hydrogen, nitrogen, helium, or a composite pellet.

For the nitrogen pulse, a scenario with `Es = 1 GV/m`, `Ea = 10 MV/m`, `L = 10⁹ m` and pellet density 1000 kg/m³ allows a radius of at most 9.25 cm at the RMS speed. At a cavity temperature of 1000 K, such a pellet needs `epsilon_P <= 0.065` to stay within the 1 MJ/kg allowance. The carbon pulse requires about 7.93 cm and `epsilon_P <= 0.060`. The faster helium factory needs about 1.34 cm and `epsilon_P <= 0.025`. These are RMS-speed requirements; the fastest nitrogen shell requires a pellet radius below 2.68 cm and absorption efficiency below 0.035, so a real launch train needs a speed-dependent pellet size or stronger acceleration. The surface Maxwell pressure in this scenario is 4.43 MPa; field emission, mechanical strength, conductivity and coating stability remain material requirements. If a separate shell carries that pressure, a 1 GPa allowable membrane stress, 2000 kg/m³ shell and 1000 kg/m³ core imply shell mass about 1.33% of core mass. A lighter conductive film requires core/adhesion support. The saved shell scenarios retain this added material explicitly; their mass is not yet included in the nuclear pulse.

Reducing the accelerator to `Ea = 1 MV/m` and `L = 10⁸ m` reduces allowed radius by a factor of 100. At 1000 K the nitrogen pellet then needs absorption efficiency below 0.00065. At 500 K the thermal constraint is sixteen times weaker. The old plant model never selected a complete collector/heat-engine temperature history, so the 300, 500 and 1000 K calculations are scenarios. The 40% thermal-conversion assumption with a 300 K sink requires a hot reservoir above 500 K in a reversible heat-engine comparison; a transient collector needs its actual temperature history and conversion cycle.

A hot cavity therefore creates a real trade between small, easily accelerated pellets and larger, thermally protected pellets. This screen neither proves a material satisfying the joint requirement nor establishes that no such material is possible. Using uncoated transparent grains, reflective shells, cold injection tunnels, a nearby sacrificial launcher, or ion beams would require separate mass, optics, charge and energy calculations.

The laboratory precedent is limited but relevant: Colorado's electrostatic dust facility accelerates micron and submicron conducting or coated particles, with measured individual charge and speed, to roughly 100 km/s. Cryogenic fusion-fueling pellets are a separate technology, usually millimetre-scale at roughly 100–1000 m/s. Neither establishes the performance of the proposed larger, much faster composite pellets. [Colorado facility](https://impact.colorado.edu/dustaccelerator); [ORNL pellet-injector history](https://www.ornl.gov/publication/pellet-injector-technology-brief-history-and-key-developments-last-25-years).

## Pointing, timing and gas

If the last trajectory correction is at the collector, limiting arrival error to 300 m requires an angular error below 1.5 × 10⁻⁸ rad. For the nitrogen pulse it corresponds to an uncorrelated transverse velocity error below 0.036 m/s, an arrival-time error below 125 microseconds, or a fractional uncorrected longitudinal velocity error below 1.5 × 10⁻⁸. A 30 m tolerance tightens each by ten. Moving the final correction to 1000 km from the focus relaxes the velocity-error requirement by 20,000, but exposes that correction system to the central pulse or requires remote force delivery.

These are free-flight sensitivity coefficients, not universal emittance limits. A lens can create position–velocity correlations; a rigid pellet's internal atomic temperature does not equal its center-of-mass velocity dispersion. Neither observation supplies the required tracking, corrections, field quality or survival automatically.

For intercepted stationary gas column `Sigma_g`, complete deposition of relative kinetic energy gives

`e_gas = 3 Sigma_g v² / (8 rho a)`.

The saved scenarios record the corresponding maximum gas column and equivalent uniformly distributed cavity-gas mass. The residual gas need not be uniform, and fast particles may penetrate small pellets, so a real collision/charging calculation is required. Recovering catalyst onto the collector and clearing gas from pellet trajectories are related but distinct requirements.

## Charge compensation before convergence

A pellet leaves an electrostatic accelerator carrying charge. The inward pulse cannot retain that charge without a charge-compensation system. In the small open-range example, 9.25 mm pellets at a 1 GV/m surface field carry about 0.00287 C/kg. Charging the complete 1.13 × 10¹³ kg pulse that way gives about 3.25 × 10¹⁰ C.

For an illustrative unscreened uniform sphere, `U_C = 3 Q_charge²/(20 pi epsilon0 R)`. Its collective Coulomb energy equals the prescribed incoming kinetic energy at about 175 km radius, long before the 30 m compressed reference. At 30 m the ratio would exceed 5800. To limit that particular geometry's electrostatic penalty to 1% of incoming kinetic energy, the residual net charge fraction must be below about 0.00131, corresponding to at least 99.87% charge compensation. Actual shell and beam geometry changes the coefficient; this bound demonstrates the missing operation rather than describing a finished neutralizer.

The compensating electrons have only about 0.185 kg of rest mass. Material supply is therefore easy compared with electron delivery, currents, field control and maintaining the desired pellet phase space. The isolated grains' stored surface-field energy is approximately 13.3 kJ/kg, which is also a discharge-heating budget to track. A driver could neutralize or discharge pellets after acceleration; then its final tracking/correction system must operate on that state. Plasma neutralization techniques in laboratory ion-beam drift compression establish a mechanism to investigate, not performance at this pulse charge or scale. [Seidl et al., 2007](https://nonneutral.pppl.gov/pdfpapers2007/SeidlHIF06.pdf).

## Where the pellets must become a fluid

A constant-velocity backward map of the nitrogen input reaches a maximum cell density of 2000 kg/m³ about 3.69 ms before the saved initial state. Its outer radius is 18.2 km and its mean density 448 kg/m³; the outer material is still sparsely packed when the densest inner cells first touch. Carbon gives a similar result. This map is a useful first-contact geometry, not an adiabatic backward solution: pressure acceleration has been omitted.

The original fully ionized EOS, continued isentropically from the central incoming cell, gives the following diagnostic contact states at 2000 kg/m³:

| Stage | Temperature | Specific internal energy |
| --- | ---: | ---: |
| N14 + p | 0.256 MK | 2.43 GJ/kg |
| C12 + 2p | 0.252 MK | 2.82 GJ/kg |
| C13 + p | 0.255 MK | 2.31 GJ/kg |
| N15 + p | 9.08 MK | 73.5 GJ/kg |

The helium-factory continuation hits the existing 0.100 MK EOS floor at 5481 kg/m³; it is not extrapolated to 2000 kg/m³. These values are **not usable cold-material specifications**. The parallel microphysics audit finds strong ion coupling and incomplete ionization in the low-temperature nitrogen/carbon contact regime. Their fully stripped ideal-gas entropy and compressibility are unreliable there. Complete ionization can cost several times the tabulated contact internal energy while still costing below roughly 0.5% of the original whole-pulse input. A small change in total driver energy can therefore accompany a major change in its adiabat.

Simple head-on collision is not a substitute. Thermalizing all of the nitrogen RMS kinetic energy at this density would deposit about 2.87 TJ/kg, over a thousand times the ideal diagnostic contact internal energy. For two symmetric streams meeting at half-angle `theta`, the kinetic energy lost on merging is `0.5 v² sin²(theta)`. Supplying only the diagnostic contact heat corresponds to about 1.67 degrees for nitrogen; allowing only one tenth of it as additional merger heat reduces this to about 0.53 degrees. These angle scales illustrate the entropy sensitivity. They are not validated tolerances for a spherical pellet assembly.

Experiments with merging plasma jets have already found shock-produced nonuniformities and dependence on jet balance and impurities. The 36-jet spherical PLX experiment uses tens of km/s jets and explicitly studies liner formation; its geometry and operating regime are much smaller and slower than this candidate. [Yates et al. 2020](https://arxiv.org/abs/2002.03006); [LaJoie et al.](https://arxiv.org/html/2401.11066v1).

## External preheating cannot be presumed uniform

For opaque spherical grains, the sum of projected cross-sections is `3 M/(4 rho a)`. Dividing by the projected area of the cloud gives a mean grain optical depth. In the nitrogen pulse with 10 cm grains at 1000 kg/m³, that depth reaches unity at about 5200 km radius, long before the approximately 18 km first-contact envelope.

If grains are vaporized at that larger radius and each vapor cloud expands independently at the diagnostic contact sound speed, the expansion distance during the remaining flight is about 110 km. Keeping that simple expansion distance below 300 m would require vaporization within roughly 14 km of the center. Thus “heat every grain externally, then let it coast into the original fluid state” has a serious shadowing and expansion mismatch. Collective hydrodynamic compression, internal energy delivery, radiation cooling or a shock sequence could change the outcome; each must be modeled. The mismatch does not disprove every pellet driver.

The useful next calculation starts cold at first contact and follows compression and shock formation forward. It asks whether the delivered mass and velocity field can create a usable dense burning region on its own, including its changed entropy, rather than imposing the desired compressed state.

## Reproduction and validation

The accompanying `delivery-requirements.png`/`.svg` figure displays the thermal/acceleration trade and launch quantiles; the rendered figure was inspected.

Run `delivery_screen.py` and `validate_delivery.py` with the existing `/tmp/protium-python` dependencies. The 151 checks include independent sphere-capacitance/work, blackbody angular-integral, intercepted-gas-energy, uniform homologous-flow and optical-depth chord-integral tests, plus source kinetic-energy reconstruction, root checks and 128/512-cell input comparisons. They verify the formulas and numerical accounting. They do not validate the material or driver assumptions.

`cold_driver_probe.py` extends the same electron/positron EOS source to 100 K while preserving every existing hot-table node. Its fully ionized cold material is deliberately an unvalidated compression surrogate. The source EOS required increasing its documented degeneracy-search upper bracket from 10⁷ to 10¹⁰ for the coldest, densest nodes. The completed table preserves every original hot node exactly. Twenty-one independent source comparisons give maximum pressure and energy errors of 7.7 × 10⁻⁷, heat-capacity error 1.9 × 10⁻⁵, and density-derivative error 5.7 × 10⁻⁶; all six interpolation checks pass. The completed bounded runs, guard stops and numerical failures are recorded in [COLD-COMPRESSION.md](COLD-COMPRESSION.md). No physically delivered final gain is established. That report also contains the 30 m prescribed-hot scale diagnostic and the explicitly chosen fictional open-range arithmetic in `prototype-range.json`. `MANIFEST.json` records packet and source hashes.
