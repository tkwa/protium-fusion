# Conditional plant allocation ✨

29 October 2026. This calculation tests whether a set of explicit engineering requirements can fit below one Jupiter mass at Thomas's preferred 1 mW/kg. It does not establish that those requirements can all be built together. The driver and processing masses are allowances, not outputs of detailed component designs.

A later [economic and replacement ledger](ECONOMICS.md) adds explicit replacement manufacture and spares. At replacement of 0.01%/year and 100 MJ/kg, it retains 1 mW/kg with 31.75% more pulse throughput and 0.035765 Jupiter masses. The original allocation below remains the zero-explicit-replacement reference.

## Closed nuclear and material cycle

`cycle_ledger.py` uses spatial convergent-flow results for N14+p, N15+p, C12+2p, C13+p and helium-to-carbon manufacture. Unburned catalyst returns to its appropriate stage. N13, O14 and O15 decay after expansion; their neutrino energy is deducted. Oxygen and heavier leakage becomes ash, with replacement C12 manufactured from helium. Unburned helium in the factory is recycled. The stationary isotope balance conserves baryons and agrees with an independent atomic-mass-defect energy balance.

With A at 128 cells, B1/B2 at 128 cells and auxiliary stages at 64 cells, the ideal-recovery ledger gives approximately 24.68 MeV of available nuclear heat and 1.78 MeV of prepared fuel energy per net helium atom. Its ordinary-hydrogen consumption is about 4.002 atoms per helium atom, including the small fraction becoming replacement catalyst and heavier ash. Those numerical results remain subject to the pulse-model limitations and incomplete resolution convergence, especially of rare channels and surface radiation.

A separate component ledger distinguishes final bulk kinetic energy, remaining hot internal energy and cold beta-decay heat. High direct-conversion efficiencies must apply only to the energy component and geometry that can realize them. Applying an 85–90% efficiency indiscriminately to all decay photons would overstate the support for such an estimate.

## Base requirements allocation

The base case uses the smaller 30 m helium factory to reduce maximum pulse size. It allows a catalyst collection loss of 1e-4 per processing stage, captures 99% of otherwise available energy, assumes 40% thermal-to-electric conversion and 30% electrical-to-prepared-fuel efficiency, and reserves 10% of the remaining cycle electricity for other auxiliaries. Under these assumptions the pre-auxiliary margin is about 3.56 MeV per net helium atom.

The collector has radius 2e10 m (20 million km), an absorber areal mass of 1e4 kg/m2, and active orbital support. The total radius is much greater than Jupiter's; the user's comparison is by mass.

| Component | Allocated mass |
| --- | ---: |
| Absorber | 5.03e25 kg |
| Driver facility | 1.00e25 kg |
| Processing and inventories | 5.03e24 kg |
| Pulse electrical storage | 1.36e24 kg |
| Insulation | 5.03e23 kg |
| Support guides, controls and redundancy | 5.03e23 kg |
| Circulating support streams | 7.67e22 kg |
| Conversion, conditioning and radiators | 4.85e22 kg |
| **Total** | **6.78e25 kg, approximately 0.036 Jupiter masses** |

The corresponding net useful output is 6.78e22 W. Gross electric generation is about 2.26e23 W; the driver consumes about 1.50e23 W. Support losses consume about 3.05e20 W, and the other auxiliary reserve is about 7.57e21 W. Waste heat is about 4.51e23 W. A 300 K radiator with emissivity 0.9 and areal mass 10 kg/m2 requires about 1.09e21 m2, approximately 22% of the collector's surface area.

The five kinds of pulses recur individually every roughly 32–40 days in this allocation, giving about one pulse per week in aggregate. Hydrogen consumption is approximately 8.8e8 kg/s. These are steady-average rates, not claims about a demonstrated schedule or processing throughput.

## What is counted, and what is still assumed

The calculation solves total mass self-consistently. Power-dependent hardware is rated against gross generation and driver demand. Radiators reject waste heat after recirculation and support losses. Pulse storage is charged at an installed specific energy of 0.1 MJ/kg, based on the largest electrical pulse requirement, approximately 1.36e29 J. Circulating support mass and its dissipation depend on the final total mass rather than on absorber mass alone.

The support streams move at 10 km/s and lose an assumed 0.1% of their kinetic energy per orbit. The guide mass, geometry and loss figure need engineering substantiation. The collector is assumed isolated from strong external tides. Its insulation, surface recovery, precision timing, local stress, fault tolerance and manufacturing process are not designed.

The whole 1e25 kg driver allocation is explicit, but no current calculation shows how to use it to create the modeled low-entropy inflow at 30% coupling. This is the largest gap between the numerical pulse and an actual machine. Multidimensional implosion stability is another unresolved gap.

## Sensitivities

- Increasing the fixed driver allowance to 1e26 kg still gives an arithmetic mass closure below 0.1 Jupiter masses, with more frequent pulses. This does not validate the driver.
- A 20% driver leaves very little cycle margin in the conservative base recovery/factory case. Its high throughput requires much more radiator area, and small model errors can reverse the energy balance.
- A catalyst collection loss of 1e-3 per stage makes this small-factory thermal case energy-negative. The larger factory or higher recovery efficiencies change that result; good collection is a real requirement.
- Increasing prepared fuel energy by 50% leaves a small positive margin; doubling it makes the cycle energy-negative at 30% driver efficiency.
- Increasing support dissipation tenfold raises its load to roughly 4.5% of net useful output in the base mass allocation.

These sensitivities identify engineering requirements. They are not a probability estimate, an optimized design, or evidence for a particular completion date.

## Economic interpretation

At 1 mW/kg, an embodied construction energy of 1 MJ/kg requires about 32 years of full-output energy repayment; 100 MJ/kg requires about 3,200 years, before maintenance and financing. A plant assembled from several Earth masses of manufactured absorber is therefore a very different economic proposition from one using mostly existing natural gas or rubble as ballast. The budget should be compared with lighter higher-power architectures, solar collection and imported or manufactured black-hole converters before selecting the article's future role.
