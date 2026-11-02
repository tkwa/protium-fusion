# Construction, replacement and competing energy sources ✨

At 1 mW/kg, a gigawatt requires a trillion kilograms of installed plant. Hydrogen abundance supplies a large fuel reserve, but construction energy, replacement rate and access to easier fuels decide whether such a plant is useful. `economic_comparison.py` evaluates those conditions for the two existing requirements allocations. Efficiencies, equipment masses and discount rates remain explicit scenarios.

## Energy invested in the plant

For net specific output `s` and construction energy `e` per installed kilogram, simple energy payback is `e/s`. At 1 mW/kg, construction energies of 1, 10 and 100 MJ/kg require about 32, 317 and 3169 years of full output. These figures describe the required manufacturing performance and service horizon; they do not assign a future material-production cost.

The conditional nonneutral plant exports 2.160 TW from 1.466e15 kg, or 1.474 mW/kg. Its 100 MJ/kg construction scenario repays in 2150 years. The stored magnetic energy is 1.261e21 J. Charging it at an assumed 95% efficiency adds 19.47 years of net output, bringing the combined energy payback to 2170 years. Recoverable magnetic energy at eventual shutdown can be treated as salvage; it still requires initial financing and charging capacity.

This reference already devotes 464.4 GW to replacement manufacture, 17.70% of its pre-replacement bus output. Its assumptions are replacement of 0.01% of installed mass per year and 100 MJ/kg of replacement material. That load is already deducted from the quoted 2.160 TW.

Local construction is valuable at these specific powers. Accelerating installed material to 10 km/s alone stores 50 MJ/kg of kinetic energy; at 100 km/s it stores 5 GJ/kg. Real transport adds propulsion and handling costs, while braking recovery can return some energy. Freight and installation work must be included in the chosen construction-energy figure or added separately.

The energy investment also bounds expansion. If every unit of net energy makes more identical plant, with construction continuously divisible and other resources unconstrained, the maximum exponential capacity-growth rate is `s/e`. At 100 MJ/kg it is 0.0465% per year, giving an ideal doubling time of 1491 years before initial field charging. Building an entire second unit in a single batch takes the 2150-year simple-payback interval instead.

To express opportunity cost without forecasting prices, hold the value of a joule fixed and discount future delivered energy at rate `r`. A plant of lifetime `L` returns

`s (1-exp(-r L))/r`

per installed kilogram, with consistent time units. Infinite life cannot repay construction if `r >= s/e`. For the 100 MJ/kg ion-plant scenario, 0.01% annual discount lengthens payback to 2421 years; 0.03% gives 3453 years. At 0.05% it never repays in this energy-value model. These rates are scenario choices. A civilization can accept a long energy repayment for strategic fuel access or resilience, but those benefits should be stated separately from cheap electricity.

## Replacement closes a gap in the opaque allocation

The original [opaque budget](PLANT-BUDGET.md) includes fuel preparation, support dissipation and an auxiliary reserve, but does not explicitly allocate replacement manufacturing and spare stock. Its quoted 1 mW/kg therefore requires this additional accounting. At 100 MJ/kg and replacement of 0.01% per year, manufacture alone consumes 0.3169 mW per installed kilogram. The original auxiliary reserve is smaller than this load.

The new calculation preserves that reserve and adds replacement separately. It solves total mass again, including a factory at an assigned 1 kg/W of manufacturing power, ten years of spares, the increased generator/driver ratings, radiator area and the self-gravitating support load. Fixed collector and pulse-driver assumptions remain the same. The zero-replacement limit exactly reproduces the old allocation.

| Quantity | Original allocation | With 0.01% annual replacement |
| --- | ---: | ---: |
| Installed mass | 6.778e25 kg | 6.789e25 kg |
| Jupiter masses | 0.035709 | 0.035765 |
| Net useful output | 6.778e22 W | 6.789e22 W |
| Net specific output | 1 mW/kg | 1 mW/kg |
| Replacement manufacturing | Unallocated explicitly | 2.151e22 W |
| Pulse throughput relative to original | 1 | 1.3175 |
| Mean interval across all pulse types | 7.194 days | 5.460 days |
| Radiator area / collector surface | 21.71% | 29.64% |

Thus this replacement scenario can be accommodated arithmetically by about 32% more pulse throughput and a modest increase in installed mass. The driver, collector recovery, thermal cycling and processing must actually sustain that cadence. The extra power is obtained by consuming more fuel, not by recovering the manufacturing energy twice.

At 0.1% annual replacement the arithmetic allocation requires 4.22 times the original throughput and radiator area equal to 102.1% of the collector surface. That case needs a different radiator arrangement and associated supports; the current spherical layout does not establish it. A stable mass fixed point alone does not settle this geometric requirement.

## Deuterium extraction

The [Galileo reanalysis by Mahaffy et al. (1998)](https://websites.umich.edu/~atreya/Articles/1998_Galileo_Probe_Measurements.pdf) measured atmospheric D/H of `(2.6 ± 0.7)e-5` by number. The comparison uses that historical measurement and its range as scenarios, plus a strongly depleted case. It does not assume one universal abundance. All feed quantities below refer to the hydrogen component; helium and other constituents increase the amount of total gas processed.

The conventional catalyzed D–D sequence, including subsequent burning of tritium and helium-3, has net reaction

`6 D -> 2 He4 + 2 p + 2 n`.

Its mass-defect energy is 43.24 MeV, or 3.453e14 J/kg of deuterium. Capturing the two neutrons on the product protons returns two deuterons and closes `4 D -> 2 He4`, with a higher available energy per net deuteron consumed. The model carries these two fuel cycles separately. It does not assign an existing commercial reactor their assumed conversion performance.

At D/H = 2.6e-5 and an assumed 30% net conversion of the first cycle's nuclear energy, each kilogram of hydrogen feed contains deuterium capable of delivering 5.382 GJ. Matching the ion-plant reference's 2.160 TW would consume 20.85 grams of deuterium per second and require processing 401 kg/s of hydrogen at complete isotope recovery. The protium reference instead consumes 17.86 grams of ordinary hydrogen per second and delivers 121.0 TJ per kilogram consumed. It accesses about 22,500 times as much useful energy from the same initial hydrogen reservoir in this comparison.

The ideal atomic-mixture separation benchmark at 300 K is

`W/N = -kT [x ln(x) + (1-x) ln(1-x)]`.

It gives 744 J/kg of hydrogen feed, or 14.31 MJ/kg of extracted D in the central abundance scenario: 1.38e-7 of the assigned deuterium output. Molecular chemistry, finite recovery, cooling, pumping and real separation equipment require additional costs. This benchmark shows why rarity alone is an inadequate argument for replacing deuterium reactors with a much heavier protium plant.

Fuel location matters. Using [IAU nominal Jupiter parameters](https://arxiv.org/pdf/1510.07674), the surface-to-infinity gravitational work is approximately 1.772 GJ/kg, ignoring rotation and atmosphere. Lifting the entire hydrogen feed would consume 32.9% of its assigned deuterium output. Separating locally and lifting only the D reduces the ideal lifting term to 1.71e-5 of the D output. Neither estimate designs an atmospheric processing facility.

A reservoir containing 1e27 kg of hydrogen at the central D/H would sustain the 2.160 TW reference load from its deuterium for about 7.9e16 years under these assumptions. At a useful power equal to the nominal solar luminosity, the same reserve lasts only 445 years. Aggregate demand changes the scarcity argument dramatically. Stellar hydrogen is also a distinct feedstock: early stellar evolution destroys deuterium, so unprocessed planetary gas and starlifted gas cannot automatically be assigned the same isotope supply.

## Solar collection

For a local collector of complete effective areal mass `Sigma`, net conversion efficiency `eta`, and distance `r` in AU from a Sun-like star,

`s_solar = eta * 1361/(Sigma r²) W/kg`.

The irradiance normalization comes from [IAU Resolution B3](https://arxiv.org/pdf/1510.07674). Here the chosen areal mass must cover the hardware needed to deliver local useful electricity. Remote transmission, storage, transport and additional support require their own additions.

At 30% net conversion, collectors of 1 and 10 kg/m² reach 1 mW/kg at approximately 639 and 202 AU. They match the conditional 1.474 mW/kg ion plant at 526 and 166 AU. These are mass comparisons under specified equipment assumptions, rather than a prediction of an actual collector's cost or lifetime. They make distant, dark or isolated locations a more plausible first market for a heavy protium plant than ordinary locations near luminous stars.

## Accretion, fuel export and coproducts

[Thorne's thin-disc calculation](https://articles.adsabs.harvard.edu/pdf/1974ApJ...191..507T), table 1, gives outgoing-radiation efficiencies of 5.70% for a Schwarzschild hole and 30.2% for its canonical isotropic-emission rotating case. With an assumed 50% capture and useful conversion of that radiation, their fuel energy yields are about 21 and 112 times the present ion-plant reference's. These are model-based fuel yields; a comparison of complete installations must also count the seed, acquisition, feeding system, converter and heat rejection. Natural stellar black holes exceed the user's one-Jupiter-mass reactor limit, while a small engineered seed introduces a separate construction problem.

Protium infrastructure could serve as a stationary fuel refinery, exporting deuterium or other intermediate fuels to lighter end-use reactors. That role requires a new product and energy ledger; exported fuel cannot also be counted as electricity generated on site. The conditional ion reference produces about 17.73 grams of helium per second. Helium can have local uses, but its value depends on the alternative supply and transport costs. No arbitrary coproduct revenue is used to make the reference profitable.

The economic case is therefore conditional on long-lived construction, low replacement demand, exhausted or inaccessible easier fuels, weak local starlight, or a valuable need for independent fuel production. The calculations do not supply a scientific breakthrough date. They do identify the kind of civilizational infrastructure and time horizon that a credible future history must explain.

## Validation and reproduction

Run `economic_comparison.py`, `validate_economic_comparison.py` and `plot_economic_comparison.py` in this directory. All 46 checks pass, including the original opaque limit, independent quadratic and numerical mass closures, electrical/thermal conservation, isotope mass-defect closure, integration of the reversible chemical-potential separation work, and discounted-output quadrature. The figure was visually inspected. The manifest records all input and output hashes. These checks verify calculations and bookkeeping; the source plant's engineering assumptions remain explicit in its own reports.
