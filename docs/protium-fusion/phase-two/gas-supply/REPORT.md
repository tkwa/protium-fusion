# Fresh hydrogen and distributed pumping ✨

## Result

About **30% of fresh 300 K hydrogen molecules emitted diffusely from the primary titanium channel wall undergo proton-impact first ionization before reaching a wall again**. The independent deterministic integral gives 0.299798; the 131,072-flight Monte Carlo gives 0.302415 ± 0.001101 standard error. These results concern the saved 137.149 T, 1.6 MeV radial source state, not an experimentally operated reactor.

This changes the interpretation of pumping coatings. With identical diffuse thermal flights after each unsuccessful wall encounter, ionization competes with capture at the next wall. For sticking probability α and first-flight survival S, the eventual first-ionization probability is (1−S)/[1−S(1−α)]. An α of 0.01 leaves about 98% of emitted molecules ionized before pumping. Even a perfect next-wall pump cannot prevent the first flight. Preventing gas release and removing electrons are separate engineering requirements.

The median first-birth radius is approximately 5.1 mm, well inside the 20.73 mm wall radius. Independent chronological trajectory integration gives a birth-weighted electrostatic potential of **10.8744 MV** relative to the wall. Positive molecular ions and subsequent dissociation products must therefore enter the charge, energy, rotation and collection accounts alongside their electrons. The co-born positive and negative electrostatic energies initially cancel; 10.9 MeV per unit charge is a gross circuit scale, not an automatic net energy debit.

## Model and verification

The actual radial proton density is recovered from the saved two-species equilibrium. A one-dimensional noncentral-Maxwell integral averages the Rudd proton–H₂ ionization fit over the thermal velocity distribution and rigid rotation. The fitted cross section at 1 MeV is 3.7414 × 10⁻²¹ m². The volume-averaged ionization frequency is 44,057 s⁻¹. The original Cartesian Gauss–Hermite velocity quadrature was biased by up to 0.87%; independent review identified and removed that approximation from the final calculation.

Diffuse wall emission uses cosine-law directions and the flux-weighted Maxwell speed distribution. Along each cylindrical chord, integration of the local ionization frequency gives survival and first-birth location. Independent angular, speed and chord quadratures converge to the first-flight probability quoted above. The sampled mean chord agrees with the kinetic-theory result 2R. Mean wall-flight duration is approximately 23 μs. The source density integral recovers the saved proton inventory to 1.4 × 10⁻⁷ relative error.

A separate uniform Knudsen-diffusion illustration has D=2R〈v〉/3. With perfectly evacuated ends, uniform wall outgassing q and no distributed sink, mean pressure is qL²/(4R²〈v〉). For the 20.24 km active length and q=10⁻¹⁰ Pa m³ s⁻¹ m⁻², this gives 0.0134 Pa and a mean diffusive residence near 16 days. It demonstrates the severe end-conductance restriction; it is not the pressure of the operating plasma, which ionizes the gas. Gas-collision rarefaction and wall accommodation must be checked before applying that reference to a different gas state.

## Supply budget and limitations

The Ti primary-channel wall area is **2.532 × 10¹² m²**. Auxiliary processor walls are additional. At the example q=10⁻¹⁰ Pa m³ s⁻¹ m⁻², complete stripping of all fresh H₂ would create at most 0.0293 electrons per primary Ti reaction. Assigning a complete electrical cost of 0.1 MeV per electron would debit 0.00293 MeV per primary; q=10⁻⁸ makes this 0.293 MeV. The assumed cost is inherited from the earlier parameterized supply budget. This model does not establish that cost, a reactor outgassing rate, or the full-plant gas load.

Room-temperature measurements on baked metal vacuum chambers provide context for q, not qualifications for irradiated carbon channels. Nuclear radiation, molecular fragments, secondary electrons, charge exchange, titanium-ion impact, wall recycling and injection can change the supply. Getter capacity, reactivation, saturation and reflection law also matter. The conditional two-electron count following the modeled first-ionization event excludes competing chemistry; the unconditional fresh-H₂ ceiling is two electrons per molecule, before recycling or nuclear pair creation.

Using measured H₂ polarizability, the saved radial electric field changes the induced-dipole potential by only about 0.013 kT at 300 K. Ordinary neutral polarization is a small correction to this flight screen.

## Sources and reproduction

- Rudd et al. (1985), [Electron production in proton collisions: total cross sections](https://doi.org/10.1103/RevModPhys.57.965). The accessible author repository abstract confirms the fitted-data method; its PDF download returned 403 and was not retried.
- Padovani, Galli and Glassgold (2009), [Cosmic-ray ionization of molecular clouds](https://doi.org/10.1051/0004-6361/200911794), equations 5–6, supplies the fit and constants used here.
- Fedchak et al. (2021), [Outgassing rate comparison of seven geometrically similar vacuum chambers](https://doi.org/10.1116/6.0000657), and its [NIST manuscript](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=930674), provide measured metal-chamber context. Conversion: 1 Pa L s⁻¹ cm⁻² = 10 Pa m³ s⁻¹ m⁻².
- NIST, [Reference measurements of hydrogen’s dielectric permittivity](https://www.nist.gov/publications/reference-measurements-hydrogens-dielectric-permittivity), supplies the measured polarizability used in the independent check.

Run `gas_supply.py`, followed by `independent_gas_check.py`, with the first-phase NumPy/SciPy environment. The latter records input hashes and verifies the flight integral in independent coordinates. No full plasma or material qualification follows from these checks.
