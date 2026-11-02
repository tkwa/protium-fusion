# ✨ Pulse mixture comparison methods

The mixture comparison follows a complete pulse and then closes the isotope, hydrogen, energy and installed-mass accounts. The quantity being optimized is exported electrical power divided by total plant mass. A high initial reaction rate, a high isolated pulse gain and a favorable plant operating point are different results.

## What the comparisons hold fixed

Let x be the number of hydrogen nuclei supplied per catalyst nucleus. The reference values are x=1 for nitrogen-14, nitrogen-15 and carbon-13, and x=2 for carbon-12. Every new pulse records its actual initial abundance of every isotope. A stage name supplies no implicit proton count to the material ledger.

The main comparison holds total pulse mass, the enclosed-mass density profile, initial temperature profile and physical end time fixed. Radial velocities are scaled together to preserve total prepared kinetic plus internal energy. A separate fixed-flow comparison preserves the velocity profile as well; its prepared energy changes with composition. The two families have separate labels and denominators.

The principal pulses contain about 1.131×10¹³ kg. The nitrogen-15 pulse contains 1.131×10¹¹ kg. Their physical integrations last approximately 63–82 microseconds. The nominal pulse radius parameter is 30 m; the saved spatial inflow determines actual cell positions and densities. Later source comparisons change the initial temperature and velocity profiles at fixed total prepared energy and retain their own explicit provenance.

An initial comparison based only on a nominal temperature, density or sound-crossing-time multiple would change several constraints at once. The runner instead fixes the physical duration against the corresponding reference and records the reconstructed energy, pressure and composition.

## What an initial-rate estimate can suggest

For a fully ionized, homogeneous mixture at fixed mass density ρ and temperature, a catalyst of mass number A and charge Z has

`n_cat = ρ / [m_u (A+x)]`, and `n_p = x n_cat`.

Its initial proton-capture rate per volume is proportional to `x/(A+x)²`, with a maximum at x=A if the reaction coefficient remains unchanged. The classical electron–ion bremsstrahlung factor is proportional to `(Z+x)(Z²+x)/(A+x)²`. Dividing the two gives `x/[(Z+x)(Z²+x)]`, whose maximum is at x=Z^(3/2). These are analytic composition factors under fixed thermodynamic conditions. Electron–electron emission, degeneracy, relativistic corrections and reaction-network evolution are absent from this estimate.

A different constraint gives a different suggestion. At fixed ideal-gas pressure and temperature, the initial capture rate per volume is proportional to `x/(Z+1+2x)²`, with a maximum at x=(Z+1)/2. Radiation and degeneracy pressure change that relation.

Hydrogen enrichment therefore deserves testing, but none of these maxima predicts pulse gain. At fixed prepared energy, changed particle numbers alter pressure, heat capacity, compression and ignition. Hydrogen also displaces catalyst in a fixed pulse mass. A composition can improve an initial rate factor while preventing the main pulse from reaching its burning state. In an opaque pulse, bremsstrahlung transfers energy into trapped radiation; transport and expansion determine how much becomes an actual cooling loss.

## Pulse model and numerical resolution

The comparisons use the preserved spherical Lagrangian hydrodynamics solver with its compression-capped artificial viscosity, composition-dependent electron/positron equation of state, ideal ions and equilibrium radiation energy. The network follows 27 isotopes and 194 reaction rates. Radiation uses grey flux-limited diffusion with a spatially constant opacity. The principal comparison uses 0.001 m²/kg, with distinct 0.01 and 0.03 m²/kg sensitivities.

The systematic sweep begins at 64 radial zones. Refinements at 128 and 256 zones preserve the initial density and temperature fields as functions of enclosed mass and renormalize boundary velocities to preserve initial kinetic energy. Temperature is piecewise constant within each original cell. The source-family refinements likewise retain the saved source prescription instead of retuning it on the finer mesh.

The ordinary time-integrator tolerances are relative 10⁻⁷ and absolute 10⁻¹². Four hundred saved states support diagnostic quadratures; the integration uses its own adaptive timesteps. Endpoint isotope and energy results use the completed integration, rather than an integral over the sparse saved diagnostics. The 2.8 GK stop is a numerical/EOS-domain guard. Any case reaching it is incomplete and is excluded from operating-point comparisons.

Energy conservation and mesh convergence answer different questions. A small energy residual verifies the numerical account while spatially unresolved compression can still change yield. The updated carbon-12 return pulse at x=2 gives gains 10.2719, 9.9176 and 9.8671 at 64, 128 and 256 zones. Its last change is −0.51%. Other stages have their own mesh trends. No common error percentage or Richardson extrapolation is assigned from this single sequence.

The archived reaction-rate family and updated carbon evaluations remain separate. The updated calculation replaces the carbon-12 and carbon-13 proton-capture curves and applies the same multiplicative change to each paired reverse rate, preserving the original detailed-balance ratio. The [nuclear review](nuclear-review.md) identifies the primary evaluations and interpolation. Published carbon uncertainty curves, alternative nitrogen-15 evaluations and independent carbon-curve corners are explicit scenarios. Their output spread is not automatically a joint 80% credible interval.

## Generalized material balance

Each stage supplies an input vector and a completed output vector normalized to one initial circulating catalyst nucleus. The output is propagated through its cold radioactive decays. The four circulating species are carbon-12, carbon-13, nitrogen-14 and nitrogen-15.

If F contains the actual stage input vectors, M contains their cold output vectors and a gives stage attempts, their stationary balance is

`(F − r_C M) a = f_C`.

Here f_C is the fresh-carbon isotope vector and r_C is eventual return to the circulating catalyst pool. The stage attempts are subsequently normalized to one net helium nucleus. Hydrogen inputs use each stage's actual x and its actual unburned output. Fresh hydrogen equals gross hydrogen input minus returned hydrogen. Nuclear heat, delayed heat, neutrino energy, prepared energy, radioactive inventories and pulse frequencies all use the same stage weights.

The primary closed-sphere comparison sets eventual material return to one. First-pass extraction inefficiency, a slow ballistic tail, implanted carbon, a retained processing reject and a retired tile are inventory or rework. Actual escape, deliberate export and conversion of a circulating isotope into a different retained nuclear product require separate entries. Carbon converted into oxygen-16 remains inside the plant but leaves the four-isotope catalytic loop; fresh carbon can replace it without asserting that matter escaped the sphere.

The 0.007 fresh C:H constraint is a mass ratio of net external feed. The 0.0056 fuel and 0.0014 structural split is a proposed allocation, with no requirement to consume either allowance. Nuclear conversion, actual structural exports and construction stock must be counted separately. Initial construction carbon is a capital resource, not a fraction of the first year's hydrogen burn.

## Coupled plant comparison

For every candidate quartet, the plant account rebuilds pulse cadence, source peak rating, two-pulse energy storage, average conditioning, circulating stocks and feed buffers. The receiver reduction uses the actual initial composition and endpoint velocity and energy profiles. A reference main arrival window includes 99.9% of pulse mass. The remaining tail is retained and routed separately; this timing convention implies no material loss.

The common comparison uses a 1.1×10⁹ m sphere, a 200 K radiator, 60% thermal conversion, 60% electrical preparation efficiency and explicit mass allowances for source hardware, receiver, insulation, services, heat exchangers, radiators, inventory, support and reserve. The allowed main-window occupancy and radiator-area fraction are dispatch and design choices. For each dispatched family, a separate chronological schedule checks that its main windows can actually be placed.

The 5,000 kg/m² receiver provides a common mixture-comparison denominator. Lighter receivers form a separate design study. Their heat capture is allowed to vary with thickness: 99% is an assumption to compare, not a minimum imposed by physics. Heat missed by the hot receiver is either assigned to escape beyond all equipment or charged to colder equipment and its radiator. The latter prevents a thin wall from gaining apparent performance by omitting downstream heating.

Named pumping, processing, cryocooling and surface-refurbishment loads are compared with retained service allowances before adding another debit. Excess hydrogen carries recovery throughput and equipment even when its separation energy is small. Surface material returned for reworking retains its carbon inventory. The [economic account](../component-design/economics.md) separates thermodynamic upgrading work, manufacturing efficiency, construction throughput and capital resources.

## Interpretation

The coarse analytic factors suggest directions to investigate. The spherical calculations establish conditional pulse yields and material flows. The component studies specify apparatus candidates and acceptance conditions. A complete engineered operating point requires the source, collection, heat transport, structural control and processing hardware to meet those conditions together.

The major remaining source question is whether a spatially distributed cold-grain array can assemble the density, momentum and entropy fields assumed by the successful pulses. The contact-state calculations constrain that task and expose failed prescriptions. They do not yet evolve the entire three-dimensional assembly. Radiation spectra, nonlocal secondary-particle transport and the actual hot/cold collection geometry likewise determine the useful capture fraction of a thinner receiver. These are consequential physical uncertainties; arbitrary per-shot replacement, immediate global clearing or perfect first-pass extraction are unnecessary substitutes for them.
