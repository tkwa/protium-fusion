# Contact-to-hot source assembly audit ✨

The carbon-return stage now has a constructive **one-dimensional pressure-work connection** between its chosen postshock thermodynamic state and the prescribed hot inflow. Nitrogen-15 has a similar connection for 82.315% of its mass, while an outer 17.603% needs a separate expansion/contact sequence. Preparation radiation changes the tested nitrogen pulse's subsequent nuclear yield by less than 1%. The remaining substantial uncertainty is the delivery and joining of finite three-dimensional streams, including the outer nitrogen material and central seeds.

This audit concerns the source128 family A H/N=3, B1 H/N=2, B2 H/C=2 and C13 H/C=4. The H/N=1.5 nitrogen alternative has its own source study and is outside these dynamical tests. [Source128 packet](../source-family-carbon-median128/source128-stage-map.json), [summary](summary.json), [exact archived records](records/).

## The missing equation was pressure work

A parcel's contact entropy specifies its temperature at another density. It does not assign the forces, velocities, contact time or pressure work that carry it there. Individual parcels exchange mechanical energy with their neighbors, so local `hot internal energy + hot kinetic energy` cannot be treated as an advected scalar conserved from the launcher.

The spherical test uses fixed cell masses and moving shell boundaries. Each parcel follows its FEOS isentrope `P(ρ,s)`. The pressure force on a boundary and the internal-energy work use the same pressure:

`dU_i/dt = −P_i dV_i/dt`,

`m_node dv_node/dt = 4πr_node²(P_inside − P_outside)`.

These exchanges cancel in the total energy of the closed fluid. Reversing the hot velocities and following an outward solution tests a possible reversible fluid segment before the hot snapshot. The flow equations, acoustic speeds `v ± c`, and entropy changes at real shocks supply the relevant distinction between smooth evolution and additional shock heating. [Clawpack authors' Euler-equation treatment](https://www.clawpack.org/riemann_book/html/Euler.html).

The supplied hot profiles have actual outer radii **170–191 m**. Their 30 m `R0` parameter sets a reference mass and time normalization; it is not the outer radius of these supplied inflow snapshots. In all four snapshots, the forward flow is compressing and the reversed velocity increases outward. A purely ballistic reversal therefore has no shell crossing. Its missing compression work is nevertheless significant: approximately **4.5% of prepared energy for the matched B1 material and 6.9% for B2**. The pressure-driven calculation accounts for that exchange.

## B2: a resolved pressure-work connection

The selected CH₂ contact has density **3,761 kg/m³**, temperature **439,707 K**, pressure **8.70×10¹² Pa** and deposited normal-collision energy **3.451 GJ/kg**. The independent primary-table Hugoniot gives nearly the same density and pressure and a 9.6% higher temperature, establishing the present contact-EOS accuracy scale. [Primary EOS comparison](../b2-source/contact-primary-check.json).

A fully ionized FastEOS continuation at that contact would give **38.6% higher pressure** than FEOS. This audit therefore follows the partially ionized FEOS adiabat for the backward segment. At the actual hot snapshot, FEOS versus FastEOS changes total internal-plus-kinetic energy by only −0.114%; the contact and hot comparisons must keep their respective states explicit.

Starting from the exact source128 hot radii, temperatures and radial velocities, the reversed flow reaches the common contact state for **all of its 99.9176% nonseed mass**. The remaining 0.0824% is the source's existing seed reserve.

| Mass-weighted contact-event quantile | 10th | Median | 90th |
|---|---:|---:|---:|
| Time before hot snapshot, 128 cells | 123 μs | 236 μs | 249 μs |
| Contact radius, 128 cells | 385 m | 631 m | 695 m |
| Time before hot snapshot, nested 256 cells | 129 μs | 235 μs | 249 μs |
| Contact radius, nested 256 cells | 381 m | 630 m | 695 m |

These are quantiles of distinct parcel contact events, not the trajectory of one spherical shock front. The events occur at different times. Beyond each event, continuing that parcel as fluid is an auxiliary boundary construction; the incoming grains and oblique collision interfaces still need their own solution.

Over the 0.4 ms reversed continuation, internal energy transfers **2.80×10²⁴ J** into radial kinetic energy. The total-energy residual stays below **1.6×10⁻¹¹**. Shell widths remain positive and grow. Refining the EOS tables and halving the maximum timestep changes contact-time quantiles by at most 4×10⁻⁷ relative. Nested spatial refinement changes the median time by −0.31% and the 10th percentile by +4.44%, so these geometric scales should be quoted to two significant figures.

The 128-cell acoustic gradients remain expansive. The nested 256-cell profile develops brief compressive acoustic gradients around the piecewise input layers; the largest compressive jump is **0.029 times the local sound speed**. The sampled global compressive-gradient envelope integrates to **0.0594**, and no shell compression appears in this B2 continuation. These are finite-resolution diagnostics. They support a modest fluid correction and do not constitute a proof that every continuum reconstruction is shock-free. The input profile is represented by piecewise cell values; physically smoothing and assembling it remains a source-design task. [Characteristic refinement](characteristic-refinement.json).

On this fixed FEOS history, the κ=0.001 m²/kg radiation screen gives a maximum contact-to-hot equivalent temperature perturbation of **0.0109%**. Subsequent nuclear burning during preparation is negligible in the evaluated fixed-composition rate screen. Radiation therefore does not invalidate this B2 pressure-work connection within the current grey model.

## B1: a separate outer-material phase and a small measured radiation effect

The nitrogen contact has density **4,453 kg/m³**, temperature **12.842 MK** and normal-collision energy **123.94 GJ/kg**. Its hot profile has lower density than that contact in the outer **17.603% of total mass**, while its velocity field is presently compressing everywhere. Those parcels cannot reach the denser contact by monotonically decompressing backward along this whole-fluid branch.

A physical sequence can accommodate this through an earlier rarefaction and density turning point, or by expansion of shocked grains/contact patches before they form the later volume-averaged fluid. The local material density during grain collision and the eventual mean spherical density need not coincide. The existing constant-entropy map leaves that geometry unresolved; it supplies no physical reason to lose the outer material.

The other **82.315%** reaches the selected contact in the reversible calculation. Contact-event mass deciles are approximately **5.0, 23.5 and 29.6 μs** before the hot snapshot, at **80, 123 and 129 m** radius. The seed reserve remains separate. Energy closes to 4.3×10⁻¹⁴ under timestep/table refinement. Small seed-region acoustic changes remain visible.

Nitrogen's lower density gives radiation a larger role in shaping temperature. A first-order integral of absolute internal radiative transfers is a loose perturbation bound: heat can move both directions during the changing flow. It must not be counted as a permanent plant energy loss or used directly as a predicted temperature change.

A nonlinear comparison therefore evolves the **same whole-fluid precursor for 30 μs**, with and without grey radiation. All three trajectories use FastEOS, the same composition and initial prepared energy, no nuclear burning and zero artificial viscosity. They carry the actual pressure work and resulting density/velocity changes. The earlier FEOS-to-FastEOS handoff at this precursor changes modeled total energy by −0.1844%, without rescaling temperatures or velocities; the no-radiation control isolates that difference. [EOS bridge](B1-precursor-EOS-bridge.json).

| Preparation opacity | Temperature change versus no-radiation control, mass deciles | Radiation escaping preparation / prepared energy | Subsequent nuclear-output change versus control |
|---|---|---:|---:|
| 0.001 m²/kg | −1.97%, −0.108%, +7.28% | 0.00229% | −0.911% |
| 0.01 m²/kg | −1.11%, −0.0415%, +1.74% | 0.000328% | −0.0819% |

The nuclear-output column comes from independent full-network continuations of the actual delivered endpoints with identical hot-run settings. The no-radiation control itself has 0.615% lower nuclear output than the original source128 reference, separating the round-trip/EOS change from the effect of preparation radiation. [Delivered-profile burn audit](../assembly-burn/summary.json).

Those runs establish that the tested radiation-induced profile changes retain useful burn performance. They keep the complete isotope outputs and preparation-escaped energy available for whole-plant integration. Their auxiliary fluid precursor still leaves the outer-material contact sequence and three-dimensional input boundary unresolved. The preparation-opacity rows are separate model choices, not a probability interval.

## A and C13: no consequential preparation burning found

A bounded thermal screen starts from each exact source128 field and follows its reversed ballistic geometry with FEOS constant-entropy temperatures for 0.5 ms. It provides a transport estimate without evolving pressure forces or closing kinetic-plus-internal energy, so it is not a realized trajectory.

At κ=0.001 m²/kg, A's largest equivalent hot-temperature radiation perturbation is **0.0200%**. C13's nonseed mass deciles are below **0.025%**; one **0.0224%-mass cell** adjoining the seed reaches a 5.34% absolute perturbation bound, falling to 0.534% at κ=0.01. Net radiation escaping either preparation screen stays below 3×10⁻⁸ of prepared energy. This does not justify another full coupled preparation calculation on radiation grounds under the current assumptions.

Holding every parcel for a full millisecond at the largest actual hot-profile temperature and density produces initial-rate nuclear-energy fractions of **3.6×10⁻¹¹ for A and 5.4×10⁻⁹ for C13**. These deliberately hot/dense envelopes are conditional on the proposed prehistory remaining below those maxima. They give substantial margin for neglecting preparation burning along the proposed contact-to-hot thermodynamic paths. An unmodeled hotter seed-production process would have its own nuclear account.

C13 has a marked seed-to-mapped-material entropy/pressure transition. Its initial outward acoustic slope has a frozen-field convergence timescale of about 2.7 μs. That diagnostic concerns the seed's momentum and joining geometry; it is not a calculated shock time or evidence that the entire C13 pulse fails. The finite contact-patch and seed-joining calculation remains the appropriate next test. [A/C13 thermal and nuclear screens](other-stages-screen.json).

## What remains most consequential

The useful remaining source test is a finite three-dimensional contact patch with the chosen cold material, normal collision velocities, packet spacing, local contact density and adjoining pressure boundaries. It should determine the spatial density/entropy/velocity distribution and how it joins neighboring patches. The nitrogen outer-material sequence and each stage's seed interface belong in that model. Actual collision loci and timing must be obtained from that trajectory; the derived hundreds-of-meters and tens-to-hundreds-of-microseconds scales above are reference geometries rather than manufacturing tolerances.

A hidden additional heater is not required by the B2 pressure-work calculation: its energy exchange can be supplied by the already large incoming kinetic inventory. Likewise, the tested B1 radiation effect does not force a new power architecture. Both conclusions are conditional fluid results. Launcher, capsule, neutralization, final packet delivery and collective assembly still have to supply the boundary conditions together.

No subjective 80% interval is assigned to an actual construction tolerance or launch reliability: those parameters do not yet have a specified contact geometry or acceptance test. Numerical refinements, opacity scenarios and EOS comparisons are reported separately so they can later inform such an estimate without being mislabeled as probabilities.

## Numerical scope, verification and reproduction

The backward model evolves two variables per shell boundary, radius and velocity, at fixed parcel entropy. FEOS pressure is interpolated on 601-point log-density tables, checked against 1,201 points. Internal-energy changes integrate the same pressure through `du = P dρ/ρ²`; exact initial FEOS energies set each parcel's reference. This preserves pressure-work cancellation. The largest table/FEOS energy-difference mismatch is 6.4×10⁻⁷ on the original B2 tables, improving with refinement. DOP853 uses relative tolerance 10⁻⁹, maximum timestep 0.2 μs, then 3×10⁻¹⁰ and 0.1 μs for refinement. The tested density and temperature domains remain recorded in each result. Radiation, nuclear reactions and viscosity are omitted from this reversible model and assessed separately.

A numerical reverse replay returns the initial dimensionless state within 9×10⁻⁹ in the original B2 calculation. This checks reversibility of the implemented ODE, while physical shocks, radiative diffusion and grain breakup remain irreversible processes outside that test. The forward radiation comparisons use the frozen radial solver with 128 cells, tolerance 10⁻⁸, actual composition, and 30 μs physical duration. Energy residuals stay below 5×10⁻¹¹. The preparation nuclear screen uses the median-carbon network on fixed histories; cold tails are evaluated with a stated 7 MK floor inside the updated carbon table, providing a negligible-reactivity diagnostic rather than a precision low-temperature prediction.

`initial_screen.py`, `backward_isentrope.py`, `omission_screen.py`, `forward_radiation.py`, `other_stages_screen.py` and `characteristics.py` reproduce the respective layers. `package.py` verifies and archives the compact records. The preceding backward runner is preserved under `runner-versions/` with its content hash. [Manifest](manifest.json) records inputs and heavy raw traces retained in the persistent research output directory. Existing source profiles and earlier research snapshots were not modified.
