# Converging-pulse instability study ✨

Started 30 October 2026 (Pacific time). Thomas authorized simulations using fewer
than eight physical cores for five hours and clarified that the question is
**growth of disturbances during the pulse**. The run launcher confines workers
to physical cores 0–6 and enforces an end time of 31 October 2026, 08:22:56 UTC
(01:22:56 Pacific). Runs and intermediate fields are in
`/var/tmp/tkwa-chatgpt/protium-stability/`.

## Question and diagnostic

Thomas selected **at least 90% of the corresponding smooth pulse's nuclear
output** as the acceptance criterion. The requested end products are central
estimates and **subjective 80% credible intervals** for (1) the maximum
mass-weighted standard deviation of arrival-time errors relative to the intended
schedule, and (2) the minimum number of independently aimed, shaped packets.
Each packet has one bulk flight direction; its outline, density profile and
longitudinal pulse can be shaped. Independently tailored directions within one
packet are excluded by Thomas's clarification. Timing and finite-packet effects
share the same 10% loss budget when combined.

Long launcher distances make arrival-time and focal-position errors the useful
delivery coordinates. At distance D and speed v, a small independent speed error
contributes sigma_t approximately D*sigma_v/v^2. This relation is a conversion
of an error budget; correlated launch-time, path-length and speed corrections
must be combined with their covariances. At v=2.4e6 m/s, one microsecond
corresponds to 5.76 m/s of uncorrected speed error over 1e6 m, or 5.24 mm/s
over the existing 1.1e9 m launcher-radius reference. The longitudinal position
error at delivery is 2.4 m in either case. A common time shift of every packet
does not spoil spherical assembly.

Does a small, specified departure from a spherical incoming flow grow enough
during convergence and burning to reduce the nuclear output? Compare every
perturbed run with an unperturbed run using the same resolution, integration,
equation of state and nuclear source. Follow angular density variation,
distortion of enclosed-mass surfaces, transverse kinetic energy, and burning.
Check small-seed scaling and grid refinement before interpreting amplification.

The background is time dependent. Compression can amplify a disturbance
geometrically; shock interactions and shear can amplify it dynamically. A
growing absolute density difference by itself does not identify an exponential
instability or its mechanism. Diagnostics normalize differences to the evolving
spherical background and keep their timing relative to burning.

## First model

- **Stage A:** the existing 64-shell prescribed hot nitrogen/proton inflow,
  reference radius 30 m, actual initial radius 170.855 m, mass 1.13097e13 kg.
  Exact input identity is recorded in `prepare.py`'s manifest.
- **Hydrodynamics:** Athena++, commit
  `053bc738d5c85dac717326b4338b34c7db4e715d`, spherical polar coordinates,
  axisymmetry, second-order piecewise-linear reconstruction and RK2. Early
  velocity and packet-sector screens use HLLC; the main timing and Cartesian
  3D branches use HLLE. Every retained-output comparison has a matched control.
- **Thermodynamics:** the existing electron/positron EOS, analytic ion and LTE
  radiation contributions. Initially hold ion number per baryon at 2/15;
  electron number per baryon remains 8/15 through the selected capture.
  Radiation is advected and supplies pressure; photon diffusion is omitted.
- **Burn:** N14(p,gamma)O15 and its reverse using the existing ReacLib
  coefficients and nuclear masses. O15 mass fraction is advected; local nuclear
  heat is coupled to total energy at every RK stage. A second scalar tracks
  original fuel separately from a low-density atmosphere.
- **Numerics:** source time steps bound fractional reaction and heating changes.
  The conserved combination `total energy - Q * product mass` independently
  audits the heat/source coupling. Density/energy floor effects and outflow
  appear in the total conservation residual; completed results require both a
  normal exit and reaching the requested physical end time.

The fixed ion number and trapped-radiation approximations are explicit model
choices. Independent radial comparisons test their effect on the smooth pulse.
The eventual conditional tolerance curve applies to this supplied hot flow and
these closures. The Cartesian cross-check below resolves three-dimensional packet edges.
Discrete cold-parcel assembly remains an extrapolation.

## Perturbations

The initial velocity is

`v_r(r,theta) = v_0(r) [1 + a F_l(theta)] / sqrt(1+a^2)`.

`F_l` is the angular-cell average of a Legendre polynomial, adjusted to zero
spherical mean and unit spherical RMS on the actual grid. Thus `a` is the RMS
fractional modulation before the common energy correction. Every velocity case
has exactly the same initial mass, internal energy and integrated kinetic energy
as its matched spherical control. The local peak error exceeds the RMS error;
the angular degree is part of every result. Even degrees avoid a net translation.

The initial radial profile is integrated over each cell, including the exact
integrals of its piecewise-linear velocity and squared velocity. Its continuum
kinetic quadrature differs slightly from the original staggered 64-zone sum;
the model records the actual initial energy rather than concealing that change.

The initial timing screen applies a ballistic phase shift `r_new = r_old + v(r_old)*delta_t`
with angular RMS `sigma_t`. Integrating every cell in the original mass
coordinate preserves mass, specific parcel internal energy and kinetic energy.
The hypothetical ballistic crossing time of any fixed radius changes by
`-delta_t`; the signs therefore have no effect on the RMS convention. Subsequent
pressure forces and burning are evolved normally. This prescribes a displaced
prepared hot flow, so the uncertainty assessment must include the thermal
history of actual packet assembly. Initial EOS range checks rejected large
polar peaks in high-order Legendre seeds. Bounded cosine patterns, also adjusted
to zero mean and unit angular RMS, avoid these artificial peaks in the main
timing sweep. Rejected initial states remain in the run archive.

The evolved-reference timing study uses a stronger construction: a 1,536-cell spherical
reference is saved every 0.075 microseconds from 0 to 22.5 microseconds. At the
seed time of 11 microseconds, each angular sector receives the conservative
state from reference time `11 microseconds + delta_t`. This supplies the physical
compression, heating and velocity evolution associated with the timing error.
Conservative radial integration and linear time interpolation preserve each
sector's mass and energy. Zero-jitter resumption is compared with the original
smooth run. The principal limitation is that angular coupling before the seed
time is omitted; the earlier ballistic construction supplies an initialization
sensitivity check. Independent uncoupled sectors would simply follow the same
reference trajectory at shifted times, which helps distinguish angular fluid
interactions from a changed radial compression state.

The earlier-seeded branch extends that smooth reference to -22.5 microseconds.
It reverses the original velocities, evolves a nonreacting expansion, and then
reverses the time labels and momenta of its saved states. Forward replay from
that prehistory changes final nuclear output by -0.0313% relative to a matched
radial-spacing HLLC reference. Sharp initial interfaces show larger volume-L1
differences: 2.19% in density and 3.27% in energy at a saved time 0.0043 microseconds
after the original start. These errors are retained in the audit. The main 3D
timing calculations use the extended trajectory to impose phase errors at the
prepared starting time, giving angular disturbances the full subsequent pulse
in which to evolve. This prehistory represents a smooth fluid preparation.

The packet screen partitions a meridional slice into contiguous shaped sectors,
each moving along its central direction. In axisymmetry the sectors represent
rings. Their radial and transverse direction moments are integrated over each
cell, with exactly the same total kinetic energy as the smooth flow. Matching
their RMS angular-direction error to a circular equal-area patch gives an
**equivalent** three-dimensional packet count. For N small-angle patches,
`<sin^2(angle)> = 2/N - 4/(3*N^2)`. `packet_geometry.py` checks this scale against
actual spherical Voronoi tilings. The conversion is geometric; the ring
hydrodynamics require a substantial three-dimensional model uncertainty.

## Verification and provenance

`prepare.py` exports the EOS and nuclear coefficients. `validate_adapter.py`
compiles the actual C++ interpolation header and compares it with the source EOS
at 10,000 sampled states; it also checks exported forward and reverse rates.
`reference_1d.py` runs the corrected Lagrangian solver with the same explicit
thermodynamic simplifications for an independent numerical comparison.

Controlled low-order velocity seeds check finite-time disturbance growth.
At 27 microseconds, RMS seeds of 0.1% and 1% produce angular
density variations, averaged over radius with fuel-mass weights, of 0.445% and 4.34%. The amplification reaches roughly seven
near 30 microseconds and then decreases during expansion. The corresponding
half-mass-surface distortions remain below the seed amplitude. These diagnostics
follow the changing spherical background; timing errors are the primary delivery
tolerance used in the engineering estimate. [Growth figure](figures/disturbance-growth.png).

The matched late-seed order-8 timing result at 5.5 microseconds retains
89.90% at 384 radial by 64 angular cells, 90.02% after doubling the angular
resolution, and 90.18% after doubling the radial resolution. Absolute smooth
burn fractions still change with refinement: approximately 0.588 at 384 radial
cells and 0.599 at 768 for HLLC, tending toward 0.608 at 4,096 cells. The close
agreement of retained-output ratios is therefore more informative for this
threshold than agreement of unrefined absolute gains.

At 256 cells the independent radial calculation changes output by about 0.47%
when ion number is allowed to evolve with reaction progress. Adding its photon
diffusion closure changes output by approximately -0.0014%. These checks
constrain the smooth-flow approximations; short angular wavelengths can have
different transport sensitivity.

The face reconstruction is limited toward the cell center when an interpolated
state falls below the EOS cold curve. This leaves conserved cell averages
unchanged. Corrections are counted, including a separate counter for denser
cells in later builds. In the completed 48-base-grid 3D packet runs, the maximum
energy-ledger residual is below 2.4e-6 of initial energy, and less than 0.031% of
final output is produced in the last two microseconds. Positivity enforcement
and numerical diffusion can still suppress unresolved shear and mixing. The
mesh checks constrain resolved output; they do not establish converged
high-Reynolds-number turbulence.

The spherical Noh shock test checks geometric compression and shock heating.
For gamma=5/3 its downstream density is 64 and its shock radius is t/3. The
128-, 256- and 512-zone tests gave relative volume-weighted density errors 4.19%,
2.41% and 1.37% over 0.06 < r < 0.5 at t=0.6, comparing with exact cell averages
through the shock. The downstream density approaches 64 (58.69, 61.42, 62.63
over 0.06 < r < 0.15); angular differences remain at float-output roundoff.
`validate_noh.py` reproduces these comparisons. Earlier point-sampled analytic
values gave a nonmonotonic L1 sequence because of shock/grid alignment; the
cell-average comparison is the appropriate finite-volume diagnostic.
[Analytic benchmark and verification discussion](https://cococubed.com/research_pages/noh.shtml).

Setup failures are retained: the first 2D attempt selected Athena's 3D `polar`
boundary rather than its `polar_wedge` boundary; an ideal-gas test called an
unavailable general-EOS getter; and the original HLLC general-EOS branch evaluated
an unused negative PVRS middle pressure during a rarefaction. The HLLC change
guards those unused calls on the existing shock/rarefaction condition. The
physical shock flux formula is unchanged. The source patch is retained with the
study. Initial small-domain runs also identified energy leaving the outer
boundary after about 38 microseconds; subsequent longer runs use a larger domain.

## Focused literature findings

| Evidence | Consequence for this study |
| --- | --- |
| [Yates et al. (2020)](https://arxiv.org/html/2002.03006v3) measured merging-jet density structure and improved symmetry after balancing the guns. Their experiments also show that ionization and impurities can change collisionality during merger. | Specify incoming nonuniformity and distinguish a collisional fluid from interpenetrating streams before modeling the driver. Laboratory mass-balance percentages do not supply a tolerance for this pulse. |
| [Cassibry et al. (2012)](https://arxiv.org/abs/1201.1879) report that a particular 30-jet hydrodynamic calculation became smoother while imploding into vacuum, with little interparticle mixing before peak compression. The abstract was accessible; the full PDF retrieval timed out. | Convergence alone does not establish catastrophic growth. The actual flow and timing matter. This paper supplies no yield prediction for the burning nitrogen/proton inflow. |
| [Samulyak, Zhang and Kim (2015)](https://arxiv.org/html/1512.01518v1) simulated a 90-jet argon liner compressing a separate plasma target. Oblique merger shocks produced uneven pressure, growing bubbles and spikes, and target breakup; their uniform comparison attained much higher target pressure. | Early shock heating and interface structure can dominate even with many directions. Their separate liner/target and atomic physics differ from this single-composition burning flow. |
| [Langendorf et al. (2018)](https://arxiv.org/abs/1805.09933) distinguish shock-forming collisional mergers from interpenetrating mergers and measure ion/electron temperature separation. | A hot single-fluid perturbation study answers a narrower question than a cold delivery/merger calculation. |
| [Haines et al. (2014)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.053302) found that simulations begun in 2D could underpredict later mixing and turbulent kinetic energy relative to fully 3D reference calculations. | Axisymmetric results are an initial screen. Azimuthal breakup and three-dimensional mixing remain open. |
| [LaJoie et al. (2024)](https://arxiv.org/html/2401.11066v2) measured a full spherical liner assembled from 36 argon jets, including a transition from interpenetration toward collisional stagnation. | A small jet count can form a liner under some conditions. The experiment does not determine the packet count for this much denser burning pulse; merger regime and subsequent compression must be considered together. |
| [Thompson and Cassibry (2019)](https://doi.org/10.1063/1.5074100) modeled oblique jet mergers and associated simulated interface undulations with Kelvin–Helmholtz instability. | Merging angle and interface shear belong in the packet study, alongside the amplification of imposed timing errors. |

Yates et al. also discuss timing requirements of order 100 ns for their PJMIF
application. That is useful evidence that timing matters, while the numerical
tolerance depends on the liner, target, speed, size and acceptance criterion.
The present pulse therefore needs its own timing sweep.

## Early merger heat and the engineering extrapolation

A supplementary radial screen keeps every cell's mass and total energy fixed
while converting a fraction f of its resolved inward kinetic energy into
internal energy. This estimates the consequence of irreversible merger heating.
At 768 radial cells, f=0.01 retains 95.73% of smooth output, f=0.025 retains
90.32%, and f=0.03 retains 88.68%. The 90% crossing is near f=0.026. A second
screen credits existing planned preheat: the initial internal energy becomes
the larger of its prescribed value and f times the original kinetic energy,
with the increase taken from kinetic energy. Its crossing is near f=0.073.
These are sensitivity experiments at the prepared state, not simulated cold
merger histories. The latter construction assumes the external heating profile
can be adjusted locally to credit merger heat.

For approximately equal angular patches, the transverse kinetic fraction is
about 2/N. If excess heat is deposited before additional effective linear compression
by a factor lambda **up to the prepared initial profile**, an ideal monatomic adiabat amplifies that heat by roughly
lambda squared. A useful conditional screen is therefore
`N >= 2 * lambda^p / f_90`, with p approximately 1.5–2 for the thermal responses
being considered, alongside the directly simulated hot-flow packet requirement. This radial
heat screen measures entropy sensitivity; it supplies no independent evidence
of instability growth.
The heating location, compression history and preheat credit are assumptions in
this extrapolation. They are not outputs of the hot-flow calculation.

For example, lambda=1 gives an entropy constraint of roughly 27–77 directions;
lambda=4 gives roughly 220–1,230 when both exponent and preheat-credit choices
are varied; lambda=10 gives roughly 870–7,700. The direct angular-flow threshold
also applies. Compression during the simulated hot pulse is already included
in f_90 and must not be counted again in lambda. Shaped, narrow incoming packets can change when their edges meet,
so the earlier spherical packing-density backmap does not determine lambda.
This dependence is the reason to retain a broad engineering credible interval.

## Three-dimensional cross-check

The Cartesian branch uses the same EOS, nuclear source and HLLE fluxes with
static nested refinement. The initial base grids have 32 or 48 cells per axis
across a 1,440 m cube. Four nested central refinement levels give minimum cell
widths of 2.8125 or 1.875 m. The main timing runs double the cube width to
2,880 m and add a fifth outer refinement level, preserving those central cell
widths. They evolve for 60 microseconds so that the delayed Gaussian tails can
finish burning. The packet-count runs evolve for 42 microseconds; their final
burn rate is checked separately. All six outer faces permit outflow. Each calculation
is paired with a smooth calculation on the identical mesh. Two-point Gaussian
quadrature in each coordinate initializes cell-averaged mass, vector momentum,
total energy and composition; its small global quadrature error is recorded.

For the packet-count comparison, antipodally paired Fibonacci directions define
spherical Voronoi patches. Every patch has a single incoming direction and the
same prescribed radial speed magnitude, density and thermal profile as the
smooth pulse. The resulting three-dimensional edges and junctions evolve under
the fluid equations. The velocity magnitude can vary along the shaped pulse.
This is the hot-stage merger of complete shaped angular portions, each carrying
a longitudinal profile; microscopic grain count and repeated subpackets are
separate quantities.

The three-dimensional timing comparison instead keeps the ideal radial velocity
direction and assigns independent patch phases on 128 Voronoi regions. Gaussian
quantiles are randomly permuted among the regions and normalized to zero mean
and unit RMS using their measured solid-angle areas. The sampled coefficient
range is -2.692 to +2.684. This fixes a light-tailed empirical error distribution;
RMS alone cannot characterize arbitrary outliers or spatial correlations.
This is one spatial permutation; variation across realizations remains part of
the subjective uncertainty. The
same evolved-reference trajectory supplies each phase. Timing and direction
errors are initially tested separately, each against its own smooth control.

In the earlier-seeded 2D branch, a two-lobed 4.5-microsecond RMS error retains
95.32% of output, while an order-8 pattern at 5.5 microseconds retains 89.22%.
The error spectrum belongs in the definition of a tolerance.

The two-dimensional large low-order **late** phase seed exposed an inadmissible internal
energy in the first radial cell. Those aborted results are excluded. A separate
sensitivity makes phase anisotropy vanish quadratically within a small central
radius and records its changed initial mass and energy. This modification tests
the origin treatment; it does not count as a numerical repair validated solely
by successful termination.

## What the instability tests establish

The multidimensional calculations evolve pressure gradients, transverse flow,
shock interactions and the burning response after a specified disturbance.
Small perturbations produce finite-time amplification, and larger timing or
direction errors measurably lower integrated output. The completed 3D hot-flow
packet tests retain 87.05% with 64 directions and 92.85% with 128; interpolation
in inverse count places their 90% crossing near 86 directions. The pictures can
remain fairly round while output falls below the acceptance criterion.

These histories do not isolate a unique exponential growth rate or assign each
lost joule to Kelvin–Helmholtz growth, geometric dephasing, or merger heating.
They answer the operational question of how the perturbed flow performs over
its finite burning interval. The separate radial heat screen estimates an
additional preparation uncertainty and is not evidence that a particular
instability occurs.

The largest next improvement would be a 3D calculation that starts while shaped
packets are still distinct and predicts contact, stopping, heat deposition and
angular mixing before the prepared profile. That would replace the dominant
packet-count judgment with a computed assembly history. Multiple spatial jitter
realizations and resolved joint timing/count tests would then tighten the
conditional engineering interval. This study covers one Stage A pulse; the
other reaction stages and sustained plant operation require their own checks.

## Reanalysis and reproduction

`results/run-inventory.json` records every retained Athena attempt, including
failures, exact arguments, executable hashes, source hashes and execution
status. Successful histories and diagnostic summaries are in `results/cases/`;
compressed logs and exact input decks are in `results/logs/`. Small independent
Lagrangian checks and EOS/Noh/limiter/reference-replay checks are retained beside
them. `results/source-versions/` preserves the exact problem, EOS and framework
patch for each executable generation. No aborted case supplies a completed
threshold point.

`results/reproduction-inputs.tar.gz` contains the exported EOS, reaction rates,
initial profile, signed reference trajectories and timing coefficients.
`results/reproduction-inputs.json` records each file's size and SHA-256 hash.
The wide reference table copies the verified uniform atmosphere beyond the
original radial boundary; `pad_phase.py` rejects a nonuniform exterior.

`infer_thresholds.py` reproduces deterministic brackets from retained summaries;
`subjective_intervals.py` propagates the explicitly stated judgment distributions
from separately supplied numerical anchors. Its draws are uncertainty
bookkeeping, not additional fluid simulations. `plot_uncertainty.py` plots the
result and assembly-prior sensitivity. Field-slice plots require the original
scratch VTK files; the rendered PNG and SVG figures are retained in `figures/`.

To rebuild a model, use a fresh checkout of the pinned Athena commit and
`build.py --source-dir results/source-versions/<version>` with the appropriate
`--geometry` and `--flux`. Exact initial tables can be unpacked into a new
scratch directory, and archived input paths relocated there. The run launcher
has the original authorization deadline hard-coded and will reject new runs
after it; a future simulation campaign needs a new authorized time window.

## Status

The study is running. Pilot numbers are provisional until the refinement,
source-integration and physical-approximation checks are complete. No article or
accepted Infinite Wiki chronology has been changed by this investigation.
