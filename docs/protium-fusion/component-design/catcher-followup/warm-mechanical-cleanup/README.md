# ✨ Warm mechanical chamber cleanup

A warm molecular-pump train is a physically credible route from the approximately 1,000 K reactor cavity to a small, isolated downstream gas system. Existing pumps establish the required helium compression scale. The unresolved engineering is a hot rotating assembly with adequate light-gas capture, controlled clearances and low idle power. Broad cold surfaces and superconducting bearings are not required by the pumping mechanism.

The bounded calculation in [calculate.py](calculate.py) and [results.json](results.json) specifies a testable configuration. It does not qualify a working 1,000 K pump. The [primary pump benchmarks](benchmarks/README.md) and [rotor/idle-power review](rotor/README.md) contain the underlying manufacturer data and inspected experimental papers. This packet preserves the existing reactor and catcher studies.

## A concrete flow and heat topology

Warm graphite louvers admit cavity gas into a warm antechamber while intercepting the approximately radial fast debris. An axial moving-blade stage supplies inlet speed; further blade and/or molecular-drag stages supply compression. The rotor, stators and broad inlet surfaces remain near 1,000 K. Gas leaves through a small sealed conduit after compression, crosses a 300 K thermal intercept, and reaches isolated downstream getter/cryo or backing pumps. The cleanup foreline is separate from the direct catcher's pulsed collection plenum.

The rotor lies behind particle-intercepting turns. These protect it from direct ions and large debris while retaining a measured gas transmission. Baffle faces collect heat through the hot receiver circuit. Their geometry must preserve the receiver's radiation shielding; an average carbon mass alone cannot certify an open gas path. The already analyzed graphite receiver supplies a useful local heat-buffer reference. A conservative additional 350 kg/m² column over a 10% intake footprint would cost 35 kg per receiver m². Reusing or reshaping receiver material may reduce that extra allocation; the study does not require two complete walls.

The heat engine's working fluid, the cleanup gas and the material-regeneration streams remain separately sealed. Collected N-bearing molecules can enter the same mechanical cleanup train; actual speed/compression and subsequent chemical routing remain species-specific. A rotor coated by condensed material is cleaned or isolated as a local module. Particle interception and coating/balance monitoring require an actual maintenance design; the present incident mass does not imply that an entire rotor must be replaced after each pulse.

## Compression with useful throughput

The [official HiPace 700 H data](https://www.pfeiffervacuum.com/global/assets/134975/1/Product%20Leaflet%20HiPace%20700%20H%20Global%20EN.pdf) give H2/He speeds of 555/655 L/s and maximum compression ratios of 2×10⁷/10⁹ under ordinary conditions. These are current hardware benchmarks, not hot ratings.

For each gas species, the molecular-flow relation is

`Q = S0 (p_in − p_fore/K0)`.

At p_in=10⁻¹⁰ Pa and p_fore=10⁻³ Pa, K0=10⁷ gives zero useful flow. K0≥10⁸ retains at least 90% of negligible-backpressure speed; K0=10⁹ retains 99%. These must be **partial pressures**, not one total pressure assigned independently to every constituent. A lower chamber partial pressure requires a larger compression ratio at the same foreline partial pressure.

[Pfeiffer's primary design account](https://www.pfeiffervacuum.com/no/en/knowledge/vacuum-technology/knowledge-book/4-vacuum-generation/4.9-turbomolecular-pumps.html) relates logarithmic compression to blade speed divided by molecular thermal speed. As a similarity screen, heating every relevant stage from 300 K to 1,000 K requires 1.826 times the speed or effective logarithmic compression length. Starting from a 400 m/s reference scale, a 700 m/s hot train needs approximately 1.286 times the reference H2 compression length to reach K0=10⁹. Helium needs approximately 1.043 times in that comparison. These factors specify requirements; they do not preserve capture or clearances automatically.

Ideal Holweck channels give an exponential in `u L/(vbar g_geom h)`, with geometry factor of order 1–3. At 700 m/s, H2 at 1,000 K requires an ideal L/h of approximately 96–288 for K0=10⁹. Real interchannel bypass can greatly reduce that result. An unrestricted gap joining the entire high- and low-pressure ends cannot be ignored; seals and staged labyrinths must prevent that bypass. A broad axial inlet followed by narrower drag channels avoids trying to obtain the whole inlet speed from a small annular Holweck opening.

Atomic H is a separate case. A simple same-speed extrapolation predicts less capture and compression than H2. It needs a longer stage train or demonstrated surface recombination before the reference H2 account applies. The dataset includes an explicit atomic-H sensitivity rather than assuming that all hydrogen is molecular.

The direct catcher reaches several pascals in its immediate-release pulse model. Sharing that pressure with a cleanup pump would change the required K0 for 90% speed from 10⁸ to roughly 4×10¹¹ at a 10⁻¹⁰ Pa inlet. Separate pumping/valving of the cleanup foreline is therefore a concrete integration choice, not a negligible bookkeeping detail.

## Required warm area

Normalize speed to the one-sided wall flux `vbar/4`. The nominal DN160 opening and measured speeds imply effective room-condition inlet captures of approximately 0.0622 for H2 and 0.1034 for He. Preserving those values at the selected hot speed is an apparatus requirement.

For a warm antechamber, baffle conductance and pump speed combine in series:

`C_b/A_wall = f_b τ_b vbar/4`

`S_p/A_wall = f_p α_p vbar/4`

`g_eff = [1/(f_b τ_b) + 1/(f_p α_p)]⁻¹ [1 − p_fore/(K0 p_in)]`.

This network includes repeated returns inside the antechamber. A one-pass probability product describes a different geometry and cannot be substituted indiscriminately.

A reference with intake fraction f_b=0.10, baffle transmission τ_b=0.30 and K0=10⁹ needs **f_p=0.244834 of receiver area as H2 pump aperture** for g=0.01. The same aperture and He benchmark capture give g≈0.01360. The baffle transmission is a specified target for the particle-protecting geometry, not a ray-traced result. Stronger baffling can be accommodated with more intake/pump area; some narrow-aperture combinations have insufficient conductance even with an arbitrarily large downstream pump. The saved parameter family makes this boundary explicit.

The quantity g is an isotropic encounter-equivalent removal probability. Directional plume capture and protected launcher corridors need the actual layout. The gas study's chamber-average inventory and pressure do not replace that local transport calculation.

## Surface and structural requirements

A smooth basal carbon sleeve is not automatically a good drag surface. [Magnico and To's primary MD study](https://doi.org/10.1016/j.ijheatmasstransfer.2023.124371) reports mostly specular, short-residence He reflections at the high-temperature end of its explored range. The accessible indexed abstract/conclusion supports this qualitative warning; no 1,000 K numerical accommodation coefficient is assigned here. Angled moving blades transfer normal momentum even for specular reflection. Corrugated or rough drag surfaces offer another engineering variable. The actual gas-surface kernel must enter the moving-wall transport model.

At 700 m/s and density 1,800 kg/m³, the thin-ring stress scale is **0.882 GPa**. An ideal uniform isotropic disk with Poisson ratio 0.2 has maximum stress **0.3528 GPa**. At 500/1,000/2,000 m/s the corresponding ring scales are 0.45/1.8/7.2 GPa. These are reference loads, not strength margins for a slotted anisotropic rotor. Blade roots, hubs, interlaminar transfer, thermal strain and fatigue need the actual component design.

Existing room-temperature composite [spin tests reached 1,310 m/s](https://doi.org/10.1163/156855102753613273). That does not qualify their matrix or hub at 1,000 K. Conversely, [measured C/C interlaminar properties](https://doi.org/10.1016/0022-3115(96)80018-3) show why high hoop-fiber strength alone is insufficient: one unidirectional material had approximately 9.3 MPa interlaminar shear strength at 1,000°C. The proposed carbon-rich rotor needs a suitable three-dimensional load path. The study does not require the 2,000 m/s endpoint.

Warm matched rotor/stator temperatures reduce differential expansion. Bearings and electronics may be thermally isolated behind the hot stages. Noncontact electrostatic or magnetic support is physically available without superconductivity; its actual force, field return, controller and loss account remains to be designed. Shielded electrostatic drive avoids imposing a new magnetic field throughout the primary debris flight. A microgravity installation reduces static bearing load, while imbalance, startup and transient forces remain.

## Power and thermal isolation

At the stated chamber pressure, gas compression itself is tiny: the separate H2 case requires **2.61×10⁻⁸ W/m² receiver** at a declared 50% isothermal-work efficiency. Rotor and controller losses can be many orders larger, so the efficient-manufacturing premise does not by itself determine operating power.

A [primary magnetic-bearing pump experiment](https://www.magneticbearings.org/publication/detail/ismb4_82/) measured 48 W total high-vacuum input for a 260 mm-diameter rotor. That is approximately 904 W per rotor-face m². For the new reference, a **chosen 1 W/m² receiver idle-power budget** requires 4.084 W per pump-aperture m², about 221 times lower in areal terms than that small-module experiment. This is an engineering target, not a thermodynamic minimum or an imposed plant failure threshold.

Larger, slower modules are relevant. If a 48 W module input could be held fixed, a radius of 1.934 m would meet that areal target. At 700 m/s such a rotor turns at approximately 3,456 rpm. Bearing force, mass, clearance and losses do not remain constant automatically with size; the calculation identifies a component scale worth examining. It does not multiply a small pump into a large one while holding every property fixed.

For a 1 m-radius example with rotor mass 100 kg per pump-aperture m², the disk-model stored energy is 3.00 MJ per receiver m². A 1 W/m² loss corresponds to a 34.7-day energy-decay time. That mass and loss are explicit requirements. The dataset provides radius, mass and loss sensitivities. Warm dissipation recovered as 1,000 K heat can reduce its net electrical penalty; if all of it is usefully recovered at 60% conversion, a 1 W input has a 0.4 W net penalty. No such recovery is silently granted.

Equally hot rotor and stator surfaces have no steady net radiative transfer. Cooled bearings or motors exposed to the warm train do. To keep their combined warm-to-cold transfer below a chosen 1 W/m² receiver, the reference area permits:

- effective emissivity × view factor × cooled-area ratio below **7.26×10⁻⁵** per pump-aperture area, for a 1,000→300 K view;
- total conductive G below **0.00583 W/(K·m² pump aperture)** across that temperature difference.

These are component acceptance curves. A whole cold stator facing the broad warm rotor would exceed them. Small isolated bearings, long thin support paths and radiation shields can be tested against the same account. Their conduction and radiation are heat bypasses; motor input is a separate electrical term.

## Small downstream cold interface

Compression makes a much smaller cold admission path possible. Use a 1,000 K warm foreline at 10⁻³ Pa, a 300 K downstream reservoir at 10⁻⁴ Pa, and explicitly include reverse molecular flux. The reverse/forward equilibrium flux ratio is `(p_cold/p_hot)√(T_hot/T_cold)=0.1826`; simply comparing the pressures would miss thermal transpiration.

For the target g flow, the ideal aperture requirement at transmission 0.3 is only 4.08×10⁻⁹ per receiver area. A more concrete 1 m-radius module example uses a 1 mm-radius, 0.1 m-long gas conduit. Including a specified factor 0.5 for its bend, its diffuse long-tube estimate has 2.67 times the required conductance. At 10 µm wall thickness and k=10 W/(m K), the total conduit-wall leak is **0.000343 W/m² receiver**. Assigning all hot radiation entering the conduit mouth to the cold system gives **0.01377 W/m²**. A first 300 K intercept precedes any colder He collector.

This conduit is a design example, not a completed baffle calculation. Its heat is much smaller than exposing the broad chamber intake directly to a cryogenic collector. Bearing and motor thermal paths remain separately charged. The global cleanup gas rate scales with the chosen cavity partial pressure; it is distinct from the much larger local catcher throughput.

## Mass, repair and next discriminator

The HiPace mass divided by nominal aperture is approximately 587 kg/m². At the required aperture fraction, this scale gives **143.7 kg/m² receiver** of pump package. Adding the conservative 35 kg/m² baffle example gives 178.7 kg/m². That fits within the frozen 350 kg/m² reference's approximately 678.6 kg/m² unassigned service allowance, leaving about 500 kg/m² before other follow-up equipment is assigned there. This is a comparison with an existing whole-pump mass scale, not a manufactured hot-pump bill of materials. Extra compression stages, isolation hardware and altered rotor construction must fit the actual allocation; the existing services must not be charged a second time.

The rotor example stores 12.25 MJ per pump-aperture m², or 38.5 MJ for a 1 m-radius module. Its local casing and stopping geometry must handle rotor fragments; the small normal impulse of the fusion debris does not qualify a catcher for a tangential rotor burst. This is a concrete enclosure function to allocate inside the package, with local module isolation. No separate containment mass is silently assumed beyond the package comparison.

Continuous local repair and reuse can address deposits and damaged surfaces under the accepted manufacturing premise. A declared 2 MJ/kg rework scenario for the benchmark package would cost approximately 0.00091 W/m² at 10⁻⁴ of its mass per year, or 9.11 W/m² for complete annual reworking. Neither rate is imposed. Hot rearrangement and selective repair can have lower dissipation; chemical work, thermal cycling and gross handling remain distinct.

The most useful next computation is a **moving-wall test-particle unit cell** with H/H2/He, actual blade angles, hot surface scattering, forward and reverse launches, and clearances. It should return inlet capture and per-stage compression for the same warm baffle geometry. A separate bearing/drive model should then close the roughly 4 W/m² aperture target, or show a better radius/loss allocation. These are large enough to change the design and can be investigated computationally before a physical experiment.

No current evidence supports a narrow subjective 80% credible interval for the complete hot pump. This packet instead distinguishes measured benchmarks, exact network algebra, geometric reference loads and explicit future performance requirements. It finds a plausible mechanical architecture without claiming that the remaining component design is complete.
