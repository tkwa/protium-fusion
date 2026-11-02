# ✨ A concrete support-guide candidate

A cold graphite carrier stream in a superconducting magnetic quadrupole gives a calculable implementation of the support-ring idea. The present screen is conditional on carrier susceptibility, field quality and the operating conductor specification. It requires modest bending acceleration despite the 10 km/s longitudinal speed. It also shows why a short failed guide segment need not cause material loss.

`guide_component.py` reads the new support topology and writes `guide-component.json`. It does not assume direct electrical recovery of fusion products, and it does not prescribe a launcher.

## Transverse confinement

Graphite is diamagnetic. Published levitation models use pyrolytic graphite susceptibility components approximately −85, −85 and −450 parts per million; experiments also show that material quality changes the effective susceptibility. Use the weak magnitude **85×10⁻⁶** as a declared scalar lower scale, density 2200 kg/m³, and a quadrupole `B²=G²(x²+y²)`. The resulting small-displacement acceleration is

`a_mag = −ωt² x`, with `ωt²=|χ|G²/(μ0ρ)`.

At `G=8 T/m`, a 0.25 m bore radius and a 10 km/s equatorial stream, the field at the bore edge is 2 T. The required centripetal acceleration is **0.0909 m/s²**. The transverse frequency is **1.403 rad/s** and the equilibrium offset is **46.2 mm**. The heaviest ring family carries 104.8 kg/m in each direction. At 50% longitudinal fill, each graphite carrier can have radius **174 mm**, leaving approximately **29.6 mm** transverse clearance in this conservative geometry.

This is a transverse potential calculation. Orientation-dependent susceptibility, carrier rotation, inter-carrier interactions, longitudinal spacing and fringe fields need explicit tracking. The weak-direction estimate avoids treating graphite as an isotropic material with its strongest response in every direction. Passive diamagnetic confinement supplies a route around the instability of freely levitated high-field-seeking ferromagnets. [Graphite force model and experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC10304365/).

## Coil mass, finite winding and refrigeration

An ideal cylindrical current sheet `Kz=K0 cos(2θ)` at radius a gives `G=μ0K0/(2a)`. Here its summed absolute axial current is 3.18 MA per bore. This is distributed winding current, not one exposed 3 MA wire. Opposite carrier directions use separate bores; the three ring families together require approximately 0.006 m of bore length per receiver m² at 1 km rail pitch.

A thin-current-sheet calculation can understate mass at low engineering current density. The code therefore also fills a finite annulus with angular winding fraction `|cos(2θ)|`, giving

`G=(μ0 Je/2) ln(b/a)`.

For operating engineering current densities of **10, 30 and 100 A/mm²**, the annular outer radii and conductor masses are recorded explicitly. The 30 A/mm² case needs about **8.1 kg/m² of complete conductor**; the deliberately low 10 A/mm² case needs about **70.5 kg/m²**. End turns, insulation, field-peak margin, cryostat, support and redundancy remain additional components within the existing 1000 kg/m² services. The finite winding calculation checks its own peak field rather than assuming every tape sees the bore-centre field.

The candidate uses 65 K conductors. Commercial REBCO data demonstrate relevant conductors and substantial current, but their self-field 77 K rating is not the operating rating for this magnet. The selected engineering current must be checked against the actual field direction, peak field, temperature and strain. [Manufacturer specification](https://superpower-inc.com/specification.aspx); [primary field/temperature measurements](https://arxiv.org/abs/1512.01930).

With both inner and outer coil surfaces behind 200 K shields, effective emissivity 0.003–0.03 and a refrigerator operating at 15% of Carnot, the code gives approximately **0.10–1.0 W/m²** of refrigerator electricity in the 30 A/mm² geometry. Supports, current leads, joints and alternating losses are additional. These are meaningful fractions of the old guide-dissipation allowance, but small fractions of station export. They must be assigned explicitly within the ordinary auxiliaries or added to them once. No room-temperature superconductor is assumed.

A complete material recipe is still needed. The conductor includes substrate and stabilizer, and a full-sphere installation can require a substantial copper inventory even when its total mass is a small plant fraction. Higher operating current density, conductor geometry and stabilizer choice affect that inventory. Electrical efficiency does not manufacture scarce elements without a material source.

## Segment faults and material retention

For a completely missing section of guide of length L, transit time is L/v. The initial tangent departure is `a L²/(2v²)`. Upon entering the next intact ideal quadrupole, the combination of displacement and transverse velocity gives oscillation amplitude

`sqrt[δx²+(δv/ωt)²]`.

| Contiguous missing guide | Tangent departure | Ideal recapture oscillation |
|---|---:|---:|
| 100 m | 4.55 µm | 0.65 mm |
| 1 km | 0.455 mm | 6.50 mm |
| 3 km | 4.09 mm | 19.9 mm |
| 10 km | 45.5 mm | 79.2 mm |

The first three fit the reference 29.6 mm clearance; the 10 km case does not. The actual maximum depends on entry phase, field ramps and adjacent structures. This makes segmented redundancy a concrete alternative to ejecting the ring whenever one element fails. The earlier seconds-to-boundary estimate describes a long continuously unsupported path; it is not the response to every local electrical fault.

Carriers ordinarily remain in circulation. Only actual unrecovered departures, nuclear conversion or unrecovered process losses require fresh material. Any outer recovery system must receive its own apparatus and energy account; it cannot be assumed free. A closed or locally redundant guide should be evaluated before imposing a routine carbon purge.

## Loss mechanism to calculate next

Graphite is electrically conductive, so longitudinal field ripple produces eddy-current heating. Patterned graphite and insulating composites containing graphite particles have been demonstrated specifically to reduce this damping. Their room-temperature levitation experiments do not demonstrate a 10 km/s guide, but identify a physical manufacturing route. [Patterned-graphite damping experiment](https://arxiv.org/abs/2211.08764); [insulating diamagnetic composite experiment](https://arxiv.org/abs/2312.01660).

The script includes the low-skin-effect conducting-sphere estimate `P/V=σ ω² Bamp² r²/20`. For a 100 m ripple wavelength, 1 mT amplitude and conductivity 10⁵ S/m, 10 µm conducting domains dissipate about 9×10⁻¹¹ W/kg; 1 cm domains give about 9×10⁻⁵ W/kg. These are small compared with the old 0.1%-per-orbit allowance of order 0.05 W/kg. A monolithic conducting carrier in appreciable fringe fields can behave differently. The calculation motivates segmented conducting domains and smooth field transitions rather than claiming that all actual guide losses are negligible.

The next useful calculation is a finite-coil unit cell with fringe fields, carrier orientation and a failed-cell transition. A static three-dimensional field solution plus neutral-particle tracking can test confinement and ripple without another fusion simulation. Global low angular modes remain a separate structural/control problem identified in [support.md](support.md).
