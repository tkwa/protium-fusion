# ✨ Support topology, pulse motion and a relevant guide fault

The existing thin-shell support arithmetic can be given a more specific candidate topology. Three radially separated families of guided latitude rings, about three orthogonal axes, avoid the repeated overpasses of an isotropic great-circle mesh. Each ring carries equal opposite streams. The new continuum calculation closes the vector force balance including the linear Solar tidal field. It leaves guide implementation, discrete attachment loads and nonradial stability for a subsequent model.

`support_topology.py` reads `sector-account.json`; its results and hashes are in `support-topology.json`. The original station files remain unchanged.

## Why a passive thin shell is a different design

The total installed areal mass is about 8014 kg/m². The old conservative self-gravity plus maximum tide produces **0.03532 Pa** inward pressure. Those loads are modest locally. At a radius of 1.1×10⁹ m, however, an ideal homogeneous 2.5 m shell with E=10 GPa and Poisson ratio 0.2 has classical buckling pressure only **6.09×10⁻⁸ Pa**. Its membrane compression would be about 7.77 MPa. Strength alone therefore gives a misleading answer: long-wavelength buckling determines whether an unbraced thin shell can stand.

The ideal equal-face sandwich scaling at the same face-material amount requires approximately 855 km face separation with that modulus, before core shear, imperfections and attachments. This is a comparison showing that passive bracing would be a substantial architectural change. It does not rule out a lighter high-modulus deep truss optimized separately. A tensile skin supplies inward force and cannot by itself balance inward self-gravity.

Classical spherical-shell buckling and its strong sensitivity to imperfections are well-established; the present numerical application uses declared future geometry. [Sieber, Hutchinson and Thompson, 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6451975/). The 10 GPa modulus is an illustrative graphite scale, not a qualified hot irradiated material property.

## Three noncrossing ring families

Let x point radially away from the Sun, y along the station's orbit and z normal to that orbit. A latitude-ring family about unit axis e_i, rotating with angular speed ω, delivers reaction per receiver area

`F_i = σ_i ω² [r − (r·e_i)e_i]`.

Here σ_i is the combined surface mass density of its two counterstreams. Their velocities have equal magnitudes and opposite signs. Within each family, rings lie at different latitudes and do not cross. Put the three families on separated cold guide decks. Local connectors carry the components needed to sum their forces; their tangential components need not vanish individually.

For an isotropic load, equal family densities give `Σ_i F_i = 2σ_i ω²r`. At a fixed equatorial speed, this requires 1.5 times the mass of a hypothetical uniformly maximum-speed stream, but has the same kinetic energy because the latitude-average squared speed is 2/3 of the equatorial value.

The actual Hill-frame calculation improves the scalar tide allowance. Let `a=GM/(2R³)` and `n` be the station's solar orbital angular frequency. A stationary shell needs support

`Σwall diag(a−3n², a, a+n²) r`.

It is supplied by

`σx = Σwall(a+4n²)/(2ω²)`;
`σy = Σwall(a−2n²)/(2ω²)`;
`σz = Σwall(a−4n²)/(2ω²)`.

All three densities are positive for this station. The required condition `a>4n²` holds. These are the leading circular-orbit Hill equations; moving stream elements also have Coriolis forces. Opposite streams cancel their leading net Coriolis reaction, while individual guide channels still carry it. At the equator it is approximately 0.39% of their centripetal force. [NASA's linear relative-motion equations](https://ntrs.nasa.gov/api/citations/19800019595/downloads/19800019595.pdf?attachment=true), §3.1; the paper uses a different in-plane axis ordering.

At a 10 km/s equatorial speed, the surface stream densities are **0.2096, 0.1174 and 0.08663 kg/m²**. Total circulating mass is **6.290×10¹⁸ kg**, only **0.0251 kg/m² above** the old scalar allocation. Total stream kinetic energy is **2.097×10²⁶ J**. Retaining the assumed 0.1% loss per revolution gives **3.034×10¹⁷ W**, below the old conservative 4.274×10¹⁷ W. No additional guide mass is declared free: magnetic structures, sensors, local attachments and repair equipment must be itemized within the existing services.

With 1 km meridional pitch, each family's combined rail line density is about **210, 117 and 87 kg/m**. Maximum centripetal reactions are **19.1, 10.7 and 7.88 N/m**, respectively. Approximately 10.4 million rings cover the three families at that pitch. Finite pole caps have small area and can be bridged by the local frame; their exact layout remains to be drawn. The continuum construction neglects the small radial separation and the streams' own redistribution of gravity, both controlled geometric corrections. It does not establish all angular/shear modes merely by satisfying static force balance.

The fixed 0.1%-per-revolution loss is still an engineering requirement. For a circular reference it corresponds to tangential drag about **7.96×10⁻⁵ of normal guide force**. Conservative magnetic deflection can in principle do no work; hysteresis, eddy currents, sensing, actuator losses and errors determine actual dissipation. A permanent magnet alone does not supply arbitrary stable three-axis confinement. The next component calculation should specify the moving elements and guide mechanism rather than interpreting this scalar loss as measured hardware.

## Permit common breathing; control differential errors separately

Outgoing pulse impulse is 6.57 N·s/m² for C13, 5.16 for A, 6.18 for B2 and 0.0434 for B1. The source team's archived initial profiles add outward launch-reaction impulses of 1.67, 1.71, 1.84 and 0.0161 N·s/m². These are spread over different times; the bounds below add their magnitudes conservatively.

One ideal estimate gives centimetre-scale motion if the stream support could be unloaded immediately after each pulse. That is **not a selected control law**: switching off a curved guide releases its moving mass. Instantaneously cancelling every pulse by changing stream speed would impose a separate large reversible-energy throughput. It is unnecessary to insist on rigid radius throughout every event.

For a periodically repeated impulse Δv with its mean force balanced, the free-mass periodic displacement range is `Δv T/8`. Summing the four stage bounds gives approximately **2.25 km** common-radius range. Individually conserved stream angular momenta give `v∝R⁻¹`, outward acceleration proportional to R⁻³ and self-gravity proportional to R⁻². The uniform radial mode is then restoring, with `ω0²=GM/(2R³)` and period approximately **3.60 years**. Holding speed strictly fixed would change that stability sign, so loss-compensation control should preserve the intended angular-momentum response.

A conservative Fourier triangle bound for the periodically forced **tide-free** linear self-gravity oscillator gives **3.04 km** summed range. It neglects damping and includes both launch and capture impulses. This is known common ephemeris motion of order 3×10⁻⁶ of the radius. It changes launch delays and aiming coordinates; it is not automatically a kilometre-scale random convergence error. The independent algebra audit also finds that fixed Hill tides matter for stability. A directional radial screen gives restoring coefficients `(a−12n², a, a+4n²)`, and here `a/n²≈9.6`: the first is negative. A constrained common spherical average instead has positive `a−8n²/3`. These are reduced-mode screens; deformation of an entire ellipsoid changes its gravity and guide geometry. They establish that static force closure cannot replace tidal-mode control. A sparse nonradial model must separately determine differential sector deformation, phase errors, tides and the response to localized source footprints. The common-mode bound is not a guarantee of those tolerances.

## Put the fast streams outboard and distinguish escape from impact

At the equator, loss of guide curvature causes a stream to depart tangentially with relative outward acceleration near 0.091 m/s². It crosses a hypothetical 1 m overhead boundary in **4.69 s**, or a 10 m boundary in 14.8 s. A 1 km rail interval contains order 10¹³ J of kinetic energy. Therefore an enclosing fragile roof close to the stream creates an unnecessary impact problem.

A preferable candidate places the guide decks **outside the radiator/receiver stack**, with open outward escape directions. A released tangent moves away from the main receiver rather than into it. Such placement removes a requirement to stop the entire ring within five seconds. It transfers the design problem to spare-guide takeover, controlled local diversion and avoiding the other guide layers and launcher structures. Their actual interception geometry must be calculated. Dumped mass is a separate material loss; recapturing it is optional recovery equipment, not an uncharged assumption.

The receiver itself moves much more slowly after loss of local support: the gravity/tide envelope produces 0.1 m differential displacement in roughly 213 s and 1 m in roughly 674 s. Neighbouring frame load transfer and local redundant guides therefore have a substantially longer response window than an enclosed-stream impact. A sector fault model should use these two different clocks. It should not substitute a generic instantaneous global collapse scenario.

## Next discriminators

1. Draw one cold guide channel and its moving material, including field gradients, confinement, drag, sensing and a spare route. Determine whether its complete linear mass fits the service allocation.
2. Assemble a low-order network of receiver sectors and the three guide families; distinguish uniform breathing, low angular modes and local loss of one support line.
3. Use source outlet positions and external accelerator tubes to check guide-layer crossings, open escape directions and actual radiator view factors.

These tasks can resolve geometry and control requirements computationally before physical qualification. Existing calculations establish a candidate force topology and tolerable common motion, not a completed magnetic bearing or irradiated structure.

## Follow-on global-mode calculation: a frame is necessary and lightweight

The next bounded calculation is now saved in `support_modes.py` and `support-modes.json`. It uses nine components of an affine deformation of the sphere, so it includes uniform size, ellipsoidal shape, shear and orientation. The receiver becomes a thin homoeoid under that map. Its exact gravitational energy is

`Ug = −GM²/4 ∫₀∞ du / sqrt[(Ax²+u)(Ay²+u)(Az²+u)]`.

For each guide family, slow deformation turns its circles into ellipses. Conserving circulation action makes its kinetic energy inversely proportional to the squared ellipse perimeter. This approximation concerns global changes slow compared with the roughly eight-day stream revolution. It does not apply to a fast local guide fault. The script includes the orbital Coriolis coupling and the Hill tidal potential.

The important result is that **rings alone do not supply adequate shear stiffness**. With no tides, an off-diagonal shear has squared frequency `−19a/20`, even though uniform breathing has `+a` and diagonal shape deformation has `+7a/40`. The corresponding unsupported growth time is approximately 215 days. The full linear-tide/Coriolis case also has growing modes. A separate numerical quadrature of the nonlinear homoeoid energy and actual mapped-circle perimeter reproduces the quadratic second variation to 4.7×10⁻⁸ in units of a. This check does not reuse the polynomial energy expression.

A cold triangulated tension network supplies a concrete missing component. At an isotropic extensional stiffness **E·h=10⁸ N/m**, its ideal material equivalent is **2 kg/m²** for density 2000 kg/m³ and modulus 100 GPa. A real lattice requires orientation factors, joints, local load spread and redundancy; a provisional **20 kg/m² service subdivision** gives ten times that material equivalent. The required stiffness is a specification, not a fabricated-network measurement. Pretension keeps cable members in tension through small shape changes. For example, 10,000 N/m isotropic pretension adds only about 1.8×10⁻⁵ Pa of inward pressure, which the stream setting must also balance.

Shape stiffness alone leaves weak attitude/shape coupling in the tidally loaded model. Adding an explicit orientation controller with 30, 100 or 300 day response and damping ratio 0.7 removes positive growth in all nine computed affine modes at the above stiffness. Those are control requirements, not an assumed passive effect. The two counterstreams provide an angular-momentum reservoir: a small difference in their speed can exchange angular momentum with the supported structure, analogous to a reaction wheel. A detailed actuator and momentum-management account remains necessary, as do local and higher angular modes.

The result supports a **cross-braced cold guide deck with slow active orientation control**. It does not require a massive rigid shell or fast cancellation of each fusion impulse. The remaining calculation is a discrete-sector model with the actual mesh, guide attachments and a localized fault. The global result and the [candidate guide calculation](guide-component.md) now supply numerical inputs for that model.

## Local guide modes after the global affine screen

The global affine result is now supplemented by [guide-local-control.md](guide-local-control.md). A flexible rail following its streams has a local growing mode near 45.7 km wavelength with a 38.9-second e-folding time in the representative reduced model. Full stream-and-guide feedback can stabilize the tested wavelength/parameter family with explicit actuator reaction. A 0.1-second first-order actuator meets the examined 20–1000 km control band; a one-second response fails some cases. This adds a concrete distributed-control requirement. It does not alter the distinction between predictable common-radius motion and differential guide error. Real-space sensing, actuator saturation, field dissipation, discrete carriers and network faults remain the next design steps.
