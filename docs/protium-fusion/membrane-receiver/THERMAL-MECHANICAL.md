# ✨ Membrane receiver: bounded thermal and mechanical assessment

The inner separator can be a supported, renewable component without its bare mass determining station performance. The strongest current candidate is a set of **locally closed, cooled cassettes**, with recoverable skins and separate launcher penetrations. The calculations establish required cooling spans, pressure-force paths and recycling duties. They do not establish a surviving irradiated skin or a complete receiver design. The controlling uncertainties are the actual soft photon spectrum, secondary-electron energy retention and return, ion-track damage, thermal contact and the operation of the braking gap.

This packet uses the frozen 64-cell open-loop station at radius **1.1×10⁹ m**, area **1.52053×10¹⁹ m²**. It preserves the earlier monochromatic screen as `thermal_reference.py` / `thermal-reference.json`, adds `spectral_proxy.py` / `spectral-proxy.json`, and calculates support/cooling requirements in `mechanics.py` / `mechanics.json`. The older reference's stopping-waveform rescaling is explicitly a historical diagnostic. All cooling numbers below use the **actual time-resolved stopping peaks** in `stopping-results.json`. Input and script SHA256 values are in the results.

## Prompt radiation: use the saved temperatures, not an artificial soft flash

The hydrodynamic calculation supplies grey escaped energy and outer-cell temperature. It does not supply a photon spectrum. The new proxy assigns each positive escaped-energy increment a normalized Planck energy distribution at that bin's midpoint outer temperature. Elemental NIST coefficients are interpolated logarithmically from **1 keV to 1 MeV**, preserving silicon's 1.8389 keV K edge; compound coefficients use mass weighting. The result is an explicit spectral-shape assumption, rather than a new radiation-transport result.

| Stage | Escaped prompt fluence (J/m²) | Escape-weighted kT (keV) | 1 µm silica average sensible rise from 1–1000 keV, full local retention (K) | Additional average rise if every sub-keV photon is absorbed (K) |
|---|---:|---:|---:|---:|
| A | 17,818 | 14.33 | 22.5 | 1.94 |
| B1 | 8,979 | 4.03 | 106.7 | 20.2 |
| B2 | 23,819 | 15.70 | 31.9 | 2.79 |
| C13 | 48,850 | 13.09 | 106.5 | 10.4 |

These are heat-capacity equivalents using constant **1,000 J/(kg K)**, before phase change, electron escape or cooling. They describe volume-averaged energy deposition. For the resolved C13 spectrum, the no-conduction surface upper screen is approximately **114 K**; finite-thickness average absorption gives **106.5 K**. The unresolved low-energy addition in the table is an average-energy bound, not a surface-temperature bound. Through-thickness diffusion across 1 µm silica has characteristic time about **1.7 µs**, versus the roughly **70 µs** C13 flash, which helps spread deposition; transient soft-photon absorption still needs resolving.

B1 is materially softer than the other stages. The all-1-keV C13 stress case that gave an approximately 15,500 K thin-limit sensible rise is therefore unsuitable as the selected spectrum. It remains a deliberately extreme sensitivity. The Planck proxy does not make soft photons irrelevant: a 10 nm silica film has resolved C13/B1 rises of **113.7/116.8 K**, but full absorption of the unresolved sub-keV energy would add **1,039/2,017 K**, respectively. The corresponding sub-keV energy fractions are only **0.0468%/0.494%**. This very thin film needs a sub-keV transfer calculation and electron-escape model; an average photon energy cannot settle it.

For C13 at 1 µm, the resolved full-retention rises are silica **106.5 K**, carbon **22.7 K**, BN **25.4 K**, and SiC **121.3 K**, under the same heat-capacity convention. These are a comparison of prompt absorption, not a ranking of complete receivers. Carbon introduces an additional carbon-return obligation; all four materials have different damage and chemical behavior.

Temporal sampling matters: assigning each C13 bin its starting or ending temperature changes the silica result from **91.8 to 126.0 K**, around the midpoint value 106.5 K. Halving/doubling all temperatures changes it to **405.8/22.5 K**. Those are sensitivity tests, not confidence limits. The Planck normalization and energy-partition checks pass **980/980**. A subjective 80% survival interval would obscure the presently unmeasured spectrum and retention assumptions, so none is assigned.

NIST's energy-absorption coefficient accounts for radiative losses during charged-particle slowing in bulk matter. It does not supply the retained fraction in a 10 nm or 1 µm film. Escaping secondary electrons can reduce film heating; electrostatic charging, sheath return and electrons arriving from the plasma or absorber can increase it. `spectral-proxy.json` also records the more conservative energy associated with all single-pass photon interactions. Multiple scattering, incoming particles and nonlocal return must be treated in a coupled transport model. Sources: [NIST coefficient definitions](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html), [C](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html), [O](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z08.html), [Si](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z14.html), [B](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z05.html), [N](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z07.html).

## Particle heat needs a physical path to a sink

For C13, 10 nm silica removes **68.43 kJ/m²**, with a timed peak **741.6 W/m²**; 1 µm removes **6.549 MJ/m²**, with peak **74.96 kW/m²**. A has the larger silica peak: **1.213 kW/m²** at 10 nm and **121.44 kW/m²** at 1 µm. Projectile energy removal is not automatically local heat. Define a retained fraction η and calculate the local load as η times the stopping result, then add incident-electron, sheath and infrared exchange separately.

An isolated film cannot store this particle energy. The C13 no-cooling sensible equivalents are approximately **3 million K** at either thickness. This is a statement that a heat path is required long before that temperature; it is not a predicted material temperature. Particle arrival extends over tens to hundreds of seconds, allowing active heat removal.

For a strip between isothermal support lines, the steady maximum rise is

`ΔT = q L² / (8 k t)`.

Using silica k=1.3 W/(m K), full retained stopping heat, and a 100 K allowable rise, C13 requires **L≈0.118 mm** at either 10 nm or 1 µm. A requires approximately **0.093 mm**. Permitting 1,000 K instead gives about 0.37/0.29 mm, but does not establish silica survival. Reducing η to 0.1 or 0.01 increases the allowable span by √10 or 10. The near cancellation of thickness occurs because thin-film stopping and in-plane conductance both scale with thickness. For the 100 K case, lateral diffusion time is a few milliseconds, short compared with the ion pulse. This calculation assumes isothermal edges, intact conductivity and negligible contact resistance.

A direct backing changes the heat path. For 1 µm silica at the C13 peak, `q t/k` is only **0.058 K**. Keeping the contact rise below 100 K requires contact conductance at least **750 W/(m² K)**, or about **1,215 W/(m² K)** for A. A supported catcher skin therefore has a credible thermal advantage if contact remains intact. It must provide vapor-separation/removal channels, compatible backing, recovered implanted material and a receiver-side replacement path. A solid backing immediately behind the skin receives the transmitted ions; it cannot also be treated as an empty downstream electromagnetic braking gap.

Three different architectures remain worth developing:

- **Free-span separator plus downstream braking gap:** fast outgoing ions cross the film, are slowed in the back gap, and reach the carbon receiver at the intended low energy. Fine support lines, electrodes, charge/current management and cooling must coexist in the gap. Frame occlusion and secondary-electron return are part of the budget.
- **Backed catcher skin:** a thin renewable surface is thermally attached to a cooled collection structure. Braking occurs upstream or in the backing, and deposited fuel is recovered there. This uses the backing's heat capacity and conduction, with a different particle-recovery design.
- **Directional renewable catcher:** swept liquid/solid curtains and nearby cold traps deliberately receive and recycle the outflow. Vapor is removed directionally or into local cells. It avoids relying on a single material being simultaneously an ion-transparent film and a robust hot vapor barrier, while introducing its own jet stability, sputter, condensation and particle-recovery requirements.

Liquid or sacrificial surfaces are replacements for a specific function, rather than automatic cures for the whole receiver. A liquid curtain that stops ions must close its momentum, isotope-recovery and nuclear-reaction ledger. An open curtain also releases vapor into the chamber unless geometry, cooling or differential pumping confines it.

## Closed-cavity infrared and local cooling hardware

Inward thermal radiation from a uniform spherical membrane returns to other membrane tiles. The diffuse-flight kernel is `K(t)=2t/tmax²`, with **tmax=2R/c=7.338 s**, mean delay **4.892 s**, and unit integral. There is **zero steady net cooling into the central cavity** in the uniform opaque case. A radiating back side exchanges energy with the real absorber/condenser and can heat the film when that surface is hotter.

For a symmetric grey sheet with absorptivity a, transmissivity t and a black back bath, the self-consistent steady outward coefficient is `a(a+2t)/(a+t)` multiplying `σ(Tfilm⁴−Tbath⁴)`. It is not generally equal to a single-face emissivity. Real materials require spectral optical properties and actual view factors. At film 1,000 K and bath 2,000 K, an effective coefficient of just **0.01** produces **8.51 kW/m²** of incoming heat. That exceeds the 10 nm C13 ion loss peak by an order of magnitude. Bath temperature, screens and optical coupling are therefore coupled design variables.

A 1% opaque support grid intercepts about 1% of the much larger complete particle flux. In the illustrative C13 10 nm case it raises the conservative cooling-peak sum from 0.742 to **7.08 kW/m²**. The 1 µm value becomes **80.56 kW/m²**; the largest all-stage value is **124.92 kW/m²**. These sums combine separately occurring maxima and are conservative for that specified grid. They omit uncomputed electron and infrared loads.

The explicit hardware requirements are:

| Quantity | 10 nm silica, 1% grid | 1 µm silica, 1% grid |
|---|---:|---:|
| Film areal mass (kg/m²) | 0.000022 | 0.0022 |
| Idealized 1 mm-deep grid material floor, ρ=3,000 (kg/m²) | 0.0299 | 0.0299 |
| Maximum peak cooling hardware if rated 1 kW/kg (kg/m²) | 7.08 | 124.92 |
| Mean captured film/grid particle heat over actual stage cadence (W/m²) | 1.052 | 25.118 |
| Pump work if it costs 1% of that heat (W/m²) | 0.0105 | 0.2512 |

A 1% grid at the required 0.1 mm pitch has line widths of order **0.5 µm**. The 1 mm depth and very small mass above are a stress-screen floor with extreme aspect ratio, not a manufactured channel design. Increasing line width, routing coolant, sealing and integrating electrodes can dominate it. `mechanics.json` separately provides **1, 10 and 100 kg/m²** reusable-frame/seal/cooling allocations and shutter allocations of `100 kg per m² of opening × opening fraction`. They are declared design budgets, not derived apparatus masses. Latest 1,000-launcher geometry is assessed separately.

A local coolant loop of round-trip length 20 m, speed 10 m/s, cp=1,000 J/(kg K) and temperature swing 100 K needs `2 q ℓ/(v cp ΔT)` of fluid per projected area: **0.142 kg/m²** for the 10 nm C13 example or **1.61 kg/m²** for 1 µm. This does not include pipes, pumps or heat exchangers. Transporting heat across a billion-metre shell instead would invalidate this small inventory. Heat conversion/rejection must remain local.

The mean values in the table are heat routed through cooling equipment, not additional station electrical consumption. The heat was already counted in the reactor energy balance. Cooling it to a less useful temperature can reduce electrical conversion; pumps, refrigeration and heat-exchanger mass must be added independently. No new inward-radiation credit or second recovery of the same energy is taken.

## Pressure closes locally, or it becomes a global shell load

Saved direct-particle momentum flux peaks are **0.104 Pa (A)** and **0.095 Pa (C13)**. The equilibrium graphite vapor pressure at 2,000 K in `results.json` is only **7.33×10⁻⁶ Pa**. Larger backside pressures in `mechanics.json` are parametric tests for a hotter or trapped plume, another material, or a failure, rather than the selected graphite equilibrium pressure. Ion sputtering and dust are separately unresolved.

For a circular, preformed cap of radius a, thickness t and sag w=δa, membrane stress is `p a(1+δ²)/(4tδ)`. A 1 µm cap of radius 0.5 m and δ=0.05 needs about **0.25 MPa** under 0.1 Pa. An initially flat film has to acquire the sag through strain; the script solves a spherical-cap elastic approximation rather than granting free slack. At sub-mm cooling spans, direct-pressure stress is small compared with thermal, defect and handling constraints.

A whole unsupported spherical film would need stress `pR/(2t)`, approximately **5.5×10¹³ Pa** for 0.1 Pa and 1 µm. Even at a 100 MPa allowed stress in a 3,000 kg/m³ support material, an **unbalanced global tension structure** has a mass floor of **1,650 kg/m²**. Tiling alone does not remove this force.

In a locally closed cassette, opposing vapor forces on the membrane and outer condenser close through radial ties. Their ideal tensile-material floor is `ρ p D/σallow`: for a 1 m gap and the same 100 MPa allowance, **3×10⁻⁶ kg/m² at 0.1 Pa**, or **0.003 kg/m² at 100 Pa**. Attachments, bending, thermal expansion, buckling and sealing are extra. Gas cannot be allowed to communicate into a geometry that recreates an unbalanced whole-sphere pressure. Direct particle momentum still reaches the receiver and its existing global support; the old gravity/tidal support calculation does not itself certify these transient local loads.

A fully restrained silica temperature change of 500 K produces an elastic thermal-stress scale **21.6 MPa**, well above the small pressure stress of a cooled microtile. Compliant, sliding or preformed mounts are required. Bulk room-temperature tensile/flexural values are not hot irradiated-film allowables. [Corning silica](https://www.corning.com/media/worldwide/csm/documents/1e624b3034474193b5b00ea6f558dd3d2.pdf) has strain/annealing points about **1,166/1,315 K** and softening point about **1,858 K**, making creep relevant well before a nominal melting limit. [Saint-Gobain bulk SiC](https://www.ceramicsrefractories.saint-gobain.com/files/5681/download) supplies hot conductivity/strength data, but a thin defective SiC film needs its own allowable. [Momentive PBN](https://www.momentivetech.com/materials/ceramic-materials/pyrolytic-boron-nitride-pbn) and [Pyroid graphite](https://www.mineralstech.com/docs/default-source/refractories-documents/pyrogenics/technical-data-sheets/physical-properties-for-pyroid-ht.pdf) are strongly anisotropic. Grade, temperature, atmosphere and radiation history matter.

## Ports, punctures and renewal

Prepared inward streams should leave **chamber-facing launcher heads through vacuum-isolated bulkhead penetrations**. Recovered slow material is supplied through sealed valves on the receiver side. The prepared incoming stream then does not traverse the separator film. If a design instead opens launch holes through the film during preparation, vapor leakage, tail overlap and exposed-source heating must be included; the outgoing stopping calculation does not apply to an incoming beam by symmetry without its own trajectory calculation.

A puncture is a localized gas-leak and recovery event. A cassette limits immediate gas inventory and supports isolation, shutters and receiver-side replacement; it does not stop ongoing vapor generation automatically. The gas packet assesses actual plume columns and packet heating. Irradiation-produced pinholes, high-permeability species, seal damage and dust require a leakage criterion rather than a perfect-seal assumption. Spent film must be captured on the **receiver side**, not released into the central chamber.

Renewal is now evaluated in `renewal_energy.py` / `renewal-energy.json` under the user's high manufacturing-efficiency premise. The old universal **100 MJ/kg** renewal cost is a sensitivity, not the selected conclusion. A 1 µm silica film weighs **3.345×10¹⁶ kg**. At all **4,736.9 shots/year**, gross reprocessing is **5.021×10¹² kg/s**; at **2 MJ/kg**, energy is **0.6605 W/m²**, giving the parent's fixed-hardware energy-only result about **2.683 mW/kg** instead of 2.766. That is compatible with the requested specific-power scale before additional qualified hardware and auxiliary costs. This example does not prove that the selected effective process cost, turnaround or film lifetime is achieved.

The same installed replacement factory would have to process routine plus film throughput equivalent to its own mass in roughly **82.4 days**. Local parallel processing, collection, hot-feed reuse and ready spare skins must satisfy both average throughput and the roughly **6,774 s** frequent-stage interval. The renewal packet includes 1/2/3.1/10/100 MJ/kg and different duty assumptions, while keeping its declared installed hardware convention explicit.

Silica renewal has a **Si/O material-return ledger separate from the carbon:H budget**. Gross recycled mass is not fresh import. A net unreturned fraction of 10⁻⁸ in the all-shot 1 µm example would require about **50,200 kg/s of Si/O makeup**, even though almost all material is internally reused. Carbon skins instead consume the already limited structural-carbon reserve in proportion to irreversible escape. Neither material choice is free of inventory and recovery requirements.

The next decisive calculation is a coupled spectral-photon/secondary-electron deposition and transient cooled-tile model, followed by an irradiation/material test or a defensible atomistic damage model. Current calculations favor continued development of supported renewable cells and backed catcher skins; they supply no fundamental prohibition and no qualified survival lifetime. The pressure/cooling algebra passes **252/252** internal checks; the preserved monochromatic reference passes **600** checks. The uncertainty is in the specified physical properties, geometry and operating process, rather than the bare film mass or an assumed universal fabrication cost.

Reproduce from the repository root with `python3 docs/protium-fusion/membrane-receiver/spectral_proxy.py`, `python3 docs/protium-fusion/membrane-receiver/mechanics.py`, and `python3 docs/protium-fusion/membrane-receiver/thermal_reference.py`. These bounded scripts perform no new hydrodynamic simulation and do not modify station or article files.
