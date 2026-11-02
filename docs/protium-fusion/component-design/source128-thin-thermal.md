# ✨ Finite backing at 500 and 750 kg/m²

**The actual source128 pulse family supports 500 kg/m² as a conditional receiver candidate.** Relative to 1000 kg/m², it adds about 57 K to the conservative material/debris temperature envelope. Its early isolated surface pulse is unchanged. The important engineering changes are a larger bulk storage swing and smaller, independently compliant receiver panels. Photon spectrum, surface contact and fatigue remain component qualifications; a deliberately arbitrary soft spectrum does not establish a required receiver temperature.

This packet uses the coherent A-H3/B1-H2/B2-H2/C13-H4 **source128 carbon-median** records in `source-family-carbon-median128/integration.json`, including their actual endpoint composition, mass, velocity, residual internal energy and prompt-energy histories. It does not relabel nominal128 or the older64 spectra. Radius is 1.1×10⁹ m. For a controlled comparison, all thicknesses retain the source128 2000/.99 capacity point's helium rate, 4.48764×10³²/s, and **B1 mean period 7230.449 s**. These are not new optimized thin-receiver plant points.

`source128_thin_thermal.py` and `source128-thin-thermal.json` reconstruct the arrival histories and backing response. `source128_prompt_proxy.py/json` replay the existing finite-depth photon calculation on these traces. `source128_prompt_radiance.py/json` test the thermal scope of its assigned colors. `source128_thermal_combine.py` and `source128-thermal-combined.json` retain the separated combination rules. The independent identities and source checks are in `source128_thermal_verify.py` and `source128-thin-thermal-validation.json`.

## Finite heat storage and pickup

The constant-property comparison uses graphite density 2000 kg/m³, heat capacity 1800 J/kg/K and through-thickness conductivity 50 W/m/K. Conductivities 25 and 100 are saved sensitivities. All source heat is retained in the backing for this upper thermal screen, including heat that a thinner plant may actually route into cold equipment or allow to escape. The model consequently receives a mean **125.881 W/m²**; it does not silently grant 99% hot capture to a thin receiver.

| Quantity | 500 kg/m² | 750 kg/m² | 1000 kg/m² |
|---|---:|---:|---:|
| Physical thickness | 0.25 m | 0.375 m | 0.5 m |
| Through-thickness diffusion time L²/α | 1.25 h | 2.81 h | 5.0 h |
| Conservative sum of bulk storage swings | 152.6 K | 101.7 K | 76.3 K |
| Isolated C13 slow surface rise, including its own bulk mode | 418.1 K | 418.1 K | 418.1 K |
| All-phase slow surface envelope, bulk floor 1000 K | 1534 K | 1496 K | 1477 K |
| Extra bulk floor needed to keep pickup ≥1000 K | 0.315 K | 0.472 K | 0.629 K |
| Extra storage allowance for ≤1 day dispatch lateness | 12.1 K | 8.1 K | 6.0 K |

The single-pulse equality follows because the early thermal disturbance has not reached the rear even in the 0.25 m slab. Thinning changes the accumulated bulk temperature, rather than exposing a new immediate surface-temperature divergence. The actual source128 maximum K+advected-U flux is about **0.955 MW/m²**. Each pulse's reconstructed slow energy is more than 99.99999% complete by the 20,000 s integration cut; the remaining kinetic, internal-energy and decay tails are still bounded separately.

The calculation uses a finite slab with a uniform energy-storage mode and exact nonuniform conduction modes. Heat leaves through a **constant mean-power back-face pickup**. The per-stage periodic bulk-energy range is bounded by one full pulse energy divided by slab heat capacity. For old nonuniform heat, large events have a one-day minimum separation; B1's actual 5269.82 s minimum main-window slot bounds possible closer spacing under dispatch jitter, while its mean power uses the 7230.449 s cadence. The large-pulse one-day guard consumes 21.6% of available time. A completed event dispatcher still has to satisfy the stated guard and any chosen lateness bound.

There is no fixed 1000 K bath. The back temperature is a positive delayed average of the finite bulk temperature, offset by −q_mean L/(2k). Raising the chosen bulk floor by the small amount in the table therefore keeps pickup at least 1000 K while allowing it to vary. The positive kernel and its mean delay L²/(6α) are independently checked. Mean-rated conversion equipment can use this stored heat; instantaneous megawatt-per-square-metre heat extraction is not assumed.

## Prompt radiation and its physical scope

The prompt flash occurs tens of seconds before the material peak. Its source128 fluences are 20.09/14.22/5.59/68.12 kJ/m² for A/B1/B2/C13. The finite-depth replay preserves those exact energies and time bins. For C13 at the bulk reference properties, the outer-temperature Planck proxy gives a 54.2 K isolated rise; the flux-equivalent temperature gives 248.8 K. Their corresponding B1 rises are 38.1 and 54.5 K. They remain assigned spectra on a grey radiation history, rather than a frequency-dependent transfer solution.

A thermal source with assigned color T_c, luminosity L and saved outer radius r needs a dilution factor

`W = L / [4πr² σ_SB T_c⁴] = (T_eff/T_c)⁴ ≤ 1`.

The outer-temperature proxy meets W≤1 in every saved emitting bin. The flux-equivalent proxy has W=1 by definition. At the same area and luminosity, **quartering the outer temperature does not generally define a possible diluted thermal photosphere**: C13's energy-weighted W is 8.18, and 64.0% of its prompt energy requires W>1. B1 requires W>1 for all its prompt energy, with weighted mean W=72.2. That alternative can be retained as an unconstrained spectral-shape sensitivity, or would require different emitting area/reprocessing/non-LTE physics. Its much hotter surface result is not imposed as a physical receiver requirement or called a credible worst case.

The following conservative envelopes use the physically scoped flux-equivalent proxy, taking the larger of prompt and slow peaks for each pulse, then allowing independent big-event and B1 contributions plus stored/background heat. They include the pickup-floor correction and the explicitly bounded one-day dispatch-lateness allowance:

| Surface requirement | 500 kg/m² | 750 kg/m² | 1000 kg/m² |
|---|---:|---:|---:|
| Uniform bulk-property reference | **1600 K** | **1557 K** | **1537 K** |
| 100 μm skin at k=4 W/m/K, bulk retained at k=50 | **1632 K** | **1590 K** | **1569 K** |

The skin row adds a conservative finite resistance for material plus delayed surface heat and recalculates only the soft prompt tail with the lower conductivity; it does not assign damaged-skin conductivity to the whole backing. Thermal contact impedance remains separate. The quarter-color sensitivity and the older all-energy-at-the-surface flash are explicitly labeled diagnostics in the JSON. No narrow credible interval for an uncomputed spectrum is inferred from those alternatives.

## Panels, strain and joints

Thinner backing needs more freedom to bend. For a freely expanding homogeneous plate, the one-dimensional temperature profile gives

`κ = 12 α_T/L³ × |∫(x−L/2) T(x) dx|`.

At the selected α_T=5×10⁻⁶/K, the C13 pulse gives curvature magnitudes 0.00509/0.00255/0.00153 m⁻¹ for 500/750/1000 kg/m². The independently checked small-slope midpoint bow κs²/8 for a 5 m panel is **15.9/8.0/4.8 mm**. If the mount allows 10 mm one-axis bow, the corresponding calculated spans are approximately **4.0/5.6/7.2 m**; allowing 100 mm gives 12.5/17.7/22.9 m. These are conditional acceptance calculations, not measured optimum panel sizes. At 100 m spans, the linear extrapolation already enters large-deflection behavior and cannot qualify a rigid panel merely because its pressure loading was acceptable.

A metre-scale tiled hot face with sliding or flexural attachments is therefore a concrete candidate. Cold global prestress remains in the separate supporting frame. Surface expansion, small panel bow and common shell motion can be accommodated without assigning all of them to random launch error. The illustrative free-plate front stress scale at E=10 GPa and ν=0.2 is about 18–22 MPa for the slow pulse; strength, fracture, transient photon gradients and irradiated joints still require their own model. A present graphite datasheet lists 12.5 GPa dynamic modulus, 4.1×10⁻⁶/K expansion over 20–200°C and 65 MPa flexural strength; these provide scale, not qualification at 1000–1600 K. [Manufacturer data](https://www.sglcarbon.com/pdf/SGL-Datasheet-SIGRAFINE-R6650-EN.pdf).

Sliding joints need overlapping or staggered shielding and insulation. For scale, open 5 mm seams on a 5 m square grid expose about 0.2% area and would leak roughly 113 W/m² from an unobstructed 1000 K black surface. That is incompatible with the retained 1 W/m² hot-to-cold bypass; preserving line-of-sight closure is part of the panel design. Neither gross panel replacement nor ablation is automatically fresh-material loss: recovered material stays in the recycling account, with processing throughput and energy charged explicitly.

## Validation and decision

The independent finite-volume operator converges to the analytic nonuniform kernel within 0.0062% at 1024 cells. Halving the actual slow-pulse time step changes the peaks by at most 0.037%. All 104 source128 photon-proxy checks pass. The new source hashes, finite pickup identity and curvature moment identities pass. Independent reconstruction of both archived source128 and nominal128 integrated plant families closes mass, electricity, nuclear heat, radiator and pulse-period identities to 2.3×10⁻¹⁶ relative.

This screen does **not** make 750 or 1000 kg/m² a necessary thermal floor. It supports carrying 500 kg/m² with finite storage, qualified bulk conductivity, compliant small panels, shielded sliding joints and a shallow-skin/contact specification. Coupled photon/electron transport, the irradiated property/contact law, and the final dispatch/pickup implementation remain the meaningful qualifications. The larger source/nuclear uncertainties and the complete plant's auxiliary allocations are unchanged by this thermal comparison.
