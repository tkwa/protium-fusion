# ✨ Conditional carbon-recapture honeycomb

A radial graphite honeycomb can recapture most carbon emitted from the catcher, with little added mass. **It does not suppress the pellet's carbon column by the same factor as the escaped atom count.** Its surviving floor emission points toward the reactor centre. For the proposed 10 mm-wide, 100 mm-deep cells, this angular effect limits even perfect local capture to about a **13-fold steady radial-column reduction**, including equally emitting 10 μm wall tips.

The live conditional candidate is a **rough or reactive, renewable inner surface**, with a separate gas-collection path for the directly irradiated tips and walls. It remains optional while the actual MeV irradiation-induced carbon release is unresolved. This packet does not resize or update the selected plant.

## Geometry and local recapture

Square cells have 10 mm clear width, 100 mm axial depth and walls aligned toward the fusion point. Wall thicknesses of 5, 10 and 100 μm add **0.200, 0.399 and 3.941 kg/m²** respectively. Their exposed tip fractions are **0.0999%, 0.1997% and 1.9704%**. The internal side area is approximately **40 m² per receiver m²**. The longitudinal graphite column is 200 kg/m² along a wall tip; the full catcher backing remains behind the open floor.

[calculate.py](calculate.py) traces uniformly distributed Lambertian floor emission in an ideal square cell. At each solid encounter an assigned fraction *s* is permanently incorporated; the remainder reflects diffusely, specularly, or as a specified mixture. Here *s* means capture over the relevant operating interval. An atom that briefly adsorbs and then desorbs is still part of the transport problem. The ray calculation closes every atom as captured or escaped.

The direct floor-to-opening view factor is **0.00316206**, independently recovered by square-to-square quadrature. For 10 μm walls, the following includes bare outward-emitting tips with the same yield per incident ion as the floor:

| Permanent capture per encounter | Reflection of survivors | Escaped carbon number / bare floor | Number suppression |
|---:|---|---:|---:|
| 0.10 | Diffuse | 0.01795 | 56× |
| 0.25 | Diffuse | 0.00768 | 130× |
| 0.25 | Half diffuse | 0.02332 | 43× |
| 0.25 | Specular | 0.09898 | 10× |
| 0.50 | Diffuse | 0.00589 | 170× |
| 1.00 | Any surviving-path law | 0.00516 | 194× |

Each local row uses one million histories. The saved binomial intervals describe sampling uncertainty, not an engineering credible interval. An independent unfolded-path calculation checks the fully specular cases in [specular-check.json](specular-check.json).

There is relevant material evidence, but it does not fix *s* for this irradiated receiver. Bossion et al. calculated thermal-carbon chemisorption of about **0.25–0.30 at 1000–1500 K** on small amorphous-carbon structures. That supports investigating a rough, reactive surface; it does not establish its long-term capture under reactor irradiation. Pristine basal-plane carbon adatoms have much weaker binding in the Lehtinen calculation, so atom incorporation and prompt return must be distinguished. [Bossion et al., 2024](https://www.aanda.org/articles/aa/pdf/2024/12/aa52362-24.pdf); [Lehtinen et al., 2003](https://doi.org/10.1103/PhysRevLett.91.017202).

The companion [getter calculation](../../secondary-ejecta/getter-results.json) shows that chemical stock and positive work are small at the examined carbon flux: a yield of one C atom per incident CNO nucleus corresponds to roughly 64 mg C/(m²·year), and dividing a 1:1 Si or Ti capture inventory across the 40-fold internal area requires only nanometres of annual consumption. Kinetics, irradiation damage, N/O poisoning and regeneration determine whether that chemical capacity becomes useful capture. No measured permanent sticking probability is assigned to a proposed Si/Ti coating.

## Spherical flights and the angular limit

Let θ be a carbon atom's departure angle from the inward radial direction. Its chord has impact parameter *b=R sinθ*. For a steady isotropic distribution of emitting wall locations, the radial-column contribution contains

\[
\int_b^R\frac{dr}{r^2\sqrt{1-b^2/r^2}}
=\frac{\arccos(b/R)}b.
\]

Relative to one bare cosine-emitted spherical flight, the angular weight is

\[
A(\theta)=\frac{\pi/2-\theta}{2\sin\theta}.
\]

Its cosine-emission average is one, but nearly radial atoms have much larger weight. The direct survivors of the 10:1 cell have **22.945 times** the average weight per escaped atom. Their unavoidable first-flight column contribution is therefore **0.072554 per floor-emitted atom**, despite their number fraction of only 0.003162.

Direct-path angular quadrature gives these lower bounds before returning flights or tip emission:

| Depth / clear width | Direct number escape | Direct radial-column factor |
|---:|---:|---:|
| 10 | 0.0031621 | 0.072554 |
| 30 | 0.0003534 | 0.024593 |
| 100 | 0.00003183 | 0.007417 |

These additional aspect ratios are geometric bounds only. Their wall alignment, gas routing and thermal design have not been substituted into the selected candidate.

[cavity.py](cavity.py) follows successive flights through a spherical cavity and back into another cell. It preserves the incidence angle on each spherical return, randomizes microscopic cell phase and azimuth, and applies the same local sticking/reflection model. It therefore does not assume every returning carbon atom is captured, or that the cavity instantly isotropizes the outgoing plume. The table uses isothermal carbon, full hood coverage, 10 μm tips and equal tip/floor yield:

| Permanent capture | Survivor reflection | Steady radial column / bare cosine source with unit return capture | Column reduction |
|---:|---|---:|---:|
| 0.10 | Diffuse | 0.1721 | 5.8× |
| 0.25 | Diffuse | 0.1051 | 9.5× |
| 0.50 | Diffuse | 0.08257 | 12.1× |
| 1.00 | Perfect capture | 0.07441 | 13.4× |
| 0.25 | Half diffuse | 0.1936 | 5.2× |
| 0.25 | Specular | 0.5431 | 1.8× |

An exact direct-flight control variate reduces sampling noise. The s=0.25 diffuse row has sampling standard error about 0.0011. These are steady transport factors; actual pulse arrival peaks, changing emission temperature, molecular species and dust need their corresponding time and speed kernels. The earlier adverse radiation-enhanced-sublimation diagnostic cannot simply be multiplied by the escaped-number fraction.

## Primary-ion interception and gas routing

The nominal radial beam strikes the tips as well as the floor. Saved stopping ranges put maximum CNO path lengths near **43 μm**, H near **88 μm**, and B1 He near **11 μm**. A 100 mm-long wall has ample longitudinal range, while a 10 μm wall does not guarantee lateral containment of recoils or charged secondaries. The thermal screen retains all intercepted nominal energy and assigns no unsupported escaping fraction. Material that exits a wall laterally may still stop in the retained backing or another wall.

For beam-axis offsets θx and θy, the geometric transmitted fraction is

\[
\frac{[p-H|\tan\theta_x|]_+[p-H|\tan\theta_y|]_+}{(p+t)^2}.
\]

For a total tilt shared equally between the two axes, the intercepted fraction is:

| Tilt | 10 μm walls | 100 μm walls |
|---:|---:|---:|
| 0 | 0.1997% | 1.9704% |
| 10⁻⁵ rad | 0.2138% | 1.9843% |
| 10⁻⁴ rad | 0.3408% | 2.1090% |
| 10⁻³ rad | 1.6061% | 3.3518% |

The point-source angular size is small compared with these sensitivities, but wall bow, local tilt and alignment must be included. The intercepted fraction is also the upper added H/He cavity-source fraction if every intercepted light ion re-emits and the wall has no collection path. Actual contribution is that fraction times release probability times failure to route gas into the closed return line. The original 150 μm floor laminate cannot silently be applied inside a 10 μm wall; a thinner wall-specific collection laminate, enclosed edge channels, or a separately demonstrated gas-routing solution is needed before claiming its benefit.

Global warm-pump inlets are part of this interface. Placing a deep honeycomb over a mouth adds molecular resistance. Leaving mouths untreated creates exposed carbon and H/He sources. The incoming primary-ion baffle faces must carry the existing local collection laminate or an equivalent closed extraction path. Carbon sources on those baffles must be recaptured in their actual recess/collector geometry, or retained explicitly as an exposed area fraction. For example, an untreated 10% area with the same radial-cosine yield contributes 0.1 of the bare first-flight column by itself. The H/He collection, carbon interception, optical path and pumping conductance must be solved jointly for that mouth. No extra global-pump performance is credited here.

## Axial heat path and temperature

[thermal.py](thermal.py) solves a finite **100 mm** graphite fin with k=50 W/(m·K), density 2000 kg/m³ and cp=1800 J/(kg·K). Its bottom follows the saved actual receiver-front history. Its tip receives the full nominal material kinetic/internal energy and delayed heat per projected tip area, plus a separate prompt-energy input. There is no inverse-tip-area amplification of the incident flux. Side radiation is omitted from this axial-conduction reference; both cooling and heating by nonlocal radiation exchange require a coupled enclosure calculation.

The calculation starts 20,000 s before each selected event, removing uncertainty about a guessed hot initial fin. The constant-mean receiver pickup and actual dispatch are retained. The highest material-phase tip temperature is **1475.66 K** for C13, or **1486.81 K** with an additional 100 μm damaged-tip layer at k=4. The corresponding A/B2 values are about **1350.47/1393.39 K** before that damaged layer. The existing finite-depth Teff prompt proxy gives microsecond tip peaks of approximately **1210–1281 K** across the large-pulse cases; its photon/electron deposition qualifications remain unchanged. No prompt peak is held for the duration of a material pulse.

Energy residuals are below 5×10⁻⁹, and refining 100 to 200 axial cells changes the major material peaks by at most 0.11 K. The B1 post-prompt tail changes by 0.26 K. Saved time histories are linked in [thermal.json](thermal.json), so temperature-dependent yields can be integrated over actual flux rather than held at a peak.

Tilt can shift appreciable energy into the long thin sidewalls. A separate [full-retained-side-energy sensitivity](side-thermal.json) uses the actual hottest C13 waveform and uniform axial deposition over the intercepted wall family. A beam tilted along one grid axis loads only one wall family; dividing that energy across both would understate its temperature. This is an energy-allocation screen, not a prediction of lateral ion escape or radiation exchange. Its results are reported as temperatures, not as a universal failure threshold.

For the 10 μm walls, the conditional side-loading peaks are **1484–1487 K at 10⁻⁵ rad**, **1557–1590 K at 10⁻⁴ rad**, and **2340–2731 K at 10⁻³ rad**, spanning diagonal and single-axis tilts. With 100 μm walls, even the 10⁻³ rad screen gives about **1556–1590 K**, while the nominal exposed-tip H/He fraction is ten times larger. These conduction-only high-temperature values do not establish actual overheating: thin-wall secondary escape and radiative redistribution between walls, floor and cavity can materially reduce the retained energy and temperature gradient. A coupled radiation/conduction calculation is the next numerical refinement if this geometry is pursued.

The basic tip-only thermal path is plausible under the stated conductivity. The stronger conditions for this optional mitigation are **durable carbon incorporation, controlled reflected angular distributions, wall-specific H/He extraction, and the actual integrated pump-mouth geometry**. All recaptured carbon remains internal material inventory. Gross redeposition and renewal do not become fresh-carbon consumption.
