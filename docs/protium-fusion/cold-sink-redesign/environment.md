# Shaded reactor radiative environment ✨

**Direct solar heating is zero for this redesign**, following Thomas's verified shading premise. The remaining site term is infrared irradiation from the Dyson swarm, local shields, machinery and neighbouring radiators. A 197–200 K radiator is physically compatible with some of these geometries; shade alone does not determine its useful heat rejection.

## What the canon fixes

The illustrative Dyson-swarm environment assumes independently orbiting collectors in layered families, roughly 90% solar-energy capture, and radiator orientations chosen to limit heating of neighbours. Collector orbital radii remain unspecified; a 0.4 AU example is an illustrative calculation rather than an adopted radius.

Consequently, there is no canonical warm enclosing shell at the 5 AU reactor, nor a specified angular infrared sky. An inner source viewed from outside, a more extensive swarm surrounding the site, and local protected radiator directions should be treated separately. The accepted zero-sunlight condition does not require changing that canon.

## Equations and area convention

For one open emitting face,

\[
q_{\rm net}=\epsilon_r\sigma T_r^4-q_{\rm abs,IR},\qquad
q_{\rm abs,IR}=\sum_j\bar\alpha_r(T_j)F_{rj}J_j.
\]

Here F is the **cosine-weighted radiative view factor**, not the fraction of all 4π steradians covered by an object; J is its emitted-plus-reflected radiosity. A black source has J=σT⁴. The examples use gray ε=αIR=0.9. A spectral design can have different Planck-mean absorptivities for hot and cold sources, with absorptivity and emissivity still related at each wavelength and direction. NASA's thermal-control treatment explicitly includes incoming infrared radiation, view factors and the radiator's own emission in the heat balance. [NASA thermal control](https://www.nasa.gov/smallsat-institute/sst-soa/thermal-control/)

For identical-temperature gray panels which see other such panels over Fself, and fixed external irradiation Hext over the remaining directions, the scalar radiosity solution is

\[
q_{\rm net}=
\frac{\epsilon_r[(1-F_{\rm self})\sigma T_r^4-H_{\rm ext}]}
 {1-(1-\epsilon_r)F_{\rm self}}.
\]

This accounts for reflected radiation as well as direct self-view. A finite shield at an unknown temperature requires another thermal balance rather than an assigned cold radiosity.

All capacities below are **per square metre of emitting face**. A two-sided panel supplies twice its planform area only when both faces are included with their actual views. For plant accounting,

\[
Q_{\rm cold}+Q_{\rm IR,collector}+Q_{\rm conducted}+Q_{\rm other}
=\sum_i A_{\rm emit,i}\,q_{\rm net,i}.
\]

If external IR is intercepted by a separate collector of area Acol, its contribution is qabs,col Acol; dividing that power by Aemit gives its contribution per radiator square metre. It must not be multiplied by Aemit as though every radiator square metre had intercepted the collector's mean flux, or charged again after already appearing in the hot-side power balance.

## External Solar-system source benchmark

For an isotropic IR luminosity fIR L☉ produced inside the reactor's orbit, the inward-normal irradiance at 5 AU is **54.4467 fIR W/m²**. In the small-angular-source limit, a spherical receiver intercepts through πR² and emits through 4πR², giving a mean **13.6117 fIR W/m²**. Thus

\[
\bar q_{\rm abs}=13.6117 f_{\rm IR}\alpha_{\rm IR}
\quad\mathrm{W/m^2}.
\]

At fIR=0.9 this is **12.2505 W/m² incident**, or **11.0254 W/m² absorbed** for αIR=0.9. This is a useful nonzero comparison point, conditional on most captured power being reradiated inside 5 AU and approximately isotropically. Energy delivered and dissipated outside that orbit or directional radiator arrangements change it. Global energy conservation makes waste-heat routing necessary; it does not require an isotropic local sky. [Wright, radiative thermodynamics of Dyson spheres](https://arxiv.org/abs/2309.06564)

| Radiator temperature | Gross emission, W/m² | Net with external spherical mean, W/m² | Net on source-facing flat face, W/m² | Net inside 300 K black enclosure, W/m² |
|---:|---:|---:|---:|---:|
| 197 K | 76.86 | 65.84 | 32.76 | −336.51 |
| 200 K | 81.65 | 70.63 | 37.55 | −331.72 |
| 250 K | 199.35 | 188.32 | 155.25 | −214.02 |
| 300 K | 413.37 | 402.34 | 369.27 | 0 |
| 450 K | 2,092.69 | 2,081.66 | 2,048.59 | 1,679.32 |

The external columns use fIR=0.9 and ε=αIR=0.9; the mean column presumes the same spherical receiving/emitting surface. The black-enclosure column is a distinct boundary-condition stress case. A negative number means environmental heat must be removed before the radiator can carry a plant load. At 197 K the spherical-mean comparison increases emitting area by about **17%** relative to zero IR; a facing flat radiator increases it by about **135%**.

The CMB contribution at 2.725 K and α=0.9 is **2.81×10⁻⁶ W/m²**. Unspecified local structures dominate any practical uncertainty at these radiator temperatures.

## Direction and source angular size

For a uniform bright source sphere with angular radius θ=asin(Rs/D), a face aimed at its centre has F=sin²θ. A tangent face has

\[
F_{\perp}=(\theta-\sin\theta\cos\theta)/\pi.
\]

A face aimed directly away has zero view of that source, provided it lies wholly behind the panel plane. A spherical receiver's exact mean F is (1−cosθ)/2, approaching sin²θ/4 only for a small source. The following are geometric examples at the same 0.9 L☉ IR luminosity, not selected swarm radii:

| Source radius | Effective black outward temperature | Angular radius at 5 AU | Exact spherical mean incident IR | Tangent-face incident IR |
|---:|---:|---:|---:|---:|
| 0.4 AU | 606.2 K | 4.59° | 12.27 W/m² | 0.83 W/m² |
| 1 AU | 383.4 K | 11.54° | 12.38 W/m² | 2.11 W/m² |
| 2 AU | 271.1 K | 23.58° | 12.78 W/m² | 4.38 W/m² |
| 4 AU | 191.7 K | 53.13° | 15.31 W/m² | 10.90 W/m² |

Directional panels can therefore give low coupling to an inner swarm without requiring its luminosity to disappear. Their support geometry, thermal transport, two-face accounting and mutual occlusion must support the stated views. Blocking the solar disk covers a much smaller angular region than blocking an extended collector/radiator swarm.

## Enclosing shells, shields and self-view

A black enclosure fills every outward hemisphere. Its irradiation is σTb⁴, independent of distance to the wall; the net exchange is εσ(Tr⁴−Tb⁴). Passive 197 K operation requires Tb<197 K. A 300 K all-sky boundary would prevent such operation.

That 300 K boundary must itself have an energy-consistent construction. A complete black shell radiating only 0.9 L☉ outward, with equal inner and outer temperatures, reaches **171.46 K at 5 AU**, **156.52 K at 6 AU**, and **121.24 K at 10 AU**. The corresponding 197 K net capacities are **32.76, 46.24 and 65.84 W/m²**. The 5 AU value is a limiting shell radius; an actual enclosing shell must lie beyond the reactor. A warmer all-sky interior would require unequal inner/outer temperatures, restricted radiative escape, additional power, or some other explicit thermal structure. It should not be inferred solely from the word “Dyson.”

A shield removes the view of whatever lies behind it and replaces that view with its own radiosity. For a 197 K radiator to retain at least 80% of clear-sky capacity, the allowed view of a black warm surface is:

| Warm surface | Maximum cosine-weighted view factor |
|---:|---:|
| 300 K | 0.0372 |
| 450 K | 0.00735 |
| 1,000 K | 0.000301 |

For example, a 1,000 K collector covering F=0.001 supplies **51.03 W/m² absorbed**, already two thirds of a 197 K radiator's gross capacity. A 450 K shield covering F=0.01 supplies 20.93 W/m². Low-emissivity surfaces can reduce these terms, but their reflected irradiation and eventual heat rejection remain in the balance. An enclosing opaque shade cannot be treated as a window to empty space.

Mutual radiator view also consumes useful escape area. At 197 K, ε=0.9 and otherwise cold sky, Fself=0.1 gives **69.88 W/m²**, and Fself=0.5 gives **40.45 W/m²**, compared with 76.86 W/m² for an unobstructed face. A protected cold-facing radiator field is the geometry that supports low IR most directly: the central swarm and hot collector remain behind its emitting faces, while its visible supports and shields are sufficiently cool or occupy small view factors.

## Bounded result and reproduction

Set direct sunlight to zero in the revised plant calculation. Retain an explicit IR term. The external-source spherical mean is an appropriate comparison benchmark when no detailed directional layout has been selected; low-IR cases require the stated view geometry. An enclosing warm structure is a separate stress case. The canon supplies neither a compulsory warm bath nor a zero-IR sky, so no probability distribution over these unknown geometries is assigned here.

[environment.py](environment.py) writes [environment.json](environment.json), including 0–400 W/m² absorbed-heat cases at all five radiator temperatures, external source orientations, enclosing-shell balances, shield limits and source hashes. Run with Python and NumPy:

```sh
python environment.py
```

The bounded calculation passes 38 checks: normal, tangent and away-facing angular integrals for 12 source cases, zero exchange at equal black-bath temperature, and heating by a warmer enclosure. No old plant packet or canon file is changed.
