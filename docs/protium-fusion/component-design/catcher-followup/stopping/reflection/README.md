# ✨ Prompt hydrogen and helium reflection from the graphite catcher

The normal-incidence screen does **not identify a reflection-driven redesign** of the delivered-B1/C13H2 catcher. It does identify a correction to the gas budget and a surface-angle qualification. With primary low-energy transport tables and a factor-two high-energy escape-model sensitivity, first-reflection fractions are **0.0380% of incident H and 0.0127% of He**. The respective mean reflected energies are **29.4 and 26.5 keV**. These particles initially form an energetic return population inside the closed sphere. Their eventual contribution to warm cavity gas depends on their next impacts and subsequent material transport.

These are calculations for a smooth, normal-incidence carbon surface. They are neither a measured hot-catcher retention efficiency nor a statistical confidence interval. The current laminate's surface orientation, defects, stored gas and fields still need the relevant material/transport qualification. No frozen stopping or plant record was changed.

## Inputs and evidence

The calculation uses every H and He packet in the archived [surface spectrum](../surface-spectra.json), at its actual stage cadence, for the coherent delivered-B1/C13H2 family at radius 1.1×10⁹ m. The graphite backing is 350 kg/m² at 2,000 kg/m³. The first reflection occurs within the surface stopping region; backing thickness establishes transmission containment but does not establish retention.

The primary low-energy comparison is the official [EIRENE TRIM database](https://www.eirene.de/old_eirene/html/body_data_eirene.html), with its [data-format documentation](https://www.eirene.de/old_eirene/trim.pdf). The newer H-on-C table uses a 1 eV projectile surface-binding parameter and He-on-C uses 0 eV. Both cover 1–5,000 eV and incidence angles 0–85°. Their numerical reflection coefficients and conditional exit-energy/cosine quantiles are extracted in [primary-tables.json](sources/primary-tables.json), alongside the older tables. They are binary-collision transport calculations for randomized targets, not measurements of the proposed irradiated laminate. Energy quantiles in the raw tables are **eV**, converted to fractions only in our extract.

For He, the normal-incidence table on printed page 51 of [Eckstein, IPP 9/132 (2002)](https://pure.mpg.de/rest/items/item_2138250_2/component/file_3587319/content) extends the comparison to 20 keV. Its factual values are in [eckstein-He-normal.json](sources/eckstein-He-normal.json). The two primary calculations differ modestly at their 5 keV overlap. Source URLs and full raw-file hashes are preserved; copyrighted reports and page images remain only in the scratch download directory.

Above these table domains, [single_collision.py](single_collision.py) implements a screened Coulomb single-hard-collision model using the [Geant4 Wentzel–Molière prescription](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/elastic_scattering/ion_coulomb.html), exact two-body recoil and archived NIST stopping. The prescription's quoted lower energy window is about 50–100 keV/u. Our low-energy results therefore receive an explicit primary-table comparison; the model alone is not a validation of that region. At 1–5 keV, the primary normal-reflection values exceed this model by factors 1.39–1.69 for H and 2.04–2.70 for He.

## Escape calculation and numerical scope

Let x be traversed mass column, R(E) the CSDA mass range, r=m/M, and μ the center-of-mass scattering cosine. A backward laboratory trajectory has μ<−r, with

```
k = (1 + r² + 2rμ)/(1+r)²                 outgoing/incoming collision energy
u = −(μ+r)/sqrt(1+r²+2rμ)                  outward laboratory cosine
Ei = R⁻¹[R(E0)−x]
escape requires R(k Ei) > x/u
Eexit = R⁻¹[R(k Ei)−x/u]
dP = [dx/(12.011 mu)] 2π(dσ/dΩcm) dμ       x expressed in kg/m²
```

The code integrates only the column and angles that satisfy escape, rather than counting all elastic collisions as reflections. It includes average electronic and nuclear stopping on both legs. That mean-loss treatment does not resolve the correlation between the selected hard collision and the rest of the collision cascade. Multiple small deflections, straggling, channeling, porous microgeometry and atomic surface dynamics are additional limitations. Doubling angular and depth quadrature gives agreement within 0.2% for fourteen H/He benchmark cases; this checks integration, not the physical approximation.

[compare.py](compare.py) uses the primary coefficients only inside their tabulated incident-energy domain. Above 5 keV H / 20 keV He it multiplies the escape model by explicit factors 1, 2, 3 and 10. The factors diagnose the consequence of a model error; they are not source-provided uncertainty curves or an 80% interval. No low-energy primary rate is silently extended to MeV energies.

| High-energy model multiplier | H first-reflection fraction | He first-reflection fraction |
|---|---:|---:|
| 1 | 0.02265% | 0.007966% |
| 2 | 0.03800% | 0.01272% |
| 3 | 0.05334% | 0.01748% |
| 10, stress case | 0.16076% | 0.05076% |

These are cadence-weighted number fractions. Contributions from the primary table domains alone are 0.007310% H and 0.003211% He. The unresolved incident population below 1 eV is at most another 9.14×10⁻⁸ H / 7.04×10⁻⁸ He of the total. The older primary tables change the factor-two result to 0.03647% H / 0.01257% He. All stage and energy-resolved values are in [comparison.json](comparison.json).

The factor-two per-stage H fractions are A 0.03559%, B1 0.05216%, B2 0.009069%, C13 0.01038%. Almost all of the cycle's He number flux is in B1; its factor-two value is 0.01272%. The rare He populations of the other stages remain explicitly included.

## From reflection to the gas account

A first-reflected atom crosses a spherical chord of length 2Ru and strikes the wall again at the same macroscopic incidence cosine. The single-collision model gives first-flight times averaging roughly 2–5 ks, with B1 He at 5.10 ks. It predicts tens of keV mean exit energies, not immediate thermal accommodation at 1,000 K. Its first-flight H/He inventory is only 1.53×10⁻⁵ / 5.67×10⁻⁶ of the existing warm-inventory comparison, using the adopted H₂/He residence times and baseline leak fractions. This inventory comparison does not establish a spatially uniform energetic-particle density or its pellet-heating rate.

Sparse primary exit-energy quantiles plus the high-energy model give a factor-two **≤1 eV exit-population bracket** of 4.67×10⁻⁷–8.30×10⁻⁶ H and 8.04×10⁻⁸–1.96×10⁻⁵ He per incident atom. For ≤100 eV the brackets are 1.96×10⁻⁵–3.02×10⁻⁵ and 2.95×10⁻⁶–2.28×10⁻⁵. These are interpolation/model brackets, not confidence bounds; the broad He upper bracket assigns its entire 5–20 keV primary contribution to the low-energy exit bin because that table lacks quantiles. Even 1 eV is above the characteristic energy of 1,000 K gas. Sub-eV surface binding, thermal target motion and the incident charge state require an atomic treatment before these proxies become a prompt thermal-gas prediction.

The useful robust bookkeeping statement is conditional on a given first-reflection fraction: **the fraction of incident atom lineages that ever reaches warm cavity gas through first reflection cannot exceed that fraction**, regardless of subsequent bounces. Repeated desorption of the same atom affects residence and net pump capture, rather than creating additional atoms. Gas-model g must mean net removal from this reservoir, as it already does in the warm-cleanup acceptance condition. Sputtered pre-existing H/He and carbon are separate source terms.

Consequently, assigning the entire factor-two reflected branch to eventual warm gas gives conservative additive allowances 0.0003801 H / 0.0001273 He, including the sub-eV incident uncertainty. Adding these to the rounded existing mobile-route leaks gives approximately **0.001410 H / 0.001002 He**. This addition deliberately overcounts any overlap with the original unreflected mobile-route prediction. It increases that comparison by about 37% / 15%; it does not multiply the total fuel inventory or imply physical escape from the sphere.

The gas follow-up's additional common thermal-source allowance is about 0.00591–0.00712 under its specified 2 cm fast-packet / g=0.01 acceptance conditions. The whole reflected branch is well inside that allowance in the normal-incidence factor-one-to-three cases, and also in the factor-ten stress case. That is a conditional interface comparison, not a new validation of g or packet survival.

A diagnostic in [single-collision.json](single-collision.json) assumes the next hit stops at u times the mean projected range and applies the existing mobile side-channel PDE. It gives eventual front-release contributions 1.03×10⁻⁴ H / 3.81×10⁻⁵ He from the single-model first-reflected branch. This is a one-rehit sensitivity: it does not follow second reflection or range straggling, and is not used to claim the smaller gas source.

There is no universal finite bounce-count bound for a perfectly smooth, grazing trajectory. If each later encounter has an independently established permanent-capture probability c, the mean number of later encounters is at most 1/c and the survival after n encounters at most (1−c)ⁿ. A chosen c is a requirement, not a result here. Particle energy is dissipated at each real encounter, but macroscopic grazing reflection can retain much more energy than a normal backward collision.

## Remaining high-impact qualifications

**Local incidence angle:** at 5 keV the primary H / He coefficients increase from 0.0197 / 0.0194 at normal incidence to 0.162 / 0.175 at 60° and 0.473 / 0.484 at 80°. These energy-specific values cannot be applied to the whole MeV spectrum. They do show why the normal result cannot qualify arbitrary rough facets, exposed channel walls, oblique baffles or magnetic deflection. For an unknown fraction a of the incident flux on such surfaces, the conservative additional first-reflection allowance is ≤a; a more useful allowance needs its actual angle and energy distribution. A blanket sub-per-mille reflection claim is unwarranted without that geometry.

**Nuclear elastic scattering:** measured proton–carbon backscattering has resonances and departures from Coulomb scattering; see the [0.7–2.5 MeV measurements](https://www.sciencedirect.com/science/article/pii/0168583X85904987) and the [Paneta et al. 2.7–7 MeV measurement](https://www.sciencedirect.com/science/article/abs/pii/S0168583X11005787). This packet does not substitute a Rutherford prediction for their evaluated differential data. As an intentionally generous importance screen, the actual full stopping column gives 4.61×10⁻⁵ H / 2.14×10⁻⁵ He events per incident particle for each barn of an assumed column-averaged backward cross section, counting every event as escape. Adding 0.1% by this mechanism alone would require 21.7 / 46.6 barn under that overcount. These are thresholds to compare against data, not asserted cross-section caps. Coulomb–nuclear interference requires a consistent full differential kernel rather than casually adding probabilities.

**Surface state and secondary release:** the idealized cap's crystallographic orientation, implantation damage, pores and H loading can change reflection and release. The primary tables support the scale of the present screen; they do not establish hot-laminate performance. Secondary ejecta of stored atoms are addressed by the separate source study. A validated BCA calculation with the actual surface geometry and an evaluated MeV elastic kernel, followed by repeated-impact transport or a relevant measurement, would close the main reflection uncertainty. No requirement for a new whole-sphere barrier or full EM braking follows from the present normal-incidence result.

## Reproduction

Use the project's Python 3.10 numerical environment. Run `single_collision.py`, then `compare.py`, then `verify.py` in this directory. `extract_tables.py /path/to/raw-download-cache` rebuilds the factual source extract from the four named official downloads. The extraction is not needed for the numerical rerun because the compact factual data are archived. [manifest.json](manifest.json) hashes this new subpacket only. Raw reports are in `output/catcher-followup/reflection/` outside the worktree.
