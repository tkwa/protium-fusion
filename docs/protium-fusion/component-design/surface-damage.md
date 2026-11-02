# ✨ Graphite stopping-skin damage and restoration

The charged-particle damage is shallow enough that it does **not** justify replacing the metre-scale receiver's 50 W/m/K thermal conductivity with a badly damaged surface value. It does justify treating the first 50–100 μm as a separately qualified, renewable surface with measured thermal contact. The present data do not establish an annual lifetime. They also do not establish a need to rebuild the surface after every pulse.

This packet uses actual wall-arrival species and velocities, the existing measured electronic stopping plus ZBL elastic stopping, and actual cycle cadence. It keeps two families separate:

- `surface-damage.json`: selected64 A3/B1-H2/B2-H2.5/C13-H4 from `plant-scan-003.json`.
- `surface-damage-current128.json`: archived median-carbon128 A3/B1-H2/B2-H2/C13-H4, using `mixture-design/current-rate128-sigma2000.json` and its receiver caches. This is the 99% capture comparison family. Surface exposure depends on its pulse rate and spectra, rather than on the metre-scale absorber's assigned areal mass.

The script `surface_damage.py` is parameterized by plant and receiver cache paths. Hashes preserve the exact source records, traces, stopping tables and EOS. No hydrodynamics was rerun.

## Exposure and damage proxy

For each projectile, the existing stopping table gives

\[
D(E)=\int_0^E {S_n(E')\over S_e(E')+S_n(E')}\,dE'.
\]

The elastic energy transferred inside a finite layer is `D(E_in) − D(E_exit)`, with the exit energy obtained from the remaining CSDA range. The new depth bins identify path length with depth in a planar normal-incidence approximation. Range straggling, angular scattering, roughness and the spatial extent of recoil cascades remain uncomputed.

The reported NRT-style proxy is `0.8 E_elastic/(2 E_d N_target)`. Using **all** elastic recoil transfer as damage energy omits the recoil's subsequent electronic loss, and no recombination or thermal annealing is subtracted. This deliberately simple exposure measure is not a count of surviving defects, a conductivity prediction, or a lifetime law. Electronic excitation damage is a separate mechanism. The original [Norgett–Robinson–Torrens method](https://doi.org/10.1016/0029-5493(75)90035-7) estimates displacement production from damage energy; its use here is explicitly a screening approximation.

| Quantity | Selected64 | Current median-carbon128 |
|---|---:|---:|
| CNO incident fluence per year | 3.456×10²¹ /m² | 3.394×10²¹ /m² |
| Total elastic transfer per year | 42.91 MJ/m² | 41.72 MJ/m² |
| NRT proxy per year, average first 50 μm, E_d=25 eV | 0.855 | 0.831 |
| Average first 100 μm | 0.427 | 0.415 |
| Average 1–5 μm band | 7.85 | 7.37 |
| Average 50–100 μm band | 3.79×10⁻⁵ | 1.78×10⁻⁴ |

The apparent low average across 50 or 100 μm should not hide the first few micrometres. For the current family, the depth-bin proxy is 5.09/year at 0–1 μm, 7.37 at 1–5 μm, 1.15 at 5–10 μm and 0.0759 at 10–25 μm. The much smaller dose farther in supports a thin-layer treatment. It does not establish a zero-damage backing: penetrating photons, secondary electrons, rare nuclear products and long-term material evolution require their own qualification.

Changing the electronic stopping multiplier from 1 to 0.5 or 2 changes the selected64 50 μm average to 1.48 or 0.480/year. Changing E_d from 25 to 12 or 60 eV rescales it by 25/12 or 25/60. These are explicit sensitivity choices, not 80% credible intervals. Molecular-dynamics work gives about 24 eV for pristine graphite, declining with accumulated structural disorder; that supports testing a range rather than treating E_d as a material constant for life. [Vuković et al., 2018](https://doi.org/10.1103/PhysRevApplied.10.064040)

Published heavy-ion experiments demonstrate that a shallow damaged layer can have much lower conductivity than the underlying material: Prosvetov et al. inferred a decrease from approximately 95 to 4 W/m/K in the irradiated region after intense Au/U exposure. Their ions and irradiation conditions differ from this reactor's hot CNO exposure. We use 4 W/m/K as a useful severe surface scenario, not as an empirical prediction for this plant. [Prosvetov et al., 2020](https://doi.org/10.1016/j.actamat.2019.11.037)

High operating temperature helps some defects anneal, but an automatic full reset between pulses is unsupported. A neutron-irradiated graphite experiment found substantial, dose-dependent recovery on heating up to 1700°C; its specimens were originally irradiated below 200°C. That establishes annealing as a real mechanism without specifying the equilibrium condition during continuous 1000 K operation here. [Maruyama and Li, 2019](https://doi.org/10.1080/00223131.2019.1633966)

## Thermal resistance of the damaged layer

The current128 maximum ballistic kinetic-plus-advected-internal-energy flux is 1.054 MW/m², below selected64's 1.173 MW/m². The maximum remains the C13 pulse. Applying the entire incoming material heat at the front of the skin conservatively overstates the resistance of volumetric stopping within the layer. For a well-contacted skin the added temperature scale is

\[
\Delta T_{skin}\le q_{material,peak}\,d/k_{skin}.
\]

At density 2000 kg/m³ and heat capacity 1800 J/kg/K:

| Skin | Conductivity | Current128 full skin drop | Skin diffusion time |
|---|---:|---:|---:|
| 50 μm | 4 W/m/K | 13.2 K | 2.25 ms |
| 100 μm | 4 W/m/K | 26.3 K | 9 ms |
| 50 μm | 1 W/m/K | 52.7 K | 9 ms |
| 100 μm | 1 W/m/K | 105 K | 36 ms |

The material pulse lasts many seconds, so this local resistance approximation is appropriate. Delayed radioactive heat and prior-pulse background remain in the larger receiver thermal model; they are not silently included in the above K+U peak. Their slow, distributed heating has to be added if it is routed through the same local skin.

A concrete engineering acceptance criterion is **no more than 50 K added material-pulse temperature across the skin and its interface**. This is a selected margin, not a fundamental threshold. For the current128 100 μm skin it requires k_skin≥2.11 W/m/K if contact resistance is negligible. At k_skin=4, contact resistance must remain below 2.24×10⁻⁵ m² K/W, corresponding to contact conductance above about 44.5 kW/m²/K. Thermal cracking or delamination can violate this requirement before the homogeneous conductivity does. A continuous graded surface or reworked monolithic graphite avoids inventing a perfect detachable-film bond, but its actual irradiated interface still needs a design.

Thus the useful maintenance signal is surface thermal impedance, together with cracking, roughness and collection performance. The present nuclear-energy proxy supplies an exposure coordinate for that qualification. It cannot turn 0.83 average dpa/year into a justified annual replacement rule.

## Prompt heating is a different calculation

The photon flash precedes the material pulse. Low conductivity must not be applied to the old all-energy-at-the-surface flash estimate as if the full receiver had become insulating. Real photon absorption has finite depth.

The saved selected64 C13-H4 finite-depth packet permits a bounded sensitivity without a new radiation model. With cp=1800 J/kg/K and k_skin=4, its resolved-energy no-conduction term plus its separately maximized sub-keV surface term gives:

- outer-temperature Planck proxy: about 25 K;
- flux-equivalent-temperature proxy: about 216 K;
- quarter-outer-temperature proxy: about 465 K.

Only the sub-keV surface-conduction term scales with inverse square-root effusivity. The resolved-depth no-conduction temperature scales with inverse heat capacity and does not depend on conductivity. These remain the earlier homogeneous first-interaction model alternatives, with ideal contact to a better-conducting backing. They are not certified spectra or an exact layered transient. Their different weighted mean temperatures do not fully order their soft tails, which explains why the quarter-temperature proxy is more demanding than the flux-equivalent proxy. See [mixture-radiation.md](mixture-radiation.md) and its NIST sources.

The current128 C13 grey photon fluence is 54.27 kJ/m² versus selected64's 61.19 kJ/m², a ratio of 0.887. Current128 B1's fluence rises from 10.65 to 12.80 kJ/m², while A and B2 decrease. A lower integrated fluence alone does not prove a lower transient peak or a harder spectrum. The saved low-k prompt rows remain labeled selected64 reference; they are **not relabeled as current128 radiation results**. Updated source timing and spectrum belong in the separate radiation qualification.

## Restoration is recycled processing

A 100 μm carbon skin contains 0.2 kg/m². Returning that material to a useful surface can involve local smoothing, annealing, moving a layer, or full reworking; it does not inherently require atomization or fresh carbon. The table deliberately prices the more demanding full-layer reworking interpretation at the existing **chosen complete process cost of 2 MJ/kg**. The user's manufacturing premise permits efficient conversion of the required free-energy increase; there is no extra blanket dissociation charge. Higher 3.1 and 10 MJ/kg sensitivities are saved separately.

| Full 100 μm reworking schedule | Recycled mass /m²/year | Average electrical work | Factory allocation at existing throughput coefficient |
|---|---:|---:|---:|
| Annual | 0.20 kg | 0.0127 W/m² | 0.634 kg/m² |
| Quarterly | 0.80 kg | 0.0507 W/m² | 2.54 kg/m² |
| Monthly | 2.40 kg | 0.152 W/m² | 7.61 kg/m² |
| Weekly | 10.44 kg | 0.661 W/m² | 33.1 kg/m² |
| Daily | 73.05 kg | 4.63 W/m² | 231 kg/m² |

A 50 μm layer halves every entry. The factory coefficient is the existing 10⁸ kg of installed factory per kg/s of processing capacity, not a fundamental manufacturing limit. Installed factory mass does not shrink merely because the process is efficient. Neither are the above allocations all additional by default: the current128 2000 kg/m² plant already carries routine reworking capacity of 0.426 kg/m²/year and 1.351 kg/m² of factory. Annual 100 μm renewal would use about 47% of that capacity; quarterly would exceed it. Other component maintenance still needs its share. A monthly surface schedule requires an explicitly reassigned or enlarged processing allocation even though its 0.152 W/m² is only about 0.29% of the current 53.34 W/m² export.

All such material remains inside the collection/recycling system in this accounting. Actual escape or irreversible nuclear conversion is a separate material-flow term. The much larger catalyst-release and isotope-processing flow must also retain its existing account; this surface-maintenance calculation does not make that process free.

The concrete revision is therefore to retain **50 W/m/K as the bulk-property requirement, add a 50–100 μm degraded-skin/contact requirement, and describe annual restoration as a budgeted scenario rather than an established lifetime**. Monthly or weekly reworking is an affordable contingency in the present energy ledger, with finite capacity and machinery allocations. None of these calculations forces per-shot surface replacement. A measured or simulated hot-CNO defect/contact evolution law is the next discriminator for cadence; a full stopping cascade and layered photon/thermal calculation would resolve the remaining deposition qualifications.

## Verification and scope

`surface_damage_verify.py` checks all source hashes, annual accumulation, and twelve independent comparisons of finite-skin transfer with the previously existing membrane routine. Elastic energy and NRT40 agree to better than 7×10⁻¹⁶ relative. The incident baryon and smooth-energy integrations close to better than 8×10⁻¹⁵. This verifies bookkeeping and units, not the omitted material physics. No narrow credible interval for restoration lifetime is justified by the available evidence; the conductivity, stopping and cadence tables are conditional engineering scenarios.
