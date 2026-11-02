# ✨ Closer collection channels relax the anisotropy requirement

**A 5 µm channel pitch with diffusivity ratio η=625 exactly matches the existing 20 µm, η=10,000 mobile-routing solution**, including its finite 100 nm cap. At η=1,000, a refined 5 µm candidate improves on that comparison: cycle-weighted front escape is **0.17467% H and 0.14508% He**. Lower anisotropy alone therefore does not force a different catcher architecture. Channel density, seals, floor collection and effective pumping speed trade against it.

The existing 20 µm, η=100,000 candidate remains the reference. No frozen record was changed. These calculations concern already-mobile particles and ideal collection boundaries, with the same evaluated stopping-data-derived mean implantation depths and the same qualifications as the [parent finite-channel model](../README.md). They do not establish an irradiated material's diffusion tensor, liberation rate or gas-tightness.

## Exact scaling and two bounded refinements

Write the backward first-passage equation as

```
Dnormal u_zz + Dparallel u_yy = 0.
```

On replacing lateral coordinate y by y/p, the coefficient is η/p². The cap thickness, source depth, rear depth and boundary conditions are unchanged. Thus η′=η(p′/p)² preserves both the collection probability and the exit time at fixed Dnormal. It also preserves the discrete finite-volume solution: all matrix and right-hand-side entries acquire the same lateral scale factor. A direct coarse-grid replay verifies this identity for 5 and 10 µm cells.

| Pitch | η matching 20 µm / 10,000 | η matching 20 µm / 100,000 |
|---|---:|---:|
| 20 µm | 10,000 | 100,000 |
| 10 µm | 2,500 | 25,000 |
| 5 µm | 625 | 6,250 |

The first reference has cycle-weighted front fractions 0.221279% H / 0.181904% He; the stronger reference has 0.103072% / 0.087548%. Thus η≈625–1,000 at 5 µm reaches the **10,000 comparison**, while matching the stronger reference still requires η=6,250 in this ideal geometry.

Only two new candidates were solved, each at depth-resolution parameters 320/640 and half-pitch resolutions 128/256:

| New candidate | H front fraction | He front fraction | Largest stage/species change on refinement |
|---|---:|---:|---:|
| 5 µm, η=1,000 | 0.174670% | 0.145084% | 7.93×10⁻⁶ absolute fraction |
| 10 µm, η=1,000 | 0.402724% | 0.325289% | 9.98×10⁻⁶ absolute fraction |

The 10 µm cell at η=1,000 misses the requested reference's front-loss level; it still has a continuous, moderate performance penalty rather than a numerical failure. For B1, a mean mobile first-exit time below 1 s requires Dparallel above about 2.09×10⁻¹² m²/s at 5 µm, or 8.31×10⁻¹² m²/s at 10 µm. Trapping, temperature dependence and the full temporal distribution remain separate. These numbers are deterministic conditional solutions and numerical refinement differences, not 80% material-property intervals.

## Slits, caps and floor impacts

Keeping slit width at 100 nm raises its projected area from 0.5% at 20 µm pitch to 1% at 10 µm and 2% at 5 µm. The unsupported span remains 100 nm, so its local span-conduction requirement is unchanged. Unsupported cap area, potential through-connected defects above slits, and directly transmitted primary flux toward the enclosed stopping floor scale by two or four. The full continuous cap remains 0.0002 kg/m²; its total exposed area and irradiation do not increase fourfold.

The mathematical channel has zero width. For finite slit projection fslot, the corresponding physical partition is

```
Ffront = (1−fslot) Fsolid + fslot Fslot,
```

where Fslot includes particles stopping in the roof and any return from the enclosed floor stream. [results.json](results.json) records the actual-spectrum straight-CSDA roof-stop fraction and sensitivities to floor return. At 5 µm pitch, a 0.1% / 1% / 10% front return from the transmitted floor stream adds at most 0.002% / 0.02% / 0.2% of the incident flux. A physically sealed, pumped floor can provide this collection path; the zero-width PDE does not qualify that hardware for free. Prompt reflection is an independent component, now screened in the [reflection packet](../../stopping/reflection/README.md).

Replacing only unsupported caps after every event, at an illustrative input of 2 MJ/kg, costs 0.000264 / 0.000529 / 0.001058 W/m² for 20 / 10 / 5 µm pitch. This is a selected maintenance scenario, not a necessary repair method or a manufacturing floor. Replacing the whole face cap has the original area and cost. Holding the total graphite inventory at 350 kg/m² requires only another 2.27 µm of physical backing depth when moving from 20 to 5 µm, beyond the existing via/plenum compensation.

## Conductance and connected-leak requirements

The comparison keeps the inherited 150 µm-deep slits, 50 µm plenum with 25% posts, 0.5 mm-radius vias at 1 cm pitch, and 175 mm backing. It retains the prior approximate diffuse-molecular formulas:

```
fslot = width/pitch
beta_slot = fslot × 8 width/(3 depth)
beta_via  = fvia × 8 rvia/(3 Lbacking)
Cslot = vbar beta_slot/4; Cvia = vbar beta_via/4
1/Seffective = 1/Cslot + 1/Cvia + 1/Spump.
```

The slit conductance rises with channel density, but the fixed vias and pump limit the benefit. Additional lateral-plenum and dogleg resistance remains an apparatus qualification.

| Pitch | Downstream Sp, m³/(s·m² receiver) | Effective H₂ / He speed, m/s | Mean total source pressure, Pa |
|---|---:|---:|---:|
| 20 µm | 0.001 | 0.0008625 / 0.0008165 | 0.001477 |
| 10 µm | 0.001 | 0.0009174 / 0.0008874 | 0.001379 |
| 5 µm | 0.001 | 0.0009476 / 0.0009277 | 0.001330 |
| 20 µm | 0.004 | 0.002442 / 0.002106 | 0.0005374 |
| 10 µm | 0.004 | 0.002941 / 0.002653 | 0.0004397 |
| 5 µm | 0.004 | 0.003275 / 0.003049 | 0.0003909 |

The full H/He throughput is retained for this conservative pressure comparison. Routing changes of a few tenths of a percent have correspondingly small effects on pump throughput and storage inventory. In the immediate-release waveform, the 5 µm / Sp=0.004 source-pressure peaks are A 1.206 Pa, B1 0.01542 Pa, B2 0.8611 Pa and C13 0.6163 Pa. Concentrated later desorption could produce a different peak.

For the same **additional 0.1% plenum leakage budget**, a unit-transmission equivalent connected area must satisfy

```
phi_max = f_leak/(1−f_leak) × 4 Seffective/vbar.
```

At 5 µm / Sp=0.004 this is 4.05×10⁻⁹ H₂ / 5.31×10⁻⁹ He of receiver area. Dividing by the larger slit projection gives a per-cap-area condition **0.949 / 0.934 times the original 20 µm / Sp=0.001 allowance**. Therefore the fourfold warm pumping-speed variant almost preserves the original connected-defect production-times-open-lifetime requirement. At unchanged Sp=0.001, the same condition tightens to 0.275 / 0.284 of reference. These are equivalent leak-area requirements, not a predicted hole yield from NRT displacement counts.

## Supplying the extra pumping speed

A conditional warm booster can supply Sp=0.004 at its inlet while feeding the existing Sp=0.001 downstream line at four times the inlet pressure. Its gas throughput and final downstream load are unchanged. Isothermal compression at 1,000 K has mean reversible work 1.74 µW/m², or 3.47 µW/m² at an assumed 50% free-energy efficiency. The archived immediate-release peak work and the minimum ideal admission areas are recorded in the JSON. This small gas-work floor does **not** establish booster mass or parasitic power: rarefied-gas admission, one-way momentum transfer, seals, drive losses, hot compatibility and peak throughput still need an apparatus. Buffer capacity need not scale fourfold merely because pumping speed does.

A different implementation—directly enlarging the existing cold-sorber interface—has a measurable accounting consequence. With its shared 10⁻⁴-area mouth and molecular transmission 0.1 held fixed, Sp=0.004 needs H-getter / He-sorber active fractions 0.000433 / 0.000528, compared with 0.0000625 / 0.0000487. Under the same full 200 K shield-view and half-Carnot assumptions, He refrigeration rises **0.344→3.739 W/m²**. The 1,000 K optical bypass remains 5.67 W/m² because that mouth is unchanged. This option cannot inherit the warm booster's tiny gas-work floor as its total electrical cost.

The bounded result favors retaining closer channels plus an explicit warm-pump requirement as a live option. It removes η=100,000 as a universal material requirement while preserving the actual qualifications: finite cap/floor collection, hot barrier leakage, diffusion speed, connected-defect repair and delivered pumping conductance.

## Reproduction and preservation

[calculate.py](calculate.py) imports the frozen PDE solver, replays the exact similarity, runs only the two new candidate refinements, and writes [results.json](results.json). It records all source and solver hashes. [verify.py](verify.py) checks the resulting algebra and preservation hashes without rerunning transport. The new [manifest.json](manifest.json) covers this subpacket only; the archived `fine-candidates.json` and all original geometry/gas records are unchanged.
