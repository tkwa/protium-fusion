# ✨ B1 needs its own collision-entropy program

The consequential new B1 sensitivity is a change in a **prescribed hot-flow profile**, not a universal temperature tolerance. The parent's full-network fixed-energy runs found that multiplying all starting temperatures by 1.2 reduced gain from 7.45–7.55 to about 0.36–0.45 for the H1.5/H2 contenders. Multiplying by 0.8 retained gains near 7.2. Because B1 supplies most pulses, its preparation path deserved a separate discriminator.

This packet supplies that discriminator. A controlled constant normal-collision energy maps the H2 feed to a productive hot profile at unchanged total prepared energy. The original-rate calculation completes at gain **7.52115075**, 99.6757% of the 7.54562254 reference. Its median-carbon-rate rerun completes at **7.52074959**, with C12 yield **0.753445692 per initial N15**. The calculations remain spherical hot-flow responses to a thermodynamically constructed input; they do not establish a delivered three-dimensional grain assembly.

## Why the A/C13 angular prescription does not transfer

The exact B1 normalization is 1e6 kg/m³, reference radius 30 m and pulse mass 1.130973355e11 kg. The prescribed inflow's actual radii and temperatures vary by shell; the initial outer radius is about 191 m. The fixed physical interval is 82.2779579 microseconds. Using the A/C13 mass normalization here would be a factor-100 error.

The chosen cold reference is isotope-enriched N15:H2 material corresponding to hydrazine's atom ratio. Density 1062.5 kg/m³ is the illustrative 1000 kg/m³ natural-material density multiplied by 34/32 under equal isotope molecular volume; 250 K and a 2 GPa reference bulk modulus are chosen inputs. They are not measured properties of an irradiated sealed N15-hydrazine capsule. Density alternatives 850 and 1200 kg/m³, and a homogenized H1.5 comparison, are retained in `entropy-map.json`. FEOS16.7 supplies the nonideal free energy; its Thomas–Fermi/Cowan model does not predict molecular fragmentation or capsule strength.

Using the same parallel-sector N=2000 proxy as the earlier A/C13 work would give B1 H2 mapped T/reference mass 10th/50th/90th percentiles **0.0144/0.0263/0.0435**. N=50 instead gives **0.303/0.691/1.210**, still a very broad entropy distribution. A count of angular regions is therefore a poor standalone specification for this stage.

Matching the reference entropy requires approximately **1.40–1.46e11 J/kg** shock heat across the nonseed parcels. Its near constancy suggests a simpler control variable: use approximately the same **transverse collision velocity** for each mass parcel, while allowing the radial velocity and crossing angle to vary. This is intentional kinetic heating during assembly. A physical design can pair streams from different source locations; the collision angle near the small assembly region is not the angle through which a final launcher must steer away from the reactor centre.

## Constructive constant-contact case

One contact-shock energy was selected from the average energy needed to approach 0.9 of the reference starting temperature. The actual resulting temperature profile is computed from common entropy, rather than being clipped to a temperature floor or ceiling:

| Quantity | Constructed H2 value |
|---|---:|
| Contact energy deposited | 1.2394264e11 J/kg |
| Fraction of whole prepared pulse energy | 4.3837% |
| Total collision heat | 1.40060e22 J |
| Normal/transverse speed of each equal opposing stream | 497.9 km/s |
| Relative transverse speed of the pair | 995.8 km/s |
| Contact density | 4452.7 kg/m³ |
| Contact temperature | 12.8424 MK |
| Contact pressure | 3.4592e14 Pa |
| Mapped T/reference, mass 10/50/90% | 0.8632 / 0.9014 / 0.9246 |
| Compensating hot radial-speed multiplier | 1.0063605 |
| Mass undergoing post-shock rarefaction to its prescribed hot density | 17.60% |

The inner approximately 0.1% seed keeps its previous prescribed temperature and remains separately assembled. At the mapped density profile, reduced internal energy is paid back to convergent kinetic energy, preserving the exact global input. The remaining parcels must compress or rarefy from the contact state as required; describing all of them as undergoing further compression would be wrong.

For a symmetric pair, `q = v_normal²/2`. An illustrative local energy allocation equal to each parcel's eventual hot internal plus radial kinetic energy gives collision half-angles with mass 10/50/90% **9.38°/12.71°/17.96°**. Real pressure work exchanges energy between parcels, so these angles are a geometric/energy scale, not a solved inverse launch trajectory or a required pointing precision.

The controlled original-rate run completes the full interval with energy residual −5.99e−7 of prepared input and maximum temperature 0.52666 GK. The updated-rate run has residual +5.24e−7 and essentially the same temperature. No numerical guard is reached. The original-profile comparison and the separate failed high-temperature sensitivities remain preserved.

## Independent warm-state check and early burning

`contact_primary_check.py` compares the contact state against an additive-volume mixture of the archived primary first-principles N and H tables at equal pressure/temperature, using N15 inertial mass. FEOS pressure is **4.77% higher**, while internal energy above the respective elemental cold references differs by **0.170%**. This supports the scale of the warm contact state; it is not an independent molecular hydrazine EOS or a validation of its entropy/phase trajectory. Elemental tables and cold chemical reference differences remain explicit in the saved result.

At the contact state, the original bare N15(p,alpha) plus N15(p,gamma) rates give an N15 consumption rate about **6.32e−16 s⁻1** in the chosen mixture. A homogeneous 1 ms residence therefore consumes only about **6.32e−19** of N15. Contact burning is negligible in that thermal rate screen. Nonthermal colliding ions, later compression and actual plasma radiation transport are separate parts of assembly. FEOS includes equilibrium radiation; the packet does not calculate escape from an initially porous collision region.

## Pulse timing and source rating

The existing prepared-hot-flow launch-bin convention gives a peak electrical rating of **1.889e21 W** for the reference and **1.914e21 W** for the mapped profile at 60% preparation efficiency: **+1.33%**. Both have the same full electrical pulse energy, **5.3250e23 J**.

A second deliberately explicit timing proxy assigns each nonseed parcel its hot internal-plus-radial energy as cold kinetic energy and uses `launch_time = −R / sqrt(2E/m)`. This produces a monotone speed sequence, launch times from about **1201.5 to 294.5 seconds before contact**, and peak midpoint-bin electrical power **2.172e21 W**, about **15.0% above** the reference hot-flow convention. The omitted seed retains its separately budgeted source. This is a useful source-power scale, not a circuit, grain arrival distribution or complete longitudinal phase-space solution.

The collision heat itself is focused stream energy. Depositing it over chosen 10/100/1000 microsecond assembly intervals would correspond to local conversion powers 1.40e27/1.40e26/1.40e25 W. Those examples are not an additional external heater or demonstrated merger duration. A physical collision/rarefaction trajectory must determine when and where that existing energy is thermalized.

## Evidence boundary and files

- `entropy_map.py` / `entropy-map.json`: stage-specific cold-reference and fixed-angular-count comparisons, exact reference profiles and FEOS hashes.
- `constant_collision.py`: generates the constant-shock-entropy hot profile, restores global prepared energy, and runs the full 27-isotope model.
- `B1-H2-constant-contact-T0p9.{json,initial.json,inputs.json}`: unchanged original-rate result and inputs; `source-path.json` points to the raw trace and records hashes.
- `contact_primary_check.py` / `contact-primary-check.json`: independent elemental-table mixture screen.
- `audit.py` / `audit.json`: 14 stored-trace checks plus six preparation/timing/contact checks, all passing.
- Median-rate raw result: `/home/tkwa/src/protium-fusion/output/b1-source/updated/B1-H2-constant-contact-T0p9-rates-updated-median.json`; it will be included in the separate consistent four-stage source-family map.

The unresolved component is the finite three-dimensional collision geometry: grain occupancy, first contact, shock junctions, pressure work, post-shock expansion/compression, radiation retention and seed coupling. The result shows a useful entropy program with explicit kinetic heating and bounded source-power arithmetic. It does not establish a universal 20% thermal allowance, a fixed number of independent directions, or a qualified capsule/launcher.
