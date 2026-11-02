# ✨ Receiver250 component extension

The coherent source128 family supports **250 kg/m² as a conditional receiver candidate**. This calculation supplies no physical reason to impose500 kg/m² as a floor. The thinner plate needs greater thermal motion accommodation and a larger local guide shield/buffer. Its nominal20.6247 versus20.0544mW/kg advantage over500 belongs to the coupled transport/plant account; this component packet keeps the existing1000kg/m² service allowance explicit.

## Preserved loading and finite backing

The new wrapper reruns the frozen finite-backing calculation at a0.125m graphite thickness, with no new hydrodynamics or source tuning. The exact family is A-H3/B1-H2/B2-H2/C13-H4,128 cells,current carbon median,constructed-source profiles. The faster source2000/.99 cadence is retained: B1 interval7230.449s, mean all-source heat125.881W/m². Current coupled250 and500 rows have slower cadence. All source heat is retained for this thermal-storage screen, independently of the thinner receiver's actual hot-capture fraction.

The model has constant mean backside extraction, finite bulk heat capacity and a nonuniform conduction response. It retains the one-day large-pulse separation assumption and gives an additional explicit bound for up to one day of dispatch lateness. It neither fixes the backface to1000K nor assigns the engine instantaneous pulse power.

| Quantity, k50W/(mK), cp1800J/(kgK), density2000kg/m³ |250kg/m²|500kg/m²|
|---|---:|---:|
|Backing thickness|0.125m|0.250m|
|Bulk stored-temperature swing bound|305.10K|152.55K|
|Isolated C13 full-front rise|418.135K|418.135K|
|Teff-proxy all-phase front envelope, including pickup floor and lateness|1728.25K|1599.70K|
|Same, with100µm k4 stopping skin|1759.97K|1631.42K|
|Free-plate C13 curvature atα=5×10⁻⁶/K|0.014609m⁻¹|0.005095m⁻¹|
|One-axis bow over5m|45.65mm|15.92mm|
|Span giving10mm one-axis bow|2.34m|3.96m|
|Illustrative free-plate residual front stress, E10GPa/ν.2|12.85MPa|18.20MPa|

The combined temperature is a conservative envelope over independent stage phases, stored backgrounds and separate prompt/slow maxima. It is not a predicted simultaneous temperature history. The thermal Teff photon shape obeys the blackbody-radiance relation at the saved emitting area. The older quarter-outer-temperature spectrum remains an unconstrained non-LTE shape diagnostic and is not imposed as a physical flash requirement. The separate k25/k100 cases test backing conductivity; a shallow damaged skin is not assigned to the whole plate.

The main short pulse has diffusion depth smaller than either backing, explaining the equal isolated early peaks. Thickness matters more to stored energy and curvature. At250, the source-specific0.25→0.125s time refinement changes the maximum slow response by at most0.01143%; the inherited1024-cell finite-volume kernel check remains within0.0062%.

## Compliant panels and material scope

Three-to-five-metre sliding or compliant panels are compatible with this screen. Their one-axis peak bows are16.4–45.7mm and edge slopes about1.26–2.09degrees. These values set mount geometry; a1cm allowable is an example, not a requirement. A simply supported5m strip under0.2Pa has only about240Pa bending stress and1µm pressure deflection at the chosen modulus. Thermal deformation dominates that local pressure load.

The thinner free plate relieves more of its thermal gradient by bending, which lowers the reported residual stress while increasing bow. A mounting system cannot simultaneously claim free-bending stress relief and enforce zero motion. Nonlinear biaxial tile compatibility, irradiation, creep and cyclic joints remain unmodeled. At5m the slope-squared geometric strain scale is roughly6.7×10⁻⁴, sufficiently close to thermal strains to merit a nonlinear plate/joint calculation. Smaller panels improve this continuously; no arbitrary rejection threshold is introduced. Overlapping or staggered sliding seams must preserve shielding and the hot/cold radiation barrier.

The [SGL R6650 datasheet](https://www.sglcarbon.com/pdf/SGL-Datasheet-SIGRAFINE-R6650-EN.pdf) gives room-temperature modulus12.5GPa and flexural strength65MPa, expansion4.1×10⁻⁶/K over20–200°C and room-temperature conductivity95W/(mK). These values establish useful scales. They do not certify irradiated graphite at the computed transient temperatures. The selected k50,cp1800,α5×10⁻⁶ assumptions remain explicit material requirements.

## Prestress and guide protection

At the current coupled250 EM20+.005 point, M=3.14760×10²²kg. The retained20kg/m² cold frame requires0.75919MN/m pretension,0.006074 biaxial prestrain and189.80MPa stress in its assumed8kg/m² equivalent load-bearing material. This frame is distinct from the hot graphite tile. All three stream families remain positive. Their total0.044483kg/m² exceeds the scalar station inventory by only0.000886kg/m²; that small correction is explicit. The30/100/300-day global affine attitude screens remain stable. Local flexible-guide control and actual joints are separate qualification conditions.

The existing local2000kg/m² radial graphite guide shield, with the measured65–70K polyethylene heat capacity and the same geometry convention, gives:

|Per receiver area|250|500|
|---|---:|---:|
|Local graphite shielding|136.113kg/m²|121.303kg/m²|
|Polyethylene buffer|24.253kg/m²|16.542kg/m²|
|Combined shield and buffer|160.365kg/m²|137.845kg/m²|
|Additional refrigeration, EM20+.005|0.987W/m²|0.592W/m²|
|Additional refrigeration, all-prompt-cold+.005|1.010W/m²|0.633W/m²|

Guide pulse sizing conservatively sends the entire saved prompt pulse into cold surroundings, including the EM20 comparison. Whole delayed heat receives the existing D1.1/local-cascade attenuation screen. It retains the measured calorimetric5% reserve. Heat accessibility, conductor quench, secondary transport and radiation lifetime remain acceptance conditions. Half-accessible heat capacity and a larger angular-interception factor are saved sensitivities: jointly they require452.77kg/m² shield+buffer and3.35–3.43W/m² additional refrigeration at250. Those are scenario results, not80% credible bounds.

## What the1000kg/m² services actually contain

Services now account for48.31% of the250 plant's installed mass. The following named allocations fit inside that existing allowance exactly once. The older helium loop screen carried a larger hot load than the current point and the same cold load, so its hardware is conservatively retained.

|Named service requirement|kg/m²|
|---|---:|
|Hot/cold helium headers and gas inventory|99.116|
|Extra loop-interface exchangers above separately retained station exchangers|1.778|
|Cold frame package|20|
|Guide coils, cryostats, trim/control and mechanical package|50|
|Local graphite guide shield and polyethylene buffer|160.365|
|Named subtotal|331.259|
|Unassigned remainder at1000|668.741|

Source acceleration/storage, station conversion and internal exchangers, the100kg/m² insulation, radiator panels, factory capacity and installed reserve retain their separate budgets. They are not charged again here. Conversely, the unassigned remainder is not evidence that recovery machinery, tile mounts, joints, valves or control hardware has been designed. Existing named packages themselves retain material and performance requirements.

A400/700/1000/1500kg/m² service sensitivity leaves68.74/368.74/668.74/1168.74kg/m² after this named subtotal. A200kg/m² alternative guide package adds150 to the subtotal. There is no completed bill of materials supporting a unique lower service allocation, and no defensible probability model for a central estimate or80% credible interval. These are conditional whole-plant sensitivities with explicit residual hardware budgets, not four equally demonstrated designs.

At the250 EM20 point, retained ordinary auxiliaries provide4.74565W/m². Conservative external pumps0.74368, baseline guide cryocooling0.31946 and new shielding/buffer refrigeration0.98747 use2.05061W/m², leaving2.69504 for other processing/control loads. A second subtraction of those same itemized loads would double-count them. If the remaining equipment exceeds its allowance, the plant must be recalculated.

## Reproduction and limits

`source128_receiver250_thermal.py` records its exact substitutions into the preserved thermal solver. `source128_receiver250_components.py` combines the separate photon proxy and runs the existing pure support/buffer equations. New JSON outputs preserve source hashes; prior packets are unchanged. The saved validation checks algebraic closure, curvature, capacity and support results. The actual forthcoming delivered-B1/C13-CH2 family must receive a separate pulse-envelope comparison before inheriting these bounds.
