# ✨ Molecular recirculation work in the hot pump

The uniform-area thermalized blade train needs **0.33–0.38 W of electrical input per receiver m²** for its moving blades and fully diffuse flat tip faces, at an assumed 90% motor efficiency. This uses the conservative global benchmark of **0.001 Pa forepressure per species**, H₂ and He simultaneously, giving 0.002 Pa total. The actual local primary-buffer mean requires **0.014–0.016 W/m²** under the same conventions. These molecular loads belong inside the proposed total motor allocations of 1 W/m² and 0.13 W/m². Bearings, drives, exposed disk surfaces and physical interstage reset hardware consume the remainder.

[calculate.py](calculate.py) reduces the existing million-particle blade trajectories; [results.json](results.json) preserves every pressure node, row area, work contribution, assumptions and input hashes. The code passes 349 checks, including independent dense pressure-node solutions. No new Monte Carlo trajectories or selected-plant changes were made.

## What is being counted

For a row of total projected area Aᵢ, open fraction f, transmission probabilities Tf and Tr, and equal-temperature reservoirs pᵢ and pᵢ₊₁,

    Q/(vbar/4) = Aᵢ f (Tf pᵢ − Tr pᵢ₊₁).

Q has pressure-volume-flow units, equivalent to watts, but is not the pump's shaft power. With K=Tf/Tr and N rows, the variable-area resistance is

    Z = sum_i [K^(−i)/(Aᵢ f Tf)],  i=0…N−1;
    p₀ = pN/K^N + [Q/(vbar/4)] Z.

A warm entrance baffle adds its own series resistance. Each row's shaft power is instead

    Aᵢ f vbar/(4 kT) × (pᵢ Wforward + pᵢ₊₁ Wreverse).

The saved W values include every collision within that row, including molecules that return to their entry reservoir. Multiplying by reflection counts again would double-count them. Forward mean work can be negative: some incident thermal molecules return work to the rotor. The pressure-weighted sum for the operating train is positive.

A stationary reset surface receives the exiting gas energy without shaft work. The saved molecular ledger separates moving-wall mechanical work and wall heat exchange. Open-row reset heating plus moving-wall heat deposition equals the shaft input. This supports thermal recovery at the actual surface temperatures; the calculation does not assign a plant-wide recovery fraction.

For a planar axial blade-tip face moving tangentially at speed u, specular reflection transfers no tangential work. Fully diffuse reflection has mean shaft work m u² per incident molecule. The separate tip term charges both exposed faces over projected fraction (1−f). Inclined edges, extra rotating annuli and roots require their own surface areas and kernels.

## Uniform-area numerical result

All rows use 45° tilt, f=0.65, axial length/pitch=1, T=1000 K and u=700 m/s. The common train length is selected by the H₂ compression: at least 10⁹ for global cleanup, or 10⁶ for the local primary booster. Helium uses the same physical row count.

| Blade diffuse probability | Global rows | Global open-blade shaft W/m² | Added fully diffuse flat tips W/m² | Global molecular electrical W/m² at ηmotor=0.9 | Local rows | Local molecular electrical W/m² at ηmotor=0.9 |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 68 | 0.2330 | 0.1084 | 0.3794 | 46 | 0.01554 |
| 0.5 | 56 | 0.2132 | 0.08983 | 0.3367 | 38 | 0.01380 |
| 1 | 57 | 0.2070 | 0.09200 | 0.3322 | 38 | 0.01367 |

Global face area is 0.244834 m² per receiver m². Both species have chamber partial pressure 10⁻¹⁰ Pa and forepressure 0.001 Pa, with the existing baffle conductance included. The resulting H₂ global uptake g is 0.0119–0.0157 and He uptake is 0.0145–0.0184; this ideal row network retains margin over g=0.01 before additional stator/clearance losses.

The local face is 1.984255×10⁻⁵ m² per receiver m². Its exact archived record-mean buffer partial pressures are 0.626985 Pa for H₂ and 0.377652 Pa for He. Actual mean booster throughputs are recovered from cold-outlet flow plus the change in buffer inventory. These give the mean work directly because the fixed-geometry pressure network and work are linear. The resulting source pressure is recomputed with the new row capture and existing channel conductance; the old empirical capture is not silently imposed on these rows.

Net isothermal gas-compression work is only about 3.0–3.9×10⁻⁸ W/m² globally and 8.7–9.4×10⁻⁶ W/m² locally. Those small numbers do not include repeatedly incident gas. Approximately 95–99% of open-row molecular shaft power occurs in the final ten rows.

At 90% motor efficiency, the uniform-area molecular result leaves about 0.62–0.67 W/m² inside the 1 W/m² global total-motor target and 0.114–0.116 W/m² inside the 0.13 W/m² local target. These are remaining allowances. If the local 0.13 W/m² were instead interpreted as a fixed idle draw inferred from a present pump, the molecular load would have to be added to that idle draw. The saved result also evaluates 50% and 100% motor-efficiency scenarios.

## Smaller high-pressure stages

Keeping the first 16 rows at full area and reducing each subsequent row's exposed moving area by 10% produces a constructive ideal-network comparison. Global molecular electrical demand falls to 0.0021–0.0063 W/m²; local demand falls to 0.0010–0.0019 W/m². H₂ global uptake remains 0.0119–0.0157. The last global stage has 0.0042–0.0148 times the first-stage area. A stronger 20% taper is included as a sensitivity, without selecting it as hardware.

This gain requires physically smaller exposed moving areas and controlled pressure compartments. A restriction drilled into a full rotor exposed to the high-pressure reservoir retains that rotor's drag. Keeping the same local tip speed can use a smaller annular span at comparable radius, or appropriately driven separate rotors. Interstage stator transmission, sealing, bypass and rotor geometry remain to be designed. The full-area global train should end at its 0.001 Pa foreline: extending the same exposed face into a roughly 1 Pa hot reservoir would raise its dominant molecular-work term by roughly three orders of magnitude. A separate small-area booster avoids that mistake.

## Collision and uncertainty scope

For a stipulated effective collision area σ, the diagnostic mean free path is kT/(√2 pσ). At σ=10⁻¹⁹ m², it is about 48.8 m at the global total forepressure 0.002 Pa, and 0.0976 m at 1 Pa. The packet varies σ from 10⁻²⁰ to 10⁻¹⁸ m². These are transport scenarios, not measured H₂–He cross sections.

Combining each row's archived mean flight time with its adjacent total pressure gives roughly 0.016–0.017 expected collisions in the most exposed local row at 1 mm pitch and σ=10⁻¹⁹ m². At 1 cm pitch the estimate is ten times larger. Local source transients and larger effective cross sections can therefore require collisional transport, particularly in the high-pressure stages. The [primary Holweck DSMC study](https://arxiv.org/abs/physics/0401028) illustrates why free-molecular performance cannot simply be carried into a collisional foreline.

The independent 100,000-particle discovery samples change the global open-row work by −0.11% to −0.59% relative to the million-particle records. Archived work variances and covariances are unavailable, so this is a numerical repeat comparison rather than a formal credible interval. The three scattering cases form a conditional model bracket. No narrow subjective interval is assigned to a fabricated pump whose stators, clearances, surface kernels and exposed rotating areas are still unspecified.

Reproduce with Python, NumPy and the archived dependencies:

    python3 calculate.py --out results.json
