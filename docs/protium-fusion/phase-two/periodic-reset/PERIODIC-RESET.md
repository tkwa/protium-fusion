# Periodic whole-inventory reset ✨

The saved 140 T plant can be given a concrete positron-clearing cycle: lower both ion end barriers, drain each channel, remove the residual leptons after the ion density falls, and refill it. The main solenoid stays energized. This calculation asks whether the energy and converter mass of such a cycle fit the existing plant allocation. It does not establish the time-dependent drain, recovery, or refill dynamics.

With the reference 95% random-ion preparation efficiency, the specified prompt-pair supply and 10–100 MeV/s per trapped positron, **none of the tested separated-efficiency cases reaches 1 mW/kg**. At the most favorable supplied pair rate and positron power, a case with 95% random-energy extraction and 99.99% coherent transfer produces a conditional **0.661 mW/kg**, about **0.994 TW** from **1.505e15 kg**. Its optimum burn lasts **2,167 s**, with 0.1 s assigned to each channel's reset. Higher positron supply or heating often makes even the optimized cycle energy-negative.

This can support an experimental operating mode or a lower-performance deployment in a future history, subject to the other engineering assumptions. It does not establish the original 1 mW/kg allocation, or eliminate continuous isotope processing during the long burn.

## Stored energy and charge

For each of the fifteen active processors, the saved radial equilibrium is integrated separately. The results are scaled by that processor's assigned channel length. In terms of number densities n_H and n_C, catalyst charge Z and mass M, common temperature T, angular velocity omega, and wall-relative potential phi:

- Thermal energy is 3 kT (N_H + N_C) / 2.
- Rotation energy is the integral of (m_H n_H + m_C n_C) omega² r² / 2.
- Self-electric energy is the integral of rho phi / 2. The independent integral of epsilon_0 E² / 2 agrees within 8.4e-8.
- Possible 40 MV end-transfer circulation is 40 MV times total ion charge. It is tabulated separately; a lowered-gate refill does not automatically incur this transfer.

| Inventory | Full network | Divided by primary rate and one MeV |
|---|---:|---:|
| Random thermal energy | 2.014e16 J | 30,105.7 MeV s |
| Ordered rotation | 2.001e16 J | 29,916.7 MeV s |
| Self-electric field | 7.926e16 J | 118,503.2 MeV s |
| Sum of these three | 1.194e17 J | 178,525.7 MeV s |
| Additional 40 MV charge circulation, if traversed | 5.582e17 J | 834,597.0 MeV s |

The primary Ti channels contain 66.8% of the first three energies. The complete ion charge is 1.395e10 C. The uniform confining magnetic field stores another 1.261e21 J, which is **excluded from reset energy** because it remains energized. Gate-electrode capacitance and the evolution of axial self-fields need an actual switching design.

## Loss and power accounting

Let tau be burn time, delta be reset downtime, b the trapped positron supply per primary event, and P the assigned energy transfer per positron per second. The assumed population is N_+(t) = b R t until reset. No stationary positron abundance or annihilation balance is asserted.

For an energy reservoir U with input efficiency eta_in, output efficiency eta_out, and residual heat conversion h = 0.7, its reset electrical loss is:

```
U [1/eta_in - eta_out - h(1 - eta_out)].
```

The model follows the existing plant convention: extraction losses thermalize and receive heat-conversion credit; input inefficiency receives no additional credit. Summing these losses and dividing by R times one MeV gives A in MeV s. If P describes ion energy transferred to radiation, its net electrical replacement cost uses the existing heater and heat engine:

```
B = b P (1/0.95 - 0.7).
```

If the supplied P is already an electrical cost, the factor is 1 instead. Per active primary event the extra costs are A/tau and B tau/2. Their loss-only optimum is tau = sqrt(2A/B), with minimum sqrt(2AB). The reported optimum instead maximizes **complete net specific power**, including downtime, converter mass, radiator mass and replacement manufacture.

The main recovery cases retain eta_thermal,in = 0.95, test eta_thermal,out = 0.95 or 0.99, and test coherent efficiency 0.999, 0.9999 or 0.99999 in each direction. They optimistically include ordered rotation with the coherent reservoir. An additional case also retains 95% rotation preparation and recovery. Uniform-efficiency cases at 0.95, 0.99, 0.999 and 0.9999 are separate capability counterfactuals. Near-perfect random-energy handling would require additional energy-recovery or adiabatic-expansion apparatus, phase-space acceptance and entropy export; it is not inherited from the reference plant.

Prompt electrons have a separate assumed all-in allowance of 0.1 MeV each. This adds only 0.0000253–0.0000906 MeV per primary in the stipulated prompt-pair scenarios. The earlier 4.7 ms escape and 77.8 keV radiation benchmark concerns electrons under prescribed fields; it does not provide the positron lifetime or the full residual-gas cost.

## Main operating cases

The following rows retain 95% thermal input/output, assign 99.99% coherent transfer, use lowered-gate refill, and assume perfectly shared converters with staggered channel resets and 0.1 s downtime. P is radiated ion energy, with the heater/heat-engine factor above.

| Pair supply per primary | P per positron | Optimum burn | Net output | Net specific power |
|---:|---:|---:|---:|---:|
| 0.000253 | 10 MeV/s | 2,167 s | 0.994 TW | 0.661 mW/kg |
| 0.000253 | 30 MeV/s | 1,244 s | 0.142 TW | 0.0923 mW/kg |
| 0.000253 | 100 MeV/s | 674 s | -1.639 TW | -1.029 mW/kg |
| 0.000906 | 10 MeV/s | 1,138 s | -0.0453 TW | -0.0294 mW/kg |
| 0.000906 | 30 MeV/s | 650 s | -1.789 TW | -1.119 mW/kg |
| 0.000906 | 100 MeV/s | 349 s | -5.294 TW | -3.087 mW/kg |

Negative outputs are external power requirements for maintaining the prescribed state. If P is already electrical power, the favorable first row drops to 0.129 mW/kg. Increasing thermal extraction efficiency to 99%, while retaining 95% reinjection and 99.99% coherent transfer, improves its radiative case to 0.732 mW/kg. Retaining 95% rotation handling as well as 95% random-energy handling gives 0.326 mW/kg at the favorable supply and power, with 99.9% self-electric transfer.

For the first-row recovery model, reaching 1 mW/kg would require P below **3.32 MeV/s** at b = 0.000253, or below **0.927 MeV/s** at b = 0.000906. These are thresholds to test, not predictions of the true positron radiation rate.

Uniform 99.9% recovery of every reservoir produces 1.168 mW/kg in the favorable supply/power case; uniform 99.99% gives 1.306 mW/kg. Those results quantify how consequential an additional random-energy recovery capability would be. They do not license silently replacing the reference 95% thermal machinery.

## Converter capacity and clearing time

Every reset processes both extraction and reinjection. The added converter throughput counts U plus U/eta_in, with appropriate efficiencies for each reservoir. For synchronous resets the sum is divided by delta. For perfectly staggered resets with shared converters, it is divided by tau + delta. The latter assumes enough switching and distribution capacity to route high-power channel discharges into common machinery; local pulse-current and voltage capability are not established.

At 0.1 s, a synchronous lowered-gate reset requires **2.399e18 W** of added converter rating. At the existing allocation of 0.1 kg/W, this alone weighs **2.399e17 kg**. Its replacement burden makes the favorable case energy-negative. Shared staggered resets at the favorable optimum require **1.107e14 W**, or **1.107e13 kg**, of added converter capacity. Ten-second synchronous resets still deliver only 0.0574 mW/kg in that favorable case. Ten-second staggered resets deliver 0.657 mW/kg, so practical clearing time matters much less than simultaneous converter demand when burns last tens of minutes.

The characteristic active-length/thermal-speed times are 1.63 ms for protons and 11.39 ms for Ti49. They are not complete-emptying times. In a fixed-field, free-streaming uniform column, the low-parallel-speed tail survives as approximately L/(sqrt(2 pi) v_thermal t). At 0.1 s the residual fractions are 0.652% for H and 4.54% for Ti; at 1 s they are 0.0652% and 0.455%. Evolving axial self-fields and rarefaction change these values, so this is a diagnostic rather than a bound. A specified residual fraction and post-dump axial sweep still require demonstration. Once dense ions have left, weak axial fields can act on electron and positron remnants in opposite directions; collector reflection and re-emission remain relevant.

## Catalyst handling cannot be folded into a long reset cycle

The primary Ti inventory divided by its reaction rate is **302.1 s**. The long-burn results therefore retain the baseline continuous selective isotope extraction and replenishment. A whole-inventory reset every 2,167 s cannot itself maintain the original Ti distribution.

For a simple fixed-proton-profile depletion diagnostic, 30 s consumes 9.45% of the initial Ti. Even granting the original constant reaction rate over that short burn, the favorable separated-efficiency reset case requires about **43.8 TW of external power** and gives **-13.0 mW/kg**, with almost 8 PW of extra shared converter rating. One-, ten-, and hundred-second results are included in the data. A changing-isotope batch reactor needs a separate reaction and kinetic calculation; replacing selective extraction by short resets has not rescued this design.

## Reproduction and interpretation

Run `periodic_reset.py` with the same NumPy/SciPy research environment as the source calculations. It writes `periodic-reset.json`, including processor inventories, all recovery/supply/power/downtime/layout scenarios, short-burn sensitivities and 1 mW/kg thresholds. Input hashes are recorded.

Seven checks pass: exact recovery of the original plant power and mass when cycling is disabled, agreement of both electrostatic integrals, agreement of saved particle inventories, successful bounded optimizations, energy closure and the independent 3NkT/2 thermal-energy identity. The maximum particle-inventory integration discrepancy is 1.63e-7; the largest relative energy-closure residual is 1.17e-14.

Inputs come from `network-plant.json`, `electron_wall_orbits.py`, `catalyst_penning.py`, the prompt-pair scenarios in `PROMPT-PAIR-SUPPLY.md`, and the electron benchmark in `lepton-escape-main.json`. The positron power range is an explicit diagnostic chosen for this test. The additional converter mass uses the existing conditional 0.1 kg/W allocation. Beam transport, switching, phase-space recovery, isotope handling, material life, gas recycling, and residual-lepton clearing remain required engineering work.

An independent review by the ion-dynamics research agent reconstructed the favorable separated-efficiency case from the saved implementation: A = 2055.3950473 MeV s and B = 0.0008921578947 MeV/s. At 2166.6153 s these produce 0.94866635 MeV of reset loss and 0.96648148 MeV of accumulated-positron loss per active primary, reproducing 0.66080154 mW/kg. The review checked the self-electric one-half factor, duty factors, accumulation average, two-leg converter throughput and radiator/replacement closure, and found no accounting correction required. It reiterated that the 302.1 s Ti consumption time and complete lepton reset remain physical qualifications.
