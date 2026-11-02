# Direct recovery of the outgoing pulse ✨

This checkpoint supports adding direct kinetic recovery to the cold-sink redesign as an explicit engineering branch. The saved source has **19.9105 MeV of outgoing bulk kinetic energy per completed helium cycle**, out of 27.0228 MeV of non-neutrino nuclear heat plus prepared driver energy. Its ordered expansion can perform electrical work before thermalization. A 40% heat engine is therefore neither a fundamental ceiling nor the only available conversion route. The separate heat-engine study also supports considering higher thermal efficiency when the cold sink is changed.

The strongest concrete candidate is a tiled, time-programmed electrostatic collector upstream of the sacrificial absorber, with local electron diversion, multiple ion collection or bypass routes, and a separately rated electrical buffer. Inductive pickup from a conducting plasma is a second candidate. **50%, 70%, 80% and 85% direct-recovery cases below are conditional converter-output requirements before storage and separately debited auxiliaries.** The calculations establish source scales and kinematic opportunities; they do not predict achieved efficiencies for this installation.

## Evidence and the useful source advantage

Barr and Moir tested plasma direct converters at up to 100 kV and 6 kW, including 70 W/cm² local power density. Their TMX experiment recovered 79 W from 138 W incident and consumed 12 W in its suppressor, giving 48% net efficiency. Gas ionization, secondary electrons and high-voltage breakdown were important limitations. This supports the conversion mechanism and a relevant **local power-density scale**, while leaving our voltage, pulse inventory, species mixture and distributed installation unqualified. [Barr and Moir, 1983](https://doi.org/10.13182/FST83-A20820); [author-hosted paper](https://www.ralphmoir.com/wp-content/uploads/2012/10/fusion-science-and-technology.pdf).

Mima and colleagues analyzed pickup coils around an expanding ICF plasma, subdividing coils to reduce output voltage. Their greater-than-80% result assumed cylindrical expansion and a capacitor recovery circuit. It supports an inductive design branch, with geometry and plasma coupling still to be established here. [Mima et al., 1992](https://doi.org/10.13182/FST92-A30054). Traveling-wave ion deceleration has also been studied experimentally; bunching and deceleration acceptance remain explicit design problems. [Togo et al., 2015](https://doi.org/10.1585/pfr.10.3405013).

The present plume has a useful feature: in the archived ballistic reconstruction, arrival time identifies outward velocity, approximately `v=R/t`. Therefore a retarding voltage can follow `V(t) ∝ v² ∝ t⁻²`, instead of applying one fixed voltage to the whole pulse's energy distribution. At a given velocity, fully stripped CNO isotopes mostly have `A/Z` near 2, whereas protons have 1. Thermal velocity spread, source evolution, charge-state dispersion, residual pressure acceleration and overlapping emissions broaden this relation. The source is the saved 64-zone endpoint with a mass-conserving continuous interpolation; the interpolation's kinetic energy differs from exact nodal accounting by 0.03–0.18%, retained explicitly.

A simple favorable kinematic screen sets `V=α E_u`, with `V` in MV and `E_u` in MeV per atomic mass unit. An ion can reach a collector at this potential if `A/Z ≥ α`; its recoverable work is `Z α E_u`. Lower-ratio ions require their own route. Results assume fully stripped ions and ideal arrival-time knowledge:

| Pulse | Ideal recovered fraction of all kinetic energy at α=1.75 | Kinetic fraction needing a separate route | Fully stripped energy-weighted median stopping potential |
|---|---:|---:|---:|
| C13 | 83.57% | 2.16% | 1.395 MV |
| A | 83.46% | 4.09% | 0.721 MV |
| B1 | 81.62% | 3.36% | 0.477 MV |
| B2 | 86.22% | 4.59% | 1.037 MV |

These fractions precede grid interception, electrical conditioning, neutral losses, electron suppression and storage. They make a **70% converter-output target before separately booked parasitics** physically motivated rather than arbitrary, conditional on charge preparation and successful extraction. They do not establish that target. Higher recovery needs more stages or more selective routing. For example, α=1.875 gives 89.54% for C13 but only 65.01% for B2 because its surviving oxygen-14 has `A/Z=1.75`. A single waveform for every stage is inadequate. If the retarder must transmit all protons together with the CNO, α≤1 limits this ideal recovery to roughly 50–54%.

## Charge, electrons and collection geometry

The hydro output describes a **quasi-neutral material plume**. Fully stripped charges in this screen are diagnostic assumptions; the companion charge-state study evaluates ionization and possible reionization separately. C13's fully stripped 1–99% energy-weighted stopping-potential range is 0.175–3.827 MV. With one positive charge per nucleus it becomes 0.996–22.963 MV, with a 9.498 MV median. Thus complete stripping reduces voltage substantially but is not a prerequisite for some direct recovery. Residual neutral material carries its kinetic energy to the thermal receiver unless reionized upstream. Reionization thresholds alone do not determine wall-plug cost: far-field optical depth, photon reuse, stripping column density, charge exchange, straggling and material recycling must also close.

The companion [charge-state study](charge-state.md) finds a full-stripping threshold cost of 0.02753 W/m² averaged over this cycle. A hypothetical 50%-efficient production route would require 0.05506 W/m², about 0.0694% of outgoing kinetic power. For C13 alone the corresponding full-stripping threshold is 7.311 kJ/m². These small threshold energies are encouraging; the neutral-to-ion optical depth and complete apparatus can still determine viability. This reionization expenditure belongs in `A` below, with any returned ionization heat credited only once.

A workable electrostatic circuit must locally separate or redirect electrons before the positive-ion retarding region. Positive ions deliver charge to the collector at elevated potential; the external load and the separated electron path close the circuit. Allowing the co-moving electron population to reach that collector would create a counter-current and return energy to electrons. Charging the entire shell or leaving the upstream plume uncompensated is not a substitute for this local circuit. Local suppressors, return conductors, insulators, electron collection and secondary-electron control belong in the converter allocation. Co-moving electron kinetic energy and carried ionization energy are components of the existing source accounting and cannot simply be added as new energy.

For C13, the fully stripped charge fluence is 35.79 C/m² and the sampled peak positive current is 0.367 A/m². A reaches 0.585 A/m². Summing local peak currents over the shell gives order `10¹⁹ A`; modular local return paths are essential. The entire sum does not pass through a central conductor. Current peaks use the archived 512 quadrature packets and can slightly underresolve a cell-edge maximum.

The familiar planar, zero-injection Child–Langmuir scale is

`j = (4/9) ε₀ sqrt(2q/m) V^(3/2) / d²`.

For N14⁷⁺ at 1 MV, `j d² = 38.65 A`. A current of 0.3 A/m² corresponds to an 11.35 m gap on this benchmark. A real incoming, slowing, multi-species beam requires a kinetic/Poisson solution; this formula is a scale diagnostic, not a current limit proved for the proposed retarder. In particular, density grows near stopping, charge exchange can produce uncollected neutrals, and slow ions must be collected or routed before they are accelerated back out. A conservative static field followed by a return to the same potential recovers no net energy by itself.

Inductive recovery can retain quasi-neutrality, but a static magnetic field alone does no work. Plasma-driven currents and induced electric fields must transfer energy to an external load. At the C13 kinetic peak, favorable fully stripped diagnostics give about `2.7×10¹⁰ ions/m³`, electron density `1.0×10¹¹/m³`, ion inertial length about 967 m and electron inertial length 16.5 m. A magnetic pressure matching the peak ram pressure needs only about 0.49 mT, yet the corresponding nitrogen gyroradius is about 657 m. Therefore a small field-strength estimate does not justify a compact reversible magnetic piston. Pickup geometry, several-hundred-metre kinetic scales, field penetration and current closure need explicit treatment. Stronger fields change those scales and must be included in coil, power-conditioning and fault-protection mass. Slowing the plume transfers its momentum to the converter structure and support system.

## Electrical ledger: prepared energy is charged once

Let all energies be MeV per completed helium cycle:

- `H=24.7280375`: non-neutrino nuclear energy.
- `U=2.29474654`: prepared driver energy delivered to the source.
- `K=19.9104911`: outgoing bulk kinetic energy, already containing recovered driver energy.
- `c=0.99`: captured energy fraction; `ηp=0.60`: driver preparation efficiency.
- `D₀=c κ ηd K`: gross direct electricity before output storage, where `κ` is the eligible kinetic fraction and `ηd` its converter efficiency.
- `b`: fraction of direct output passing through storage; `r`: its round-trip efficiency; `L=b(1−r)D₀`.

If remaining captured energy thermalizes at the useful heat-engine input, new auxiliary consumption is `A`, and a justified fraction `χ` of auxiliary heat returns there, then

`E_net = (D₀−L) + ηh[c(H+U)−D₀] − U/ηp − A + ηh χ A`.

Storage loss is rejected separately in this checkpoint; crediting it as useful high-temperature heat would require an additional temperature and routing model. Relative to the same heat engine, the ideal incremental gain is `(1−ηh)D₀−L−(1−ηhχ)A`. Directly routing recovered electricity to the next driver can reduce external power demand, but does **not** add a second credit. The existing `U/ηp` debit is the only driver preparation debit.

For `κ=1`, before new converter auxiliaries, station support and replacement:

| Thermal efficiency | Gross direct fraction of captured K | Direct output storage | Net MeV/He |
|---|---:|---:|---:|
| 60% | 0% | — | 12.2270 |
| 60% | 50% | lossless | 16.1692 |
| 60% | 70% | lossless | 17.7461 |
| 60% | 70% | all output, 90% round trip | 16.3663 |
| 60% | 80% | lossless | 18.5346 |
| 60% | 80% | all output, 90% round trip | 16.9577 |

The 50% case is a cautious development target; 70% is a useful conditional redesign case; 80–85% are stretch cases requiring efficient stage selection and circuit recovery. None is a measured performance estimate for this plume. Actual neutral fraction reduces `κ` unless separately processed.

## Peak hardware and steady-output storage

At `R=1.1×10⁹ m` and the frozen `mass999`, 80%-occupancy schedule:

- Mean bulk-kinetic flux is **79.350 W/m²**.
- Largest kinetic peak is **636,158 W/m²**, a peak/mean ratio of **8,017**.
- C13's exact nodal kinetic pulse energy is **5.19338×10²⁶ J**.
- At 80% gross recovery and 99% capture, direct output peaks at **503,837 W/m²**, or **7.661×10²⁴ W** across the receiver, and the largest output pulse contains **4.11316×10²⁶ J**.

A complete converter rating based on captured input peak, including collection structures, fields, return paths, local switching and conditioning, gives:

| Complete peak input rating | Converter mass per receiver area | Whole converter mass |
|---:|---:|---:|
| 10 W/kg | 62,980 kg/m² | 9.576×10²³ kg |
| 100 W/kg | 6,298 kg/m² | 9.576×10²² kg |
| 1,000 W/kg | 630 kg/m² | 9.576×10²¹ kg |
| 10,000 W/kg | 63 kg/m² | 9.576×10²⁰ kg |

These are sensitivity allocations, not proven complete-system ratings. If a vendor rating instead uses electrical output peak, convert the numerator accordingly; do not count the same conditioner twice. Sizing the braker solely from average output would understate its pulse duty by several thousandfold. Spatial concentration on a smaller area increases local flux and current; it does not remove the total pulse-rating requirement.

For steady export, a storage constraint is the maximum cumulative mismatch between converter output and permitted export or scheduled reuse. From the C13 1–99% energy-arrival window alone, the required capacity is at least **4.02344×10²⁶ J** after subtracting constant cycle-average withdrawal. For independently periodic stage streams, the sum of recovered energy in one pulse of each stage, **1.00915×10²⁷ J**, is a conservative sufficient capacity under ideal continuous average withdrawal. This is a lower-bound/sufficient-allocation bracket, not a solved joint schedule. It corresponds to **26.46–66.37 MJ/m²**, or **264.6–663.7 kg/m² at 100 kJ/kg**, before structural reserve and charge/discharge hardware. Storage usable-energy convention, efficiency, lifetime and peak charging capability must be specified together.

The existing two-driver-pulse store is **1.35813×10²⁶ J** and is kept separate. It is already assigned to driver service and is smaller than one large recovered pulse. A joint chronology could share some capacity; no such credit is taken. Exporting pulses instead is possible only if the interconnect and receiving loads accept the corresponding peak. Thermalizing excess output preserves a heat-recovery route but forfeits the claimed direct electrical gain.

A useful incremental break-even against a 60% thermal station of specific output `s₀` is

`ΔM < ΔP/s₀`, hence `r_peak > s₀ c P_K,peak / ΔP`

for converter-only added mass. At 80% gross direct recovery, `ΔP=25.138 W/m²` with lossless storage, or `18.854 W/m²` when all output traverses 90%-round-trip storage. For illustrative `s₀=4 mW/kg`, the minimum complete input-peak ratings are **100.2 W/kg** and **133.6 W/kg** respectively. Added storage and auxiliaries raise those thresholds; reduced radiator and heat-engine loads can lower them. The parent plant model must evaluate the complete mass feedback.

The final matched `cases.json` reference, `working_200K_60pct_IR11_bypass1`, has `s₀=5.394993 mW/kg` and an ordinary auxiliary load equal to 10% of the process margin. Applying that convention gives 80%-direct marginal benefits of **22.6243 W/m²** and **16.9682 W/m²**, corresponding to total added-mass allowances of **4,193.6 kg/m²** and **3,145.2 kg/m²** before new mass-dependent loads. Converter-only minimum peak ratings become **150.2 W/kg** and **200.2 W/kg** respectively. These match the thermal reference; they do not establish the direct plant's final specific output.

## Position relative to the separating film

The existing stopping screen predicts C13 projectile-energy removal of 0.2006%, 2.0000%, 19.2006% and 89.9659% in 10 nm, 100 nm, 1 μm and 10 μm silica respectively. An upstream film therefore subtracts that ordered energy from the direct-recovery budget; scattered trajectories and charge-state changes impose additional acceptance losses. The energy removed can contribute to heat if retained and routed, but cannot also remain in `K` credited to the converter.

Placing the converter first preserves the largest direct opportunity. Strongly braked residual material then stops in the film, so the film becomes a material catcher requiring a controlled recycling and electron-current route. Partial braking leaves a reduced transmitted kinetic spectrum and requires the species-resolved stopping calculation to be repeated. Neither placement establishes membrane lifetime from fluence alone. The previous film damage, sputtering, charging and fuel-return limits remain relevant.

## Reproduction and scope

Run `PYTHONPATH=/tmp/protium-python OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 direct_recovery.py --root /home/tkwa/.codex/worktrees/protium-fusion/website`. Inputs, code hashes, all stage currents/potentials, electrical cases, storage bounds, hardware sensitivities and energy-identity checks are saved in `direct_recovery.json`. Weighted nodal kinetic energy agrees with the frozen cycle to `4.19×10⁻⁸` relative; it is not silently rescaled. The electrical global-energy residual is below `10⁻¹² MeV/He`.

Remaining decisive work is a measured or validated atomic charge distribution, coupled electron/ion collector or inductive transport, a timed collection-and-return circuit, secondary emission and breakdown qualification, a complete peak-rated module design, and storage/export scheduling. These requirements support a specific future development program and conditional plant comparison without assigning an unearned efficiency or deployment date.
