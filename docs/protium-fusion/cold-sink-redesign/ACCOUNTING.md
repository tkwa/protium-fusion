# Cold-sink station accounting ✨

The working thermal case uses **60% converter-package heat-to-electric efficiency**, a 1000 K hot receiver, a 200 K radiator, zero direct sunlight, and an explicitly specified infrared environment. With 11.0254 W/m² absorbed over the full outer sphere, a 1 W/m² hot-to-cold bypass, and 2 MJ/kg effective routine refurbishment, the conditional station exports **6.5743×10²⁰ W from 1.21859×10²³ kg: 5.395 mW/kg**. Its radiator covers 79.925% of the outer area, within the selected 90% available area.

This package efficiency is 75% of the external 1000-to-200 K Carnot bound. It includes the converter's internal pump and bearing allowance. Driver electricity, ordinary station auxiliaries, refurbishment and support are then subtracted. The resulting net export is **40.55% of total nuclear energy**, or 43.87% of the non-neutrino nuclear energy. Those percentages describe different boundaries.

The old 40% assumption accompanied a 450 K radiator: its external Carnot bound was 55%. A 200 K sink raises that bound to 80%. The [explicit conditional gas-cycle study](heat-engine.md) obtains 61.4% package output with finite intercooling, reheating, recuperation, pressure losses and internal mechanical/pump allowances. This supplies a concrete basis for the 60% working case; actual components and irradiation lifetime remain engineering requirements.

## Frozen inputs and retained equipment

The pulse rate, radius, fuel recovery, 99% pulse-energy capture, 60% driver efficiency, prepared-energy ledger, source peak rating, source storage and fuel inventories are those of `high_source_rating_surface` in the frozen `open-loop-redesign` packet. No radius/cadence optimization or additional nuclear yield is used.

The full receiver, 100 kg/m² insulation allocation, 1000 kg/m² generic service allocation, source hardware and inventory remain in the denominator. The old gross-electric converter mass is **replaced** by explicit package bus power divided by a selected 10 We/kg equipment rating. Dedicated heat exchangers are **additional** to generic services, with no mass credit claimed for functions that might eventually be integrated. The baseline detailed service design is unresolved; these named new allocations should be reconciled with it before claiming an optimized plant.

Factory mass follows repaired throughput, independently of energy cost. The original factory allocation implies 10⁸ kg of machinery per kg/s repaired. That coefficient is retained as total installed mass changes. Reducing routine electricity to 2 MJ/kg therefore does not shrink factory mass or improve its implied throughput. Ten-year spares, 20% installed reserve, self-gravity/solar-tide support streams and support dissipation remain included.

Thermal-engine hardware is 9.636×10¹⁹ kg in the working case. Dedicated exchanger mass is 8.861×10¹⁸ kg at 10 W/K/kg; radiator panels are 2.431×10²⁰ kg at 20 kg/m². The factory-capacity allocation is 3.861×10¹⁹ kg. These are requirement-based allocations, not manufactured component measurements.

## Heat-engine temperatures and conductance

For gross package efficiency η and selected internal Carnot quality φ,

`Tc_internal / Th_internal = 1 − η/φ`.

`finite_engine()` minimizes the sum of the two external exchanger conductances subject to this temperature ratio. At η=0.60, φ=0.80 and a 1000/200 K external pair, it gives internal temperatures 911.696/227.924 K, with finite approaches 88.304/27.924 K. The 60%/250 K/φ=0.8 boundary has zero approach and is rejected with a numerical tolerance; floating-point rounding cannot silently create an infinite-conductance solution.

The two-exchanger result is a minimum, not a complete recuperated engine. The working case adds the independently calculated six-stage cycle's recuperator and generator-cooling conductance to the optimized external-exchanger minimum. This yields **5.828 W/K per outer m²** after the 1 W/m² bypass, or 8.861×10¹⁹ W/K total. The internal conductance allowance is scaled with actual hot input. At 65% the stronger cycle's internal allowance is used. At 70% this remains a minimum requirement; no matching complete 70% cycle has been solved. The 220 K cases likewise describe a finite-temperature and conductance requirement, rather than claiming the separately calculated 200 K gas cycle works unchanged.

The 60% package already contains the independent cycle's internal loss allowance. Separate 0.5% and 2% hot-input pump loads appear only as **additional external-loop sensitivities**. They are not silently charged a second time in the main case.

## Energy and environmental boundaries

Let H be non-neutrino nuclear power, U prepared pulse power, c=0.99, Qcap=c(H+U), and B the hot-to-cold bypass. With all captured heat eligible for the hot reservoir,

```
Qhot = Qcap − B
G = eta_heat Qhot
D = U / eta_driver
Aordinary = 0.1 (G − D)
Pexport = G − D − Aordinary − Pextra_pumps − Prepair − Psupport
Qcold = Qcap − U − Pexport
```

The final expression includes engine rejection, bypass, driver preparation losses and all dissipated electrical auxiliaries. This accounts for prepared pulse energy once. Increasing bypass reduces electrical export by η times its heat diversion before the ordinary-auxiliary response; the same energy remains a cold-sink load. The code also exposes a hot-eligible fraction. A sensitivity sends the entire delayed-decay component to the cold sink rather than crediting it at 1000 K.

Direct sunlight is zero. `environment_absorbed_W_m2` is required explicitly: it represents already absorbed external IR after spectral absorptivity and external source view factors. The principal cases multiply it by **the full 4πR²**, including non-radiating outer area whose absorbed heat must be routed to rejection. A separate case charges an identical numerical flux only on actual radiator area; this is a different environmental/geometry assumption.

For a cold-opening view factor f whose complement sees identical gray radiator faces, effective emission coefficient is

`epsilon_escape = epsilon f / [1 − (1 − epsilon)(1 − f)]`.

The corresponding radiator equation is

`Arad = [Qcold + q_full_sphere Aouter] / [epsilon_escape sigma Tr^4 − q_radiator_only]`.

The already absorbed environmental input is not multiplied by emissivity or view factor again. The [environment packet](environment.md) derives the 11.0254 W/m² benchmark from an illustrative compact Solar IR source; it is a declared environment, not an established local measurement. Zero-background cases are explicitly optimistic controls. Warm self-views and actual radiator orientation remain separate geometric requirements.

## Bounded comparison

All rows below use 60% package conversion, φ=0.8, 2 MJ/kg refurbishment, full-sphere IR absorption, and 90% available radiator area unless stated otherwise. An area failure means the displayed electrical ledger cannot be realized in that radiator footprint.

| Radiator | Absorbed IR, W/m² | Bypass, W/m² | Net mW/kg | Radiator area / outer area | Area result |
|---:|---:|---:|---:|---:|---|
| 197 K | 11.0254 | 1 | 5.394 | 0.849 | fits |
| 197 K | 11.0254 | 10 | 4.787 | 0.912 | exceeds 0.9 |
| 200 K | 11.0254 | 1 | 5.395 | 0.799 | fits |
| 200 K | 11.0254 | 10 | 4.788 | 0.859 | fits |
| 220 K | 11.0254 | 1 | 5.399 | 0.546 | fits |
| 220 K | 11.0254 | 10 | 4.792 | 0.587 | fits |
| 200 K | 0 | 1 | 5.397 | 0.664 | fits |
| 200 K | 30 | 1 | 5.391 | 1.032 | exceeds available area |
| 200 K | 100 | 1 | 5.377 | 1.889 | exceeds available area |

At the working point, the maximum full-sphere environmental input fitting 90% area is about 19.25 W/m² under its other assumptions. A cold-opening view factor of 0.9 increases area to 87.92%. Routing all delayed heat to the cold sink gives 3.884 mW/kg and 94.74% area at 200 K; at 220 K it gives 3.888 mW/kg and 64.71% area. Insulation and where decay heat is captured materially affect the result.

With the working geometry fixed, routine cost 3.1/10/100 MJ/kg gives 5.392/5.370/5.084 mW/kg. Changing the exchanger rating from 10 to 1 W/K/kg gives 5.391 mW/kg; changing thermal converter rating from 10 to 1 We/kg gives 5.347 mW/kg. Extra external pump charges of 0.5%/2% of hot input give 5.329/5.131 mW/kg. The existing source and receiver masses dominate these particular sensitivities.

Higher-efficiency requirements give 5.988 mW/kg at 65%, φ=0.85, and 6.581 mW/kg at 70%, φ=0.9. Those cases are not selected merely because their output is higher. Their component performance and regenerative heat exchange need a corresponding design. The 80%/150 K example fails radiator area despite its larger calculated electrical output.

## Historical controls and direct recovery

The original evaluator is reproduced without modification. Its 197 K/60%/100 MJ/kg case gives **5.150416 mW/kg**, with 90.804% radiator coverage while retaining the old sunlight load. Removing sunlight from that historical accounting gives **5.153446 mW/kg** and 73.096% coverage. Shading alone at 450 K barely changes specific power; it reduces heat-rejection area. A colder sink improves electrical output when the converter is redesigned to exploit it.

The direct-recovery entries are marginal **energy bounds**, separate from the integrated thermal family. If D becomes direct bus electricity from captured kinetic energy, the thermal ledger loses that same D before conversion. With otherwise hot residual energy, the gross increment is `(1−eta_heat)D`, not D on top of unchanged thermal output. Ordinary station auxiliaries and any new direct auxiliary load are also counted. Cold-dumped unrecovered kinetic energy loses its thermal conversion opportunity.

No direct-converter plant is ranked. Its peak handled power, pulse buffer/export capacity, charge-state handling, supports and new losses are unresolved. The existing two-driver-event store is already assigned. At 85% ideal direct recovery, the marginal bound allows 6.775×10²² kg of extra apparatus before erasing the *specific-power* improvement, ignoring new mass-dependent loads. This is a break-even ceiling, not a converter mass prediction or permission to ignore peak power and storage.

## Reproduction and validation

Run `PYTHONDONTWRITEBYTECODE=1 python3 docs/protium-fusion/cold-sink-redesign/cases.py`. The bounded family contains 36 cases, with every input recorded in `cases.json`. Eight numerical comparisons reproduce the historical mass, power, specific power and radiator area to relative 10⁻¹². The exact zero-temperature-margin case is rejected. All computed cases close installed mass, electrical power, total nuclear energy including neutrinos, separate cold-heat components and radiator/environment energy to relative 10⁻¹⁰; the largest observed residual is below 5×10⁻¹⁶.

The `feasible` field means these declared scalar mass/energy/temperature/area conditions pass. It does not establish a component-resolved plant, irradiation lifetime, insulation bridge performance, fluid maps or transient operating schedule. All previous packets remain historical snapshots.
