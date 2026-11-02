# ✨ Refurbishment, selective extraction and material inventory

The useful form of per-pulse nanorefurbishment is local defect repair, renewed seals and selective extraction from a shallow collection structure. Whole-skin processing is a separate option with a calculable throughput. Recycled graphite and held hydrogen/helium stay in the material inventory.

This account uses the coherent delivered-B1/C13H2 family, H/catalyst 3,2,2,2, the 350 kg/m² receiver at radius 1.1×10⁹ m, and the final 370-day dispatch. [calculate.py](calculate.py) produces [results.json](results.json), with input hashes and material/energy checks. The hydro calculation uses baryon masses A mᵤ; the process ledger uses neutral atomic masses. Exact nucleus counts connect them.

## Incoming material

| Pulse | Total baryonic material, µg/m² | H+He, µg/m², baryonic | H+He, µg/m², neutral atomic | H₂+He volume at 1000 K and 1 bar, µL/m² |
|---|---:|---:|---:|---:|
| A | 743.802 | 96.1324 | 96.8847 | 3,996 |
| B1 | 7.43802 | 1.92184 | 1.92686 | 50.8805 |
| B2 | 743.802 | 20.7432 | 20.9055 | 862.344 |
| C13 | 743.802 | 30.3319 | 30.5693 | 1,260.97 |

The larger pulses deliver **0.744 mg/m²**, and B1 delivers **7.44 µg/m²**. The rounded user checks reproduce correctly. Fully neutral returned H is 3.2270404×10⁶ kg/s, and He is 2.758883×10⁶ kg/s. Their enormous whole-sphere rates correspond to small local loads over 1.52053×10¹⁹ m².

These gas volumes are ideal-gas storage requirements if all H becomes H₂. Atomic H would contribute twice the hydrogen pressure; hydrocarbons would change both particle counts and carbon routing. A 50 µm-equivalent enclosed void holds one complete A gas pulse at approximately 7.99 kPa and one B1 pulse at 102 Pa at 1000 K. The pressure vessel capacity does not establish that the irradiated front can retain that pressure. The [barrier analysis](../barrier-review/) distinguishes storage behind an isolated valve from pressure applied directly to the catcher skin.

## Whole-layer handling

A 50 µm graphite layer at 2000 kg/m³ is **0.1 kg/m²**. The selected mean rates give 4,173.02 stage events/year; the actual 370-day dispatch has 4,228 events, including 75 larger pulses. The steady cadence is used for annual plant capacity; the event list supplies pulse counts and transient checks.

| Layer and service schedule | Handled kg/m²/year | Input at 10 kJ/kg, W/m² | Input at 2 MJ/kg, W/m² | Hardware at the existing 10⁸ kg/(kg/s) factory coefficient, kg/m² |
|---|---:|---:|---:|---:|
| 50 µm, every pulse | 417.30 | 0.1322 | 26.447 | 1,322.35 |
| 100 µm, every pulse | 834.60 | 0.2645 | 52.894 | 2,644.70 |
| 150 µm, every pulse | 1,251.91 | 0.3967 | 79.341 | 3,967.05 |
| 100 µm, larger pulses | 14.734 | 0.00467 | 0.9338 | 46.689 |
| 100 µm, monthly | 2.400 | 0.000761 | 0.1521 | 7.605 |

Both the specific input energy and machinery coefficient are scenario parameters. A layer translated, inspected or selectively repaired need not undergo the same chemical transformation as a new high-energy product. Moving a 1 m-long strip through a service station once per mean stage interval requires only 0.132 mm/s. For 100 µm every-pulse handling, a 1 kg/m² local machine must process 2.285 kg of graphite per kg of machine per day. A proposed machine should demonstrate that throughput, its contact/rebonding work, its vacuum seals and its service lifetime.

Cooling 100 µm of graphite from 1000 to 300 K takes 252 kJ/m² per pass with the present constant cₚ=1800 J/(kg K). Doing this after every pulse diverts 33.32 W/m² of sensible heat. At 95% or 99% heat recuperation, the residual diversion is 1.666 or 0.333 W/m². Hot in-situ repair avoids this particular cooling cycle. The lost high-temperature heat and refrigeration loads are separately assigned in the plant account.

A moving or exchanged skin also needs an interface that restores heat conduction. A spare skin adds inventory, while gaps, contact resistance and local supports set the thermal performance. The [collection geometry](../geometry/) specifies this interface separately from the sealed heat-engine working-fluid circuit.

## Targeted repair

The current all-species elastic-transfer proxy is 37.358 MJ/m²/year. The NRT conversion gives 3.7307×10²⁴ proxy displacement events/m²/year at E_d=25 eV, or 2.3317×10²⁴ at 40 eV. Prompt cascade recombination, hot annealing, defect clustering and repeated visits to the same site determine how many surviving defects need service.

Let f be serviced events divided by the NRT proxy, and let w be the **net electrical input in eV per serviced event**. At the nominal stopping calculation and E_d=25 eV,

\[
P_{\rm repair}/A = 0.01894085\, f w\quad{\rm W/m^2}.
\]

For f=1, input costs of 1, 10 and 100 eV/event give 0.0189, 0.189 and 1.894 W/m². The one-carbon-atom equivalent of all proxy events is only 2.356×10⁻⁹ kg/(m² s). This small mass does not establish machinery throughput: reaching each defect, identifying a connected leak, moving neighbors and restoring the lattice can require many operations. It establishes the distinction between local service and processing every atom in a 100 µm layer.

The formation energy of an irradiation defect was already supplied by the incident particle. Recombination can release stored energy; an activation barrier is a kinetic constraint. [Ewels et al.](https://doi.org/10.1103/PhysRevLett.91.025505) identify a metastable graphite Frenkel-pair mechanism associated with stored Wigner energy. The repair ledger therefore uses explicit net work per completed service operation. It does not charge the formation energy again as a compulsory manufacturing input.

The connected-leak requirement can be stricter than a bulk conductivity requirement. A small area of open pores can matter even when most lattice sites and most heat conduction remain intact. Seal inspection and repair should operate during irradiation where necessary; restoring the skin only after the pulse does not recover gas that already entered the cavity. The [hydrogen](../hydrogen-evidence/) and [helium](../helium-evidence/) packets give the relevant mobile and trapped population bounds.

## Selective gas processing

If a positive free-energy upgrade of Δg eV is required per H or He nucleus, the current gas rate requires

\[
P_{\rm useful}/A=2.46920\times10^{-5}\,\Delta g\quad{\rm W/m^2}.
\]

At the project's 50% conversion efficiency, 1, 10, 100 and 1000 eV/nucleus cost 0.0000494, 0.000494, 0.00494 and 0.0494 W/m². Higher efficiency reduces these inputs. These scenarios cover a specified useful transformation; they do not turn a desorption activation energy into a thermodynamic work requirement. Pumps still need conductance, thermal isolation, finite capacity and a regeneration path. Their standing and refrigeration loads may exceed the atom-handling work.

Graphite regeneration and isotope separation return recovered carbon to the process. Fresh carbon is required only for matter exported from the accessible loop, deliberately stored outside it, or made unavailable on the relevant plant timescale. Gas circulating in the sealed cavity, graphite held in a service station, captured dust and gas in an enclosed extraction line all remain inventory. Availability and return delay can force additional inventory without changing elemental conservation.

## Plant-account interface

The selected point exports 45.6381 W/m² at 2195.862 kg/m², or 20.7837 mW/kg. It retains an ordinary auxiliary allowance of 5.07283 W/m². Previously named pumps and guide cooling occupy about 1.82798 W/m², leaving **3.24485 W/m²** for unassigned services. Of the 1000 kg/m² general-service allowance, **678.606 kg/m²** remains unassigned after existing named components. These are allocation limits awaiting machinery, rather than measured capacities.

The existing generic routine repair entry is 0.01392 W/m² and 0.6958 kg/m² of factory capacity. It cannot silently stand for every-pulse whole-skin processing. A new explicit repair allocation must identify whether it replaces part of this entry, occupies another retained allowance, or adds a new load.

The selected 200 K radiator already uses its specified 90% shell area. An additional electrical load outside the retained allowance creates both a net-power debit and extra cold rejection. Reassigning an already budgeted auxiliary does neither twice. Diverting receiver heat to a colder collection system also changes heat-engine input and must be carried through the same ledger. [results.json](results.json) contains fixed-rate sensitivities; their required radiator temperature is a heat-rejection calculation, and changing it requires a new heat-engine solve.

## Evidence and uncertainty

Mass, event counts and the stated energy/throughput sensitivities are deterministic calculations from the current files. NRT values are damage proxies. The efficiency premise is supplied by the project. Net work per repair, the service fraction, seal lifetime, damaged-material conductivity and machinery specific throughput remain engineering inputs. No credible interval for an actual repair machine is assigned before a design specifies those operations. The sensitivity tables make the required performance explicit without presenting a chosen operation cost as a measurement or a forecast.

Reproduce with `python3 calculate.py` from this directory; source paths are resolved relative to the repository. Numerical results retain more digits than warranted by material-property uncertainty so material/energy closure can be checked.
