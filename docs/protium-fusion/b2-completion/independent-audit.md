# Independent trace and process audit

The arithmetic audit independently reconstructs physical radii, cell densities,
Lagrangian cell and node masses, velocity/energy units, kinetic and EOS internal
energy, atomic mass defect, weak energy and radiation losses from each saved
trace. It does not import either hydrodynamic runner. Nuclear energy is computed
using isotope mass excesses, which avoids the large cancelling common baryon
rest mass. The newer traces also carry separately integrated nuclear heat and
signed weak nuclear-charge change. Both are compared with abundance changes.

The constitutive EOS and isotope mass library remain shared inputs. These
checks test code accounting and the saved trajectory; they do not independently
validate the physical EOS, cross sections, opacity or driver.

## Legacy defect

All five legacy small-stage traces reproduce their reported energies, isotope
fractions and baryon counts. Nevertheless, all have negative scalar artificial
viscous heating in expanding spherical cells. The legacy switch activates
positive isotropic pressure from a negative radial velocity difference, even
when spherical volume is increasing. The resulting term `−q dV/dt` is negative.

The B2 trace contains 317 saved cell-samples of this defect, with most-negative
summed power about −1.42e30 W. A sparse trace trapezoid gives −3.16e24 J,
7.75% of preparation energy. Corresponding diagnostic ratios are 4.45% for A,
3.07% for B1, 6.97% for C13 and 0.544% for He3. These sparse integrals are
neither accurate shock integrals nor predicted changes in yield.

Force work and internal pressure work cancel to roughly 3e−16 in relative
arithmetic. Thus excellent global conservation alone does not detect the
entropy defect. The legacy B2 endpoint additionally stops at 2.23 of its
requested four sound times; its original `final_time` label is misleading.

The first corrected `volume_gate` trace respects the heating sign and closes
independently integrated nuclear heat to 1.8e−14 of input energy and weak
charge to 2.4e−16 per baryon, but stalls at 1.30 sound times before substantial
burning. An abrupt gate is discontinuous at zero volume compression when the
radial velocity difference remains negative. The continuous correction is
being tested separately; this failed trace is not a replacement yield.

## Completed corrected trajectories

The continuous `compression_cap` correction completed all five stages to the
requested four sound times. Each passes all ten independent trace checks:
reported quantities, reconstructed units and masses, baryons, trace positivity,
energy closure, pressure-force/work cancellation, viscous heating sign,
completion, separately integrated nuclear heat, and signed weak-charge change.

| Stage | Largest saved-sample energy residual / input |
| --- | ---: |
| A | 1.41e−7 |
| B1 | 1.32e−6 |
| B2 | 6.18e−7 |
| C13 | 3.56e−6 |
| He3 | 2.38e−8 |

The completed B2 state has positive, outward-increasing boundary velocities,
minimum temperature 2.20 MK and minimum outward velocity 150 km/s. Its
nuclear heat accumulator agrees with mass-excess energy to 1.29e−13 of input,
and signed weak charge agrees to 4.92e−16 per baryon. Thus a cold EOS or
cold-rate replacement was not required to complete this corrected calculation.
These figures qualify the accounting of this reduced spherical calculation;
mesh convergence and multidimensional physical fidelity remain separate tests.

## Independent reaction accounting

`network_stoichiometry_audit.py` reconstructs baryon and nuclear-charge changes
from reaction reactant/product metadata. All 194 reactions conserve baryon
number. All 182 strong reactions conserve nuclear charge; the other 12 are
weak reactions and require compensating leptonic charge. Mass-excess Q agrees
with the runner's atomic-mass Q to 5.46e−12 MeV.

Reaclib embedded Q differs from the newer isotope mass table by up to 4.07 keV.
The runner's heating uses the mass table, so this does not produce an energy
bookkeeping discrepancy. The comparison does not quantify possible reverse-fit
detailed-balance effects or rate uncertainty.

## Independent process balance

`independent_cycle_audit.py` uses a five-variable linear balance normalized to
one net helium atom, rather than the production ledger's factory-normalized
four-variable fixed-point solve. The independent mass-excess energy and
preparation values agree with `cycle_ledger.py` to 3.2e−11 MeV per helium over
collection losses from zero to 0.001. This validates the algebra with the
legacy inputs; it does not restore their physical validity.

The same script supplies a six-variable open-loop balance: four catalytic
stage attempt rates, fresh carbon and fresh hydrogen, normalized to one
helium product. It separately recovers C12/C13/N14/N15 and hydrogen, keeps
all unreturned isotopes in the outgoing inventory, and compares external
feed minus outgoing-product mass defect with the stage-summed Q. Fresh carbon
has a selectable C13 atom fraction; the default 0.0107 is an assumption pending
the chosen source composition. Fresh C:H uses actual neutral-atom masses.
There is no helium-to-carbon factory in this open-loop calculation.

With the legacy trajectories, the diagnostic minimum equal CNO/H recovery for
fresh C:H mass ratio 0.007 is 0.9995292242. Heat before collection is
24.73274274 MeV per helium and preparation is 2.30074752 MeV; a selected
60% driver, 40% conversion and 99% collection gives 6.87068294 MeV before
other auxiliaries. These are diagnostic bookkeeping values from invalid
trajectories and must be replaced before selecting a plant scenario.

All code, reports and source hashes are in this directory. Full trace paths are
recorded in each report. The reusable verifier requires the corresponding
`.npz` and `.json`; its completion check independently requires the requested
elapsed sound times as well as solver success and the final-time label.

## Corrected process comparison

The corrected external-carbon report has been checked at all 21 reported
operating points using the independent six-variable solve. All 84 checks pass:
reported physical fields, external mass-defect energy, baryons and stationarity.
The largest field differences are about 3.3e−10 relative, dominated by the
neutral-atom mass constants' rounding. `open_process_comparison.py` reproduces
this comparison against the selected `process.json`.

For the corrected trajectory set, a fresh C:H mass ratio of 0.007 with equal
CNO and hydrogen recovery requires 0.9995280479 return. Closed-station heat
before collection is 24.73328718 MeV per helium; preparation is 2.29511307 MeV.
The selected 60% driver, 40% conversion and 99% capture arithmetic yields
6.87805805 MeV per helium before other auxiliaries. These efficiency and
recovery values remain engineering assumptions, rather than outputs of the
pulse calculation.

The consistently corrected factory-inclusive cycle at collection loss 1e−5
per stage gives heat 24.72265857 MeV, preparation 3.58130479 MeV and selected
preaux electricity 5.23952817 MeV per helium. Relative to the old conditional
cycle, heat changes by +0.00049734 MeV and preparation by −0.00837056 MeV.
The completed expansion shifts 1.97 MeV per helium into bulk kinetic energy
and removes 2.05 MeV from the final internal-energy component. This partition
change includes the longer B2 evolution, so it cannot be attributed solely to
the viscosity correction.

## Additional numerical controls and collector interface

The tighter-tolerance 64-zone B2 run, eight-sound-time continuation, 128-zone
nested refinement, energy-matched 128-zone refinement, and alternate
divergence-viscosity run each pass the same ten independent trace checks.
The 64-zone tight run lowers the final energy residual to 1.92e−8 of input.
The eight-time run preserves the nuclear yield while continuing conversion of
internal energy into ordered expansion.

A 128-zone mesh raises B2 yield by about 5.9% relative to 64 zones. Matching
the initial kinetic energy changes that refined yield only about 0.019%.
This is material spatial uncertainty, which excellent conservation does not
remove. Substituting the unscaled 128-zone B2 into the corrected external-carbon
cycle raises preaux electricity from 6.87806 to 6.91312 MeV per helium,
about 0.51%. The alternate entropy-safe viscosity changes the single-pulse
yield by about 1.1%. These controls do not establish grid convergence.

`collector_interface_audit.py` independently checks the reduced collector
inputs against the raw traces and exact preparation energy. All 44 checks
pass: initial kinetic and internal energy, source-bin energy integral and peak,
prompt cumulative energy directly from raw units, trace hashes, the three
incoming/outgoing timing definitions, summed occupancy, and cold/received heat.
The total prepared source peak is 3.122783331e23 W at a 1e9 m radius and 60%
driver efficiency. The summed 99.9%-mass occupancy is
1.922865168e−33 s per helium, limiting throughput to 4.160458119e32 helium/s
at 80% occupancy. The source peak includes initial EOS internal energy.

The collector interface check does not revalidate thermal convolution,
photon deposition, stopping-induced nuclear transmutation or active removal
of the zero-speed central tail. A 99.9%-mass schedule leaves a larger slow
fraction than the allowable unrecovered CNO fraction; it requires protected
recovery, not disposal of the tail.


## Final refinement and whole-station closure

The energy-matched 256-zone B2 trace also passes all ten independent checks.
The eleven completed corrected trajectories therefore pass 110 trace checks
in total. The 256-zone calculation reaches four sound times with Q=11.43272617
and target yield 0.8985603883. Its largest saved-sample energy residual is
4.33e−6 of input, baryon error 4.09e−14, and separately integrated nuclear heat
agrees with abundance mass defect to 2.59e−13 of input. The all-trace minimum
mass fraction is −9.54e−12; the final minimum is −4.35e−19. There is no negative
viscous heating at the saved samples. Gain rises another 2.98% and yield 3.25%
from matched 128 to 256 zones, so spatial convergence remains incomplete.

The coherent 64-zone family remains the selected station baseline. The refined
B2 results are sensitivity cases; replacing only one stage does not produce a
uniformly refined system. With 0.0056 fuel C:H and 0.0014 net material reserve,
selected process electricity before station auxiliaries is 6.87644, 6.91226 and
6.93291 MeV per helium for the 64-, 128- and 256-zone B2 families respectively.

`station_independent_audit.py` independently solves the component mass sum with
a numerical root, instead of using the production quadratic. All 18 checks
pass for the rounded balanced and aggressive points and the finite-scan winner.
The selected rounded aggressive point gives 1.215203326e23 kg and
3.360979011e20 W, or 2.765775 mW/kg. The balanced point gives 2.296788532e23 kg
and 3.009007883e20 W, or 1.310094 mW/kg. One-day process inventory, full solar
heat at 5 AU, replacement work, stored energy, installed reserve and the
conservative maximum solar tidal load are included.

The support check independently starts from g_self=GM/(2R²) and maximum
rotating Hill-frame acceleration g_tide=3GM_sun R/a³. Equating their summed
scalar load to stream centripetal force gives stream mass
M(g_self+g_tide)R/v². Charging the selected 0.1% kinetic loss per stream circuit
then reproduces support dissipation. The maximum tide is 31.2% of the
self-load in the rounded aggressive point. The bound does not supply a guide
topology or orbital stability analysis.

At the rounded radius the mean captured hot heat flux is 106.618 W/m². A
2.5-metre carbon layer with conductivity 100 W/(m K) has an illustrative steady
through-thickness gradient of only 2.665 K. The finite engine approaches are
47.72 K on the hot side and 26.14 K on the cold side. Large single-pulse surface
excursions and their periodic accumulation are distinct constraints. The
finite-slab kernel in `collector-periodic.py` has the correct positive mixed
Neumann/Dirichlet modes and integrates to the steady resistance L/k; its
1965.3 K bound additionally requires the declared daily large-shot guard,
constant properties, maintained back temperature and active slow-tail handling.

The source rating and service mass are explicit apparatus requirements. The
independent closure does not validate their construction. A component-resolved
electromagnetic brake, particle-spectrum transport, material survival and
recycling model remain necessary to qualify the carbon budget.
[Wall stopping requirements](wall-stopping-requirement.md) quantify the
remaining irreversible carbon allowance and a candidate deeply sub-barrier
entry-energy screen without asserting a universal nuclear threshold.
