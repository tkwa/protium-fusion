# Mixed-stage dispatch and retained slow tails ✨

The current provisional mixture point has a separate executable 370-day dispatch in `mixture-schedule-003-events.csv`. It contains **4,802 shots: 4,721 B1, 46 A, 33 B2 and two C13**, using A-H3/B1-H2/B2-H2.5/C13-H4 at 4.4734×10³² net He/s. All main delivery/capture windows are disjoint. Their finite-period occupancy is **79.007%**, versus **78.988%** in the continuous sizing calculation.

The largest dispatch delay is **5,577 seconds**, or 93 minutes. B1 intervals range from **5,271 to 11,625 seconds** around a nominal **6,772 seconds**; the longer intervals accommodate rare large shots, and later slack restores phase. This verifies one selected chronological dispatch, rather than requiring every B1 interval to equal its long-run mean. `mixture-schedule-003-input.json` archives actual isotope recipes, outputs, timing and plant data, so this calculation can be regenerated without raw hydro traces.

Its ready-stock upper bounds are **3.80×10¹³, 1.31×10¹⁴ and 1.2235×10¹⁵ kg** for one-day, thirty-day and one-year solid holds, with a separate **one-hour local hydrogen hold**. For each species the calculation withholds up to the entire omitted 0.1% pulse mass, with a 1% atomic/baryonic mass margin, throughout the finite horizon. These independently conservative species bounds retain slow tails even if their composition is concentrated. They do not charge fresh replacement for those tails. A steady rotating-skin plant additionally arranges its resident age distribution. No station-power penalty for these tiny additional stock masses is implied beyond their explicit ledger entries.

The original-family seed and tail calculations below remain controls for the earlier compositions. The new mixture dispatch establishes timing and stock arithmetic; its angular assembly, central seed and interaction with outgoing tails still require the selected source/collector implementation.

## Original-mixture control

A chronological representative schedule can accommodate the rare large pulses and the frequent B1 pulses without waiting for every slow particle to reach the outer wall. The saved `mixed-stage-events.csv` contains **4,799 shots over 370 days**: 4,719 B1, 50 A, 28 B2 and two C13. This preserves the integer-shot demonstration of the historical 64-zone rate family. The whole-period counts differ from their continuous expectations by less than one shot each.

## Dispatch rule and observed delays

Nominal periodic jobs are sorted by requested focus time. A new main delivery begins after the preceding main capture window ends; its focus is delayed when necessary. The main windows contain the frozen 99.9% incoming/outgoing mass quantiles. All windows are explicitly checked to be disjoint. Their occupied fraction is **80.008%** over this finite history, matching the old average 80% screen.

Rare pulses cause local delays, not indefinite accumulation of jobs. The longest focus delay is **7,716 seconds**, about 2.14 hours. B1-to-B1 intervals span **5,341–14,490 seconds**, or 89.0–241.5 minutes, while the long-run nominal B1 period remains 6,774 seconds. The longest interval follows a cluster of rare jobs; the queue subsequently uses its spare time to recover the nominal phase. The script saves a four-day excerpt around the first C13 shot as well as the complete CSV.

These nominal phases are a selected demonstration, not an optimum. A steady-dispatch controller can choose other phases or postpone individual rare shots, provided isotope stock and the eventual material ledger are respected. New hydrogen ratios alter stage attempts and pulse counts; `material_cycle_compare.py` supplies the inputs for a new dispatch run.

## The incoming core is an assembly problem, rather than a forced long flight

A literal backward extension of the initial radial flow gives an arbitrarily long launch tail as central speed approaches zero. That extrapolation need not define the launcher architecture. The example assigns the slow inner **0.1% of pulse mass** to faster opposing seed packets assembled near focus. Their kinetic energy supplies the core's already-budgeted hot internal energy.

For the original large stages, this core weighs **1.13×10¹⁰ kg**, has hot-snapshot radius **9.44–9.46 m**, and contains **3.55–4.35×10²¹ J internal energy**. B1's core is 100 times lighter and contains 4.30×10¹⁹ J. Packet speeds of approximately **0.79–0.88 million m/s** carry that energy, giving 21–23 minute flights from the outer sphere. They can leave after the previous shot. Core kinetic-plus-internal energy is only **0.0117% of the large-pulse preparation budget**, or 0.0144% for B1.

This is an energy-and-timing construction. A multidimensional merger must still produce the required density, entropy and velocity distribution. The fixed-grey optical-depth diagnostic is large for the intended compressed seed, but does not establish how cold packets get there. The construction provides a concrete target for the driver model without adding its kinetic heat twice to the existing prepared energy. See `central-seed-energy.json`.

## Outgoing tails can remain in the chamber and in the inventory

The remaining old material is carried in the material account. It may reach the wall later or be entrained into a later pulse. No mass is discarded at the 99.9% timing cutoff.

`tail_overlap.py` calculates the old tail density from each saved cell's uniform spatial mass density and linear radial velocity, then integrates the density along incoming main-packet trajectories. For the worst adjacent stage pairs in the schedule, the intercepted relative column is at most approximately **9.3×10⁻¹⁰ kg/m²**. If all intercepted relative kinetic energy heats a spherical grain of density 1,000 kg/m³, the largest heat input is approximately **331 J/kg for a 1 mm radius grain**, or 3.31 kJ/kg at 0.1 mm. An illustrative 100 kJ/kg allowance corresponds to about **3.3 µm radius**. Actual grain composition, fields and surface response determine its usable heating budget.

This screen explains why mere kinematic intersection with a dilute tail does not establish a destructive collision. It does not cover ambient warm gas accumulated from a collector, which can have a much larger column; prompt H/He handling is assessed separately. It also excludes the separately assembled central incoming seed. Doubling the quadrature from 12 to 24 points changes the sampled column, heat and impulse integrals by less than 7×10⁻¹⁶ relatively; this confirms integration accuracy within the adopted ballistic reconstruction, not the accuracy of the underlying hydro or geometry.

At the next focus, the largest old tail still inside the sphere is **0.171% of the next small B1 pulse's mass**, following a large A event. The corresponding B1-to-B1 case is 0.00241%. Only a fraction necessarily joins the new compressed fuel. Treating the entire old tail as an additional incoming contaminant is a conservative material-loading sensitivity; it should be passed to the finite-purity model if it affects the selected driver layout. The eventual-return ledger continues to count it.

## Working stocks and startup

A one-day average-flow stock is not enough to stage all rare large isotope batches. The explicit historical-boundary sequence requires roughly **3.54×10¹³ kg** initial ready stock across C12/C13/N14/N15/H when every returned batch is delayed one day, including about 9.88×10¹² kg C13 and 1.36×10¹³ kg N14. This exceeds the old generic 3.87×10¹² kg one-day inventory. Its station-mass effect remains negligible.

Thirty-day and one-year uniform batch delays raise the finite-run initial-stock bounds to about **1.40×10¹⁴ and 1.35×10¹⁵ kg**. Those two script cases deliberately retain the old common delay for every species as an inventory upper-bound comparison. The current architecture assigns long residence primarily to C/N-bearing solids and gives H/He a prompt gas loop. Exact future stocks must use the selected mixture, eventual-return fractions and separate material delays.

Startup isotope stock is separate from steady fresh-carbon consumption. Building a pool of enriched N15/C13 and arranging tile ages needs a commissioning sequence; spreading ordinary recovery over months does not make that stock free. The mass is nevertheless around 10⁻⁸ of the station for the annual case. Source storage, launch energy, retained fuel and receiver material are separate inventories.

## Scope

The CSV establishes chronological main-window compatibility and finite material-buffer arithmetic. The tail calculation supplies a bounded collision-column discriminator. The fast seed supplies a physically explicit energy/timing target. The remaining qualification is assembly of that seed and the main angular flow, together with local gas capture; a universal ten-hour shot spacing is not implied by the archived radial tails.
