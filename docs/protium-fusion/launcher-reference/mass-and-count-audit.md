# Launcher count and mass allocation audit ✨

A 1,000-complex reference fits the existing source account as a partition of total capacity. The present calculation supplies no optimum count and no component design proving that count. Each complex may deliver many independently controlled packets through many internal channels. The number of complexes, independently controlled apertures, simultaneous beam directions and physical packets must remain separate design quantities.

The audited point is `high_source_rating_surface` in [operating-cases.json](../open-loop-redesign/operating-cases.json): radius 1.1×10⁹ m, installed mass 1.21520×10²³ kg. Its source model assigns 100 W/kg peak electrical handling, 0.1 kg/W mean electrical conditioning, and 100 kJ/kg stored electricity. These are selected engineering requirements. They give 2.83889×10²¹ kg peak-source equipment, 2.31763×10¹⁹ kg mean conditioning and 1.35813×10²¹ kg two-event storage: **4.22020×10²¹ kg altogether**. Generic station services, heat conversion, reserve and support are additional existing allocations. None has been resolved into launcher-specific hardware.

Equal sharing of every event gives:

| Number of source complexes | Peak electrical W per complex | Mean electrical W per complex | Assigned source + conditioning + storage kg per complex | Stored electrical J per complex |
|---:|---:|---:|---:|---:|
| 100 | 2.839×10²¹ | 2.318×10¹⁸ | 4.220×10¹⁹ | 1.358×10²⁴ |
| 1,000 | 2.839×10²⁰ | 2.318×10¹⁷ | 4.220×10¹⁸ | 1.358×10²³ |
| 10,000 | 2.839×10¹⁹ | 2.318×10¹⁶ | 4.220×10¹⁷ | 1.358×10²² |

The total assigned mass and power remain equal across these rows because the same capacity is divided among N complexes. This does **not** hold the individual complex mass fixed. At N=1,000, each complex supplies 1.131×10¹⁰ kg to a large event or 1.131×10⁸ kg to a B1 event; its time-averaged circulating fuel share is 44,821 kg/s. Its largest-event electrical allocation is 6.791×10²² J. Changing N scales these quantities by 1/N. Event cadence remains common to all complexes in this sharing model; the B1 interval remains approximately 6,774 s.

The mean/peak source-power ratio is 0.0008164. The largest event's total electrical input divided by booked total peak power is 239.2 s, independent of N. That is a lower bound on its source-energy delivery duration at the booked peak; the converged fusion event has a much shorter duration. A proposal to deliver the full source electricity on the fusion timescale would require a different peak-power and storage-discharge account.

For explicit fixed overhead h per complex, additional shared overhead H, and N complexes, use

`M_source(N) = P_peak / r_peak + 0.1 P_mean + E_store / 100000 + H + N h`.

Here h and H represent hardware not already included in the allocated terms: mounts, independent services, local shielding, recoil handling, diagnostics, feed routing or other actual functions. A documented reassignment from existing service or reserve can cover a function once. Increasing N cannot make h disappear. Source-control and handling power likewise needs a stated additional load or an identified share of existing auxiliaries.

As a simple illustration, adding an identical 10¹⁹ kg of previously unallocated fixed hardware to every complex adds 10²¹, 10²² or 10²³ kg for N=100, 1,000 or 10,000. At the corrected 2 MJ/kg manufacturing sensitivity, a mass-only diagnostic gives 3.051, 2.842 or 1.688 mW/kg respectively. Added support, refurbishment, reserve and local operating power lower those figures further. The 10¹⁹ kg value is a transparent sensitivity, not an estimated necessary overhead.

The manufacturing correction is applied through the existing [fixed-hardware sensitivity](../manufacturing-efficiency/sensitivity.json). At a selected effective refurbishment cost of 2 MJ/kg, P₀=3.73835×10²⁰ W; this cost remains a process-allocation sensitivity under the accepted ≥50% efficiency premise. The old 100 MJ/kg value is not treated as necessary. For all-in additional installed mass x and additional electrical load S, the 1 mW/kg boundary in the existing support model is

`P0 − S − r x − k x (2 M0 + Mt + x) = 0.001 (M0 + x)`.

Here r is the repair power per installed kilogram, k=0.001 G v/(8πR²), v=10 km/s, and Mt=3.79457×10²² kg is the station model's equivalent solar tidal load. The quadratic term accounts for incremental support dissipation as mass increases.

With S=0, the naive mass-only headroom is 2.523×10²³ kg. Including additional maintenance and support power reduces it to **2.479×10²³ kg of all-in extra installed mass**. Reserving 20% of that addition, ten-year spares, support streams, extra radiators and the baseline manufacturing-machine mass per repaired throughput leaves **1.979×10²³ kg for bare extra apparatus**. At 1,000 complexes, that is 1.979×10²⁰ kg per complex if every remaining allowance goes to them. Other unresolved plant apparatus competes for the same headroom.

The corresponding bare headroom is 1.969×10²³ kg at 3.1 MJ/kg refurbishment, and 1.908×10²³ kg at 10 MJ/kg. Adding a new standing load equal to 1% or 10% of the original mean driver power reduces the 2 MJ/kg bare headroom to 1.961×10²³ or 1.798×10²³ kg. These are conditional ceilings within the old support geometry. Concentrated launcher loads or changed structural topology can require additional hardware beyond this scalar model.

There is no ledger evidence that 1,000 is too high or too low. A uniform sphere partition at this radius has area-equivalent cap half-angles 11.48°, 3.624° and 1.146° for 100, 1,000 and 10,000 complexes, respectively. These describe allocation areas, not beam convergence angles. A single complex can serve a structured family of trajectories. Conversely, calling a large family one complex does not eliminate its apertures or control hardware.

Count must follow a specific requirement: independently controllable angular structure; source emittance and allowable merger heating; minimum aperture and channel scale; charge compensation; local electrical distribution and recoil; dose, failure isolation and maintenance access; or the acceptable loss of an angular sector when a complex fails. The existing radial calculations and generic source W/kg rating do not solve those requirements. More complexes reduce individual event duty and sector size while increasing repeated fixed hardware. With fixed positive h and no other constraint, the mass equation simply favors fewer complexes; that observation cannot select a physically sufficient count.

**Recommended status:** retain 1,000 as a provisional organizational reference and 100/10,000 as explicit alternatives. An 80% credible interval would require a distribution for aperture, fidelity, reliability and fixed-overhead requirements; none has been established. These three scenarios should not be labeled a credible interval. The immediate concrete design task is to specify one complex's channels and charged/neutral fuel preparation, then assign every function to source, service, reserve or genuinely additional mass and power.

`mass_count.py` reproduces the equal-share allocations and solves the overhead quadratics, saving provenance and checks in `mass-count.json`. Share sums and threshold power/mass equalities close to relative 10⁻¹². Original station inputs and article are unchanged.
