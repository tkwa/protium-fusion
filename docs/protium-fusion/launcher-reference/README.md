# 1,000-launcher reference ✨

Use **1,000 launcher complexes** as the provisional reactor reference. A launcher can emit many packets; participation, channel count and firing sequence remain design variables. The count has not been qualified by a discrete-source implosion simulation.

The existing receiver radius is 1.1×10⁹ m. An equal allocation gives each complex about 1.13×10¹⁰ kg of fuel and catalyst per large event, 2.84×10²⁰ W peak electrical preparation power and 4.22×10¹⁸ kg of the existing source, conditioning and storage equipment. These are divisions of the total station account. Ports, shielding, mounts and independent controls must receive explicit allocations from existing services or additional mass and power.

An equal-area layout has a characteristic neighboring-direction scale of 6.42° and a surface spacing scale of 123,310 km. The equivalent service-cell radius is 69,570 km. Those dimensions describe the region served by a complex. Actual outlets follow the packet geometry and source optics.

The reference distinguishes two installation arrangements:

1. **Chamber-facing sources supplied through protected feed systems.** Slow material can pass through sealed, staggered feed locks before acceleration. The high-speed inward path can begin on the cavity side of the receiver. Local source shielding and closure still have to withstand outgoing radiation and fuel.
2. **Sources behind the receiver firing through it.** Each high-speed path needs a clear aperture through the receiver and any membrane. Its opening interval, closure mechanism, radiation path and material escape enter the account explicitly.

For 1,000 identical small radial circular ports of radius a, the isotropic open-area fraction is approximately 1,000a²/(4R²). At a=100 km this is 2.07×10⁻⁶. If the aggregate open area stays fixed, dividing it among more ports does not change this fraction. Local shielding dose and outgoing-material interception can still dominate the hardware design. The whole fuel-loss allowance is already shared with collection and nuclear losses; it cannot be reassigned to holes without revising that account.

The comparison set is 100, 1,000 and 10,000 launchers, with intermediate counts where useful. Total pulse inventory, prepared energy and intended radial profile remain fixed in those comparisons. More packets from the same compact site refine radial and temporal sampling while leaving nearly the same incoming direction. Additional spatial coverage requires distributed outlets, more sites or an explicit steering/merger system.

## Evidence and acceptance

The frozen pulse models are spherical. Their 64, 128 and 256 radial zones are spatial resolutions, not launcher counts. They cannot infer a minimum number of discrete directions. Related plasma-liner experiments and simulations motivate multidimensional merger and stability calculations; they do not supply a transferable minimum count. See the [primary literature assessment](literature.md).

The proposed initial comparison retains at least 95% of the spherical reference nuclear yield at the same prepared energy, then checks the full CNO/H return and whole-station power account. This is an explicit screening criterion. The accepted supply budget remains 0.007 kg fresh carbon per kg fresh hydrogen, and 1 mW/kg remains the useful minimum station target. Manufacturing follows the [corrected efficiency premise](../manufacturing-efficiency/README.md).

Report central estimates and subjective 80% credible intervals for uncertain tolerances where a physical acceptance criterion makes them meaningful. State the assumed packet, source geometry and error correlation. Exact geometry, selected sensitivity values and inferred probability intervals have separate roles. The 100–10,000 comparison range is a design sweep, not a credible interval for a minimum count.

## Calculations

- [Source mass and count audit](mass-and-count-audit.md), [script](mass_count.py) and [results](mass-count.json): equal shares, explicit fixed overhead, and conditional mass/power headroom.
- [Ports and shielding](ports-and-shielding.md), [geometry script](geometry.py) and [results](geometry.json): angular allocation, aperture losses and protection geometry.
- [Incoming timing](incoming-timing.md): packet count, distinct source directions, preparation and arrival tolerances.

The receiver comparison uses this same reference in the adjacent [membrane assessment](../membrane-receiver/README.md). It does not grant a smaller receiver, a higher source rating or faster pulse cadence merely by choosing a launcher count.

[Driver physics and cold-compression studies](driver-physics/README.md) retain the first-contact, neutralization, source-timing and numerical-failure research. Their earlier pulse sizes and material assumptions are recorded alongside the current launcher reference.
