# ✨ Coherent source-profile family with median carbon rates

The completed source tests now form one explicit A H/N=3, B1 H/N=2, B2 H/C=2, C13 H/C=4 family. All four use 64 radial cells, the full 27-species network, the same median LUNA carbon-capture update, FastEOS, 0.001 m²/kg grey opacity, compression-cap viscosity and the original physical integration durations. Their total pulse masses, compositions and prepared energies equal the corresponding nominal64 reference cases. `stage-map.json` names every record, exact initial profile, input record, raw trace and hash. Earlier source maps and the unsuccessful B2 geometry control remain separate.

The resulting four-stage material ledger needs **0.38475% more preparation energy per net helium nucleus** than the same-recipe nominal64 family. Under the stated common plant assumptions, capacity specific power falls **0.1225%**. This small whole-cycle change is more informative than interpreting a fixed percentage of single-pulse yield as a pass/fail boundary.

| Stage | Nominal Q | Source-profile Q | Change | Source construction represented at the hot snapshot |
|---|---:|---:|---:|---|
| A, H/N=3 | 10.288815 | 10.252517 | −0.3528% | N=2000 local collision-entropy mapping plus controlled initial preheat floor, paid from radial kinetic energy |
| B1, H/N=2 | 7.545426 | 7.520750 | −0.3270% | Nearly constant contact-shock entropy from controlled opposing-stream normal velocity |
| B2, H/C=2 | 10.271886 | 10.275337 | +0.0336% | Common CH₂ contact-shock entropy; original density profile and total energy retained |
| C13, H/C=4 | 23.862778 | 22.787949 | −4.5042% | N=2000 local methane collision-entropy mapping |

Here Q is deposited prompt nuclear energy divided by prepared pulse energy. The C13 desired-product statistic counts N14 at the endpoint; it rises 21.4% in the mapped case even as Q falls, because competing further capture changes N14/O15 routing. `integration.json` uses the full isotope outputs and subsequent cold decays. It does not apply a scalar Q multiplier to the nominal material loop.

## Controlled plant comparison

These are conditional accounting points at receiver radius 1.1×10⁹ m, 200 K radiators, 60% thermal-package conversion, 60% electrical-to-prepared efficiency, 99% energy capture, complete eventual ordinary-material return, maximum 80% main-flight occupancy and 90% radiator area. The frozen hardware/maintenance allocations are retained. The thinner receiver is an explicitly separate mass sensitivity; its stopping, photon capture, lifetime and support must be checked in the receiver packets.

| Family / receiver | Net power | Installed mass | Net specific power | Main-flight occupancy |
|---|---:|---:|---:|---:|
| Nominal64 / 5000 kg/m² | 7.994314×10²⁰ W | 1.219200×10²³ kg | 6.557016 mW/kg | 0.788891 |
| Source64 / 5000 kg/m² | 7.984557×10²⁰ W | 1.219205×10²³ kg | 6.548984 mW/kg | 0.798214 |
| Nominal64 / 2000 kg/m² | 8.006287×10²⁰ W | 6.480133×10²² kg | 12.355126 mW/kg | 0.789431 |
| Source64 / 2000 kg/m² | 7.996523×10²⁰ W | 6.480186×10²² kg | 12.339958 mW/kg | 0.798760 |

The source family has U=1.98313135 MeV prepared per net helium, versus 1.97553040 MeV nominal, and H=24.70866590 MeV non-neutrino nuclear release per helium. Gross fuel circulation rises 0.40277%. At the 5000 kg/m² point, B1 fires once per 6763.2 s; its main incoming/outgoing window is 5319.5 s. The combined 0.7982 occupancy nearly reaches the selected 0.8 scheduling allowance. This is a capacity constraint in this accounting scheme, rather than a universal maximum reactor duty cycle.

The fuel-network fresh carbon/hydrogen mass ratio is 0.00032417, versus nominal 0.00035864. These values include network-held material outside the four-stage routing and distinguish retained ordinary material from permanent escape. They leave room below the user's 0.007 total ongoing-supply target for actual receiver, capsule and port losses; those apparatus losses have not been derived by this integration.

## What source formation has established

The simulations establish that each prescribed hot profile reached by these **conditional contact-entropy constructions** can burn at roughly the required gain in the one-dimensional model. B1 needs a different controlled contact geometry from A/C13: its pulse density is 100 times lower, and transferring the same N=2000 parallel-sector prescription produced much lower entropy. The selected B1 contact heat is 1.2394×10¹¹ J/kg, 4.3837% of the total prepared pulse energy, supplied by conversion of already budgeted stream kinetic energy. Each symmetric stream has normal speed 4.9788×10⁵ m/s. B2 uses 3.4510×10⁹ J/kg, 0.095715% of its prepared energy, and normal speed 8.3079×10⁴ m/s. These are collision-energy choices, not added external flash-heater loads. Their conditional local parcel energy allocations correspond to B1 half-angle mass quantiles of about 9.4°, 12.7°, 18.0° and B2 1.38°, 1.86°, 2.63° near contact. They do not prescribe how much each remote launcher points away from the centre.

The completed B2 N=2000 control yielded Q=0.17636. Only 0.0824% of its mass reached 0.5 GK, despite its 1.40 GK peak temperature. The common-contact-entropy profile heats 95.4% of the mass above 0.5 GK and yields Q=10.27534. This is a profile-shape failure and a constructive model-level remedy. It demonstrates why neither peak temperature nor pulse-mean entropy alone specifies the driver.

The A controlled-preheat case replaces an earlier numerical temperature-guard termination. The 2.8 GK guard was an EOS/numerical operating boundary, not a demonstrated physical reactor temperature limit. Both records remain archived; no accepted trajectory was dynamically clipped at the guard.

## Remaining source timing and geometry

The integrated peak source rating, 2.83927×10²³ W, and stored electrical energy, 1.35813×10²⁶ J, come from the existing exact-hot-profile source adapter. That adapter accounts for the actual input energy and radial velocity bins. It does not solve the cold acceleration, neutralization, finite-aperture delivery, collision sequence and energy-bearing central seed together. B1's separate cold point-contact proxy gives about 2.17×10²¹ W, 15% above its nominal hot-profile estimate and still well below the stagewise B2 peak; this proxy has a stated local-energy allocation and is not a delivered waveform. Replacing slow central parcels with faster colliding seeds must also use their already budgeted energy.

The main unresolved source question is three-dimensional assembly: finite packets and gaps must create the contact density, pressure, entropy and predominantly radial velocity field assumed here without uncontrolled mixing or early expansion. A useful next computation is a finite contact patch with the selected composition, three or more incident directions, controlled normal velocities, actual spacing/packet material, and a nonideal contact EOS. It should return spatial entropy, density and velocity covariance and their correlations, then propagate those distributions through burn and material/plant accounting. A specified 95% single-pulse yield retention or a fixed outlet count is not required to make that comparison useful. The completed spatial128 refinement is archived separately under `component-design/source-family-carbon-median128/`, with no contact-energy retuning.

## Reproduction and verification

Run `integrate.py` with the pinned research Python environment. It uses `mixture-design/variant_integrate.py` for the unchanged adapters and archives reduced inputs sufficient for `--replay integration.json` without raw hydro. `nominal64-stage-map.json` explicitly records the comparison controls. The integration passes 76 checks covering composition, mass, energy, physical duration, receiver reconstruction and plant closure. `trace-audit.json` independently passes 78 checks on the four selected trajectories and excluded B2 control. An independent read-only plant review reproduced all four rows' mass, electrical output, global energy, radiator rejection and stage pulse intervals within 2.3×10⁻¹⁶ relative numerical error.

The nominal64 rate source hash 56e07d… and source64 hash 9e3252… differ only by initialization and an assertion that both carbon tables use identical temperature nodes. The exact archived code diff changes no rate construction for these identical tables. Each hash remains recorded; the two carbon-pair updates do not imply that other nuclear-rate uncertainties have been resolved. Surface thermal proxies are diagnostics and do not qualify real capture or material lifetime.
