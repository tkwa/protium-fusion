# ✨ Warm return gas and cold-packet exposure

The circulating gas requires a local return interface even in a perfectly closed, mass-conserving station. The useful criterion is encountered gas column and the resulting grain heating; gross circulation is not fresh-feed loss. Helium product belongs in this calculation alongside unburned hydrogen.

`return_gas.py` evaluates explicit ledger scenarios at 3.78224727×10³² He/s and radius 1.1×10⁹ m. The reference returns 1.85904845 H atoms/He; the A-H2 plus B1-H2 candidate returns 3.817. These source strengths are inputs supplied by the nuclear-ledger work, not a new network calculation. Hydrogen source rises from 1.177 to 2.415 million kg/s. Helium product contributes another 2.514 million kg/s.

For uniform warm gas at temperature T, mean molecular speed v̄, surface capture fraction g, and fraction f of the source entering the central cavity,

- Σ = 4 Q R f/(A v̄ g), radial mass column;
- residence time = 4R/(3v̄g);
- P = Σ kT/(mR).

The capture fraction means removal into a return/product line, rather than merely hitting a 1000 K wall. The equations follow directly from Maxwellian wall flux nv̄/4 and cavity volume AR/3. A uniform Maxwellian is also a stationary collisional equilibrium. Gas-gas collisions do not automatically add a diffusion penalty to this uniform steady solution; localized sources, heterogeneous traps and directional plumes need transport. Cross-section sensitivities in the JSON flag that the completely injected cavity need not be collisionless.

At 1000 K and unit capture, residence is approximately 5.24 days for H₂ and 7.38 days for He. Complete source injection gives:

| Gas scenario | Radial column | Total pressure | Full-retention heat into a 1 cm-radius, 1000 kg/m³ grain at 2.4 Mm/s |
|---|---:|---:|---:|
| Reference H₂ + He | 4.213×10⁻⁷ kg/m² | 9.91×10⁻¹⁰ Pa | 91.0 kJ/kg |
| Richer-mixture H₂ + He | 5.318×10⁻⁷ kg/m² | 1.406×10⁻⁹ Pa | 114.9 kJ/kg |

Atomic hydrogen reduces these respective columns to 3.905×10⁻⁷ and 4.687×10⁻⁷ kg/m². The very low pressure does not itself establish acceptable packet heating: the flight is a gigametre long and the relative speed is millions of metres per second.

The conservative energy-deposition criterion is Δu = 3Σv_packet²/(8ρ_grain a). Actual retention, surface ablation and ion penetration need the grain material model. The following are three **chosen acceptance allocations**, not measured tolerances or credible intervals:

| Gas heat allocation, 1 cm grain | Maximum common f/g, reference H₂ + He | Maximum common f/g, richer mixture |
|---|---:|---:|
| 1 kJ/kg | 0.01099 | 0.008705 |
| 10 kJ/kg | 0.1099 | 0.08705 |
| 100 kJ/kg | 1.099 | 0.8705 |

A value above one means that complete injection with unit capture passes that particular bulk-energy test. For 1 mm grains the permissible f/g is ten times smaller. All limits scale linearly with grain density and radius, linearly with the allocated heat, and inversely with the square of packet speed. The calculations include 300, 1000 and 2000 K gas cases: the steady column scales as T⁻¹ᐟ² while pressure scales as T¹ᐟ².

At the 10 kJ/kg benchmark and g=0.01, 1 cm grains require local return to keep approximately 99.890% of reference gas or 99.913% of the richer-mixture gas from entering the main cavity. This is a specification to test against a design. A source directly connected to a closed backplenum can maintain gas inventory there for much longer without creating this chamber column. A front-facing hot oxide alone does not supply that connection or retain H/He indefinitely. A combination of local capture channels, selective membranes, pumping or colder collectors can supply removal; global cryocondensation is one option rather than a requirement.

The source change raises the combined H₂+He column by 26%, despite roughly doubling the unburned hydrogen. A baseline local recovery system with adequate margin may accommodate this without a different architecture. The substantial-change event is failure of its actual gas-routing/packet-heating criterion, not the mere presence of excess hydrogen. Current data establish a finite throughput and isolation target; they do not establish the mass or power of the recovery hardware. Eventual recovery, temporary holdup, central-cavity contamination and actual escape through openings remain separate accounting categories.

Reproduce with standard-library Python: `python return_gas.py --output return-gas.json`. Twenty-seven independent inventory/column identities pass. The JSON records formulas, source assumptions, temperatures, both hydrogen states, both grain sizes and all three heat allocations. This model assigns no probability distribution to unspecified trap efficiency or sputter production.

The later **provisional** A-H3/B1-H2/B2-H2.5/C13-H4 integration uses 4.473×10³² He/s,3.997 million kg/s returned H and 2.973 million kg/s He product. At1000 K its fully injected H₂+He column is 7.308×10⁻⁷ kg/m² and pressure 2.045×10⁻⁹ Pa, giving 157.9 kJ/kg to the 1 cm reference grain. The1/10/100 kJ/kg allocations permit f/g=0.00633/0.06335/0.6335. At the 10 kJ/kg allocation and g=0.01, local routing must prevent 99.9367% from entering the cavity. This is a parameterized consequence of the provisional ledger; it is not a claim that the integrated point has an engineered gas-return system.

This screen concerns thermalized neutral H/H₂/He after collection. Prompt charged outflow, sheath acceleration, recombination, charged-gas confinement and schedule overlap are separate transport questions. A long residence in a closed local plenum does not equal a long residence in the central delivery cavity.
