# Manufacturing efficiency and graphite reuse ✨

Thomas's premise is a **minimum 50% efficiency, often higher, for a declared manufacturing conversion that raises useful product free energy using electricity or accessible chemical energy**. For a positive useful increase ΔG, use η = ΔG / E_supplied and E_supplied = ΔG/η. Count consumed chemical free energy alongside electricity. The conversion's dissipative losses are already included: adding a second generic manufacturing penalty to ΔG/η would count those losses twice.

For a graphite-to-graphite refurbishment loop, the bulk chemical ΔG can be near zero. Its cost follows the actual disassembly, purification, repair, handling and thermal operations. The premise does not assign this zero-net transformation a universal energy per kilogram. Bond-breaking barriers can be traversed with recoverable work; the sum of all broken bond energies is not automatically dissipated.

Drexler's *Nanosystems* §13.3.8 estimates >99% conversion for a mechanochemical H₂/O₂ power apparatus and suggests similar performance for oxidation of small organic molecules. This is a theoretical apparatus estimate; the enthalpy-scale molecular energy quoted there should not be read as 99% of complete ambient chemical enthalpy becoming external work. §13.3.7 separately estimates manufacturing dissipation, including an illustrative 1.5 MJ/kg mean moietal-processing cost. [Primary book text, chapter 13](https://nanosystems.contact.ms/molecular_sorting_processing_and_assembly/).

Chapter 14 adds 0.5 MJ/kg block assembly and discusses a possible 0.1 MJ/kg computation cost. Its later acetone example states 3.1 MJ/kg free-energy dissipation; the visible preceding component values do not reproduce that total, so this audit retains 2 and 3.1 MJ/kg as distinct sensitivity values. The example consumes oxygen and accounts for entropy; its feedstock chemistry supplies useful energy. §14.5.2 also anticipates inexpensive reuse by disassembling and reassembling micron-scale modules. These are theoretical process examples, not measured costs for the proposed receiver. [Primary book text, §§14.4.8 and 14.5.2](https://nanosystems.contact.ms/molecular_manufacturing_systems/).

New carbon from the canonical methane source can be represented by CH₄(g) → C(graphite) + 2H₂(g). At 298.15 K and standard 1 bar states, JANAF gives ΔG = +50.768 kJ/mol and ΔH = +74.873 kJ/mol. Per kilogram of carbon, this is 4.227 MJ useful free-energy increase and 6.234 MJ enthalpy increase. At 50% efficiency the required usable input is 8.454 MJ/kg C; at 99%, 4.270 MJ/kg C. Alternatively, an explicit 2 MJ/kg dissipation allocation gives 6.227 MJ/kg and 67.9% efficiency; 3.1 MJ/kg gives 7.327 MJ/kg and 57.7%. These alternatives satisfy the premise without adding another efficiency penalty. Both carbon and hydrogen are retained as products. Heat-driven chemistry at other temperatures and partial pressures needs its own balance. [NIST-JANAF methane data](https://janaf.nist.gov/tables/C-067.html).

Useful physical scales for graphite reuse are:

| Operation and stated endpoint | Calculated energy scale |
|---|---:|
| Graphite → monatomic carbon gas, 298 K, standard pressure | 55.89 MJ/kg free energy; 59.67 MJ/kg enthalpy |
| Warm graphite from 298.15 to 1000 K | 0.982 MJ/kg heat; 0.517 MJ/kg thermal availability relative to 298.15 K |
| Retain two basal surfaces per 1 μm graphite lamella, density 2000 kg/m³ | 195 J/kg surface energy |
| Same lamella model at 10 nm | 19.5 kJ/kg surface energy |
| Ideal complete C12/C13 separation, 1.07% C13, 1000 K | 41.0 kJ/kg reversible mixing work |

The monatomic endpoint is an optional energy-intensive route; ordinary reuse can preserve bonds and recover adhesion energy. The heating cost applies to material actually cooled and reheated; hot recirculation and heat recuperation change that duty. Surface numbers use 0.39 J/m² basal cleavage energy, comprising two new surfaces. They describe aligned basal faces, with the material value inferred from measurements near room temperature; edges, radiation damage and hot interfaces require separate properties. The isotope calculation uses −RT Σx ln(x)/molar_mass and is an illustration: the bulk absorber need not undergo complete isotope separation. [NIST carbon gas](https://janaf.nist.gov/tables/C-003.html), [NIST graphite](https://janaf.nist.gov/tables/C-002.html), [Wang et al. primary cleavage measurement](https://arxiv.org/abs/1506.00536).

The selected station's old allowance processes 10⁻⁴ of installed mass per year at 100 MJ/kg. The primary [fixed-hardware sensitivity](sensitivity.json) retains all installed equipment and changes only recurring manufacturing electricity:

| Refurbishment allowance for all components | Net output | Specific power | Increase over old net |
|---|---:|---:|---:|
| 100 MJ/kg | 3.361 × 10²⁰ W | 2.766 mW/kg | — |
| 10 MJ/kg | 3.708 × 10²⁰ W | 3.051 mW/kg | 10.31% |
| 3.1 MJ/kg | 3.734 × 10²⁰ W | 3.073 mW/kg | 11.10% |
| 2 MJ/kg | 3.738 × 10²⁰ W | 3.076 mW/kg | 11.23% |
| Zero-cost arithmetic limit | 3.746 × 10²⁰ W | 3.083 mW/kg | 11.46% |

These are conditional bookkeeping sensitivities. A whole-plant 2 MJ/kg allocation still needs a material-and-process recipe. The attached `energy_examples.py` separately evaluates a secondary mass-reallocation sensitivity, using the old 1 kg/W manufacturing mass allowance and resized radiators. That produces 3.0775 mW/kg at 2 MJ/kg, versus the preferred fixed-hardware 3.0763 mW/kg. Lower energy consumption alone does not establish a smaller factory, so this secondary result is not used for the selected operating point. Applying 2 MJ/kg only to returned graphite, keeping other costs and installed mass fixed, saves 2.361 × 10¹⁹ W, about 7.0% of previous net output.

The distinction between gross handling and fresh synthesis is especially large here: the assumed receiver refurbishment flow is 2.409 × 10¹¹ kg/s, whereas the entire net structural-carbon makeup allowance is only 3548.6 kg/s. Producing that makeup from methane at 50% conversion efficiency costs 3.00 × 10¹⁰ W. Applying fresh methane chemistry to every kilogram of circulating graphite would incorrectly charge for carbon already present. Conversely, initial construction from methane must pay that upgrade once for new carbon. The 2 MJ/kg reuse sensitivity cannot establish a 2 MJ/kg initial construction cost from methane.

Transport should retain its own ledger. At 30 km/s, kinetic energy is 450 MJ/kg. A separately selected 50%-efficient acceleration process would draw 900 MJ/kg; braking recovery, orbital trajectory and which material is transported determine the net journey cost. Efficient local recycling can greatly improve the fabrication term while long-range bulk delivery remains significant.

Recommended bookkeeping: separate fresh-feed chemical upgrade, retained-material refurbishment, net escape/replacement, local handling, and external delivery. Select concrete processing dissipation consistent with the setting and avoid a compulsory atomization path. Retain material-throughput, damage, collection and mass constraints independently of energy efficiency.

Reproduction: run `PYTHONDONTWRITEBYTECODE=1 python3 docs/protium-fusion/manufacturing-efficiency/energy_examples.py` from the dedicated worktree. The JSON records source rows and hashes of the frozen station inputs and code. The original 100 MJ/kg station point is reproduced to relative 10⁻¹²; every sensitivity passes the existing mass and global/electrical energy closure checks to relative 10⁻¹⁰. No original station packet or article was changed.
