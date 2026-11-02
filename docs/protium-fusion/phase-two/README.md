# Second-phase protium investigation ✨

Thomas authorized a further five-million-token research allowance on 29 October 2026, followed by an article revision giving concrete fictional late-2040s deployments. The allowance is a ceiling. He explicitly chose to let the research determine outcomes, with pilots and failures allowed. The first phase is preserved in [ASSESSMENT.md](../ASSESSMENT.md).

The bounded investigation is complete. [CONCLUSIONS.md](CONCLUSIONS.md) states what the calculations establish and what remains open. [DEPLOYMENT-DECISIONS.md](DEPLOYMENT-DECISIONS.md) records the selected fictional history: the Oriel pilot opens in 2047, the Moravec range in 2048, and both consume net electricity through 2049. Actual measurements, numerical results, engineering assumptions and invented outcomes remain separate here; the article stays in universe.

## Research packets

| Packet | Main result | Reproduction record |
|---|---|---|
| [Driver delivery](driver/REPORT.md) | Charging, heating, aiming, neutralization and low-entropy merger requirements; cold-inflow gain remains unresolved | [Cold compression](driver/COLD-COMPRESSION.md), scripts, compact traces and MANIFEST.json |
| [Pulse transport](pulse-transport/README.md) | Strong hot-core absorption and scattering; cold contact EOS and interface transport remain consequential | Scripts, validation, compact results and source hashes in the packet |
| [Ion dynamics](ion-dynamics/README.md) | Passive self-consistent electrostatic response includes radial ion absorption; finite capture and acceleration remain separate requirements | Scripts, saved response/capture results and manifest |
| [Gas supply](gas-supply/REPORT.md) | About 30% of diffuse incoming H₂ is ionized on its first flight in the chosen source; wall pumping and neutron economy interact | gas_supply.py, independent_gas_check.py and both saved results |
| [Whole-inventory reset](periodic-reset/PERIODIC-RESET.md) | Explicit thermal, coherent and electrostatic restart costs; favorable separated-efficiency case reaches only 0.661 mW/kg | periodic_reset.py, parameter scan and independent reconstruction recorded in the report |
| [Plant audit](plant-audit/README.md) | Coatings and impurities consume neutron margin; heat transport requires substantial mass and power | Three calculation/validation pairs, figure and manifest.json |

Each report identifies its inputs, numerical checks and physical limitations. Hash manifests preserve provenance. Larger orbit and hydrodynamic intermediates are omitted where scripts and compact archived results reproduce the relevant calculation. No packet alone validates an integrated power plant.

The author retained article authorship. Independent agents checked the technical packets, the final physics and economic statements, and the reset calculation. [article-questions.json](article-questions.json) preserves the ten ranked questions prepared before drafting. The current editorial review is [the article review record](../../slop-reviews/articles/protium-fusion.json).

[Final verification](VERIFICATION.md) records independent package checks, calculation validation, editorial review and browser diagnostics. MANIFEST.json binds the final research packet and article inputs.
