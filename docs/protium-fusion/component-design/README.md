# ✨ Current reactor components and unresolved design work

The current component study supports a **conditional plant account near 20.8 mW/kg** using a staged, radiation-trapping CNO cycle. The representative receiver has 350 kg/m² of hot graphite at a radius of 1.1 billion metres. Its source family uses ammonia-like nitrogen-14 feed, hydrazine-like nitrogen-15 feed, and CH₂ polymer for both carbon stages. The [mixture report](../mixture-design/README.md) and [selected point](../mixture-design/selected-plant/point.json) contain the quantitative operating account.

The most consequential remaining task is **assembling cold, finite packets into the density, velocity and entropy fields used by the successful spherical burns**. The calculations now specify several parts of that path and identify where a three-dimensional model is needed. They also separate useful heat capture from recovery of ordinary material.

## Incoming fuel and pulse assembly

The launcher candidate charges condensed packets, accelerates them through independently timed electrical stages, corrects their trajectories, neutralizes them, and merges converging streams near the focus. Charge is a delivery tool. The final reacting material is locally close to charge neutrality; a macroscopic electron deficit supplies neither confinement nor a solution to bremsstrahlung.

A working architecture uses about 100,000 physical transport columns and roughly 2,000 effective angular merger regions. Those quantities describe different parts of the apparatus. Beam-envelope and staged-acceleration calculations give explicit field, bore, timing and neutralization requirements; they do not impose that outlet count as a universal minimum. The selected plant allocates source hardware for a peak electrical rating of 2.86 × 10²³ W and two-pulse storage of 1.36 × 10²⁶ J. [Launcher study](launcher.md).

Carbon-return and carbon-13 polymer contact states are particularly useful. Their selected collision-normal speeds are about 83 km/s, while most prepared energy remains in the larger converging radial motion. An independent primary EOS table agrees with the polymer contact compression and pressure within about 0.13% and 0.044%; its temperature differs by roughly 11%. The resulting carbon-13 source-mapped burn agrees with its prescribed-profile gain within 0.006%. [Polymer calculation](c13-polymer-source/README.md).

The carbon-return stage also has a constructive spherical pressure-work connection between the contact state and the hot profile. Its nonseed mass reaches contact roughly 123–249 microseconds before the hot snapshot, at radii of hundreds of metres. The nitrogen-15 calculation connects 82.3% of the mass on the corresponding simple branch; the outer 17.6% needs an expansion/contact sequence. Central seeds remain a separate joining problem. These are derived trajectory scales, rather than manufacturing tolerances. [Assembly audit](assembly-audit/README.md).

Continuing the nitrogen preparation through radiation transport and then the full burn reduces nuclear output by 0.911% in the principal opacity comparison, relative to the same delivered control without preparation radiation. Its actual precursor energy, already escaped radiation, outgoing composition and receiver timing all enter the final family. The model handoff between two EOS implementations remains separately documented. [Delivered burn](assembly-burn/README.md), [complete family](final-family/README.md).

## Capsules, surface optics and neutralization

Bulk phase changes provide useful thermal margin. The capsule study allows melting and finite pressure. Cold polymer can absorb the reference flight dose while warming toward its contact-reference temperature; nitrogen recipes can use latent heat. Methane and ethane require explicit pressure and headspace calculations, and their carbon vessels change the delivered H:C ratio. These results favor the simpler polymer carbon-13 feed when the small whole-plant power differences are included. [Capsule study](capsule-survival/README.md).

The important remaining capsule quantities are the manufactured reflector's absorption, thermal contact with the core, seams, strength, permeability and the deposition of neutralizing-electron energy. The ordinary smooth-aluminum optical reference has 1000 K-spectrum absorptivity about 0.022. A **subjective 80% design range of 0.007–0.06**, conditional on a manufactured reflective surface, expresses the uncertainty used in this research. It is separate from the underlying optical measurement. Grain size, protected flight distance and accepted phase change can be varied if the absorbed dose exceeds the selected material's enthalpy allowance. [Radiation interface](radiation-interface.md).

Mean energy alone does not determine coating survival. Concentrating the stored electric energy in a thin skin can produce a severe local heat load even when the core's average temperature is acceptable. Electron escape, slower discharge, conduction into the core and a different skin are concrete alternatives to calculate. A protected exit also needs a local collector or shield for outgoing pulses; its small global intercepted fraction can coexist with a large local fluence.

## Outgoing material and radiation

A closed, massive receiver catches fuel products without requiring every ion to be electromagnetically braked first. Implantation, slow tails, sorting rejects and worn surfaces create retained inventories and processing streams. Conversion of catalyst into oxygen changes the usable isotope pool while keeping matter inside the plant. The final fresh-carbon account follows those distinctions. [Material ledger](material-reproducibility.md).

Hydrogen and helium need local collection and return paths so that gas does not accumulate throughout the central cavity. Carbon and nitrogen can be recovered through surface reworking and separation; the study includes batch catcher processing as a candidate. Every gram does not need to be atomized or separated in the same operation. A phase-resolved receiver/recovery design still needs actual retention, residence, purity and equipment specifications.

Photon transport now follows photoelectric absorption, Compton scattering, pair conversion, electron and positron slowing, bremsstrahlung and annihilation in a closed carbon sphere. The source spectra and geometry are explicit. At the selected column the principal hard-source scenario retains about 93.57% of the hot-boundary energy; the plant charges the remainder to cold equipment. This heat-routing choice does not imply ordinary matter escapes. An alternative sends the entire saved prompt-radiation account to cold equipment and lowers specific power from 20.78 to 20.51 mW/kg. [Coupled transport](coupled-photon-transport/README.md), [column refinement](coupled-photon-extension/README.md).

The main 370-day schedule contains 4,228 pulses with 68.39% occupied main arrival windows. Nitrogen-15 pulses have a mean interval of 128.3 minutes. The longest queue delay is about 83.5 minutes. Slow material outside the main-window convention stays inventoried. Rare large pulses can be separated by 2.75 hours in the explicit schedule; thermal calculations must use that chronology rather than inventing a daily-shot restriction. [Dispatch](final-family-dispatch/README.md).

## Receiver, structure and cold equipment

Finite-depth photon deposition, conduction through the backing, surface damage and recurring pulse heat have separate calculations. The hot receiver requires compliant panels with overlapping, shielded joints. Panel bow and attachment compatibility matter more than the very small mean mechanical pressure. The exact-family 350 kg/m² thermal screen gives an all-phase front envelope of about 1,639 K, or 1,668 K with the shallow damaged-skin sensitivity, including the actual dispatch-separation bound. A shallow damaged surface layer can be reworked while retaining the carbon; it does not establish a new fresh-carbon loss for every pulse. [Surface damage](surface-damage.md), [exact-family thermal and components](final-receiver.md).

The sphere needs active support. The candidate combines three families of circulation paths, a cold braced frame, prestress and controlled guideways. It replaces the mechanically unstable passive-monolith picture with explicit support loads, electrical losses and local control timescales. Global attitude and local guide models have been checked separately. Flexible guide modes, joints, actuators and actual materials remain apparatus requirements. [Support](support.md), [prestress](receiver-prestress.md), [local control](guide-local-control.md).

Cold superconducting guideways need local radiation shielding, a temperature buffer and refrigeration. A graphite shield with an abundant polymer buffer can fit within the retained service allocation. The energy intercepted by cold equipment is counted before assessing refrigerator load. The calculation does not place all missed receiver heat at the guide's cryogenic temperature. [Guide buffers](guide-buffer-materials.md).

The 1,000 kg/m² general-service allocation remains the largest installed-mass item. Named headers, interfaces, frame, guides and shielding occupy only part of it. Lower total allocations could improve specific power materially, but recovery equipment, mounting details and other services need a complete bill before that reduction becomes a supported component choice. The [16–32 mW/kg allocation sensitivity](../mixture-design/selected-plant/services-sensitivity.json) is a set of conditional cases; no probabilistic mass bound is assigned to unspecified equipment.

## Priorities that can change the design

| Open question | Why it matters | Useful work before a physical experiment |
|---|---|---|
| Finite stream merger and seed joining | Can change delivered density, entropy, ignition and source geometry | Three-dimensional contact-patch calculation with neighboring streams, realistic boundaries and the outer nitrogen sequence |
| Capsule and neutralizer behavior | Can change grain size, coatings, free-flight distance and actual composition | Coupled electron deposition, heat transfer, phase/pressure and stress calculation; propagate optical uncertainty |
| Receiver panels, joints and cold paths | Can change segmentation, shielding and useful heat capture | Finite panel/contact mechanics and a geometric radiation/particle model using the actual dispatch |
| Recovery and service equipment mass | Can move whole-plant power per kilogram by tens of percent | Complete throughput and residence-based component bill with explicit retained inventories |
| Dense-plasma EOS, opacity and spectra | Can move local compression and heat deposition | Consistent constitutive tables and multigroup radiation calculations, with resolution and source comparisons |

Measured optical properties, material strength after irradiation, electron-surface behavior and component lifetime will eventually require experiments. The numerical work can first narrow the geometry and acceptance criteria. Current Monte Carlo intervals describe transport sampling; published rate curves and mesh differences describe their own comparison classes. They are kept separate from subjective design judgments.

## Manufacturing and economic consequences

The manufacturing premise is at least 50% efficiency for a specified useful free-energy increase, often higher. At the selected specific power, each 1 MJ/kg of actual total embodied energy requires 1.525 export-years. The receiver alone contains 5.32 × 10²¹ kg of carbon; sourcing that stock at C:H = 0.007 corresponds to access to 7.60 × 10²³ kg of associated hydrogen. That is a construction-resource scale, distinct from the operating fuel ratio.

Construction throughput and factory growth can therefore matter as much as conversion efficiency. Long service life improves the energy account; ordinary refurbishment already has its own throughput, equipment and electrical allocations. Illuminated solar power and deuterium-based systems remain economic competitors whose complete delivery and operating systems should be compared. [Selected economics](../mixture-design/selected-plant/economics/README.md).
