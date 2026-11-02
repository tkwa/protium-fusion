# ✨ Unretuned source128 family

All four source-shaped A H3 / B1 H2 / B2 H2 / C13 H4 pulses complete at 128 cells. Each starts from a nested remap of the exact completed source64 state, retaining composition, piecewise density and temperature, total mass, total prepared energy and physical integration duration. No contact energy or source entropy was retuned. The median carbon table and runtime kernel are identical across the four source runs and their nominal128 controls.

The source128 family exports **6.61734 mW/kg** in the conditional 5000 kg/m² receiver budget. The separate 2000 kg/m² mass case gives **12.46663 mW/kg**, retaining the same assumed 99% capture for comparison. Relative to the same-recipe nominal128 hot-flow family, the source profiles cost 0.3503% of net power per kilogram. These figures describe the declared source/plant model; receiver containment and three-dimensional source realization remain separate qualifications.

| Stage | Source64 Q | Source128 Q | Nominal128 Q | Source versus nominal128 |
|---|---:|---:|---:|---:|
| A H3 | 10.252517 | 10.764814 | 10.817796 | −0.490% |
| B1 H2 | 7.520750 | 8.082346 | 8.261904 | −2.173% |
| B2 H2 | 10.275337 | 9.918671 | 9.917629 | +0.0105% |
| C13 H4 | 22.787949 | 23.293987 | 24.698571 | −5.687% |

The 64→128 changes include increases and decreases. B2's decline is shared by its nominal128 control and is not a new source-specific failure. C13 remains below the 2.8 GK numerical/EOS guard, reaching 2.2557 GK. A single refinement does not establish spatial convergence; these are explicit separate families rather than a silent replacement of the 64-cell evidence.

## Integrated point

| Receiver mass assumption | Net power | Installed mass | Net specific power | Main-flight occupancy |
|---|---:|---:|---:|---:|
| 5000 kg/m² | 8.069393×10²⁰ W | 1.219431×10²³ kg | 6.617345 mW/kg | 0.739829 |
| 2000 kg/m², separate comparison | 8.081418×10²⁰ W | 6.482440×10²² kg | 12.466631 mW/kg | 0.740335 |

Both rows use radius 1.1×10⁹ m, 200 K radiators, 60% thermal-package conversion, 60% electrical-to-prepared efficiency, the frozen maintenance/hardware assumptions and 99% capture. The radiator fills the selected 90% area allowance. Reducing receiver mass does not itself prove equal capture, acceptable local temperatures or lifetime; the thinner-receiver transport, support and thermal packets determine whether that assumption can be retained.

The cycle uses 1.91758735 MeV prepared per net helium nucleus and releases 24.71568142 MeV non-neutrino nuclear energy. Source preparation raises energy per helium 1.1556% and gross fuel circulation 1.2436% relative to nominal128. The fuel-network fresh C:H mass ratio is 0.00031633, leaving most of the user's 0.007 allowance for actual apparatus losses. Ordinary material held for processing remains inventory; this ratio does not measure net erosion or escape from a receiver or capsule.

## Remap and verification

`refine_source.py` splits every original cell at the midpoint of its mass-radius coordinate. It interpolates radius cubed in enclosed mass, copies the parent temperature and interpolates velocity in physical radius. A global velocity correction around 0.9995 restores the original total internal-plus-kinetic energy. Maximum density discrepancy is about 1.2×10⁻¹⁴, temperature is unchanged, and integrated initial energies agree near machine precision. This compensates the changed staggered-node kinetic quadrature; it is not a new optimization of the source profile.

`source128-stage-map.json` and `nominal128-stage-map.json` identify unchanged compact records, exact initial arrays, input metadata, full raw trace paths and hashes. `trace-audit.json` passes 56 checks; its largest sampled energy residual is 4.03×10⁻⁶ of prepared energy. `integration.json` passes 72 checks on energy, mass, composition, physical duration, receiver reduction and plant accounting, and embeds the reduced inputs. `package.py` reproduces this archive and comparison. The two startup failures in B2/C13 were corrupt standard-library bytecode imports before simulation; the logs were preserved and the unchanged runs completed using a fresh bytecode cache prefix.

Some input free-text descriptions are inherited from the source64 contact construction. The explicit result mesh, refinement metadata and rate metadata describe these completed128 runs. The original source64 packet is unchanged.

## Remaining interpretation

These calculations establish that the chosen nonideal-contact entropy mappings still lead to high modeled burn on a finer radial grid. They do not evolve the actual finite packets from first contact through the intervening pressure work and shocks into the hot initial state. In particular, a parcel's final internal-plus-radial energy is a diagnostic allocation, not an advected invariant through pressure work. The separate assembly audit is testing the contact-to-hot connection.

The source-power and timing reducer also still projects the prescribed hot profile back to source bins; it does not prove a complete cold acceleration/neutralization/merger waveform. Finite outlets, slow-core seeds, capsule impurities and thermal protection retain their explicit component requirements. The modeled power margin means these tasks should be assessed by their actual energy, mass and burn effects, rather than an imposed 95% pulse-yield cliff.

B1's additional H1.5/H2 source-rate tests are archived under `component-design/b1-source-robustness/`. The controlled H2 source eliminates the weak-burn outcome of the nominal prescribed H2 profile under the low N15 curve. Completed high-curve controls also retain strong burn using the same profiles without retuning. The six-row matched envelope is in `b1-source-robustness/source128-envelope.json`; it does not alter this fixed median-carbon family.

The schema-only interchange files `source128-reduced-inputs.json` and `source128-sigma5000.json` / `source128-sigma2000.json` expose the same archived values to the existing photon and receiver-tradeoff tools. `interchange.py` derives only missing H/catalyst labels from the exact initial isotope abundance; no physical field is changed.
