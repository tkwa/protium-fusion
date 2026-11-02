# Coupled receiver transport at 200 and 350 kg/m² ✨

This extension adds 90,000 completed Geant4 histories to the frozen [coupled electromagnetic transport packet](../coupled-photon-transport/README.md). It uses the same Geant4 11.4.1 implementation, datasets, option4 electromagnetic physics, standard two-photon in-flight annihilation, 10 µm production cuts, natural graphite, and closed spherical geometry. At each new column, 20,000 isotropic photons are sampled at 2.5438523 and 5.657344 MeV, and 5,000 radial photons at 20 MeV. Delayed sources begin 25 µm within the shell; prompt sources enter from the cavity. This remains a conditional static response model with the source-position, angular-distribution, photonuclear and material-collection limits documented in the parent packet.

`combined-basis.json` is the composable handoff to the existing `reweight.evaluate` API. It contains the original 250, 500, 750, 1000 and 2000 kg/m² columns unchanged, with their event paths redirected to the neighboring archive, plus the new 200 and 350 columns. The two new columns support the broad delayed-source scenario; they do not contain the nominal-line or primary-positron kernels. `basis-extension.json` contains only new columns. Compact run metadata, actual-cut logs and lossless per-event energy tallies are archived under `records/`, `logs/` and `events/`. `run-provenance.json` records the raw file, code and command provenance. `reduce_extension.py` is the original raw extraction; `verify_archived.py` validates the saved packet without raw simulation files.

For the delivered-B1/C13H2 source128 family, the estimated fractions leaving the hot graphite are:

| Column, kg/m² | Broad delayed + EM20 prompt | All saved prompt assigned to cold | 80% sampling half-width |
|---:|---:|---:|---:|
| 200 | 10.5654% | 10.8727% | 0.0673 percentage points |
| 250 | 8.8485% | 9.2544% | 0.0651 percentage points |
| 350 | 6.4290% | 6.9881% | 0.0592 percentage points |
| 500 | 4.1436% | 4.8978% | 0.0499 percentage points |
| 750 | 2.1199% | 3.1097% | 0.0365 percentage points |
| 1000 | 1.1412% | 2.3198% | 0.0267 percentage points |
| 2000 | 0.1482% | 1.6813% | 0.0086 percentage points |

The denominator is the hot-stage non-neutrino nuclear plus prepared energy, H+U_hot. B1 preparation radiation is carried separately by the final-family plant adapter. Values include the explicitly inherited 6.4 ppm hot-bypass sensitivity; no new capture-flow calculation is implied. The all-prompt policy replaces the prompt transport response once. No additional 0.5% missed-heat allowance is included. Fractions may be assigned to cold equipment or external escape by a separate geometry-dependent plant scenario. Monte Carlo sampling intervals do not represent uncertainty in that geometry, source placement or omitted physics.

The new histories pass all 30 initial event-accounting checks. The archive verifier also reopens every original and new compressed event record, checks the means and standard errors, verifies complete tracking and energy closure, and confirms that all five original columns are unchanged. The response table enables a nearby-column plant comparison; it does not by itself identify a global optimum receiver thickness.

Run `python3 verify_archived.py` from any directory to repeat the archive audit. The sibling original packet must remain present because the combined basis intentionally references its event records.
