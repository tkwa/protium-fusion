# Coupled electromagnetic receiver transport ✨

This packet follows photons, electrons and positrons with **Geant4 11.4.1**, through a closed natural-carbon shell. The exact source128 inventory gives **95.880% hot capture at 500 kg/m²** in the broad hard-photon scenario. Assigning the entire saved prompt-radiation account to cold equipment gives **95.148%**. These are conditional radiation-loss calculations within the existing heat ledger; the surrounding material collection and plant remain separate calculations.

The calculation is an actual coupled electromagnetic cascade, including photoelectric absorption, Compton and Rayleigh scattering, pair conversion with recoil, electron/positron ionization and multiple scattering, bremsstrahlung, and annihilation. It uses the standard two-photon in-flight annihilation model after a reproducible conservation failure in the installed optional three-photon model was isolated. That approximation is explicit below.

## Results and source definitions

The table uses the coherent source128 quartet and its **26.633268776 MeV per net helium** total non-neutrino nuclear plus prepared energy. It includes the separately carried 6.4 ppm hot-bypass sensitivity and assigns the very small untransported nuclear lines their entire energy as loss. All percentages are fractions of that whole heat account.

| Hot carbon column, kg/m² | Broad scenario loss, % | 80% sampling half-width, percentage points | Loss if all saved prompt heat goes to cold equipment, % |
|---:|---:|---:|---:|
| 250 | 8.8179 | 0.0651 | 9.2119 |
| 500 | 4.1200 | 0.0499 | 4.8521 |
| 750 | 2.1019 | 0.0364 | 3.0628 |
| 1000 | 1.1282 | 0.0265 | 2.2723 |
| 2000 | 0.1450 | 0.0085 | 1.6334 |

`weighted-results.json` preserves both nominal128 and source128 reweightings. The source model treats the dominant 2.312593 MeV decay line as isotropic photons at 2.5438523 MeV, and all beta kinetic plus annihilation energy as isotropic 5.657344 MeV photons. The factor 1.10 increases these source energy accounts by drawing from already counted outgoing bulk kinetic energy. It does not create additional nuclear energy. All saved grey prompt energy is a separate, normally incident 20 MeV photon scenario. The nominal scenario retains the 511 keV and 2.312593 MeV lines and replaces beta kinetic energy alone with 5.14304 MeV photons.

These are specified source distributions. An upper source energy by itself does not prove an upper leakage probability for a nonlocal cascade with energy-dependent interactions. The broad result is a hard-spectrum sensitivity, rather than a theorem about every possible spectrum. Source motion, decay before wall arrival, radioactive material routed behind the shielding, angular correlations and actual prompt spectrum remain source-geometry questions. The 6.4 ppm carry is inherited from the earlier pulse-flow screen, is not recomputed here, and can overlap already saved prompt energy.

## Positrons reaching the wall

A separate response grid launches primary positrons normally from the cavity into the receiver. Every history is normalized to **kinetic energy plus 1.02199782 MeV**; the annihilated material electron is accounted for consistently with the positron annihilation energy already in the nuclear heat ledger. The native code also counts an escaping positron's remaining kinetic plus annihilation energy.

| Positron kinetic energy, MeV | Leakage from its complete calorimetric energy, Σ500 | Leakage, Σ1000 |
|---:|---:|---:|
| 0.5 | 1.6967% | 0.02836% |
| 1.0 | 1.4088% | 0.02666% |
| 1.8 | 1.0848% | 0.03116% |
| 4.12 | 1.2692% | 0.16947% |

Each row has 10,000 histories. Slowing, in-flight annihilation, annihilation position and subsequent gamma transport are included. These responses support materially better retention than replacing all positron energy with hard photons. They do not determine the late plume's fields or the true beta spectrum. The supplied energy-normalized whole-cycle comparisons vary the effective positron count with the chosen monochromatic energy, so they are explicitly response brackets, not replacements for the actual parent-decay inventory. No additional positron or recombination energy is credited to the plant.

## Geometry and numerical checks

The shell has density 2000 kg/m³, temperature 1000 K, inner radius 10⁶ m and thickness Σ/2000 m. Isotropic photons start 25 µm below its inner face. Normally incident prompt photons and positrons start 1 µm inside the cavity. The external vacuum ends 10 m beyond the shell. Inner exits traverse the actual cavity and can strike the opposite wall; they are not counted as lost. The reduced radius preserves an extremely large radius-to-thickness ratio. The geometry is static and calculates energy partition, not the actual gigametre light-travel delay or evolving wall temperatures.

Production uses verified 10 µm secondary-production cuts, a 100 eV photon production floor, and a 100 eV electron tracking floor. The actual built graphite table gives electron/positron production thresholds of 31.535/31.229 keV. A 1 µm control gives 742.7/733.4 eV. This cut concerns explicit secondary production; continuous energy loss is retained.

There are 605,000 production histories and 150,000 control histories. All final histories finish with zero residual tracked kinetic energy. The largest production event energy error is **5.65×10⁻¹¹ of incident calorimetric energy**. The archived verifier passes **414/414 checks**, reconstructing statistics from compressed per-event energy arrays. A separate agent independently reproduced all 605,000 production histories and 110 component/partition identities: maximum statistics difference 5.56×10⁻¹⁷ and ledger identity error 2.35×10⁻¹⁷.

At Σ500 the following changes in source128 broad leakage were measured against the production baseline:

| Control | Difference, percentage points | 80% sampling interval, percentage points |
|---|---:|---:|
| 1 µm instead of 10 µm cut | −0.0557 | [−0.1200, +0.0087] |
| Radius 10⁴ m instead of 10⁶ m | −0.00035 | [−0.0649, +0.0642] |
| Inner-face source instead of 25 µm depth | −0.0202 | [−0.0846, +0.0442] |

Each control has 25,000 histories at each of the two dominant broad delayed energies. The results are statistically compatible with the baseline. The cut comparison does not establish a rigorous error bound below 0.1 percentage point. Reported production intervals are Monte Carlo sampling intervals; they do not include source, geometry, cross-section or apparatus uncertainty.

## Physics and implementation qualifications

The implementation uses [Geant4 option4](https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/electromagnetic/Opt4.html), the [11.4.1 constructor source](https://github.com/Geant4/geant4/blob/v11.4.1/source/physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option4.cc), and the isolated [geant4-pybind](https://github.com/HaarigerHarald/geant4_pybind) 0.1.3 wheel. Official EMLOW8.8, ENSDFSTATE3.0, PhotonEvaporation6.1.2 and PARTICLEXS4.2 datasets were verified against published MD5s; SHA256s and exact download URLs are archived. Native `G4PSEnergyDeposit` scorers tally hot and vacuum deposition. Python callbacks tally tracks, outer escape and annihilation products; no biasing, track roulette or artificial termination is used.

Three implementation findings are preserved rather than hidden:

1. Inherited Python physics-list methods could report requested cuts while the actual material table used a different setting. Production applies cuts through Geant4's native `/run/setCut` interface and verifies `/run/dumpCouples`. Earlier misleading-cut smoke runs are excluded.
2. The installed option4 default enables in-flight three-photon annihilation. One 20 MeV history lost 2.133881 MeV at that single annihilation step. A detailed native step audit and the archived [three-photon model source](https://github.com/Geant4/geant4/blob/v11.4.1/source/processes/electromagnetic/standard/src/G4eplusTo3GammaOKVIModel.cc) identify a nonzero summed CMS photon momentum before boosting. Final production explicitly disables that branch and uses standard two-photon in-flight annihilation; the at-rest model is unchanged. Separate default-branch tests at Σ500 found zero three-photon histories in 10,000 each for 5.657 MeV photons, 20 MeV photons and 4.12 MeV positrons. Each zero count has an exact central-80% binomial upper fraction of 0.00023023. This occurrence screen does not bound all soft radiative corrections. The original failed history remains archived and is excluded from production.
3. An earlier 5.657 MeV event showed a 217 eV excess from a 1.691 keV carbon recoil crossing the inner boundary. The step audit is preserved. Every final production and control event satisfies the much tighter closure quoted above; no deposited or escaped energy was renormalized. The binding also faults during interpreter destruction after normal completion. Each isolated simulation process verifies and flushes its completed output before an explicit process exit, avoiding teardown without truncating any history.

Photonuclear reactions are omitted. The actual major decay gamma lines lie below the relevant carbon breakup thresholds. The broad 5.657 MeV sensitivity crosses natural C13(γ,n), while the 20 MeV prompt sensitivity crosses additional C12 channels. The independent [photonuclear screen](../../mixture-design/photonuclear-screen/README.md) finds these are small relative to electromagnetic interactions at the supplied energies, and distinguishes isotope conversion from permanent carbon loss. A first-interaction estimate is not a cascade calculation. The separate **all saved prompt heat to cold equipment** policy is therefore supplied as a conservative accounting sensitivity, charged once; it is not imposed as a physical requirement. Geant4's recoil ions from electromagnetic pair conversion are included and are distinct from nuclear breakup.

The shell is natural graphite. Fields, guide ports, cold equipment geometry, photonuclear secondary neutrons, changing charge state in the far-field plume and material damage are outside this radiation transport model. Escaped energy can be allocated to external escape or cold equipment by the callable ledger. Other non-photon energy retention remains the material receiver's assumption.

## Reuse

`basis.json`, the archived spectra and `reweight.py` are sufficient for downstream calculations; external raw traces and Geant4 are unnecessary for reweighting. For example:

```python
result = evaluate(basis, spectrum, 500, scenario="broad_D1p10",
                  prompt_policy="all_to_cold", cold_intercept_fraction=1.0,
                  hot_bypass_ppm=6.4)
```

Run `python verify.py` inside this packet to reconstruct all reported sufficient statistics and ledger checks. `transport_v3.py` and `final_campaign.py` reproduce transport using an isolated installation of the named Geant4 wheel and datasets. `reduce_basis.py` extracts the compact archive from those raw runs. Dataset payloads and excluded full failed-run arrays remain at the hashed raw paths; the completed event-energy arrays needed for independent verification are included here.
