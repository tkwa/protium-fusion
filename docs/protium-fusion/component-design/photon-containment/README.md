# ✨ Photon containment and receiver thickness

For the selected 64-zone A-H3/B1-H2/B2-H2.5/C13-H4 cycle, the dominant delayed photons are 511 keV annihilation radiation and the 2.313 MeV oxygen-14 decay line. The old screen that treated receiver photons as 20 MeV radiation is too restrictive to establish a minimum absorber thickness. A **2,000 kg/m² carbon receiver is a physically credible candidate** for a subsequent coupled transport calculation. A 1,000 kg/m² receiver is also worth testing; this packet does not certify its capture fraction. The original 5,000 kg/m² mixture-comparison denominator and assumed plant capture fraction remain unchanged.

These results concern radiation deposited in the designated receiver. Heat reaching cold structures, openings, radioactive processing lines, or an inadequately shielded inventory needs a separate account.

## Actual delayed spectrum

The source is the original ideal-EOS, completed 64-zone quartet used by plant-scan-003, weighted by its perfect eventual-return cycle. It releases 24.6976734 MeV of non-neutrino nuclear energy per helium nucleus and receives 1.9604524 MeV of prepared pulse energy. Their sum, **26.6581258 MeV/He**, is the heat-boundary denominator below.

The actual output inventories contain 0.0374339 nitrogen-13, 0.962685 oxygen-14, and 1.000120 oxygen-15 decays per helium nucleus. Primary delayed photons carry **4.2555195 MeV/He**, or **15.9633%** of that heat boundary:

| Photon energy | Photons/He | Energy, MeV/He |
|---|---:|---:|
| 511 keV | 3.995983 | 2.0419474 |
| 1.63520 MeV | 0.000500596 | 0.000818575 |
| 2.312593 MeV | 0.956793 | 2.2126733 |
| 3.94750 MeV | 0.0000203127 | 0.0000801842 |

Nitrogen-13 and oxygen-15 decay to their daughters' ground states. Their annihilation intensities and mean positron energies come from the [LNHB nitrogen-13 evaluation](https://www.lnhb.fr/nuclides/N-13_tables.pdf) and [LNHB oxygen-15 evaluation](https://www.lnhb.fr/nuclides/O-15_tables.pdf), both updated in March 2026. Oxygen-14's nuclear gamma intensities come from the [ENSDF nitrogen-14 decay evaluation](https://www.nndc.bnl.gov/ensnds/14/N/ec_decay.pdf), dated July 1990: 99.388%, 0.052%, and 0.00211% per parent for the three lines, respectively. Cascade intensity differs from direct beta feeding.

Mean positron kinetic energy adds **1.4984683 MeV/He**. It generally heats matter, with some converted into bremsstrahlung or in-flight annihilation photons. If annihilation occurs before slowing, both the positron kinetic energy **and the 1.022 MeV annihilation rest energy** enter the harder component. The broad bound below therefore moves all annihilation energy plus all positron kinetic energy into photons with energy no greater than 5.14304 MeV; it does not move only the kinetic part. The whole atomic Q of all other unstable endpoint species sums to only 5.3×10⁻¹⁰ MeV/He.

Neutral-atom branching is used. Suppressing electron capture in stripped ions changes the dominant annihilation counts by at most about 0.2%. The existing cycle uses approximate neutrino means; its delayed-heat total differs from this updated neutral-decay photon-plus-positron sum by 0.00134 MeV/He. The old energy ledger has not been silently revised.

## What attenuation does and does not establish

The [NIST graphite table](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html) gives total attenuation κ and mass energy absorption κ_en. The table's cm²/g becomes m²/kg by multiplying by 0.1. The uncollided fraction at normal incidence is exp(−κΣ). Scattering removes photons from the original beam without necessarily depositing their energy. In particular, **1−exp(−κ_en(E₀)Σ) is not used as a thick-shield transport prediction**.

[NIST's definition of κ_en](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html) includes radiative energy produced while secondary charged particles slow. Applying it as a local collision-energy kernel assumes that this slowing and photon regeneration occur near the interaction point. Explicit electron/positron transport, surface escape, and spatially displaced bremsstrahlung are outside that approximation.

A useful bound within that local-cascade model is obtained with

`κ* = min κ_en(E), for all degraded photon energies E ≤ E₀`.

For a photon population with energy weight `W = Σ E exp(κ*x)`, propagation contributes at most κ*W per unit mass path, while interaction removes κ_en(E)W on average. Thus the expected weight cannot grow. An inner-cavity crossing returns to another receiver patch at x=0 without increasing W. Outer escape at x=Σ then gives

`expected escaped energy / incident energy ≤ exp(−κ*Σ)`.

This includes repeated scattering and photon branching **within the specified local energy kernel**. It requires no photon up-scattering, a closed cavity, and source positions at the inner face or inside the cavity. Carbon's tabulated minimum below about 3 MeV occurs near 80 keV; the script takes a rounded-down minimum of the whole interpolated coefficient curve, rather than using the incident value. It assumes radiation below the table's 1 keV endpoint ultimately thermalizes in opaque carbon. Receiver temperatures around 1,000–2,000 K cannot regenerate keV photons thermally.

The broadest row below adds the actual saved grey-prompt energy, 0.3457872 MeV/He, placing all of it at 20 MeV. This is a deliberately hard sensitivity, not a predicted escaping spectrum.

| Carbon column, kg/m² | Primary delayed-photon loss upper bound, % of heat boundary | All annihilation + beta as hard photons, % | Previous column + all grey prompt at 20 MeV, % |
|---:|---:|---:|---:|
| 1,000 | 2.0966 | 3.5411 | 3.9357 |
| 2,000 | 0.27535 | 0.59541 | 0.71546 |
| 5,000 | 0.000624 | 0.003165 | 0.006545 |

The preceding table uses rest-frame decay energies. The fastest saved endpoint bulk flow is C13-H4 at 0.076563c, giving a maximum forward Doppler factor of 1.0797323. The script records each stage’s speed and uses both that factor and a rounded **D=1.10** sensitivity. It multiplies delayed photon energies and total source weights by D, conservatively drawing the additional laboratory photon energy from the existing kinetic-energy budget. The saved grey-prompt sensitivity remains specified at 20 MeV in the laboratory frame.

At D=1.10, the broad photon-loss bounds become **4.45654%, 0.838416%, and 0.007955%** for 1,000/2,000/5,000 kg/m². This factor also exceeds the largest endpoint bulk speed plus a thermal velocity margin of 21.3 standard deviations, even assigning every radioactive parent the maximum endpoint temperature and mass number 13. The corresponding Maxwell speed-tail probability is below 5.7×10⁻⁹⁸. This is an endpoint Maxwellian diagnostic; subsequent pressure acceleration or externally energized particles need to remain inside the selected velocity envelope.

For 2,000 kg/m², the D=1.10 broad bound leaves **0.161584% of H+U**, or **0.0430754 MeV/He**, for all remaining losses under a 99% capture target. Reserving the entire Doppler-adjusted direct hot-gamma bypass as an additional unrecovered loss reduces this to about **0.16095%**. This extra reserve can double-count part of the grey radiation, which is conservative. Neutrinos have already been excluded from H+U. Charged-particle escape, ports, hot-to-cold leakage, and radioactive material routed outside the credited shield must fit the remaining allowance or prompt a changed design.

The 2,000 kg/m² result leaves a limited remaining loss allowance if the complete plant is to capture 99%. It is not a proof that the entire plant meets that target. The script also removes 20 and 50 kg/m² from the credited column as boundary-layer sensitivities. Those are chosen length scales, not measured electron ranges or a substitute for coupled electron-photon transport: repeated displaced secondary emission is not rigorously bounded by one simple subtraction.

Decay nuclei implanted within tens of micrometres of the inner surface fit the front-source approximation relative to a 0.5–2.5 m carbon receiver. Nuclei moved behind that receiver before decaying need an enclosing local shield. Openings, cooling channels, and thin structural paths are similarly separate geometry requirements.

## Explicit Compton transport sensitivity

`containment.py` also transports 50,000 photons at each principal rest-frame line energy and column through a cold, free-electron Compton model. It uses the [Klein–Nishina sampling equations documented by Geant4](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/gamma_incident/compton/compton.html). A photon leaving the inner boundary crosses the closed spherical cavity and reaches another patch with the same absolute incidence cosine. A reflected coordinate in the thin-slab calculation represents that chord geometry, rather than a reflective graphite surface.

Recoil energy is deposited locally. Photoelectric absorption, coherent scattering, pair creation, charged-particle transport, and secondary bremsstrahlung are omitted. All residual photons below 1 keV are charged to leakage, although real graphite would strongly absorb them. The result is a controlled transport sensitivity, not a full material simulation or rigorous upper bound on real leakage.

| Column, kg/m² | 511 keV energy reaching outer surface | 2.313 MeV energy reaching outer surface | Both lines, including untracked sub-keV reserve, % of total heat |
|---:|---:|---:|---:|
| 1,000 | 0.5678% | 4.7679% | 0.4468% |
| 2,000 | 0.02298% | 0.1502% | 0.03118% |
| 5,000 | zero sampled | 0.000000904% | 0.01855% |

The final column at 5,000 kg/m² is almost entirely the deliberately discarded sub-keV residual; it is not a physical leakage floor. At 1,000 kg/m² the weighted result's sampling-only 80% interval is 0.4381–0.4556% of the heat boundary. At 2,000 it is 0.02971–0.03265%. These intervals describe Monte Carlo sampling, not uncertainty in the material physics or plant design. Exact energy balance is checked for every history, and sampled uncollided escape is checked against the analytic exponential.

This is why a thinner receiver deserves a coupled transport calculation. The difference between the analytic bound and the model result is large enough to matter; neither justifies keeping or discarding 1,000 kg/m² by itself. A useful acceptance test is total radiation lost from the hot recovery system below its assigned share of the plant's 1% heat-loss budget, with all radioactive plumbing and openings represented.

## Prompt nuclear gamma bypass

`hot_bypass.py` reconstructs every positive nonweak reaction flow from the actual 400-state saved traces. Capture and reverse cycles are retained. Their gross photon source can exceed the net nuclear release. Charged-particle channels are additionally assigned entirely to photons for a deliberately loose second source envelope.

The summed capture-cascade energy is Q plus the reacting-pair kinetic energy. For a two-body Maxwellian rate this is `Q + kT[3/2 + d ln(rate)/d ln(T)]`, using analytical derivatives of the saved Reaclib coefficients; the generalized n-body normalization is also implemented. This avoids treating barrier-selected reacting pairs as an unweighted thermal average. Photon cascade spectra remain unresolved.

Every cell's source is placed at its outer edge and directed outward. First-interaction opacities 10⁻⁴, 3×10⁻⁴, 10⁻³ and 3×10⁻³ m²/kg are explicit scenarios. Even at the smallest, the cycle-weighted direct, uncollided source is at most **0.00015313 MeV/He**, **5.74×10⁻⁶ of the heat boundary** in the local rest-energy account. Multiplying all source energies by the D=1.10 envelope gives 0.00016844 MeV/He, or 6.32×10⁻⁶ of the heat boundary, provided the stated interaction-opacity scenario applies across the boosted photon energies. Almost all of that deliberately loose result comes from assigning the outermost cell no overlying material. Per-stage fractions of initial plus net nuclear energy are at most 1.74×10⁻⁵. A uniform-in-cell source sensitivity is much smaller.

Minimum total radial columns during the saved trajectories are 2.28×10⁷, 2.45×10⁵, 2.79×10⁷, and 1.06×10⁷ kg/m² for A, B1, B2, and C13. Maximum saved temperatures are approximately 0.102, 0.0415, 0.142, and 0.130 MeV. This strongly disfavors a substantial fresh uncollided nuclear-gamma component at the receiver within the resolved spherical model. A first collision is not complete thermalization: hot electron scattering, pair cascades, frequency redistribution, the escaping photospheric spectrum, and unresolved subcell transport still need a consistent treatment. The model has no holes or asymmetries.

Time quadrature reproduces integrated net nuclear Q to 0.018%, 0.128%, −0.214%, and 1.576% respectively. Every-other-saved-state bypass integrals differ by at most 0.541%. These are sampling diagnostics, not mesh convergence. A small hard prompt tail is therefore carried as a bounded scenario rather than silently eliminated.

## Reproduction and scope

Files: [spectrum.py](spectrum.py), [spectrum.json](spectrum.json), [hot_bypass.py](hot_bypass.py), [hot-bypass.json](hot-bypass.json), [containment.py](containment.py), [containment.json](containment.json), and [integrity.json](integrity.json). The stage paths and source hashes are saved in the results. `spectrum.py` reads the existing generalized material-cycle ledger; `hot_bypass.py` additionally needs the retained original NPZ traces under `/home/tkwa/src/protium-fusion/output/mixture-runs`. Those large raw inputs are not duplicated here. There is no new hydrodynamic run.

Using the existing scientific Python environment, run these in order:

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/tmp/protium-python
export MPLCONFIGDIR=/tmp/protium-mpl
export OPENBLAS_NUM_THREADS=1
export NUMBA_CACHE_DIR=/home/tkwa/src/protium-fusion/output/nonideal-runs/.numba-cache
PY=/tmp/protium-python-interpreters/cpython-3.10.20-linux-x86_64-gnu/bin/python3.10
$PY spectrum.py
$PY hot_bypass.py
$PY containment.py
$PY verify.py
```

The large practical uncertainty is the coupled radiation-and-source geometry, not the small statistical uncertainty in this Compton calculation. The deterministic scenarios here should not be read as an 80% engineering credible interval. A credible interval for whole-plant capture requires a specified receiver, plumbing, port, and secondary-particle design.
