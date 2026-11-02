# ✨ N14(p,gamma) currency check

The December 2025 Chen et al. measurement should be cited, but it does not supply a justified 15% multiplier for our GK pulse calculation. The published abstract gives S(0)=1.93 ± 0.10 keV b. The accessible public preprint instead gives 1.92 ± 0.08; those versions are recorded separately. The new proton energies are 110–260 keV laboratory, or approximately 103–243 keV center of mass. Its R-matrix fit uses higher-energy data too, while explicitly targeting behavior below the 278 keV laboratory resonance. The preprint adopts 13.0 ± 0.4 meV for that resonance strength, compared with Imbriani 2005's 12.9 ± 0.9. We found R-matrix parameters but no numerical hot-rate table in the accessible material. The published full text and supplement require subscription and were not retrieved. [Published PRL](https://doi.org/10.1103/q756-hzmt), [public author preprint and supplement](https://arxiv.org/html/2410.16086v1).

A transparent relevance screen is saved as `n14_chen_relevance.py` and `n14-chen-relevance.json`. It adds a *specified illustrative* constant deltaS=0.25 keV b to a Maxwellian average,

`delta<sigma v> = sqrt(8/(pi mu)) (kT)^(-3/2) deltaS integral exp[-E/kT-sqrt(EG/E)] dE`,

with EG=2mu c²(pi alpha Z1 Z2)² and mu=(14/15)u. Multiplication by NA and the SI-to-cm³ conversion gives the same units as the saved REACLIB rate. One integration ends at 0.242667 MeV; another extends the *same assumed constant* to 10 MeV. Neither reproduces the authors' R-matrix fit, assigns a confidence interval, or bounds all resonance/interference changes.

| T(GK) | Saved rate(cm³mol⁻¹s⁻¹) | Added constantS below0.243MeV,relative | Same added constantS across all energies,relative |
|---:|---:|---:|---:|
|0.684|48.5019|0.03569%|0.33786%|
|1.000|113.1856|0.02648%|0.87627%|
|1.058|123.6648|0.02547%|1.02547%|
|1.147|138.7395|0.02407%|1.28801%|

The A-H3 saved reaction-flow 5th/50th/95th temperature percentiles are approximately 0.684/1.058/1.147 GK. At 1 GK the nonresonant Gamow peak is 0.4364 MeV, beyond the newly measured range. The archived rate's resonance-labeled fit terms provide 98.59% of its value at 1 GK. That last number describes the fit representation; individual REACLIB sets are not independent measured reaction channels.

The low-energy enhancement is consequently a poor proxy for the pulse-temperature rate. The calculation does **not** establish that every change to the full hot R-matrix rate is negligible. A reproduced published cross-section evaluation, including resonance contributions and covariance, would be needed for a source-supported replacement. On the evidence retrieved here, keep the present N14 rate and its documented currency limitation; do not add an unsupported uniform multiplier or infer a 15% pulse-yield error.

Provenance: `rate-sources/chen2025-public-preprint.html` is a byte-for-byte public download, with URL/SHA256 in `chen2025-provenance.json`. The rate manifest hash and calculation hash are in the output JSON. No hydro endpoint was changed by this check.
