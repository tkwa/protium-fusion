# ✨ Nuclear-rate uncertainty: completed bounded checkpoint

The updated carbon data change the preferred B2 mixture materially. The remaining tested carbon uncertainty curves produce a much smaller change in the complete cycle than in individual pulse gains. Nitrogen-15 exposes a real sensitivity to the prepared source profile: a low-rate collapse in one optimized profile disappears in the explicitly constructed source profile at the same preparation budget. Mixture, source, nuclear-rate and mesh comparisons must therefore stay paired.

All station numbers below use the unchanged receiver, isotope-routing, cadence and plant accounting, with 128-zone hydro endpoints. They remain conditional model outputs. Their many digits make comparisons reproducible; they do not represent physical accuracy. The root-owned 256-zone and apparatus studies are separate qualifications.

## Carbon data and four tested corners

The exact 194-rate manifest is `../network-provenance.json`. The nominal revision replaces C12(p,gamma)N13 and C13(p,gamma)N14 with the LUNA public evaluated tables through 10 GK. Each reverse capture is multiplied by the same factor as its forward partner. Identity and detailed-balance checks pass; the independently selectable kernel gives bit-identical co-shifted rates to the earlier kernel. Exact runtime versions, inputs and hashes remain archived.

For B2 H2, the tested carbon curves give pulse gain 9.581576–10.228667 around median 9.917629. B2 responds almost entirely to C12. For C13 H4, the corresponding range is 24.607797–24.817039 around 24.698571; its response is dominated by C13 and has the opposite sign. Changing each carbon rate independently confirms that this separation is real in the simulated network.

The actual isotope endpoints were routed through the complete cycle, with A H3 and B1 H2 held at their median-carbon profiles:

| C12 curve | C13 curve | Net mW/kg, 5000 kg/m² absorber | Net mW/kg, 2000 kg/m² absorber |
|---|---|---:|---:|
| median | median | 6.640608 | 12.510481 |
| minus1sigma | minus1sigma | 6.623351 | 12.477984 |
| plus1sigma | plus1sigma | 6.655616 | 12.538740 |
| minus1sigma | plus1sigma | 6.622287 | 12.475981 |
| plus1sigma | minus1sigma | 6.656486 | 12.540379 |

Thus the tested 5000 kg/m² range is −0.276% to +0.239% around the median. All five points remain radiator-limited. Altered incomplete burning changes the frequency of each recycled stage, which moderates the station response. Scaling plant output by one pulse-gain ratio would miss this behavior. `carbon-integrated.json` and `carbon-reduced-inputs.json` preserve the full calculation; 68 independent energy, mass, reconstruction and occupancy checks pass.

These are four conditional corners from the source's named one-sigma curves. They establish neither a joint probability interval nor a rigorous bound on interior parameter values. The calculation does not assign cross-reaction covariance or sample the full temperature-dependent rate posterior.

## Nitrogen-15 and the source profile

NACRE-II supplies adopted, low and high N15(p,alpha)C12 evaluations. Its low/high limits have no probability assigned here. At the original optimized B1 profiles, H2 has a very small central plant advantage over H1.5, but the low curve reverses that ordering:

| N15 curve | H1.5, mW/kg at 5000 kg/m² | H2, mW/kg at 5000 kg/m² |
|---|---:|---:|
| saved | 6.633048 | 6.640608 |
| adopted | 6.631935 | 6.639219 |
| low | 6.615183 | 5.509809 |
| high | 6.620179 | 6.625552 |

The low-H2 point reaches the allowed 0.8 main-window occupancy before filling the available radiator. The source profile built by the launcher study changes this result. Its existing contact-entropy allocation, with the same prepared-energy budget, gives:

| Constructed source case | Saved N15, mW/kg | Low N15, mW/kg |
|---|---:|---:|
| H1.5 | 6.620868 | 6.633942 |
| H2 | 6.626477 | 6.640488 |

These comparisons hold A H3/B2 H2/C13 H4 at their nominal 128-zone endpoints. All four constructed-source points are radiator-limited, with main-window occupancy 0.720–0.747. The H2 source low-rate result is within 0.002% of the original nominal plant number. `b1source-integrated.json` and its replay archive preserve 64 passing checks; `n15-integrated.json` preserves the original fixed-profile family and 97 checks.

The concrete source adjustment resolves the particular low-rate loss observed in the original profile. The subsequent coherent-source calculation below also completes the high-curve test. The response remains nonmonotone; these finite source-controlled cases carry no assigned probability coverage. In particular, the original collapse must not be presented as an unavoidable future reactor limit or as a probability of failed deployment.

## Completed coherent-source envelope

The subsequent [source128-envelope.json](../../component-design/b1-source-robustness/source128-envelope.json) holds the other three stages at their **actual constructed source128** endpoints. This differs from the nominal128 A/B2/C13 baseline used in the two preceding tables; those earlier records remain unchanged.

| B1 source recipe | Original N15 control | NACRE-II low | NACRE-II high |
|---|---:|---:|---:|
| H2 | 6.6173445 | 6.6313366 | 6.6024311 |
| H1.5 | 6.6117466 | 6.6248025 | 6.5979080 |

Values are net mW/kg with the 5000 kg/m² absorber. H2 narrowly leads at every tested curve. All cases remain radiator-limited; maximum main-window occupancy is 0.76848 at 5000 kg/m² and 0.76901 when the 2000 kg/m² cases are included. The archive contains 82 passing accounting/profile checks and 28 passing trace checks. An independent read-only check confirmed the reported power-to-mass ratios, radiator limit, matched ranking, check counts and source/code fingerprints.

The control retains the original REACLIB N15 rate. The low/high cases use the NACRE-II evaluation curves with median carbon rates. There is no adopted-NACRE-II source point in this family and no probability assigned to the tested envelope. These simulations close the previously outstanding source-high comparison while preserving the distinction between sampled rate sensitivity and overall design uncertainty.

## Resolution and remaining model qualifications

Updated B2 H2 gains at 64, 128 and 256 zones are 10.271886, 9.917629 and 9.867099. The last refinement changes gain by −0.5095%, much less than the preceding −3.45%. This is encouraging resolution behavior, without an established asymptotic error bound. At 128 zones, updating the carbon rates changes the same profile from 11.101894 to 9.917629; rate choice and resolution are separate effects. A small conservation residual does not measure either model accuracy or convergence.

The independent N13 source-fit substitution has a small effect on the tested 64-zone B2 profiles despite a large disagreement between the saved hot REACLIB expression and the cited paper. It remains a separately labeled sensitivity, supported by the independent STARLIB comparison. The provenance discrepancy is preserved rather than called a proven library error.

The new Chen2025 N14 measurement concerns low-energy capture. Its solar S(0) change cannot be transplanted as a uniform GK multiplier; `N14-CURRENCY.md` records the primary evidence and bounded energy-window calculation. Screening, EOS, radiation transport, source asymmetry and collection apparatus uncertainties remain distinct from this nuclear-rate packet. Their combined effect has not been converted into an 80% station-performance interval.

The useful outcome is a reproducible conditional comparison and a better operating family. A statistically defensible overall credible interval would require explicit probability models and correlations for the remaining physics and hardware, together with source retuning rules. The present measured-rate curves should retain their source conventions rather than acquire invented coverage.
