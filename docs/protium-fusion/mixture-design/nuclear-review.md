# ✨ Nuclear-data review and controlled rate revision

This packet audits the exact 27-species, 194-rate network used in the saved pulse calculations. It supplies a controlled replacement of the two carbon proton-capture rates, plus a separate nitrogen-13 source-fit sensitivity. It does not establish the full reactor design or turn model uncertainty into a measured performance interval.

## Current interpretation

The completed source-supported uncertainty studies and coupled plant comparisons are summarized in [rate-uncertainty/INTERPRETATION.md](rate-uncertainty/INTERPRETATION.md). All four independent carbon corners, six N15 evaluation cases, and four constructed B1 source controls are complete. The constructed source profiles remove the original low-N15 H2 collapse; the earlier fixed-profile table below remains a documented sensitivity rather than a universal reactor limitation. The latest N14 measurement is covered in [N14-CURRENCY.md](rate-uncertainty/N14-CURRENCY.md).

## Original 64-zone comparison

The consequential measured-data update is carbon capture. In exact-profile reruns, the LUNA public evaluated tables reduce B2 H2 gain from 10.5550 to 10.2719 and B2 H2.5 from 10.7999 to 9.3500. They raise the hydrogen-rich C13 gains slightly while reducing the residual N14 endpoint inventory. Whole-cycle routing must therefore use the saved isotope endpoints. Multiplying an old plant result by a pulse-gain ratio would miss the changed routing.

The saved N13(p,gamma)O14 REACLIB expression differs sharply above about 1 GK from Eq.17 of the paper cited by its `lg06` label. The independent official STARLIB table supports the lower high-temperature rate. Nevertheless, substituting the published Li expression changes the carbon-updated B2 H2 and H2.5 gains by only +0.0056% and +0.166%, respectively. This discrepancy deserves provenance follow-up, but it does not undo the carbon-update effect in these two profiles. The Li expression remains a separately named sensitivity; the reason for the JINA version-5 coefficients has not been established.

Classical screening changes the flow-weighted capture-rate diagnostic by roughly 1–3% on these trajectories. An EOS-consistent screening calculation is useful near a mixture optimum, but there is no basis here for generic factors of ten. A bounded omitted-channel screen finds very small direct entry exposure into the missing light species. The current priorities are the carbon update, mesh/mixture convergence, and consistent plasma thermodynamics, followed by selected source-supported rate-band reruns.

## Exact provenance

[network-provenance.json](network-provenance.json) contains every rate in saved order: reactants/products, stoichiometry, density exponent, statistical prefactor, each seven-coefficient REACLIB set, raw source record, Q values, weak-energy convention, source metadata and sampled values. The runtime is Python 3.10.20, pynucastro 3.1.0, NumPy 2.2.6 and SciPy 1.15.3. The library is `reaclib_default2_20250330`, SHA256 `24e1f37c502cf8521109b7dd91cb73d2d4acf38f94ffd1fe1fb1e4beec3ba862`. The pulse-model source SHA256 is `69596156a457eecda91c19188684d70a848197bfb7f1f07a1cada5545ad581c5`.

The reconstructed seven-term expression agrees with the library evaluator to a maximum relative difference of 1.46e-11. Baryon stoichiometry is exact. The packaged `ls09` bibliography has malformed metadata; the packet preserves it rather than treating the malformed author string as a valid citation. A library's release date does not imply that all its component measurements are recent.

| Main forward reaction | Saved label | Relevant source |
|---|---|---|
| C12(p,gamma)N13 | ls09 | Li et al., DOI 10.1007/s11433-010-0128-8 |
| C13(p,gamma)N14 | nacr | Angulo et al. 1999, DOI 10.1016/S0375-9474(99)00030-5 |
| N13(p,gamma)O14 | lg06 | Li et al. 2006, DOI 10.1103/PhysRevC.74.035801 |
| N14(p,gamma)O15 | im05 | Imbriani et al. 2005, DOI 10.1140/epja/i2005-10138-7 |
| N15(p,alpha)C12 | nacr | Angulo et al. 1999 |

## Carbon tables and implementation

The nominal revision uses the collaboration's public evaluated [C12 table](https://luna.lngs.infn.it/index.php/scientific-output/reaction-rate-repository%3Flayout%3Dedit%26id%3D178.html) and [C13 table](https://luna.lngs.infn.it/index.php/scientific-output/reaction-rate-repository%3Flayout%3Dedit%26id%3D180.html), associated with [Skowronski et al., PRL131,162701 (2023)](https://doi.org/10.1103/PhysRevLett.131.162701). The public repository supplies 54 temperature nodes from 0.007 to 10 GK and median plus/minus one-, two- and three-sigma columns. Consequently no independently invented continuation above 1 GK is necessary. The shorter table in the public manuscript ends at 1 GK; its rounded medians differ by less than 0.92% from the repository values. This packet uses the repository values exactly.

[updated-carbon-rates.json](updated-carbon-rates.json) holds both numerical tables and source hashes. `rate-sources/luna-web-tables.json` and `luna-web-tail.json` preserve the retrievals; `luna2023-source.tar` and `supplement.tex` preserve the authors' public source. The table SHA256 is `9166934c7f82868daf7ff95054c78f6078994ab1c0a8e53040a37c1d651f4b99`.

The separate [Felsenkeller C12 experiment](https://doi.org/10.1103/PhysRevC.107.L062801) and [Kettner et al. high-energy measurement](https://doi.org/10.1103/PhysRevC.108.035805) support revisiting the old normalization. The [2025 C13 evaluation](https://doi.org/10.1103/h1yt-5myk) reports agreement with the LUNA evaluation through 10 GK and approximately 7% rate uncertainty. Its exact numerical supplement was not retrieved; these runs should not be described as using that newer table. A [2025 C12 uncertainty-method study](https://arxiv.org/abs/2502.14040) also warns that multimodal fit structure matters; pointwise covariance errors need not represent the full model uncertainty.

`updated_network.install(radial, table_path=None, quantile='median', mode='updated', interpolation='pchip')` replaces only the calling process's nuclear RHS. All numerical arrays are runtime arguments to an uncached Numba kernel. No cached global constants are mutated. It multiplies both directions of each C12/N13 and C13/N14 capture pair by the same temperature-dependent factor, preserving the original equilibrium ratio, stoichiometry, Q and weak-energy treatment. Log(new/old) is interpolated versus log(T9); outside the table a constant endpoint ratio is used. Relevant burning remains inside the tabulated range.

The identity kernel is bit-identical to the archived kernel on the check inputs. Updated rate nodes agree to 1.42e-14 and paired equilibrium ratios to 2.22e-16. PCHIP versus linear interpolation of log factors differs by at most 0.4033% over 0.3–2 GK. Input guards check the common carbon temperature grid, physical-duration aliases and compression-cap viscosity. `rate-sources/runtime-v1/` preserves the original runner/kernel associated with the first nine result hashes; later guards change no numerical physics.

Source sigma bands remain source sigma bands. Coherently selecting the same quantile at every temperature and for both reactions defines a sensitivity direction. It is not a sample of the full correlated rate posterior, a joint confidence interval, or an 80% plant-performance interval.

## Full-hydro results

All initial energy, composition, radial profile, opacity, requested physical duration, EOS and tolerances are retained from each baseline. The first nine nominal runs finish at their requested durations. Full endpoint mass fractions are in [carbon-nominal64.json](carbon-nominal64.json), and individual metadata/profiles in `rate-results/carbon64/`. Raw NPZ traces remain under `/home/tkwa/src/protium-fusion/output/mixture-rate-runs/nominal64/` with hashes in the summary.

| Case | Old gain | Updated gain | Old target yield | Updated target yield |
|---|---:|---:|---:|---:|
| A H1 | 8.196339 | 8.196339 | .539367 | .539367 |
| A H3 | 10.288815 | 10.288815 | .767340 | .767340 |
| B1 H1 | 6.381970 | 6.381946 | .601869 | .601887 |
| B1 H2 | 7.545623 | 7.545426 | .755968 | .756016 |
| B2 H2 | 10.555031 | 10.271886 | .821884 | .799328 |
| B2 H2.5 | 10.799859 | 9.350042 | .879866 | .750786 |
| C13 H1 | 15.538979 | 15.373712 | .816621 | .792060 |
| C13 H3 | 23.576308 | 23.729640 | .310566 | .292113 |
| C13 H4 | 23.655667 | 23.862778 | .222027 | .202573 |

Gain is deposited nuclear energy divided by prepared input energy. Target yield is the original stage's named endpoint product per initial catalyst, not total cycle utility. The largest endpoint energy residual is 1.24e-5 of prepared energy, in C13 H4. The additional 64-zone bracket cases and ten 128-zone exact-profile cases are separate outputs in `extension64/` and `nominal128/`; their status/results must be read before treating the rate/mesh comparison as complete.

## N13 source-fit discrepancy

The exact coefficients of [Li et al. 2006 Eq.17](https://arxiv.org/pdf/0903.3094), visually checked on PDF page14, are recorded in [n13-source-fit.json](n13-source-fit.json). Its stated less-than-5% representation error over 0.01–10 GK is a fit residual, not experimental uncertainty. The JINA page identifies `lg06` version5, modified 2010-06-11; the explanation for its different hot continuation remains unresolved.

| T9 | Saved REACLIB | Li Eq.17 | Official STARLIB where checked |
|---|---:|---:|---:|
| 1.0 | 904.312 | 949.376 | — |
| 1.5 | 6289.676 | 3774.222 | 3513 |
| 2.0 | 32473.377 | 6429.797 | 6398 |
| 3.0 | 472374.514 | 9972.886 | — |

Rates are cm3 mol−1 s−1. The [official STARLIB repository](https://github.com/Starlib/Rate-Library) archive `starlib_v610_120222.dat.zip` supplies an independent `nacr` evaluation. Its exact N13 record, archive/member hashes and URL are saved in [starlib-n13-provenance.json](starlib-n13-provenance.json) and `rate-sources/starlib-n13-official.txt`. Its uncertainty factor is a STARLIB lognormal approximation to the older evaluation bounds, not a new measurement.

`updated_network_n13.py` adds the Li forward/reverse paired factor to the carbon revision. In completed tests:

| Case | Carbon-only gain | Carbon+Li gain | Carbon+Li target yield |
|---|---:|---:|---:|
| B2 H2 | 10.271886 | 10.272462 | .799009 |
| B2 H2.5 | 9.350042 | 9.365565 | .752522 |

Both complete with energy residual magnitudes below 2.6e-6 of prepared energy. [li2006-nominal64.json](li2006-nominal64.json) and `rate-results/li2006_64/` retain the full endpoints. The limited sensitivity of total gain does not imply the two rate formulas are interchangeable in every possible profile.

## Temperature coverage, screening and omitted species

[trace-network.json](trace-network.json) integrates gross reaction flows on the archived saved trajectories. Typical flow-weighted 5th/50th/95th temperature percentiles are A N14 capture .684/1.058/1.147 GK; B1 N15(p,alpha) .243/.351/.465 GK; B2 C12 capture .679/1.235/1.618 GK; B2 N13 capture .966/1.384/1.622 GK; C13 H4 C13 capture .677/1.030/1.350 GK. About 93% of the B2 N13 flow is above1 GK. Solar S(0) precision alone therefore does not describe these rate windows.

These are sampled quadrature diagnostics. Nuclear-energy quadrature residuals range from roughly 0.02% to 1.6% for these cases. Solver-integrated energies and isotope endpoints are the preferred quantitative ledger.

The screening diagnostic uses classical Debye response, kappa_D²=e² n_b(sum Z_i²Y_i+Y_e)/(epsilon0 kT), and h=Z1 Z2 e² kappa_D/(4 pi epsilon0 kT). An ion-sphere expression is retained as a second scale, not combined into an asserted interpolation. Flow-weighted classical rate multipliers are approximately 1.011–1.026. Carbon coupling parameters are usually around0.1–0.3 over the important flow. A consistent calculation should couple screening chemical potentials to the Coulomb EOS and appropriate electron response. Primary foundations are [Salpeter1954](https://doi.org/10.1071/PH540373) and [Chugunov and DeWitt2009](https://doi.org/10.1103/PhysRevC.80.014611).

[omitted-channels.json](omitted-channels.json) evaluates 132 library entry reactions from present species into omitted species on unmodified trajectories. Total direct entry exposures are about7.4e-6 per catalyst in A H3, below1e-16 in B1 H2, 1.2e-8 in B2 H2 and1.75e-5 in C13 H4, dominated where relevant by proton-unbound F15/F16 channels. This is an entry-flux screen, not an enlarged closed network. Tiny neutron-fed flows below inventory and solver-error scales are numerical diagnostics, not resolved physical yield predictions. C14 entry and deuteron production through the available library channels are smaller still.

The absent C13(p,d)C12 channel has a2.722 MeV center-of-mass threshold. An explicitly assumed constant1-barn cross section above threshold gives an exposure of8.3e-6 per catalyst for C13 H4. That is a deliberately transparent sensitivity scale, not measured cross-section data or a rigorous upper bound. The screen provides no reason to expect D/T/He3/C14 omission to reorder these particular mixtures, while leaving unrepresented nonthermal distributions and unlisted channels outside its claim.

## Remaining rate priorities

The saved N14 `im05` representation agrees with official STARLIB to about0.3% over the important window. [Frentz et al.2022](https://doi.org/10.1103/PhysRevC.106.065803) extends experimental coverage and documents remaining fit/data issues; recent Bellotti work is still developing the next evaluation. The familiar historical reduction of the solar N14 rate is already in `im05` and must not be applied a second time.

For N15(p,alpha), [NACRE II](https://doi.org/10.1016/j.nuclphysa.2013.09.007) supplies a newer adopted/low/high table and explicitly discusses the two relevant resonances. Across the B1 main burn window, its adopted rate is1.013–1.026 times the saved expression; its low/high limits are roughly0.86/1.18 times that expression. Its source bounds are evaluation limits, not an asserted Gaussian probability interval. The visually verified Table35 and comparison are stored in [nacreii-n15.json](nacreii-n15.json). Its low-energy S(0) uncertainty cannot simply be transplanted to the B1 resonance window. A targeted B1 rate-bound check can follow if the remaining flat mixture optimum requires it; a blanket tenfold change has no evidential basis.

The immediate planned calculations are the authorized composition brackets, exact128 profiles, and carbon plus/minus-one-sigma sensitivity curves on B2 H2 and C13 H4. Rate uncertainty will remain separate from opacity, EOS, unresolved multidimensional flow and apparatus uncertainty.

## Added evaluation sensitivity API

`updated_network_n15.install(..., n15_curve="adopted"|"low"|"high")` adds the paired NACRE-II N15(p,alpha) factor with the carbon tables held at their selected quantile. `run_updated_case_n15.py --record BASELINE.json --out OUTPUT --n15-curve low` preserves the exact baseline physical profile. Identity, all three table-node sets and paired equilibrium ratios pass [n15-checks.json](n15-checks.json). All six B1 H1.5/H2 adopted/low/high 128-zone runs completed. These cases are conditional evaluation sensitivities; their bounds carry no assigned 80% probability.

A final currency check examined [ETR25, published2026](https://doi.org/10.3847/1538-4365/ae2bdc), including its [primary table index](https://arxiv.org/html/2601.20059v1). It does not provide a replacement for the five main proton reactions listed above. It updates several alpha-breakout and oxygen/fluorine channels; their entry flows in the present short saved pulses are small. Its rigorous percentile conventions should not be retrospectively assigned to NACRE-II bounds.

## Completed nominal128 comparison

All ten requested128-zone profiles completed with exactly unchanged prepared energy and physical duration. [carbon-nominal128.json](carbon-nominal128.json) contains full isotope endpoints and trace hashes; `rate-results/carbon128/` contains the metadata and initial profiles. [rate-mesh-comparison.json](rate-mesh-comparison.json) joins the64/128 results without removing nonmonotone changes.

| Case | Saved-rate128 gain | Updated-rate128 gain | Updated target yield |
|---|---:|---:|---:|
| A H1 | 8.713302 | 8.713302 | 0.573386 |
| A H3 | 10.817796 | 10.817796 | 0.806791 |
| B1 H1 | 7.284268 | 7.284246 | 0.686979 |
| B1 H1.5 | 8.165561 | 8.165321 | 0.794052 |
| B1 H2 | 8.262244 | 8.261904 | 0.827739 |
| B2 H2 | 11.101894 | 9.917629 | 0.766377 |
| B2 H2.5 | 11.214630 | 8.071872 | 0.636641 |
| C13 H1 | 15.994264 | 15.863752 | 0.818893 |
| C13 H3 | 24.300363 | 24.546266 | 0.279123 |
| C13 H4 | 24.392279 | 24.698571 | 0.180265 |

The carbon-update penalty in B2 is larger at128 than64. Updated B2 H2 gain changes10.271886→9.917629, and H2.5 changes9.350042→8.071872. This is material rate–resolution interaction;128 alone does not establish asymptotic convergence. The root-owned additional B2 H2 carbon-updated 256-zone control is complete, as recorded below. The nominal128 maximum energy-residual magnitude is3.04e-6 of prepared energy; a small conservation residual does not imply spatial convergence.

The four co-shifted carbon cases, four opposite carbon corners and six N15 evaluation cases completed in `carbon-bands128/`, `carbon-corners128/` and `n15-bands128/`. Their complete conditional comparison is in [rate-band-comparison.json](rate-band-comparison.json) and the integrated results linked above. A redundant unconditional duration lookup was removed from the carbon/N15 wrappers after the nominal family; old wrappers are preserved in `rate-sources/runtime-v2/`. All nominal inputs contained `duration_s`, so their calculations are unchanged. The correction permits older source inputs carrying only `physical_duration_s`.

## First completed carbon-band checkpoint

At128 zones, B2 H2 gives gain9.581576730 /9.917629469 /10.228667338 for simultaneous carbon minus-one-sigma /median /plus-one-sigma curves. Target yields are0.73792391 /0.76637728 /0.79289842. These correspond to−3.39%/+3.14% gain sensitivity around the median, comparable to the−3.45% change of that median between64 and128 zones. The physical rate uncertainty cannot be reported independently of the visible resolution uncertainty as a complete yield interval. Both carbon curves are shifted coherently; cross-reaction and cross-temperature correlations are not sampled.

## Independent256-zone B2 control

Root completed the same median-carbon B2 H2 profile at256 zones: gain9.86709879164, O14 target0.76156353023, energy residual−1.2654e-6 and peak temperature1.57635 GK, at the full requested duration. The128→256 gain change is−0.5095%, substantially smaller than64→128. This reduces the initial convergence concern, while not proving an asymptotic error bound. It is kept as a separate mesh comparison in [rate-mesh-comparison.json](rate-mesh-comparison.json); the published-rate curves remain conditional sensitivities at128. Additional root-owned256 profiles are separate ongoing checks.

## Matched N15 evaluation and whole-cycle check

All six B1 H1.5/H2 adopted/low/high128 runs completed. Under the low evaluation curve, H2 gain falls to6.212850 versus7.945344 for H1.5. Under adopted/high, H2 retains only a small pulse-gain advantage. These profiles are not reoptimized at each rate.

The actual endpoints, source/receiver reconstruction and retained material cycle were integrated using unchanged plant routines, with nominal128 A H3/B2 H2/C13 H4 fixed. [rate-uncertainty/n15-integrated.json](rate-uncertainty/n15-integrated.json) and its reduced-input replay archive preserve all results;97 independent energy/mass/occupancy checks pass.

| N15 curve | H1.5 mW/kg,5000kg/m² | H2 mW/kg,5000kg/m² | H1.5 mW/kg,2000kg/m² | H2 mW/kg,2000kg/m² |
|---|---:|---:|---:|---:|
| saved | 6.633048 | 6.640608 | 12.496244 | 12.510481 |
| adopted | 6.631935 | 6.639219 | 12.494148 | 12.507864 |
| low | 6.615183 | 5.509809 | 12.462604 | 10.379646 |
| high | 6.620179 | 6.625552 | 12.472012 | 12.482128 |

At the adopted rate, H2 gains only0.110% in specific plant power. Under the low curve, it is16.71% below matched H1.5. The low-H2 plant reaches main-window occupancy0.8 while using radiator fraction0.792; low-H1.5 stays radiator-limited at0.9 with occupancy0.741. For these original profiles, retaining H1.5 has a small nominal accounting cost. The later constructed-source controls restore H2 performance under the same low rate, as documented in the current interpretation above. The result is a fixed-profile design sensitivity, not a prediction that future hardware cannot be retuned. The bounds carry no assigned probabilistic coverage.
