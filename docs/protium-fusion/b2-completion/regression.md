# Four-stage viscosity regression ✨

All four previously completed small stages—A, B1, C13 and the 3 m helium factory—complete their requested four sound times with the corrected continuous compression-cap viscosity. **Their deposited nuclear gains change by less than 0.4%, and their principal product yields by less than 0.43%.** Escaped-radiation fractions change more substantially and must be carried into the receiver calculation. These results replace the four corresponding legacy outputs when assembling a cycle with the corrected B2 calculation.

This is a narrow regression of the existing prescribed hot inflows. It adds no driver model, new fuel search, multidimensional stability result or complete reactor claim. `He` here denotes the helium-4 triple-alpha carbon factory at a **3 m reference radius**; it does not denote helium-3 fuel.

## Correction and exact scope

The old scalar artificial pressure used a negative difference between the two radial boundary velocities as its compression trigger. In spherical flow, that velocity difference can be negative while the cell volume is increasing. The old pressure then makes a negative contribution to internal heating. It can conserve total energy while violating the intended dissipative sign, so a small global energy residual does not detect the error.

The corrected prescription used here is

```
C_radial = max(-Delta v, 0)
C_volume = max(-(dot V / V) Delta r, 0)
C = min(C_radial, C_volume)
q = rho (a_quadratic C^2 + a_linear c_s C)
```

The original coefficients, 2.0 and 0.5, are retained. The same effective pressure remains in force and pressure work. Since q vanishes in volume expansion, `-q dot V / m` is nonnegative. The cap is continuous at zero volume compression, unlike a hard switch that can jump while the radial velocity difference remains negative. It preserves the old measure for homologous compression and approaches it in a thin planar shock. It also limits some nonhomologous compressive states, so this is a changed viscosity prescription rather than a subtraction of the previously erroneous heat alone.

The preliminary hard-switch implementation failed the BDF integration in B2 and in the helium regression; it is excluded from the results. Three other preliminary regression attempts lost their output while saving traces because `/tmp` reached its quota. Their solver status is not inferred from those failures. The completed results below were rerun in `/var/tmp/tkwa-chatgpt/protium-b2-regression`, using metadata saved before trace output.

For all four production regressions:

- The saved scaled initial radius, velocity and temperature arrays were loaded **exactly**, without reconstructing the scale transformation.
- The 64-zone grid, full 27-isotope network, original EOS/rates, radiation diffusion, constant opacity and requested four sound times were retained.
- A, B1 and C13 retain 30 m reference radii; the helium factory retains 3 m. Actual initial surfaces are larger for these prescribed convergent inflows.
- Tolerances changed from the legacy `rtol=1e-6, atol=1e-11` to `rtol=1e-7, atol=1e-12`. Two passive nuclear-energy/weak-charge diagnostic variables were added by the corrected runner. The comparison therefore includes those numerical changes; it is not an exact same-tolerance one-line experiment.
- Each corrected trace retains 120 accepted-step snapshots. This sampling changes storage, not the integration.

The exact corrected module used for the four runs has SHA-256 `3af8f36fad063bde4ac5b3491c324c9add973a5459ecb1a9e913630564eb70cc`. A pinned copy is in the scratch directory's `model/radial_completion.py`; the B2 continuation owns its canonical source archive. Input profiles, legacy results, wrappers, EOS, network and trace hashes are recorded in the per-run inputs and [comparison JSON](regression-summary.json).

## Gains, yields and energy balance

Gain means `(nuclear mass-difference energy - the assigned escaping weak energy) / initial internal-plus-kinetic energy`. It does not deduct photons that leave the pulse, because their recovery is a separate receiver calculation.

| Stage | Legacy gain | Corrected gain | Relative change | Legacy principal yield | Corrected principal yield | Yield change |
|---|---:|---:|---:|---:|---:|---:|
| A | 8.1816413 | 8.1963390 | +0.17964% | 0.5384000 | 0.5393672 | +0.17964% |
| B1 | 6.3672050 | 6.3819705 | +0.23190% | 0.6004770 | 0.6018693 | +0.23186% |
| C13 | 15.4794667 | 15.5389791 | +0.38446% | 0.8131295 | 0.8166214 | +0.42944% |
| He, 3 m | 0.0009525407 | 0.0009522103 | −0.03469% | 0.0003494121 | 0.0003492909 | −0.03469% |

The principal products are O15 for A, C12 for B1, N14 for C13 and C12 for the helium factory. The first three yields are product counts per initial catalyst nucleus. The helium entry is the carbon mass fraction, following the existing factory normalization. Full isotope vectors are retained; the cycle should use those vectors rather than just the principal yields.

| Stage | Corrected nuclear energy, J | Total-energy residual / initial energy | Completion |
|---|---:|---:|---|
| A | 2.8631152×10^26 | −1.347×10^-7 | Four sound times |
| B1 | 2.0390503×10^24 | +1.215×10^-6 | Four sound times |
| C13 | 5.1717862×10^26 | −3.473×10^-6 | Four sound times |
| He, 3 m | 2.3107423×10^22 | +2.377×10^-8 | Four sound times |

The independently accumulated reaction energy agrees with the final isotope mass-difference energy within **5.15×10^-11 of initial energy** across these cases. The maximum stored-trace baryon residual is 1.74×10^-14. The minimum stored mass fraction is −2.10×10^-13, within the declared 10^-12 trace tolerance; the largest final-state negative fraction is much smaller. Initial masses and energies are unchanged to the recorded precision.

The [comparison script](regression_compare.py) reports **44/44 passing checks**: full-duration success, energy residual, independent nuclear-energy accumulation, endpoint and sampled baryon balance, trace abundance tolerance, unchanged initial energy/mass, positive shell widths and absence of scalar artificial pressure in expanding stored cells. Nonnegative artificial heating also follows algebraically from the corrected cap. These checks do not establish grid convergence or physical validity of the assumed hot inflow.

## Radiation and receiver consequences

The radiation partition changes more than the nuclear gain:

| Stage | Legacy escaped radiation / initial energy | Corrected escaped radiation / initial energy |
|---|---:|---:|
| A | 0.00171984 | 0.00775600 |
| B1 | 0.340415 | 0.427329 |
| C13 | 0.00976042 | 0.02231726 |
| He, 3 m | 0.000046336 | 0.000066272 |

For B1, final kinetic energy falls by 1.086% despite the small increase in nuclear gain. A's and C13's final internal energies rise by approximately 3.0% and 2.6%. Escaped photons are energy leaving the pulse; whether the installation recovers them depends on its collector. The regression therefore supports using the corrected full metadata in the heat/kinetic/radiation ledger, rather than multiplying the old collector results by a gain ratio. No updated complete-cycle or receiver result is claimed by this report.

## Independent legacy sign check

The [sign audit](regression_sign_audit.py) applies the geometric tests `Delta v < 0` and `dot V > 0` to the original 160 stored snapshots for each stage. Positive dimensional scale factors cancel. Mass weighting uses the original shell masses, proportional to differences of radius cubed.

| Stage | Snapshots containing negative scalar-viscous heating | Maximum affected mass fraction |
|---|---:|---:|
| A | 64 / 160 | 7.996% |
| B1 | 71 / 160 | 4.831% |
| C13 | 61 / 160 | 4.615% |
| He | 59 / 160 | 4.615% |

All four legacy endpoints themselves have zero affected fraction. Endpoint-only diagnostics would therefore miss the earlier violation. This audit identifies its sign and occurrence; it does not integrate the erroneous heat or assign the full gain difference to it. Reproduction data and source hashes are in [regression-legacy-sign.json](regression-legacy-sign.json).

## Files and reproduction

The unchanged corrected stage metadata are archived as [A](regression-A-compression_cap.json), [B1](regression-B1-compression_cap.json), [C13](regression-C13-compression_cap.json) and [He](regression-He-compression_cap.json). Raw NPZ traces, inputs, logs and pretrace metadata remain under `/var/tmp/tkwa-chatgpt/protium-b2-regression`. The original `/tmp/protium-small-stages` traces are read-only inputs to the sign audit.

Run [regression_runner.py](regression_runner.py) once for each of `A B1 C13 He`, with `--root` pointing to `docs/protium-fusion`, `--model-directory` pointing to the directory containing the exact corrected `radial_completion.py`, `--output` pointing to a writable scratch directory and `--mode compression_cap`. The existing dependency runtime is selected with `PYTHONPATH=/tmp/protium-python`; set `NUMBA_CACHE_DIR` and `MPLCONFIGDIR` to writable scratch directories and limit BLAS/OpenMP to one thread per run. Then run the sign audit and comparison script with the corresponding trace, scratch and output paths. Their input hashes make changed inputs visible.

Only new regression files in the dedicated worktree's `docs/protium-fusion/b2-completion/` were written. Historical source/results, the article, unrelated blog edits and commits were left untouched.
