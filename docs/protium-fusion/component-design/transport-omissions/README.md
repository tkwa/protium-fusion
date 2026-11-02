# Omitted transport on the saved opaque-pulse histories ✨

Thermal neutrino cooling and electron conduction are too small on the ten tested histories to motivate a pulse or station redesign. The largest calculated thermal-neutrino loss is **4.95 × 10⁻¹² of prepared energy**. The deliberately generous electron-conduction model transfers at most **3.82 × 10⁻⁵ of prepared energy**, summed absolutely over interfaces, and produces a largest frozen-history local temperature perturbation of **5.91 × 10⁻⁵** (0.00591%). Conduction redistributes energy within the pulse; this interface sum is not an escaping-energy debit.

This is a saved-history significance screen. It uses the original corrected 64-cell four-stage family, the selected-003 mixture family A-H3/B1-H2/B2-H2.5/C13-H4, and its two modified source-profile runs. The source family shares the selected B1/B2 histories. All retain the archived ideal/pair EOS and original nuclear-rate family. Updated nuclear-rate and refined-mesh runs are outside this packet. The hottest source C13 state reaches 1.983 GK. This screening result does not qualify the source assembly, grey radiation model, or multidimensional burn.

## Results

Each dimensionless energy ratio uses that pulse's delivered prepared energy. The conductive quantity is the absolute heat crossing every internal interface, which can count repeated redistribution of the same energy. The local temperature column integrates the absolute conductive heating or cooling of each cell using its actual total heat capacity, then reports the largest cell value.

| Saved case | Peak T (GK) | Thermal neutrino / prepared | Absolute conductive transfer / prepared | Largest local ΔT/T |
|---|---:|---:|---:|---:|
| original_B2 | 1.647 | 2.6e-12 | 1.84e-05 | 4.13e-06 |
| original_C13 | 1.598 | 2.7e-12 | 2.12e-05 | 7.24e-06 |
| original_A | 1.145 | 2.75e-13 | 9.95e-06 | 2.39e-06 |
| original_B1 | 0.460 | 5.59e-19 | 3.01e-05 | 6.53e-06 |
| selected003_A | 1.186 | 3.76e-13 | 1.43e-05 | 3.01e-06 |
| selected003_B1 | 0.484 | 5.39e-19 | 3.82e-05 | 7.15e-06 |
| selected003_B2 | 1.646 | 2.36e-12 | 1.96e-05 | 4.08e-06 |
| selected003_C13 | 1.513 | 4.85e-12 | 2.6e-05 | 8.74e-06 |
| source64_A | 1.186 | 3.78e-13 | 1.42e-05 | 3.22e-06 |
| source64_C13 | 1.983 | 4.95e-12 | 2.42e-05 | 5.91e-05 |

The three stage-weighted cycle values use their already archived material balances and capacities, without changing those ledgers:

| Family | Additional thermal neutrinos (MeV/net He) | Thermal neutrino power at archived capacity (W) |
|---|---:|---:|
| original64 | 1.88e-12 | 1.14e+08 |
| selected003 | 1.76e-12 | 1.26e+08 |
| source64 | 1.76e-12 | 1.26e+08 |

These losses are separate from the nuclear beta-decay neutrino energy already included in the cycle. Even their full thermal-energy debit is negligible beside the roughly 10²¹ W plant energy flows. A coupled recalculation could change the detailed loss values without affecting this conclusion.

## Thermal-neutrino calculation

The calculation uses `pynucastro 3.1.0`'s port of Timmes's `sneut5`, implementing Itoh et al.'s pair, plasma, photo, bremsstrahlung and recombination processes. Inputs are actual saved density, temperature and composition moments. The code accepts density in g/cm³ and returns erg/g/s; the conversion is `1 erg/g/s = 10⁻⁴ W/kg`. Each component is integrated over cell mass and time. The source implementation and every input trace are hashed. [Timmes implementation description](https://cococubed.com/code_pages/nuloss.shtml), [pynucastro implementation and numerical example](https://pynucastro.github.io/pynucastro/neutrino-cooling.html), [Itoh et al. 1996](https://adsabs.harvard.edu/pdf/1996ApJS..102..411I).

The effective ion uses `1/abar = ΣY` and `zbar/abar = ΣZY`. A second evaluation preserves both electron density and `ΣZ²nᵢ` to test the mixture treatment of bremsstrahlung. Pair emission dominates the hot cases; for the hottest source C13 run it supplies 99.93% of the loss. Its effective-ion alternative changes total cooling by only 1.35 × 10⁻⁷ relatively. The composition copy follows the existing EOS convention of nonnegative, baryon-normalized abundances; saved accepted abundances remain untouched.

The installed fit returns zero below 10 MK. Evaluating these cells at 10 MK separately gives a negligible diagnostic contribution. Some cold outer states also lie below the photo-fit density range. At most 2.26 × 10⁻⁷ of a case's calculated cooling comes from `rho Ye < 1 g/cm³`; this is an applicability flag rather than a new fitted prescription. Itoh's quoted fit accuracy applies where the respective mechanism dominates, so the tiny result should not be represented with the printed numerical precision as a physical uncertainty.

Log-interpolated midpoint temperatures and densities change the time integral by 0.17–4.85% across the ten cases. This tests the saved sampling, not new hydrodynamic time resolution. A deliberately generous additional screen maximizes `neutrino power / EOS internal energy` over a 16 × 32 temperature/density grid and composition corners, then applies it to the largest saved internal energy for the entire pulse. Its maximum is **4.51 × 10⁻⁹ of prepared energy**. This is a sampled envelope over the frozen domain, rather than a certified bound on arbitrary future trajectories. A separate all-mass-at-peak-temperature corner test is retained in the JSON for transparency; it combines thermodynamically incompatible low-density/high-temperature states and is not used as a realistic prediction.

## Electron-conduction calculation

The conductivity screen takes the larger of two published limiting forms: classical Spitzer–Härm electron-ion transport in the Lorentz limit, and relativistic-degenerate electron-ion transport. It sets the Coulomb logarithm to one and omits electron-electron resistance. The degenerate collision frequency is

$$
\nu_{ei}=\frac{4\pi[e^2/(4\pi\epsilon_0)]^2\sum_i Z_i^2n_i\,\Lambda}{p_F^2v_F},\qquad
\kappa_{deg}=\frac{\pi^2k_B^2Tn_e}{3m_*\nu_{ei}}.
$$

Here `pF = ħ(3π²ne)^(1/3)`, `m* = sqrt(me² + pF²/c²)`, and `vF = pF/m*`. The classical Lorentz coefficient deliberately exceeds the ordinary hydrogen Spitzer–Härm value. An additional Maxwell-Boltzmann equilibrium pair-carrier multiplier reaches 1.577 in the hottest source trace; its scattering resistance is omitted. This is a generous interpolation sensitivity, not a calibrated conductivity table or a strict universal upper bound. [Spitzer and Härm 1953](https://doi.org/10.1103/PhysRev.89.977), [Cassisi et al. 2007, equations 2 and 8](https://arxiv.org/pdf/astro-ph/0703011), [Potekhin, Pons and Page 2015, appendix A.3](https://arxiv.org/pdf/1507.06186).

For each saved time the maximum neighboring conductivity multiplies the actual center-to-center temperature gradient. Spherical interface area gives the heat flow. Flux is zero at the symmetry center and free-vacuum surface. This measures internal collisional transport; it does not model an external plasma contact or nonthermal electron escape. All interface exchanges cancel in the total internal-energy balance and produce nonnegative entropy. Two heat capacities are retained: full ions/electrons/pairs/radiation, and electrons/pairs alone. The largest electron-only local change is **1.18 × 10⁻⁴** (0.0118%), still small.

The traces include partially and strongly degenerate cells; the source C13 `T/TF` range is 0.0347–254. The dominant-ion coupling diagnostic reaches Γ ≈ 35 in that source case, below the approximately 175 one-component freezing value; intermediate coupling remains a conductivity-model limitation. A simple impact-parameter logarithm is also saved, without treating it as an ion-correlation calculation. The selected limits explicitly include degeneracy rather than applying a classical formula throughout.

The effective heat-transport mean free path `3κ/(Cv,e v*)` stays below **1.72 × 10⁻⁶ of the sampled temperature-gradient scale**. Its largest absolute value is about **23 µm**; minimum zone widths are 0.139–0.305 m. Thus these sampled gradients do not call for collisionless saturated transport. Full-pulse, mass-99th-percentile diffusion lengths are **2.1–3.5 cm** with total heat capacity and **6.1–14.2 cm** using electrons alone. These lengths compare diffusivities over time; a length divided by a minimum width attained at another time is only a conservative geometric diagnostic. The local heating integral above is the more direct impact measure.

As a deliberately severe sensitivity, multiplying every conductivity by **100** would give at most **0.382%** summed redistribution and **0.591%** one-temperature local change on the same histories; the electron-only maximum would be **1.18%**. This multiplier is a stress test, not a subjective 80% credible interval. A coupled temperature response is not exactly linear, and unresolved gradients could differ from saved gradients. The wide margin makes a more precise conductivity interpolation a low priority relative to the existing assembly, nuclear-rate and radiation-transport questions.

## Reproduction and checks

`summary.json` gives the compact results and three cycle-weighted losses. `sources.json` explicitly identifies the ten original records/traces and three family maps. Each case JSON retains source and code hashes; its `.diagnostic.npz` contains density, temperature, composition moments, geometry, heat capacities, neutrino components and conductive face powers. These compact arrays permit independent integration and conservation checks without the full raw hydro archives. No hydro or historical output was changed.

The independent checker passes **78/78** checks: the published pynucastro numerical example, cgs/SI normalization, three degenerate-collision algebra checks, three degenerate-conductivity checks, classical hydrogen normalization, and seven checks for each pulse. Per-pulse checks independently recompute mass-weighted cooling, interface redistribution, local temperature perturbation, spherical geometry, energy-transfer cancellation and entropy sign, and verify archived array hashes.

Use Python with numpy/scipy and `pynucastro==3.1.0`. From the research root `docs/protium-fusion`, with the existing `microphysics.py` and EOS dependencies:

```sh
python component-design/transport-omissions/validate.py \
  --research-root . \
  --sources component-design/transport-omissions/sources.json \
  --results component-design/transport-omissions \
  --output /tmp/transport-validation.json
```

To rerun the complete microphysics extraction, use `screen.py` with the same `--research-root` and `--sources`, and a new `--out-dir`; `--case` can select one explicit case. That command requires the full original traces at the declared source paths, the archived FastEOS table, and the stated pynucastro version. The independent archived-array check above does not require those full raw traces. `manifest.json` records every delivered file. The repository's science environment uses CPython 3.10; system Python 3.14 is incompatible with its existing scientific binary packages.

This screen supports retaining the current omission of thermal-neutrino cooling and electron conduction for these specific pulse families. Re-evaluate if new source or rate variants create substantially hotter, longer-lived or much steeper-gradient states.
