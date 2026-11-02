# Carbon-reserve and B2-grid scenarios ✨

The frozen `process.json` remains the corrected 64-cell process family. `process-scenarios.json` compares explicit B2 spatial grids and carbon supply reserves; it does not silently replace the original family. A/B1/C13 remain at corrected 64 cells in every current scenario. All included B2 histories complete four sound times, but the 64-to-128 change means spatial convergence has not been established.

`process_scenarios.py` exposes:

```
scenario_case(stages, fuel_C_H_target=.0056, total_C_H_budget=.007, **kwargs)
```

The returned dictionary retains every flat field from `process.case`, including stages, pulse metadata, heat, prepared energy, fresh hydrogen, return fraction and electrical scale factors. It adds:

- `net_structural_carbon_allowance_kg_per_He` and `_kg_per_preaux_J`;
- `total_fresh_carbon_budget_kg_per_He` and `_kg_per_preaux_J`;
- `fuel_C_H_target`, `total_ongoing_C_H_mass_budget` and `carbon_budget_reserved_fraction`.

The JSON path is `families.<family>.cases.<target>`, with targets formatted as `"0.0056"`, `"0.0063"`, `"0.0070"`. The completed families are `b2_64`, `b2_128_energy_matched` and `b2_256_energy_matched`. The original 64-cell JSON remains unchanged; the two-family checkpoint is preserved as `process-scenarios-through128.json`. The source mapping is `b2-process-scenario-map.json`.

## Carbon budget

The user-selected ongoing budget is .007 kg **total fresh carbon** per kg fresh protium. Fuel targets .0056, .0063 and .007 reserve 20%, 10% and 0% of that budget for other ongoing net carbon losses. The reserve uses the actual fresh hydrogen denominator of each newly solved case:

```
structural_C_allowance = (total_C_H_budget - fuel_C_H_target) * fresh_H.
```

This includes receiver, coating, substrate, processing and other carbon makeup after actual recovery and reuse. It does not identify every replaced component with new imported material. Construction and startup inventory are separate one-time carbon requirements. If a global lifetime abundance budget includes those inventories too, their amortized charge must be added explicitly.

At the default 60% driver, 40% thermal conversion and 99% energy capture:

| B2 grid | Fuel C:H mass target | Required CNO/H return | Preauxiliary MeV/He | Reserved structural C, kg/preauxiliary J |
|---|---:|---:|---:|---:|
| 64 | .0056 | 99.9627607% | 6.8764449 | 8.51585e-18 |
| 64 | .0063 | 99.9577830% | 6.8772514 | 4.25798e-18 |
| 64 | .0070 | 99.9528048% | 6.8780581 | 0 |
| 128, energy matched | .0056 | 99.9622637% | 6.9122554 | 8.47174e-18 |
| 128, energy matched | .0063 | 99.9572195% | 6.9130780 | 4.23592e-18 |
| 128, energy matched | .0070 | 99.9521749% | 6.9139008 | 0 |
| 256, energy matched | .0056 | 99.9619829% | 6.9329133 | 8.44651e-18 |
| 256, energy matched | .0063 | 99.9569013% | 6.9337453 | 4.22330e-18 |
| 256, energy matched | .0070 | 99.9518190% | 6.9345775 | 0 |

For example, the 64-cell .0056 case permits about 8,516 kg/s of additional net carbon loss per 10^21 W preauxiliary output, or 2.69e11 kg/year. If a station holds 10^23 kg of structural carbon, that is only 2.69e-12 of its stock per year. Substantial component replacement is compatible with this budget only if the replaced material is almost completely recovered and reused. Compared with gross pulse-fuel throughput, the allowance is about 7.9e-5 kg extra carbon loss per kg circulated fuel. Actual receiver erosion, coating recovery, contamination handling and escape must be assessed against that allowance.

## Spatial sensitivity

The 128-cell B2 case preserves the 64-cell initial profile and matches its total injected energy, 4.07439303594e25 J. Its prompt nuclear/prepared gain is 11.1018937 and desired-product yield .87026450, compared with the 64-cell gain 10.55503 and yield .82188448. The corresponding full-process electric energy increase is about .52%, because B2 is one of several stage costs and much of the cycle's nuclear release is fixed by its overall reaction.

The completed energy-matched 256-cell B2 run gives gain 11.4327262 and desired-product yield .89856039. Its energy residual relative to injected energy is -4.33e-6. At fuel C:H=.0056, the complete process gives H=24.7185481, U=2.2473492 and W=6.9329133 MeV/He: W is .0564684 MeV/He (+.8212%) above 64 cells and .0206579 MeV/He (+.2989%) above 128 cells. Fewer B2 retries reduce preparation energy. At the fuel-only .007 target, the corresponding increases are +.8217% and +.2991%.

These comparatively small changes in full-process energy do not make the spatial solution converged. It also does not validate receiver arrival histories, multidimensional merging, or apparatus efficiencies. `comparison_with_b2_64` and the 256-family `comparison_with_b2_128_energy_matched` record energy and return-fraction changes without conflating them with a hardware improvement.

The module accepts `--b2-map` containing an explicit family-to-source-JSON mapping, including `b2_64`. The default loads the explicitly saved `b2-process-scenario-map.json`; every source hash and the original frozen process hash are recorded. All non-B2 source files remain identical across these scenarios.
