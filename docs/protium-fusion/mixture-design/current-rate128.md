# Median-carbon 128-zone material, receiver and plant integration ✨

The best tested mixture is **A: 3 H/N14, B1: 2 H/N15, B2: 2 H/C12, C13: 4 H/C13** at both receiver masses. The conditional plant calculation gives **6.6406 mW/kg at 5000 kg/m²** and **12.5105 mW/kg at 2000 kg/m²**. Both reach the 90% radiator-area allocation before the 80% principal pulse-window occupancy limit.

This packet contains ten explicitly named, completed 128-zone prescribed-flow runs with the median updated carbon-capture rates. It evaluates all 36 four-stage combinations through the actual-composition material and receiver reducers, then the existing plant ledger. Source-mapped runs, 64-zone studies, the frozen-rate family and the separate 256-zone B2 run remain separate. These results establish a comparison within this model family; they do not establish mesh convergence or a constructed launcher/receiver.

## Capacity points

| Quantity | Receiver 5000 kg/m² | Receiver 2000 kg/m² |
|---|---:|---:|
| Net export | 8.09772 × 10²⁰ W | 8.10976 × 10²⁰ W |
| Total allocated mass | 1.219424 × 10²³ kg | 6.482376 × 10²² kg |
| Specific export | 6.640608 mW/kg | 12.510481 mW/kg |
| Helium production | 4.491705 × 10³²/s | 4.494779 × 10³²/s |
| Main-window occupancy | 0.721236 | 0.721729 |
| Radiator fraction | 0.900000 | 0.900000 |
| Prepared energy | 1.895680 MeV/net He | 1.895680 MeV/net He |
| Fresh carbon/H mass ratio | 0.000343745 | 0.000343745 |
| Fresh H | 3.00715 × 10⁶ kg/s | 3.00921 × 10⁶ kg/s |
| Fresh carbon | 1033.69 kg/s | 1034.40 kg/s |
| Returning unburned H | 3.46907 × 10⁶ kg/s | 3.47145 × 10⁶ kg/s |
| Returning quartet CNO | 3.70490 × 10⁷ kg/s | 3.70744 × 10⁷ kg/s |
| Gross fuel handling | 4.34730 × 10⁷ kg/s | 4.35028 × 10⁷ kg/s |

Common assumptions are radius 1.1 × 10⁹ m, 200 K radiator, 60% thermal conversion, 60% preparation efficiency, 99% energy capture and eventual return of retained ordinary H/CNO. The existing mass ledger rebuilds driver peak power, storage, circulation inventory, reserve, replacement and support allocations. Changing the receiver areal mass changes those coupled allocations too. Receiver thermomechanical qualification remains a separate component calculation. The surface-deposition prompt-temperature proxy, approximately 1148 K rise for the most demanding pulse, is recorded as a diagnostic rather than used as a qualified apparatus limit.

## Mixture comparisons

The following changes hold A-H3 and C13-H4 fixed, change one of the two middle stages, and recompute station capacity. The percentage is relative to the best quartet at the same receiver mass.

| Changed stage | 5000 kg/m² (mW/kg) | 2000 kg/m² (mW/kg) | Change in specific power |
|---|---:|---:|---:|
| B1 H2, B2 H2 | 6.640608 | 12.510481 | Reference |
| B1 H1.5, B2 H2 | 6.633048 | 12.496244 | −0.114% |
| B1 H1, B2 H2 | 6.435438 | 12.116637 | −3.09% / −3.15% |
| B1 H2, B2 H2.5 | 6.517967 | 12.280959 | −1.85% / −1.83% |

**B1 H2 has a narrow numerical preference over H1.5.** The latter reduces returning H from 3.469 million to 3.071 million kg/s in the 5000 kg/m² comparison, lowers fresh C:H to 0.000334492, and raises prepared energy to 1.902799 MeV/He. A 0.114% plant difference is too small to treat as a robust engineering optimum given the remaining source, mesh and processing uncertainties. Both choices should remain live candidates.

**B2 H2 is materially better than H2.5 within this family.** H2.5 raises prepared energy to 2.014177 MeV/He and returning H to 4.247 million kg/s at 5000 kg/m². Fresh C:H is virtually unchanged. The conditional integrated benefit of H2 is larger than the B1 distinction; its final selection should still be checked against the separate refined B2 and source-profile results.

No probability distribution is inferred from this finite grid. These are deterministic conditional comparisons, not subjective credible intervals or a continuous mixture optimum. The updated rates are the two carbon-capture median tables described in [nuclear-review.md](nuclear-review.md); the other 190 reaction rates are unchanged.

## Retained material and processing

Fresh C:H uses only **4.91% of the available 0.007 mass-ratio budget**, leaving about 0.006656 per unit fresh H for other uses. This residual is an allowance, not a requirement to consume carbon. The exact eventual-return case treats incomplete first-pass extraction, adsorption and implanted material as retained reservoirs. The fresh-carbon term mainly balances catalyst transferred to side-product inventories outside the four-stage circulation. It is not a prediction of ordinary matter escaping through the wall.

The existing H-processing sensitivity at 1 MJ–1 GJ/kg adds 3.47 × 10¹²–3.47 × 10¹⁵ W at the best point, approximately 4.28 × 10⁻⁹–4.28 × 10⁻⁶ of export. These loads are compared with the existing ordinary-auxiliary allowance; they are not automatically debited again. The allocated ordinary auxiliary load is about 9.01 × 10¹⁹ W. Local gas routing, residence time and cold handling remain geometric and hardware questions even when chemical energy is a small fraction of output.

A fixed-capacity algebraic sensitivity gives roughly **2.32 × 10¹² J/kg of additional returned-H processing** before the H1.5 case overtakes H2 in specific power. This compares the two saved points after an extra linear processing debit; it does not re-optimize radiator loading or processing hardware. It demonstrates why ordinary chemical/ionization energy costs alone do not erase the H2 benefit, while the small numerical preference remains sensitive to other component effects.

## Cadence and retained tails

At 5000 kg/m² the mean intervals are A **8.99 days**, B1 **123.07 minutes**, B2 **10.57 days**, and C13 **107.75 days**. The 2000 kg/m² capacity changes these by less than 0.1%. Thus the old approximately 111-minute B1 cadence should not be carried unchanged into this quartet.

The 0.721 occupancy counts the incoming and outgoing principal 99.9% mass windows from the actual profiles. Their omitted slow material remains retained. The resulting slack is useful for chronological scheduling, but does not itself prove full cavity clearance or route every central slow tail. Existing protected-path/recapture options must be evaluated with the selected profile. This packet checks mean cadence and occupancy; it does not regenerate a full year-long event calendar.

## Validation and reproduction

All ten receiver reductions complete with monotone incoming/outgoing profiles and no warnings. They reconstruct source energy with actual initial composition, recover the named prompt-radiation accumulator without rescaling, and retain the nodal kinetic ledger alongside the smooth mass-conserving arrival reconstruction. All 36 combinations pass the material baryon/mass-defect/energy checks and plant mass/global-energy/occupancy checks: **238/238 checks per receiver-mass scenario**. An independent checker additionally validates 274 identities and replay/boundary checks (see `current-rate128-validation.json`).

The first output attempt encountered a NumPy integer serialization error in the new wrapper's check count. That wrapper-only issue was corrected; the final fresh extraction and reduced-input replay give identical physical results. No archived solver, material, receiver or plant implementation was changed.

Files:

- `current-rate128-source-map.json`: ten explicit record/profile/input/trace hashes.
- `current-rate128-receivers/`: ten actual-composition receiver reductions.
- `current-rate128-reduced-inputs.json`: complete material/receiver/record inputs for replay without raw hydro traces.
- `current-rate128-sigma5000.json` and `current-rate128-sigma2000.json`: all 36 rankings, ten detailed leading combinations, controlled mixture changes, checks and provenance.
- `current_rate128_integrate.py`: explicit extraction or reduced-input replay.
- `current_rate128_verify.py`, `current-rate128-validation.json`: independent accounting and boundary checks.
- `current-rate128-manifest.json`: delivered-file hashes.

From `docs/protium-fusion`, with the existing scientific Python environment:

```sh
python mixture-design/current_rate128_integrate.py \
  --research-root . \
  --replay mixture-design/current-rate128-reduced-inputs.json \
  --output-dir /tmp/current-rate128-replay
```

Use `--manifest mixture-design/current-rate128-source-map.json` instead of `--replay` to repeat raw extraction; that path requires the declared original traces. To audit archived outputs, use `current_rate128_verify.py` with `--research-root .`, `--inputs mixture-design/current-rate128-reduced-inputs.json`, `--results-dir mixture-design`, and a new `--output` path. The packet records the exact material, receiver, plant and environment-ledger implementation hashes consumed.
