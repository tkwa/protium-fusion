# ✨ B1 burn after the delivered assembly profiles

The three exact delivered profiles all sustain the modeled B1 burn. Relative to the delivered control with preparation radiation disabled, preparation at κ = 0.001 m²/kg reduces deposited nuclear output by **0.9112%**; preparation at κ = 0.01 m²/kg reduces it by **0.08185%**. These are conditional numerical responses of the stated profiles and opacity cases, rather than an uncertainty interval.

| Preparation radiation | Preparation escape / input | Q over burn-start energy | Q over common preparation input | C12 per initial N15 | Nuclear output vs delivered control |
|---|---:|---:|---:|---:|---:|
| Radiation disabled | 0.000000% | 8.02938352 | 8.02938352 | 0.80472062 | +0.000000% |
| κ = 0.001 m²/kg | 0.002285% | 7.95640027 | 7.95621847 | 0.79738236 | -0.911216% |
| κ = 0.01 m²/kg | 0.000328% | 8.02283751 | 8.02281123 | 0.80406123 | -0.081853% |

Here Q means `(nuclear release − the model's conservative weak-energy escape) / prepared energy`; it is a pulse energy ratio. All three hot runs use the **same κ = 0.001 m²/kg**, 128 cells, actual H/N = 2, median updated carbon captures with paired reverse rates, compression-cap artificial viscosity (linear 0.5, quadratic 2), and **82.2779579 μs** physical burn duration. The opacity differences in the table apply only to the preceding 30 μs preparation. No velocity, temperature, density, or energy retuning was performed. Within this diffusion model, the smaller preparation opacity permits more energy escape; these two constants are comparison cases, not measured material opacities.

The original constructed source128 hot profile had Q = 8.08234646 and C12 conversion 0.80969456. The delivered no-radiation control produces 0.6146% less nuclear energy than that original profile. This separate control includes the preceding equation-of-state handoff and spherical preparation round trip. Relative to the original hot reference, the two radiative preparation outputs are lower by 1.5202% and 0.6959%. Comparing each with the delivered control isolates the incremental preparation-radiation response.

## Energy handoff

The FEOS → FastEOS handoff occurs at the common precursor, 30 μs before these delivered profiles. Both the forward preparation and all three hot burns already use FastEOS; no further equation-of-state conversion occurs at burn start.

- Original constructed hot source FastEOS energy: **3.195016782566003 × 10²³ J**.
- Common precursor direct FEOS energy: **3.202230949940097 × 10²³ J**.
- Common precursor FastEOS energy: **3.196326437256525 × 10²³ J**.
- Shared precursor radial kinetic energy: **3.105768545281137 × 10²³ J**.
- Raw model handoff difference, FastEOS minus FEOS: **−5.904512683572 × 10²⁰ J**, or **−0.184387%** of the direct FEOS total.
- Mass-weighted FastEOS/FEOS pressure ratio at that precursor: **0.99138384**; maximum local ratio is about 1.113 in dilute outer material.

The two EOS implementations have their own internal-energy reference conventions. Their exposed difference is a model handoff discrepancy; assigning it physical heater work or recoverable energy would require a consistent reference alignment. The numerical continuation ledger uses the common FastEOS precursor energy and retains the raw FEOS quantities separately. Using the raw FEOS total as an alternative denominator would reduce each preparation-based Q by 0.184387%, while leaving the calculated nuclear energy unchanged.

Preparation escaped radiation is respectively **0**, **7.303449010114459 × 10¹⁸ J**, and **1.047115722729 × 10¹⁸ J**. Each burn-start `U + K`, plus that already escaped energy, closes to the common FastEOS preparation input within **5 × 10⁻¹¹**. The ledger does not credit receiver recovery of this radiation. Combined preparation-plus-burn energy residuals are below **1.84 × 10⁻⁶** of the precursor input.

## Validation and scope

`summary.json` records **51 passing checks**: exact delivered T/r/v retention, exact 128-cell mass grid, correct H/N = 2, starting trace/composition, duration completion, settings, baryon conservation, source/trace hashes and energy closure. `provenance-regeneration.json` adds seven byte-identity checks after replacing a correct literal reference energy with its hashed source record. Economics independently reviewed the preparation and continuation ledger.

These runs address the nuclear response to the specified delivered profiles. The precursor remains an auxiliary spherical whole-fluid construction with a fixed vacuum outer boundary. Its simple backward branch does not show contact crossing in the outer 17.6% of B1 mass. Cold launch, merger, heat transfer between components, multidimensional stability, and a consistent FEOS-to-FastEOS physical energy reference remain separate source-design questions. The hot calculation retains the established simplified opacity, EOS, and 27-species nuclear network assumptions. No plant output or full source qualification follows from this checkpoint alone.

## Reproduction and records

Run `prepare.py`, then `run.py --label k0`, `--label k0p001`, and `--label k0p01`, followed by `summarize.py`, using the recorded Python 3.10 / NumPy 2.2.6 / SciPy 1.15.3 / pynucastro 3.1.0 research runtime. All runner/table/source fingerprints are in the inputs and records. `bridge.json` contains the full per-cell EOS comparison; `inputs/` preserves delivered profiles; `records/` contains completed hot endpoints. Full traces are under `/home/tkwa/src/protium-fusion/output/assembly-burn/`, with their SHA-256 hashes in the canonical records. Earlier assembly and source records were preserved unchanged.
