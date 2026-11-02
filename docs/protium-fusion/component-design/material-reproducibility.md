# ✨ Material and schedule packet

The main findings are in [material-loop.md](material-loop.md) and [schedule.md](schedule.md). The material ledger uses actual feed composition and eventual returns. The finite scheduling examples retain slow-tail material and delayed extraction in inventories.

`material_cycle.py` is the reusable API. `load_stage()` requires actual initial composition; `solve_cycle()` accepts separate eventual CNO/H returns; `solve_supply_target()` finds the maximum-loss boundary which spends a chosen carbon allowance. `material_cycle_compare.py` accepts stage maps and reports exact input/output isotope counts, net-helium normalization, fresh H/C, energy partitions, stage attempts, pulse counts and radioactive processing loads. It does not infer a recipe from its name.

The baseline module verification passes 45 checks. These cover the old ledger, a nonidentity recipe basis, excess hydrogen return, explicit multiple-daughter decay, requested physical duration and contradictory input rejection. The archived baseline composition map was extracted from the accepted initial trace states. An optional `--verify-raw-traces` flag repeats that extraction if the original large traces are available.

The original-family schedule and tail calculations are distinct from the provisional 003 mixture schedule. `material-plant-003.json` saves the complete conditional plant/cycle point. `mixture-schedule-003-input.json` also saves actual stage inputs/outputs and receiver timing; subsequent schedule and stock calculations need no raw traces. This mixture example is a 64-zone comparison point rather than a claim of mesh convergence or completed source/receiver engineering.

## Reproduction

Run with Python, NumPy, SciPy and the research network's `pynucastro==3.1.0`; the saved runs used CPython 3.10. The research root is `docs/protium-fusion`. Set `PROTIUM_RESEARCH` to its absolute path and place generated outputs in a separate directory.

```bash
python "$PROTIUM_RESEARCH/component-design/material_cycle_verify.py" \
  --root "$PROTIUM_RESEARCH" --output /tmp/material-cycle-recheck.json

python "$PROTIUM_RESEARCH/component-design/material_cycle_compare.py" \
  --root "$PROTIUM_RESEARCH" \
  --stage-map "$PROTIUM_RESEARCH/component-design/material-cycle-baseline-map.json" \
  --output /tmp/material-cycle-comparison.json

python "$PROTIUM_RESEARCH/component-design/graphite_pool.py" \
  --plant "$PROTIUM_RESEARCH/component-design/material-plant-003.json" \
  --output /tmp/graphite-pool-recheck.json

python "$PROTIUM_RESEARCH/component-design/mixture_schedule.py" \
  --input "$PROTIUM_RESEARCH/component-design/mixture-schedule-003-input.json" \
  --output /tmp/mixture-schedule-recheck.json \
  --events-csv /tmp/mixture-schedule-recheck.csv

python "$PROTIUM_RESEARCH/component-design/isotope_tracer.py" \
  --root "$PROTIUM_RESEARCH" \
  --history-root "$PROTIUM_RESEARCH/component-design/tracer-history" \
  --substeps 8 --output /tmp/isotope-tracer-recheck.json
```

The passive tracer histories are archived compact temperature/density/composition samples. Running without `--history-root` extracts them from the large raw traces, with `--cache` choosing a new cache directory. The archived final4/final8/zero-neutron results use one implementation and one rate table. `isotope_tracer_compare.py` compares their probability changes and conservation checks. Sparse saved-time interpolation, quantified by native-isotope reconstruction, limits the screen; no full hydro response to contamination is claimed.

The `tail_overlap.py` reconstruction can be repeated with `--quadrature-order 12` or `24`, using the historical `mixed-stage-schedule.json`. It integrates actual spatial mass dilution and relative collision energy. Its inputs come from the earlier archived collector/loading packets. This tail result applies to the original compositions and is distinct from the new-mixture chronological compatibility calculation.

## Boundaries

The source JSONs record hashes and the raw source identities. Analysis can be rerun from archived compact inputs even where provenance names an external original trace. Incomplete first-pass extraction and wall implantation create residence and inventory. Ordinary matter loss requires an explicit physical escape or deliberate retirement mechanism. Nuclear routing outside the selected C12/C13/N14/N15 quartet creates a held side-product stream, dominated by O16; an optional NO return loop has been identified stoichiometrically and receives no plant credit.

The graphite processing costs and 30–365 day batch durations are declared design scans. They are not subjective credible intervals, because the proposed irradiated-skin release kinetics have not calibrated a useful probability distribution. Natural wall carbon abundance is an explicit pool assumption, with C12 and C13 separately balanced. Prompt H/He routing is constrained independently in [return-gas.md](return-gas.md).
