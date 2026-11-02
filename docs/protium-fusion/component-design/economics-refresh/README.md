# Economics refresh for an explicit final plant selection ✨

This evaluator consumes the selected plant's mass, net electricity and material ledger. Source profiles, receiver capture, preparation radiation, heat routing and operating auxiliaries stay in the upstream plant solve. It produces a current economic report without altering earlier calculations.

## Exact interface

```python
from calculate import evaluate, render
result = evaluate(plant, cycle=None, assumptions=None)
markdown = render(result)
```

`plant` needs `mass_kg`, `net_W` and the complete `mass_breakdown_kg` with a `receiver` entry. Specific power, receiver area and fresh C:H identities are checked if supplied. Flow, driver and repair fields are consumed directly; missing optional sections remain unspecified. A plant marked infeasible is rejected. `cycle` is retained as an exact reference. Assumptions are explicit and replaceable in [assumptions.json](assumptions.json).

For the usual ranked plant document:

```sh
python3 docs/protium-fusion/component-design/economics-refresh/calculate.py \
  --input /absolute/path/to/final-plant.json \
  --select /top_detailed/0 \
  --out-dir /absolute/path/to/a/new/economics-report-directory
```

If the input is itself a `{ "plant": ..., "cycle": ... }` record, omit `--select`; if a selected plant is nested elsewhere, provide its explicit JSON pointer. `--assumptions` accepts a separate scenario file. The evaluator never picks a different reactor point itself. Existing report files are protected from replacement.

The output directory receives `selected-input.json`, `summary.json` and a generated `README.md`. The summary records source, selection pointer, input/code/assumptions hashes and primary sources. [regression64](regression64/README.md) is a regression example on the original 2,000 kg/m² case; it is not the latest source or receiver selection.

## Interpreting the refreshed result

Whole-plant embodied energy is an input scenario in J/kg. Simple energy payback is `embodied energy / (net specific power × availability)`. Manufacturing efficiency acts on a specified useful free-energy change, with 50–100% allowed in this implementation. The methane-to-graphite example is one chemical upgrade at declared standard conditions. Acquiring carbon, moving it, manufacturing the selected structures and installing the plant retain their own accounts. Recycled bulk graphite does not incur the initial methane-upgrade energy on every pass.

The receiver's carbon stock defaults to its full mass, with an explicit adjustable carbon fraction. At source C:H = 0.007, the associated hydrogen is a scale for resource access. This construction inventory is separate from operating fresh C:H. Services, reserves and other hardware require an elemental bill; the calculator does not silently assign them graphite composition.

Operating fuel C consumption comes directly from the selected network/material balance. The remaining 0.007 allocation is reported before identified terminal structural makeup, alongside the original 0.0056/0.0014 split. Retained processing tails, skins, working material and repaired components remain inventory. The calculator supplies no fictitious per-shot material loss.

Repair energy scenarios report increments above the repair electricity already present in the plant. Hydrogen processing scenarios mean **additional** energy above any existing auxiliary allocation. Neither is silently applied to the plant's net power: a materially different load would need a new thermal/capacity solve.

Construction scenarios report installation throughput, carbon acquisition and the external factory stock implied **if** the chosen repair productivity applies to construction. They do not infer a deployment year from energy efficiency. The solar comparison likewise exposes its collector mass, conversion and delivery assumptions; it supplies a scale for further equal-service comparison.

If the final plant reaches exactly **20 mW/kg**, then every **1 MJ/kg** of actual embodied energy takes **1.584 years** of export to repay at full availability, and twenty years supplies **12.623 MJ/kg**. These are general scaling coefficients, not an assertion that the final point has that specific power. At the declared 30%-efficient, 1 kg/m² collector comparison, mass-specific solar parity would be about **143 AU**.

The automatic report ranks mass items from the selected record. In a 500 kg/m² receiver design with a retained 1,000 kg/m² services allowance, services have twice the receiver's mass; the earlier conclusion that receiver carbon is the largest item would no longer hold. The complete new mass ledger determines the actual ranking.

## Verification

[validation.json](validation.json) records 33 regression and physical-scaling checks against both exact earlier plant cases, payback/availability and power scaling, variable carbon composition, missing-data behavior and rejection of inconsistent records. The command-line regression adds five ledger identities. `validate.py` reruns those deterministic checks. Existing economics inputs, outputs and earlier research snapshots are preserved.

[Alternative-energy comparisons](alternative-energy.md) retain general solar, deuterium-separation and construction-energy scales for future economic comparisons.
