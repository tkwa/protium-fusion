# ✨ Whole-dispatch heat-spreading audit

The chosen 1250 K and 1200 K conductivity comparisons retain the same controlling events after a bounded audit of the full dispatch. **No hotter event was missed by the seven original representative histories.** The archived representative packet remains unchanged.

[all_event_audit.py](all_event_audit.py) resolves all **75 A/B2/C13 events** from 0.25 to 4000 s after prompt arrival at **0.25 s** spacing. It also evaluates the entire 370-day record on a 20 s grid. Outside every large-pulse window it constructs a positive upper bound from the sampled bulk mode, a sampling guard, maximum B1 response and its tails at the actual minimum spacing, maximum late tails of each large stage, and the retained 0.1 K historical-tail allowance.

| Effective k, W/(m·K) | Fine all-large-event maximum | Controller | Bound outside large-pulse windows |
|---:|---:|---|---:|
| 325.14298 | 1250.0000 K | C13 event index 75 | 1189.857 K |
| 933.64369 | 1200.0000 K | A event index 1623 | 1189.634 K |
| 1000 | 1198.8928 K | A event index 1623 | 1189.623 K |

Indices refer to the existing zero-based chronological dispatch array. The large-window maxima equal the representative values to numerical precision. Their original source-tail qualification is retained; no extra numerical digits establish an exact physical temperature ceiling.

The audit also checks the concern that the stored peak-detail index was selected before applying the material-phase time mask. The saved supplement lists the masked and unmasked peak times for every frozen table row, and separately for the three audited conductivities. Any discrepancy is exposed in the JSON rather than silently changing the original record.

This covers the resolved material/decay history and prompt tail. The separate microsecond photon-deposition proxy remains outside these comparison temperatures. [Full results](all-event-audit.json).
