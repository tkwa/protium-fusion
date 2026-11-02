# ✨ Final nominal mixture and mesh refinement

The refined 128-cell full-cycle comparison retains **B2 H/C = 2** as the best of the tested choices. Nominal C13 H/C = 2 has only a **0.24%** station-power disadvantage relative to H/C = 4, making its fuel chemistry a consequential separate design tradeoff. These values are conditional engineering accounting outputs, not measured reactor performance or confidence intervals.

All rows hold radius **1.1 × 10⁹ m**, absorber **5000 kg/m²**, energy capture **0.99**, heat conversion **0.60**, and driver efficiency **0.60** fixed. The common model limits main pulse-window occupancy to 0.8 and radiator area fraction to 0.9. Actual initial compositions, complete endpoint isotope distributions, cold-decay/recycle routing, prepared energy, source peak and inventory masses are used. Eventual material return is the stated closed-reservoir assumption; unfinished per-pass recovery is not automatically a permanent loss.

## Refined B2 mixtures at fixed A H3 / B1 H2 / C13 H4

| B2 H/C | Net mW/kg | Relative to H2 | Prepared MeV per He | Main-window occupancy |
|---|---:|---:|---:|---:|
| 1.5 | 6.56797846 | -1.0937% | 1.964820 | 0.717766 |
| 1.75 | 6.62631744 | -0.2152% | 1.909150 | 0.720471 |
| 2 | 6.64060846 | +0.0000% | 1.895680 | 0.721236 |
| 2.25 | 6.60552316 | -0.5283% | 1.929457 | 0.719804 |
| 2.5 | 6.51796696 | -1.8468% | 2.014177 | 0.715929 |

All five rows are radiator-limited. The direct new-run B2 H2 endpoint is bit-identical in gain to the earlier median-carbon128 reference. The H1.5, H1.75 and H2.25 points therefore refine the same model family rather than changing the source or nuclear-rate convention. The full 5 × 4 B2/C13 grid is archived; H2 remains best within each tested C13 choice.

## Nominal C13 chemistry comparison at fixed A H3 / B1 H2 / B2 H2

| C13 H/C | Net mW/kg | Relative to H4 | Prepared MeV per He | Main-window occupancy |
|---|---:|---:|---:|---:|
| 1 | 6.59552988 | -0.6788% | 1.938287 | 0.719431 |
| 2 | 6.62486024 | -0.2372% | 1.910535 | 0.720599 |
| 3 | 6.63640907 | -0.0632% | 1.899638 | 0.721065 |
| 4 | 6.64060846 | +0.0000% | 1.895680 | 0.721236 |

These are prescribed hot-profile comparisons. Changing to a molecular or polymer feed, adding a capsule, and mapping actual preparation entropy require their own source runs. The small nominal power differences should be weighed against those construction burdens; this packet does not automatically debit or credit them.

## Matched original-endpoint mesh results

| Stage | Q at 128 cells | Q at 256 cells | Q change | Target fraction at 128 | Target fraction at 256 |
|---|---:|---:|---:|---:|---:|
| A H3 | 10.81779596 | 11.13140574 | +2.8990% | 0.80679102 | 0.83018003 |
| B1 H2 | 8.26190402 | 8.61106356 | +4.2261% | 0.82773875 | 0.86269288 |
| B2 H2 | 9.91762947 | 9.86709879 | -0.5095% | 0.76637728 | 0.76156353 |
| C13 H4 | 24.69857103 | 25.25187614 | +2.2402% | 0.18026500 | 0.17284679 |

Q is deposited nuclear energy divided by initial prepared energy. The target species are O15 for A, C12 for B1, O14 for B2, and N14 for C13. A stage's named target fraction does not summarize all useful routing: the full isotope endpoint is used in the material cycle. At the original end times, the selected 256-cell quartet needs **1.85379883 prepared MeV per He**, versus **1.89568044** at 128 cells (-2.2093%). This matched-time material/energy cycle can be evaluated independently of a ballistic receiver projection.

The original 256-cell B2 endpoint contains four decreasing-velocity intervals, cells 31–34, involving **0.0779867%** of pulse mass. All material is moving outward, but the velocity decreases from about 973 to 912 km/s over this region; neighboring-shell ballistic crossing estimates are about **31–164 μs** beyond the stored endpoint. The receiver guard correctly rejects that endpoint. Sorting velocities or silently borrowing the 128-cell receiver profile would conceal the unresolved hydrodynamic evolution.

## Later B2 outflow and the full 256-cell ledger

The authorized B2 rerun to twice the reference end time completed and passed the ballistic receiver guard. Its Q is 9.867098792; the reported change from the original endpoint is only 1.30e-11 relative. The new quartet uses original endpoint times for A/B1/C13 and this later B2 endpoint. It gives **6.67763075 mW/kg**, **+0.5575%** relative to the selected 128-cell quartet, with occupancy 0.692872 and radiator fraction 0.900000. This row combines mesh and B2 endpoint extension; it is not a matched-time mesh convergence estimate. Full mass, power, pulse cadence, isotope-flow, energy and receiver records are in `results.json`.

## Validation and provenance

The current packet has **263/263 passing checks**. Failed ballistic eligibility is retained explicitly in `receiver_projection_exclusions`; it is not counted as a qualified plant row. All consumed model/source hashes are saved. `integrate.py` reuses the unchanged receiver, material-cycle and constrained station APIs, with cache validity tied to exact source and implementation hashes. `run_tail.py` preserves the original B2 initial state and doubles only the requested physical duration. `report.py` generates this review from the saved results.

A and B1 gains still move about 3–4% between 128 and 256 cells, so asymptotic numerical convergence has not been established. The 256-cell family tests the already selected quartet, rather than every composition at 256 cells. The C13 H2 result in this packet is at 128 cells; the H4 mesh comparison does not establish H2-specific mesh convergence. A constructed H2 source calculation at 128 cells likewise needs its own source-specific mesh qualification. The pulse models retain spherical symmetry, simplified opacity/EOS and the specified network. Source construction, assembly, receiver transport, material lifetime and hardware qualification remain separate. Mesh spreads and recipe differences are not statistical confidence intervals.
