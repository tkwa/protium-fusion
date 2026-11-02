# ✨ Helium hold requirements under the dispatched thermal histories

The actual B1 temperature histories make rapid local collection a substantial alternative to waiting for a later refurbishment cycle. At the abstract-reported Jung mobile/orientation diffusion coefficient, an unchannelled mobile helium population reaches the cavity-facing surface quickly. A separate trapped population can survive longer, but neither its abundance nor its reactor-temperature lifetime has been measured for the proposed receiver.

The new calculation is [thermal_trapping.py](thermal_trapping.py); [thermal-trapping.json](thermal-trapping.json) records the exact source hashes, 177 numerical checks, all seven representative event histories, and both thermal models. The earlier helium evidence and spectrum packets are unchanged.

## Actual histories and the hold criterion

The calculation uses the geometry agent's chronological 350 kg/m² receiver histories, the actual helium packet numbers and arrival times, and the stopping packet's mean projected depths. It asks what single first-order trap barrier would retain **95% of the number-weighted helium cohort** until a repair deadline:

\[
S_j(E)=\exp\left[-\int_{t_{\rm implant,j}}^{t_{\rm repair,j}}\nu\exp\left(-\frac{E}{kT(t)}\right)dt\right],\qquad\sum_jw_jS_j=0.95.
\]

Attempt frequencies of \(10^{12}\) and \(10^{13}\,\mathrm{s^{-1}}\) are explicit kinetic scenarios. The barrier is a required trap property. It is not an energy charge per atom, a measured barrier distribution, or a prediction of the fraction of helium initially trapped. Retrapping, irradiation detrapping, bubbles and connected pores require additional equations.

The first 0.25 s of each temperature record is replaced by the separately resolved prompt upper curve. It is not added a second time. Newly implanted helium arrives after its own prompt; helium remaining from previous shots experiences that prompt. A 0.1 K allowance covers the archived omitted old thermal tails. The front temperature is used for the micron implantation region, with the bulk-conductivity and 100 µm damaged-skin histories retained as separate scenarios.

| B1 history, 100 µm low-conductivity skin | Hold after local implantation | Required barrier for 95% survival, ν=10¹²–10¹³/s |
|---|---:|---:|
| Ordinary midyear | 100 s | 3.26–3.47 eV |
| Ordinary midyear | 1,000 s | 3.47–3.69 eV |
| Warmest selected background | 100 s | 3.61–3.85 eV |
| Warmest selected background | 1,000 s | 3.85–4.08 eV |

Complete histories cover 99.9925% and 99.9833% of B1 helium nuclei for the 100 s and 1,000 s deadlines respectively. These tabulated values refer to those covered cohorts. Including every arrival while bounding the missing durations at 1,000–1,668 K changes the warm-background requirements by less than 0.0004 eV. The bounds are conditional thermal scenarios, not credible intervals.

A 4,000 s deadline **after implantation** extends beyond the archived record, which ends 4,000 s after prompt arrival. No cohort is completely covered. For the warm B1 history, the all-cohort barrier lies between 3.98 and 5.18 eV at ν=10¹²/s, or 4.22 and 5.51 eV at ν=10¹³/s, when the missing interval is held at 1,000 or 1,668 K. A narrower value requires a longer chronology with every subsequent prompt represented.

Helium already trapped before the hottest archived C13 event has a complete 4,000 s record. Its corresponding requirement is **4.39–4.68 eV** for the low-conductivity skin. At those barriers the first 0.25 s contributes only about 4–6×10⁻⁷ of the release hazard; the later material-heating excursion dominates this example. A microsecond peak held for thousands of seconds would give a misleading requirement.

## Mobile and trapped helium remain separate populations

The mobile calculation uses \(F=\operatorname{erfc}[x/(2\sqrt{\int Ddt})]\) in an absorbing-front half-space. Two coefficient policies are retained: extrapolating the Jung Arrhenius reference above its reported range, or holding it at its 1,050 K value. Across these policies and the ordinary/warm B1 histories, **98.9–99.5% of the covered mobile population reaches the front within 100 s**; about **99.67–99.84%** does so within 1,000 s. These are mobile-population results. The total receiver release requires the unknown initial mobile/trapped partition.

For a two-population illustration, the fraction still in an unchannelled receiver is

\[
f_{\rm trapped}S_{\rm trapped}+(1-f_{\rm trapped})(1-F_{\rm mobile}).
\]

Setting a trap barrier that gives 95% survival does not establish 95% total retention when a mobile population is also present. Conversely, release from the lattice can feed a deliberately placed collecting channel; it need not become a permanent inventory loss. The channel model must provide the collecting boundary and the gas model must carry away the released atoms.

The projected depths here are fixed-energy means, with no range-straggling distribution. For the half-space erfc kernel the convex depth dependence makes the mean-depth approximation optimistic for front retention at a fixed energy. The calculation contains no prompt reflection or implantation rejection. The kinetic reference and its limited material interpretation are documented in [the evidence review](README.md), including the accessible full Li et al. article and the abstract-only Jung evidence.

## Audit of the associated accounting

Read-only inspection found no material normalization error in `refurbishment/calculate.py` or the stopping energy partitions. Refurbishment keeps NRT event proxies separate from surviving defects and actual input work, labels its storage inventory as a mean-rate scenario, and separately counts the selected dispatch events. The named service and auxiliary allocations agree with the selected sigma350/EM20 component row to rounding accuracy. The stopping calculation separates total elastic transfer, a conventional NRT proxy, and actual retention/damage qualification.

The channel-routing finite-volume equations also match their declared boundaries: absorbing front, ideal absorbing side collectors below a continuous cap, reflecting rear and lateral symmetry. The time equation has the correct cell-volume source, and all current mean projected depths lie inside the 150 µm modeled layer. Its ideal side collector is a boundary condition to be supplied by a physical channel/plenum/pump design. It is not evidence for a selected diffusion anisotropy, finite channel capture probability, damage tolerance or low-backpressure operation.
