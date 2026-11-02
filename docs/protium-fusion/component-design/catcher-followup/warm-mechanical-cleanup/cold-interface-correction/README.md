# ✨ Passive cold-panel correction to the hot buffer

The retained hydrogen and helium panels support the specified gas throughput at a higher cold-interface pressure than the original prescribed total of 0.0001 Pa. Closing their finite pumping speed changes the 0.5 m³/m² hot-buffer peak from **1.69060 to 1.69143 Pa**. The source-plenum change stays below (8.6\times10^{-10}) Pa because the warm booster isolates it with a declared zero-flow compression ratio of (10^6).

This supplement preserves [the original buffer](../BUFFER.md) and replays its actual 4,228-event dispatch, 50 s source bins and 1,024 arrival quantiles. It retains the same physical conduit, initial mean-flow precharge, hot-source reservoir and warm booster. It changes only the downstream pressure boundary. The files contain 39 replay, flux, inventory and limiting-case checks.

## Passive boundary

For each species, write (Q_h) as the throughput in Pa·m³/s per receiver m², referenced to (T_h=1000) K. At the 300 K cold interface,

\[
P_c=P_{eq}+\frac{T_c}{T_h}\frac{Q_h}{S_{panel}},\qquad
Q_h=S_c\left[P_b-\sqrt{T_h/T_c}P_c\right].
\]

Eliminating (P_c) gives

\[
S_{eff}=\frac{S_c}{1+S_c\sqrt{T_c/T_h}/S_{panel}},\qquad
Q_h=S_{eff}\left[P_b-\sqrt{T_h/T_c}P_{eq}\right].
\]

This effective speed replaces the conduit speed in the two-reservoir time-dependent solver. The cold partial pressure is then reconstructed from the actual outlet flow at every time step. No mean pressure is substituted for its transient history.

The retained active areas are (6.2483931\times10^{-5}) m²/m² for hydrogen and (4.8666936\times10^{-5}) for helium. Their reference speeds are 10 m³/(s·m² active area) at 300 K and 13 at 273.15 K, respectively, with the existing conditional square-root temperature scaling.

| Quantity | H₂ | He |
|---|---:|---:|
| Panel speed at 300 K, m³/(s·m² receiver) | 0.0006248393 | 0.0006630364 |
| Original hot conduit speed, m³/(s·m²) | 1.4065646e−6 | 9.9815038e−7 |
| Combined hot speed, m³/(s·m²) | 1.4048325e−6 | 9.9732803e−7 |
| Speed reduction | 0.123145% | 0.0823875% |
| Cold pressure, finite-record mean, Pa | 0.0004233318 | 0.0001705332 |
| Cold pressure, individual sampled peak, Pa | 0.0008869423 | 0.0001717207 |

The individual pressure peaks need not occur simultaneously. Their combined sampled peak is **0.001056819 Pa**. The combined exact finite-record mean is **0.000593865 Pa**. These are reservoir-equivalent pressures at the stated interface; the helium adsorbent retains its separate 5 K thermal design.

## Corrected buffer result

| Quantity | Frozen prescribed-pressure reference | Passive panels, (P_{eq}=0) |
|---|---:|---:|
| Total minimum buffer pressure, Pa | 0.67935787 | 0.68006624 |
| Total maximum buffer pressure, Pa | 1.69060054 | 1.69143114 |
| Total exact mean buffer pressure, Pa | 1.00463742 | 1.00553855 |
| Peak outlet / continuous mean throughput | 1.77591364 | 1.77488411 |
| Initial hot-buffer gas mass, kg/m² | 1.6636095e−7 | 1.6650021e−7 |

With zero equilibrium pressure, the reverse molecular flux fraction is constant during the passive transient: 0.123145% for H₂ and 0.0823875% for He. This follows from (1-S_{eff}/S_c). Evaluating a fixed mean cold pressure at the smallest hot pressure would instead describe an independently regulated boundary; it would overstate this passive model's backstreaming at that time.

A helium equilibrium pressure of (10^{-6}) Pa shifts its entire hot-buffer curve upward by (1.825742\times10^{-6}) Pa and its source-plenum curve by (1.825742\times10^{-12}) Pa. The precharged reference flow history is unchanged. The maximum He reverse-flux fraction then becomes 0.082879%. `results.json` also records the distinct regulated-pressure comparison.

## Heat, apparatus and scope

The shared conduit still receives 0.01377056 W/m² of admitted hot radiation and 0.000342768 W/m² of solid conduction, totaling **0.01411333 W/m²**. Neither its geometry nor the active cold-panel area changes. Refrigeration, two-bank cycling, chemical separation and panel working capacity remain separately charged in the gas handoff. The same primary gas cohort is counted once.

The panel boundary here is quasi-steady and has negligible additional pneumatic storage. An extra manifold conductance of at least ten times the panel speed limits its added pressure rise to 10% of the panel-side rise. Any substantial cold-manifold volume, valve lag, loaded-surface capture change or sorber breakthrough needs its own response. These are component requirements; the finite conductance calculation supplies no measured hot-pump or irradiated-sorber qualification.

Pressure extrema are sampled exact bin endpoints. The archived 25/50/100 s controls already measure the source-discretization sensitivity; this supplement does not repeat that sweep. The separate global-cleanup capacity envelope, N/O chemistry and deferred regeneration cohorts retain their archived interfaces.

## Atom-lineage closure

This response sends 100% of arriving H/He directly into the local buffer to screen its capacity. The original buffer packet separately reports a global-cleanup increment of about 0.1% of that mean flow. Adding those two curves describes an independent capacity envelope; it cannot describe additional recovered H/He from the same incident atoms.

For a physical coupled history, let (I) be incident atom flow, (fI=J) the portion entering the cavity, (L=(1-f)I) the directly collected portion, (C) the later cavity capture and (O) the final outlet. Then

\[
\dot G=J-C,\qquad \dot B=L+C-O,\qquad
\dot G+\dot B=I-O.
\]

If solid trapping delays local release, include the held inventory (S): (\dot S=(1-f)I-L), and the exact total is (\dot G+\dot B+\dot S=I-O). The long-term recovered H/He flow therefore remains the incident flow, adjusted for changing stored inventory and any explicitly identified external supply or escape. Recapture itself creates no new feed requirement. These identities preserve the finite-record initial/final inventory convention used by the buffer solver; no changed gas-routing history is claimed in this supplement.

Run `calculate.py` with the shared scientific Python runtime. `results.json` records source and code hashes; `curves.json` contains downsampled pressure and flow traces. The exact inventory integrals include the initial precharge and final stored gas.
