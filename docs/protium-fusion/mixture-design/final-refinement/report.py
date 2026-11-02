"""Render the compact refinement review from the validated integration packet."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent
def main():
 d=json.loads((P/'results.json').read_text());rows=d['rows'];base=next(r for r in rows if r['summary']['case']=='median128-selected');q0=base['summary']['mW_kg'];lines=[]
 for r in sorted((r['summary'] for r in rows if r['summary']['family']=='nominal128-composition' and r['summary']['H_per_catalyst']['C13']==4),key=lambda r:r['H_per_catalyst']['B2']):
  lines.append(f"| {r['H_per_catalyst']['B2']:g} | {r['mW_kg']:.8f} | {100*(r['mW_kg']/q0-1):+.4f}% | {r['prepared_MeV_per_He']:.6f} | {r['main_window_occupancy']:.6f} |")
 chem=[]
 for r in sorted((r['summary'] for r in rows if r['summary']['family']=='nominal128-composition' and r['summary']['H_per_catalyst']['B2']==2),key=lambda r:r['H_per_catalyst']['C13']):
  chem.append(f"| {r['H_per_catalyst']['C13']:g} | {r['mW_kg']:.8f} | {100*(r['mW_kg']/q0-1):+.4f}% | {r['prepared_MeV_per_He']:.6f} | {r['main_window_occupancy']:.6f} |")
 mesh=[]
 for r in d['mesh_stage_comparison']:
  mesh.append(f"| {r['stage']} H{r['H_per_catalyst']:g} | {r['gain128']:.8f} | {r['gain256']:.8f} | {100*r['gain_relative_change']:+.4f}% | {r['target128']:.8f} | {r['target256']:.8f} |")
 late=next((r for r in rows if r['summary']['family']=='mesh-with-late-outflow'),None)
 if late:
  r=late['summary'];e=json.loads((P/'cases/B2-H2-z256-tail2.json').read_text());old=json.loads((P/'cases/B2-H2-z256.json').read_text());x=e['record'];y=old['record'];late_text=f"The authorized B2 rerun to twice the reference end time completed and passed the ballistic receiver guard. Its Q is {x['deposited_nuclear_over_initial_energy']:.9f}; the reported change from the original endpoint is only {x['deposited_nuclear_over_initial_energy']/y['deposited_nuclear_over_initial_energy']-1:.2e} relative. The new quartet uses original endpoint times for A/B1/C13 and this later B2 endpoint. It gives **{r['mW_kg']:.8f} mW/kg**, **{100*(r['mW_kg']/q0-1):+.4f}%** relative to the selected 128-cell quartet, with occupancy {r['main_window_occupancy']:.6f} and radiator fraction {r['radiator_fraction']:.6f}. This row combines mesh and B2 endpoint extension; it is not a matched-time mesh convergence estimate. Full mass, power, pulse cadence, isotope-flow, energy and receiver records are in `results.json`."
 else:
  late_text="The authorized B2 rerun to twice the reference end time remains pending or unqualified. A full 256-cell plant point is withheld until its receiver projection passes; no 128-cell timing or sorted velocity profile has been substituted."
 cy=d['matched_time_256_material_cycle_only'];prepared=cy['cycle']['U_prepared_MeV_per_He']
 text='''# ✨ Final nominal mixture and mesh refinement

The refined 128-cell full-cycle comparison retains **B2 H/C = 2** as the best of the tested choices. Nominal C13 H/C = 2 has only a **0.24%** station-power disadvantage relative to H/C = 4, making its fuel chemistry a consequential separate design tradeoff. These values are conditional engineering accounting outputs, not measured reactor performance or confidence intervals.

All rows hold radius **1.1 × 10⁹ m**, absorber **5000 kg/m²**, energy capture **0.99**, heat conversion **0.60**, and driver efficiency **0.60** fixed. The common model limits main pulse-window occupancy to 0.8 and radiator area fraction to 0.9. Actual initial compositions, complete endpoint isotope distributions, cold-decay/recycle routing, prepared energy, source peak and inventory masses are used. Eventual material return is the stated closed-reservoir assumption; unfinished per-pass recovery is not automatically a permanent loss.

## Refined B2 mixtures at fixed A H3 / B1 H2 / C13 H4

| B2 H/C | Net mW/kg | Relative to H2 | Prepared MeV per He | Main-window occupancy |
|---|---:|---:|---:|---:|
'''+ '\n'.join(lines)+'''

All five rows are radiator-limited. The direct new-run B2 H2 endpoint is bit-identical in gain to the earlier median-carbon128 reference. The H1.5, H1.75 and H2.25 points therefore refine the same model family rather than changing the source or nuclear-rate convention. The full 5 × 4 B2/C13 grid is archived; H2 remains best within each tested C13 choice.

## Nominal C13 chemistry comparison at fixed A H3 / B1 H2 / B2 H2

| C13 H/C | Net mW/kg | Relative to H4 | Prepared MeV per He | Main-window occupancy |
|---|---:|---:|---:|---:|
'''+ '\n'.join(chem)+'''

These are prescribed hot-profile comparisons. Changing to a molecular or polymer feed, adding a capsule, and mapping actual preparation entropy require their own source runs. The small nominal power differences should be weighed against those construction burdens; this packet does not automatically debit or credit them.

## Matched original-endpoint mesh results

| Stage | Q at 128 cells | Q at 256 cells | Q change | Target fraction at 128 | Target fraction at 256 |
|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(mesh)+f'''

Q is deposited nuclear energy divided by initial prepared energy. The target species are O15 for A, C12 for B1, O14 for B2, and N14 for C13. A stage's named target fraction does not summarize all useful routing: the full isotope endpoint is used in the material cycle. At the original end times, the selected 256-cell quartet needs **{prepared:.8f} prepared MeV per He**, versus **{base['cycle']['U_prepared_MeV_per_He']:.8f}** at 128 cells ({100*cy['relative_prepared_MeV_per_He']:+.4f}%). This matched-time material/energy cycle can be evaluated independently of a ballistic receiver projection.

The original 256-cell B2 endpoint contains four decreasing-velocity intervals, cells 31–34, involving **0.0779867%** of pulse mass. All material is moving outward, but the velocity decreases from about 973 to 912 km/s over this region; neighboring-shell ballistic crossing estimates are about **31–164 μs** beyond the stored endpoint. The receiver guard correctly rejects that endpoint. Sorting velocities or silently borrowing the 128-cell receiver profile would conceal the unresolved hydrodynamic evolution.

## Later B2 outflow and the full 256-cell ledger

'''+late_text+f'''

## Validation and provenance

The current packet has **{d['checks']['passed']}/{d['checks']['total']} passing checks**. Failed ballistic eligibility is retained explicitly in `receiver_projection_exclusions`; it is not counted as a qualified plant row. All consumed model/source hashes are saved. `integrate.py` reuses the unchanged receiver, material-cycle and constrained station APIs, with cache validity tied to exact source and implementation hashes. `run_tail.py` preserves the original B2 initial state and doubles only the requested physical duration. `report.py` generates this review from the saved results.

A and B1 gains still move about 3–4% between 128 and 256 cells, so asymptotic numerical convergence has not been established. The 256-cell family tests the already selected quartet, rather than every composition at 256 cells. The C13 H2 result in this packet is at 128 cells; the H4 mesh comparison does not establish H2-specific mesh convergence. A constructed H2 source calculation at 128 cells likewise needs its own source-specific mesh qualification. The pulse models retain spherical symmetry, simplified opacity/EOS and the specified network. Source construction, assembly, receiver transport, material lifetime and hardware qualification remain separate. Mesh spreads and recipe differences are not statistical confidence intervals.
'''
 (P/'README.md').write_text(text)
 print('Wrote README.md; late-outflow row:',late is not None)
if __name__=='__main__':main()
