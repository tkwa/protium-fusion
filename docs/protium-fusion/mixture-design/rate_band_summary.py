"""Matched conditional nuclear-rate comparisons, without probabilistic relabeling."""
from pathlib import Path
import argparse,json,hashlib
ap=argparse.ArgumentParser();ap.add_argument('--runs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
nom={}
for p in (a.runs/'nominal128').glob('*-rates-updated-median.json'):nom[p.name.split('-fixed')[0]]=json.loads(p.read_text())
rows=[]
for family in ['carbon-bands128','n15-bands128','carbon-corners128']:
 for p in sorted((a.runs/family).glob('*.json')):
  if any(t in p.name for t in ['.inputs.','.initial.','.pretrace.']) or p.name.endswith('status.json'):continue
  d=json.loads(p.read_text())
  if 'rate_update' not in d:continue
  key=p.name.split('-fixed')[0];m=d['rate_update'];q=m['quantile'];qs=q if isinstance(q,dict) else {'C12_p_to_N13_reaclib':q,'C13_p_to_N14_reaclib':q};n15=m.get('nitrogen15_evaluation',{}).get('curve','saved')
  rows.append({'case':key,'family':family,'record':str(p.resolve()),'record_sha256':sha(p),'trace_sha256':sha(p.with_suffix('.npz')),'carbon_curves':qs,'N15_curve':n15,'gain':d['deposited_nuclear_over_initial_energy'],'gain_change_from_carbon_median':d['deposited_nuclear_over_initial_energy']/nom[key]['deposited_nuclear_over_initial_energy']-1,'target_yield':d['desired_product_per_initial_catalyst'],'energy_residual_fraction_initial':d['energy_residual_fraction_initial'],'complete':d['success'] and d['termination_reason']=='final_time','output_mass_fractions':d['output_mass_fractions']})
matched=[]
for curve in ['adopted','low','high']:
 rr={r['case']:r for r in rows if r['family']=='n15-bands128' and r['N15_curve']==curve}
 if all(k in rr for k in ['B1-H1p5','B1-H2']):matched.append({'N15_curve':curve,'H1p5_gain':rr['B1-H1p5']['gain'],'H2_gain':rr['B1-H2']['gain'],'H2_gain_relative_advantage':rr['B1-H2']['gain']/rr['B1-H1p5']['gain']-1})
carbon=[]
for key in ['B2-H2','C13-H4']:
 rr=[r for r in rows if r['case']==key and r['family'] in ['carbon-bands128','carbon-corners128']]
 if rr:
  lo=min(rr,key=lambda r:r['gain']);hi=max(rr,key=lambda r:r['gain']);carbon.append({'case':key,'tested_corner_count':len(rr),'median_gain':nom[key]['deposited_nuclear_over_initial_energy'],'minimum_tested_gain':lo['gain'],'minimum_curves':lo['carbon_curves'],'maximum_tested_gain':hi['gain'],'maximum_curves':hi['carbon_curves']})
result={'scope':'Conditional source-curve simulations at fixed128 profiles. Simultaneous and opposite carbon quantiles probe discrete directions/corners; no joint posterior or guaranteed interior bound. N15 low/high are NACRE-II evaluation limits without assigned probability. All other nuclear rates and hydro/EOS/opacity choices remain fixed.','code_sha256':sha(__file__),'rows':rows,'matched_B1_comparisons':matched,'carbon_tested_corner_ranges':carbon}
a.out.write_text(json.dumps(result,indent=2)+'\n')
for r in rows:print(r['case'],r['carbon_curves'],r['N15_curve'],f"Q{r['gain']:.9f}",f"target{r['target_yield']:.8f}")
print('Matched B1:',matched)
