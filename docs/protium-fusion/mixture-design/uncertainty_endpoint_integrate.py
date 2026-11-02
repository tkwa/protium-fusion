"""Conditional matched endpoint comparisons: carbon corners or B1 source controls.

No stochastic joint confidence interval is inferred. The existing receiver,
material-cycle and constrained plant models are held fixed and checked.
"""
from pathlib import Path
import argparse,json,math,sys
from rate_uncertainty_integrate import sha,read,write,module

RKEYS=['C12_p_to_N13_reaclib','C13_p_to_N14_reaclib']
def curve_pair(record):
 q=record['rate_update']['quantile']
 return tuple(q[k] for k in RKEYS) if isinstance(q,dict) else (q,q)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--base-reduced',type=Path,required=True);ap.add_argument('--runs',type=Path,required=True);ap.add_argument('--family',choices=['carbon','b1source'],required=True);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--replay',type=Path);a=ap.parse_args()
 r=a.research_root.resolve();o=a.output_dir.resolve();o.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(r/'component-design'));import material_cycle as mc
 receiver=module('corner_receiver',r/'mixture-design/receiver.py');pc=module('corner_plant',r/'mixture-design/plant_compare.py');vi=module('corner_checks',r/'mixture-design/variant_integrate.py')
 paths=[Path(__file__),r/'mixture-design/rate_uncertainty_integrate.py',r/'mixture-design/receiver.py',r/'mixture-design/plant_compare.py',r/'mixture-design/variant_integrate.py',r/'component-design/material_cycle.py',r/'cold-sink-redesign/station.py',r/'cold-sink-redesign/cases.json',r/'open-loop-redesign/operating-cases.json'];codes={str(p):sha(p) for p in paths};checks=[]
 base=read(a.base_reduced)['cases'];common={s:base[k] for s,k in [('A','A-H3'),('B1','B1-H2'),('B2','B2-H2'),('C13','C13-H4')]}
 if a.replay:archive=read(a.replay)['sensitivity_cases']
 else:
  if a.family=='carbon':
   sources=[]
   for folder in ['carbon-bands128','carbon-corners128']:
    for p in sorted((a.runs/folder).glob('*-rates-updated-*.json')):
     if any(x in p.name for x in ['.inputs.','.initial.','.pretrace.']):continue
     d=read(p);q=curve_pair(d);sources.append((d['stage']+'_'+q[0]+'_'+q[1],p))
   assert len(sources)==8,[(x,str(y)) for x,y in sources]
  else:
   sources=[('B1-H2-source-carbon',a.runs/'source128/B1-H2-source128-carbon-median.json'),('B1-H2-source-low',a.runs/'b1-source-robustness/h2low/B1-H2-source128-carbon-median-rates-carbon-n15-low-updated-median.json'),('B1-H1.5-source-carbon',a.runs/'b1-source-robustness/B1-H1p5-source128-carbon.json'),('B1-H1.5-source-low',a.runs/'b1-source-robustness/B1-H1p5-source128-low.json')]
  archive={}
  for label,p in sources:
   d=read(p);assert d['zones']==128 and d['success'] and d['termination_reason']=='final_time',str(p)
   duration=d['mixture_inputs'].get('duration_s',d['mixture_inputs'].get('physical_duration_s'));assert duration is not None and math.isclose(d['elapsed_s'],duration,rel_tol=1e-12)
   rr=receiver.reduce_case(p,r.parents[1],radius=1e9,thermal_points=65536);assert rr['usable_ballistic'] and rr['completed_requested_duration'];st=mc.load_stage(p)
   archive[label]={'name':p.stem,'path':str(p.resolve()),'record':d,'receiver':rr,'stage_data':st,'label':label,'source_sha256':{str(q.resolve()):sha(q) for q in [p,p.with_suffix('.initial.json'),p.with_suffix('.inputs.json'),p.with_suffix('.npz')]}}
   write(o/'receivers'/(p.stem+'.receiver.json'),rr);checks.append({'name':label+' completed endpoint and receiver','passed':True});print(json.dumps({'reduced':label,'gain':d['deposited_nuclear_over_initial_energy']}),flush=True)
 reduced={'scope':__doc__,'family':a.family,'base_reduced':str(a.base_reduced.resolve()),'base_sha256':sha(a.base_reduced),'common_cases':common,'sensitivity_cases':archive,'code_sha256':codes};write(o/(a.family+'-reduced-inputs.json'),reduced)
 variants=[('nominal128',[common[s] for s in ['A','B1','B2','C13']])]
 if a.family=='carbon':
  for q12,q13 in [('minus1sigma','minus1sigma'),('plus1sigma','plus1sigma'),('minus1sigma','plus1sigma'),('plus1sigma','minus1sigma')]:
   selected=[common['A'],common['B1'],archive['B2_'+q12+'_'+q13],archive['C13_'+q12+'_'+q13]];variants.append((q12+'_'+q13,selected))
 else:
  for label,entry in archive.items():variants.append((label,[common['A'],entry,common['B2'],common['C13']]))
 rows=[]
 for sigma in [5000.,2000.]:
  for label,selected in variants:
   cy=pc.cycle_for(selected);pl=pc.best_rate(selected,cy,absorber_sigma=sigma);ch=vi.result_checks(selected,cy,pl)
   for v in ch:v['passed']=bool(v['passed'])
   assert all(v['passed'] for v in ch);checks.extend([dict(v,case=label,sigma=sigma) for v in ch]);assert pl['feasible'],(label,pl['reasons'])
   summary={'case':label,'absorber_sigma_kg_m2':sigma,'mW_kg':1000*pl['net_W_kg'],'net_W':pl['net_W'],'mass_kg':pl['mass_kg'],'helium_per_s':pl['helium_per_s'],'main_window_occupancy':pl['sequential_occupancy_mass999'],'radiator_fraction':pl['radiator_area_fraction'],'prepared_MeV_per_He':cy['U_prepared_MeV_per_He'],'fresh_C_H':cy['fresh_C_H_mass_ratio'],'gross_fuel_kg_s':pl['gross_fuel_kg_s'],'pulse_interval_s':pl['pulse_interval_s'],'stage_attempts_per_He':cy['stage_attempts_per_He']}
   rows.append({'summary':summary,'selected_records':[s['path'] for s in selected],'cycle':cy,'plant':pl,'checks':ch});print(json.dumps(summary),flush=True)
  nominal=next(x['summary']['mW_kg'] for x in rows if x['summary']['case']=='nominal128' and x['summary']['absorber_sigma_kg_m2']==sigma)
  for row in rows:
   if row['summary']['absorber_sigma_kg_m2']==sigma:row['summary']['relative_to_nominal128']=row['summary']['mW_kg']/nominal-1
 for p,h in codes.items():assert sha(p)==h,('Implementation changed',p)
 qualification=('A3 and B1H2 held at nominal carbon-median128 endpoints. B2H2 and C13H4 use matching paired carbon quantiles, with all other physics/source assumptions fixed. These discrete corners do not establish a joint posterior interval or bound all interior values.' if a.family=='carbon' else 'A3, B2H2 and C13H4 held at nominal carbon-median128 endpoints. Only B1 source profile/composition and B1-only N15 low evaluation curve vary. The concrete source profiles retain their own prepared-energy and launch budgets; plant/source differences cannot be inferred from pulse gain alone. No probabilistic confidence attached.')
 out={'scope':__doc__,'family':a.family,'rows':rows,'checks':{'count':len(checks),'passed':sum(bool(c['passed']) for c in checks),'all_passed':all(c['passed'] for c in checks)},'code_sha256':codes,'base_reduced_sha256':sha(a.base_reduced),'qualification':qualification};write(o/(a.family+'-integrated.json'),out);print(json.dumps({'checks':out['checks']}),flush=True)
if __name__=='__main__':main()
