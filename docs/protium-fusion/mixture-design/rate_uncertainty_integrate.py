"""Integrate matched B1 N15 curves with fixed median-carbon128 A3/B2H2/C13H4.

Conditional evaluation sensitivities, not probabilistic yield intervals. Reuses
unchanged receiver/material/cycle/plant APIs; caches complete reduced inputs.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,sys

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def module(name,p):sp=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--base-reduced',type=Path,required=True);ap.add_argument('--runs',type=Path);ap.add_argument('--replay',type=Path);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();assert bool(a.runs)!=bool(a.replay)
 r=a.research_root.resolve();o=a.output_dir.resolve();o.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(r/'component-design'));import material_cycle as mc
 receiver=module('rate_unc_receiver',r/'mixture-design/receiver.py');pc=module('rate_unc_plant',r/'mixture-design/plant_compare.py');vi=module('rate_unc_checks',r/'mixture-design/variant_integrate.py')
 paths=[Path(__file__),r/'mixture-design/receiver.py',r/'mixture-design/plant_compare.py',r/'mixture-design/variant_integrate.py',r/'component-design/material_cycle.py',r/'cold-sink-redesign/station.py',r/'cold-sink-redesign/cases.json',r/'open-loop-redesign/operating-cases.json'];codes={str(p):sha(p) for p in paths}
 base=read(a.base_reduced)['cases'];common={s:base[k] for s,k in [('A','A-H3'),('B2','B2-H2'),('C13','C13-H4')]};checks=[{'name':'receiver self checks','passed':bool(receiver.self_check()['passed'])}]
 if a.replay:archive=read(a.replay)['sensitivity_cases']
 else:
  archive={}
  for p in sorted(a.runs.glob('*-rates-carbon-n15-*-updated-median.json')):
   d=read(p);u=d['rate_update'];curve=u['nitrogen15_evaluation']['curve'];h=d['mixture_inputs']['hydrogen_per_catalyst'];label=f'B1-H{h:g}-{curve}'
   assert d['stage']=='B1' and d['zones']==128 and h in [1.5,2.] and curve in ['adopted','low','high'];assert u['quantile']=='median'
   assert d['success'] and d['termination_reason']=='final_time' and math.isclose(d['elapsed_s'],d['mixture_inputs']['duration_s'],rel_tol=1e-12)
   rr=receiver.reduce_case(p,r.parents[1],radius=1e9,thermal_points=65536);assert rr['usable_ballistic'] and rr['completed_requested_duration'];st=mc.load_stage(p)
   archive[label]={'name':p.stem,'path':str(p.resolve()),'record':d,'receiver':rr,'stage_data':st,'label':label,'source_sha256':{str(q.resolve()):sha(q) for q in [p,p.with_suffix('.initial.json'),p.with_suffix('.inputs.json'),p.with_suffix('.npz')]}}
   write(o/'receivers'/(p.stem+'.receiver.json'),rr);checks.append({'name':label+' complete endpoint and receiver reconstruction','passed':True});print(json.dumps({'reduced':label,'gain':d['deposited_nuclear_over_initial_energy']}),flush=True)
  assert len(archive)==6
 reduced={'scope':__doc__,'base_reduced':str(a.base_reduced.resolve()),'base_sha256':sha(a.base_reduced),'common_cases':common,'nominal_B1':{k:base[k] for k in ['B1-H1.5','B1-H2']},'sensitivity_cases':archive,'code_sha256':codes};write(o/'n15-reduced-inputs.json',reduced)
 allrows=[]
 for sigma in [5000.,2000.]:
  group=[]
  variants=[('B1-H1.5-saved',base['B1-H1.5'],'saved'),('B1-H2-saved',base['B1-H2'],'saved')]+[(key,entry,entry['record']['rate_update']['nitrogen15_evaluation']['curve']) for key,entry in archive.items()]
  for label,b1,curve in variants:
   selected=[common['A'],b1,common['B2'],common['C13']];cy=pc.cycle_for(selected);pl=pc.best_rate(selected,cy,absorber_sigma=sigma);ch=vi.result_checks(selected,cy,pl)
   for v in ch:v['passed']=bool(v['passed'])
   assert all(v['passed'] for v in ch)
   checks.extend([dict(v,case=label,sigma=sigma) for v in ch]);assert pl['feasible'],(label,pl['reasons'])
   summary={'case':label,'N15_curve':curve,'H_per_N15':b1['record']['mixture_inputs']['hydrogen_per_catalyst'],'absorber_sigma_kg_m2':sigma,'mW_kg':1000*pl['net_W_kg'],'net_W':pl['net_W'],'mass_kg':pl['mass_kg'],'helium_per_s':pl['helium_per_s'],'main_window_occupancy':pl['sequential_occupancy_mass999'],'radiator_fraction':pl['radiator_area_fraction'],'prepared_MeV_per_He':cy['U_prepared_MeV_per_He'],'fresh_C_H':cy['fresh_C_H_mass_ratio'],'gross_fuel_kg_s':pl['gross_fuel_kg_s'],'pulse_interval_s':pl['pulse_interval_s'],'stage_attempts_per_He':cy['stage_attempts_per_He']}
   group.append({'summary':summary,'cycle':cy,'plant':pl,'checks':ch});print(json.dumps(summary),flush=True)
  baseline=next(v for v in group if v['summary']['case']=='B1-H2-saved')['summary']
  for row in group:row['summary']['relative_mWkg_to_saved_H2']=row['summary']['mW_kg']/baseline['mW_kg']-1
  allrows.extend(group)
 matched=[]
 for sigma in [5000.,2000.]:
  for curve in ['saved','adopted','low','high']:
   pair={x['summary']['H_per_N15']:x['summary'] for x in allrows if x['summary']['absorber_sigma_kg_m2']==sigma and x['summary']['N15_curve']==curve}
   matched.append({'absorber_sigma_kg_m2':sigma,'N15_curve':curve,'H1p5_mW_kg':pair[1.5]['mW_kg'],'H2_mW_kg':pair[2.]['mW_kg'],'H2_relative_advantage':pair[2.]['mW_kg']/pair[1.5]['mW_kg']-1,'preferred_of_two':1.5 if pair[1.5]['mW_kg']>pair[2.]['mW_kg'] else 2.})
 for p,h in codes.items():assert sha(p)==h,('Implementation changed',p)
 out={'scope':__doc__,'fixed_cases':{s:v['path'] for s,v in common.items()},'rows':allrows,'matched_B1_plant_comparisons':matched,'checks':{'count':len(checks),'passed':sum(bool(c['passed']) for c in checks),'all_passed':all(c['passed'] for c in checks)},'code_sha256':codes,'base_reduced_sha256':sha(a.base_reduced),'qualification':'Fixed source profiles and median128 A/B2/C13 endpoints; only B1 N15 curve and its actual receiver/material output varied. Other stages have negligible N15(p,alpha) flow in the saved diagnostic, but are explicitly held fixed here. Capacity is reoptimized within unchanged radiator/cadence constraints; source profiles are not retuned. Evaluation curves carry no assigned joint probability or80% yield interval.'};write(o/'n15-integrated.json',out);print(json.dumps({'matched':matched,'checks':out['checks']}),flush=True)
if __name__=='__main__':main()
