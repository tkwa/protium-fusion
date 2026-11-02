"""Current-rate mesh and composition refinement through unchanged plant accounting.

Only completed declared nominal hot-profile records are admitted. Receiver
reductions are cached against source/code hashes; pending records stay explicit.
"""
from pathlib import Path
import json,hashlib,importlib.util,sys,itertools,math
P=Path(__file__).resolve().parent;R=P.parent.parent;RAW=Path('/home/tkwa/src/protium-fusion/output/mixture-rate-runs')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def module(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 runner_hash=sha(__file__)
 sys.path.insert(0,str(R/'component-design'));import material_cycle as mc
 rec=module('final_refinement_receiver',R/'mixture-design/receiver.py');pc=module('final_refinement_plant',R/'mixture-design/plant_compare.py');vi=module('final_refinement_checks',R/'mixture-design/variant_integrate.py')
 paths=[R/'mixture-design/receiver.py',R/'mixture-design/plant_compare.py',R/'mixture-design/variant_integrate.py',R/'component-design/material_cycle.py',R/'cold-sink-redesign/station.py',R/'cold-sink-redesign/cases.json',R/'open-loop-redesign/operating-cases.json']
 codes={str(p):sha(p) for p in paths};basepath=R/'mixture-design/current-rate128-reduced-inputs.json';basepack=read(basepath);base=basepack['cases'];checks=[]
 def check(name,ok):
  checks.append({'name':name,'passed':bool(ok)})
  if not ok:raise AssertionError(name)
 for p,h in codes.items():check('baseline/current implementation: '+Path(p).name,basepack['code_sha256'][p]==h)
 check('receiver analytic identities',rec.self_check()['passed'])
 common={s:base[k] for s,k in [('A','A-H3'),('B1','B1-H2'),('B2','B2-H2'),('C13','C13-H4')]};archive={};pending=[];projection_exclusions=[]
 declarations={}
 for s,h in [('A',3),('B1',2),('B2',2),('C13',4)]:
  declarations[f'{s}-H{h}-z256']=RAW/'root256'/f'{s}-H{h}-fixed_energy-z256-k0p001-t1-rates-updated-median.json'
 for h in ['1p5','1p75','2','2p25']:
  declarations[f'B2-H{h}-z128']=RAW/'refine128'/f'B2-H{h}-fixed_energy-z128-k0p001-t1-rates-carbon-median.json'
 declarations['C13-H2-z128']=RAW/'source-aware128/C13-H2-fixed_energy-z128-k0p001-t1-rates-carbon-median.json'
 declarations['B2-H2-z256-tail2']=Path('/home/tkwa/src/protium-fusion/output/final-refinement/B2-H2-z256-tail2.json')
 for label,path in declarations.items():
  if not path.exists():pending.append({'case':label,'record':str(path)});continue
  d=read(path);inp=d['mixture_inputs'];u=d['rate_update'];duration=inp.get('duration_s',inp.get('physical_duration_s'));stage=d['stage']
  check(label+' completed duration',d['success'] and d['termination_reason']=='final_time' and math.isclose(d['elapsed_s'],duration,rel_tol=1e-12) and math.isclose(duration,pc.BASE_DURATION[stage]*(2 if label.endswith('tail2') else 1),rel_tol=1e-12))
  check(label+' fixed original physics',d['viscosity_mode']=='compression_cap' and d['opacity_m2_kg']==.001 and inp['constraint']=='fixed_energy' and inp.get('temperature_factor',1)==inp.get('prepared_energy_factor',1)==1)
  check(label+' same median carbon table',u['mode']=='updated' and u['quantile']=='median' and u['table_sha256']==common[stage]['record']['rate_update']['table_sha256'])
  source={str(q):sha(q) for q in [path,path.with_suffix('.inputs.json'),path.with_suffix('.initial.json'),path.with_suffix('.npz')]};cache=P/'cases'/(label+'.json')
  if cache.exists():
   e=read(cache);check(label+' cached source immutable',e['source_sha256']==source);check(label+' cached reduction unchanged',e['reduction_code_sha256']==codes)
  else:
   rr=rec.reduce_case(path,R.parents[1],radius=1e9,thermal_points=65536);st=mc.load_stage(path);e={'name':path.stem,'path':str(path),'record':d,'receiver':rr,'stage_data':st,'label':label,'source_sha256':source,'reduction_code_sha256':codes};write(cache,e)
  check(label+' receiver physical endpoint',e['receiver']['completed_requested_duration'] and e['receiver']['radius_m']==1e9)
  if not e['receiver']['usable_ballistic']:
   q=d['profile'];v=[0.]+q['velocity_m_s'];bad=[i for i in range(len(v)-1) if v[i+1]<=v[i]];projection_exclusions.append({'case':label,'reason':'Nonmonotone/outgoing velocity guard; excluded from plant projection','nonmonotone_intervals':bad,'affected_mass_fraction':sum(q['mass_fraction_of_total'][i] for i in bad),'warnings':e['receiver']['warnings']})
  check(label+' actual initial composition',e['stage_data']['initial_Y']==d['initial_Y'])
  archive[label]=e;print(json.dumps({'reduced':label,'gain':d['deposited_nuclear_over_initial_energy']}),flush=True)
 variants=[('median128-selected','mesh-reference',[common[s] for s in pc.CAT])]
 if all(f'{s}-H{h}-z256' in archive and archive[f'{s}-H{h}-z256']['receiver']['usable_ballistic'] for s,h in [('A',3),('B1',2),('B2',2),('C13',4)]):
  variants.append(('median256-selected','mesh-comparison',[archive[f'{s}-H{h}-z256'] for s,h in [('A',3),('B1',2),('B2',2),('C13',4)]]))
 tailkeys=[('A','A-H3-z256'),('B1','B1-H2-z256'),('B2','B2-H2-z256-tail2'),('C13','C13-H4-z256')]
 if all(k in archive and archive[k]['receiver']['usable_ballistic'] for _,k in tailkeys):
  variants.append(('median256-B2-late-outflow','mesh-with-late-outflow',[archive[k] for _,k in tailkeys]))
 b2={2:common['B2'],2.5:base['B2-H2.5']}
 for h,t in [(1.5,'1p5'),(1.75,'1p75'),(2.25,'2p25')]:
  if f'B2-H{t}-z128' in archive and archive[f'B2-H{t}-z128']['receiver']['usable_ballistic']:b2[h]=archive[f'B2-H{t}-z128']
 c13={1:base['C13-H1'],3:base['C13-H3'],4:common['C13']}
 if 'C13-H2-z128' in archive and archive['C13-H2-z128']['receiver']['usable_ballistic']:c13[2]=archive['C13-H2-z128']
 for h2,h13 in itertools.product(sorted(b2),sorted(c13)):
  variants.append((f'B2-H{h2:g}_C13-H{h13:g}','nominal128-composition',[common['A'],common['B1'],b2[h2],c13[h13]]))
 rows=[];rejected=[]
 for label,family,selected in variants:
  cy=pc.cycle_for(selected);pl=pc.best_rate(selected,cy,absorber_sigma=5000.);ch=vi.result_checks(selected,cy,pl)
  for v in ch:check(label+': '+v['name'],v['passed'])
  check(label+' capture .99',cy['energy_capture_fraction']==.99)
  if not pl['feasible']:rejected.append({'case':label,'family':family,'reasons':pl['reasons']});continue
  check(label+' receiver fixed',pl['receiver_areal_mass_kg_m2']==5000 and pl['radius_m']==1.1e9)
  summary={'case':label,'family':family,'H_per_catalyst':{e['record']['stage']:e['record']['mixture_inputs']['hydrogen_per_catalyst'] for e in selected},'mW_kg':pl['net_W_kg']*1000,'net_W':pl['net_W'],'mass_kg':pl['mass_kg'],'helium_per_s':pl['helium_per_s'],'main_window_occupancy':pl['sequential_occupancy_mass999'],'radiator_fraction':pl['radiator_area_fraction'],'prepared_MeV_per_He':cy['U_prepared_MeV_per_He'],'fresh_C_H':cy['fresh_C_H_mass_ratio'],'gross_fuel_kg_s':pl['gross_fuel_kg_s'],'pulse_interval_s':pl['pulse_interval_s'],'stage_attempts_per_He':cy['stage_attempts_per_He'],'surface_proxy_rise_K':pl['max_material_plus_surface_prompt_rise_K']}
  rows.append({'summary':summary,'cycle':cy,'plant':pl,'selected_records':[e['path'] for e in selected]});print(json.dumps(summary),flush=True)
 nominal=rows[0]['summary']['mW_kg']
 for row in rows:row['summary']['relative_to_median128_selected']=row['summary']['mW_kg']/nominal-1
 mesh=[]
 for s,h in [('A',3),('B1',2),('B2',2),('C13',4)]:
  k=f'{s}-H{h}-z256'
  if k not in archive:continue
  b=common[s]['record'];d=archive[k]['record'];mesh.append({'stage':s,'H_per_catalyst':h,'gain128':b['deposited_nuclear_over_initial_energy'],'gain256':d['deposited_nuclear_over_initial_energy'],'gain_relative_change':d['deposited_nuclear_over_initial_energy']/b['deposited_nuclear_over_initial_energy']-1,'target128':b['desired_product_per_initial_catalyst'],'target256':d['desired_product_per_initial_catalyst'],'energy_residual256':d['energy_residual_fraction_initial'],'record256':archive[k]['path']})
 if 'B2-H2-z128' in archive:check('direct-run B2H2 duplicates original median128 endpoint',archive['B2-H2-z128']['record']['deposited_nuclear_over_initial_energy']==common['B2']['record']['deposited_nuclear_over_initial_energy'])
 for path,h in codes.items():check('implementation unchanged during calculation '+Path(path).name,sha(path)==h)
 cycle_only=None
 meshkeys=[f'{s}-H{h}-z256' for s,h in [('A',3),('B1',2),('B2',2),('C13',4)]]
 if all(k in archive for k in meshkeys):
  cy=pc.cycle_for([archive[k] for k in meshkeys]);check('256 material-only baryon closure',abs(cy['baryon_residual_per_He'])<1e-7);check('256 material-only energy partition',abs(cy['energy_partition_residual_MeV_per_He'])<1e-7);cycle_only={'scope':'Matched original endpoint times, independent of receiver projection. This material/energy cycle is not a plant point when ballistic projection is rejected.','cycle':cy,'relative_prepared_MeV_per_He':cy['U_prepared_MeV_per_He']/rows[0]['cycle']['U_prepared_MeV_per_He']-1}
 check('runner unchanged during calculation',sha(__file__)==runner_hash)
 ranked=sorted([r['summary'] for r in rows if r['summary']['family']=='nominal128-composition'],key=lambda q:q['mW_kg'],reverse=True)
 result={'scope':__doc__,'comparison':{'receiver_sigma_kg_m2':5000.,'radius_m':1.1e9,'energy_capture_fraction':.99,'heat_conversion':.6,'driver_efficiency':.6,'maximum_main_window_occupancy':.8,'maximum_radiator_area_fraction':.9,'eventual_material_return':1,'profile_family':'prescribed hot inflow; no constructed source or transport correction'},'pending':pending,'receiver_projection_exclusions':projection_exclusions,'mesh_stage_comparison':mesh,'matched_time_256_material_cycle_only':cycle_only,'rows':rows,'ranked_nominal128_compositions':ranked,'rejected':rejected,'base_reduced_path':str(basepath),'base_reduced_sha256':sha(basepath),'new_cases':{k:str(P/'cases'/(k+'.json')) for k in archive},'code_sha256':codes,'runner_sha256':runner_hash,'checks':{'passed':sum(x['passed'] for x in checks),'total':len(checks),'details':checks},'qualification':'Conditional model comparison. The 256 mesh row tests the already selected quartet, not every possible composition at 256. A separately labeled 256 late-outflow row uses twice the original B2 end time if available, so combines mesh and endpoint-extension effects. Nominal compositions preserve prepared pulse energy and actual isotope endpoints. Construction, charge preparation, capsule additions, assembly response and receiver transport remain separate model layers; no confidence interval inferred from mesh or recipe spread.'}
 write(P/'results.json',result)
 print(json.dumps({'checks':result['checks']['passed'],'pending':pending,'best_nominal128':ranked[0]},indent=2),flush=True)
if __name__=='__main__':main()
