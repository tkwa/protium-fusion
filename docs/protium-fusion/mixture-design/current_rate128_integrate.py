"""Explicit median-carbon 128-zone integration, with two receiver mass scenarios.

Only declared stage records are loaded. Raw extraction and reduced-input replay
are separate. Uses unchanged actual-composition material, receiver and plant APIs.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,itertools,json,math,sys

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def module(name,p):
 sp=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def run():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--research-root',required=True,type=Path);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--manifest',type=Path);g.add_argument('--replay',type=Path);ap.add_argument('--output-dir',required=True,type=Path);a=ap.parse_args();r=a.research_root.resolve();a.output_dir.mkdir(parents=True,exist_ok=True)
 sys.path.insert(0,str(r/'component-design'));import material_cycle as mc
 rec=module('current_rate128_receiver',r/'mixture-design/receiver.py');pc=module('current_rate128_plant',r/'mixture-design/plant_compare.py');vi=module('current_rate128_checks',r/'mixture-design/variant_integrate.py')
 codepaths=[Path(__file__),r/'mixture-design/receiver.py',r/'mixture-design/plant_compare.py',r/'mixture-design/variant_integrate.py',r/'component-design/material_cycle.py',r/'cold-sink-redesign/station.py',r/'cold-sink-redesign/cases.json',r/'open-loop-redesign/operating-cases.json'];codehash={str(p):sha(p) for p in codepaths}
 checks=[{'name':'receiver analytical self-checks','passed':rec.self_check()['passed']}]
 if a.replay:
  old=read(a.replay);archive=old['cases'];source=old['source_manifest'];checks.append({'name':'replay explicit archived reduced inputs','passed':True})
 else:
  source=read(a.manifest);archive={};tables=set();kernels=set()
  for label,declared in source['cases'].items():
   p=Path(declared['record']);d=read(p);inp=read(p.with_suffix('.inputs.json'));i=d['mixture_inputs'];u=i['rate_update'];duration=i['duration_s']
   assert d['stage']==declared['stage'] and d['zones']==128 and i['hydrogen_per_catalyst']==declared['hydrogen_per_catalyst']
   assert d['success'] and d['termination_reason']=='final_time' and math.isclose(d['elapsed_s'],duration,rel_tol=1e-10,abs_tol=1e-15)
   assert d['viscosity_mode']=='compression_cap' and d['opacity_m2_kg']==.001 and i['constraint']=='fixed_energy'
   assert i.get('temperature_factor',1)==1 and i.get('prepared_energy_factor',1)==1 and u['mode']=='updated' and u['quantile']=='median'
   assert math.isclose(duration,pc.BASE_DURATION[d['stage']],rel_tol=1e-10,abs_tol=1e-15)
   assert inp['initial_Y']==d['initial_Y']==i['initial_Y'];tables.add(u['table_sha256']);kernels.add(u['kernel_sha256'])
   for path,h in declared['source_sha256'].items():assert sha(path)==h,(label,path)
   receiver=rec.reduce_case(p,r.parents[1],radius=1e9,thermal_points=65536)
   assert receiver['completed_requested_duration'] and receiver['usable_ballistic'],label
   stage=mc.load_stage(p)
   assert max(abs(stage['initial_Y'].get(n,0)-d['initial_Y'].get(n,0)) for n in set(stage['initial_Y'])|set(d['initial_Y']))<1e-12
   entry={'name':p.stem,'path':str(p),'record':d,'receiver':receiver,'stage_data':stage,'label':label}
   archive[label]=entry;write(a.output_dir/'receivers'/(p.stem+'.receiver.json'),receiver)
   checks.extend([{'name':label+' exact family, rates and physical endpoint','passed':True},{'name':label+' actual-composition source and receiver reconstruction','passed':True}])
   print(json.dumps({'case':label,'receiver':'complete','gain':d['deposited_nuclear_over_initial_energy'],'warnings':receiver['warnings']}),flush=True)
  assert len(tables)==len(kernels)==1
  checks.append({'name':'one common median rate table and kernel','passed':True})
 write(a.output_dir/'current-rate128-reduced-inputs.json',{'schema':1,'scope':__doc__,'source_manifest':source,'cases':archive,'code_sha256':codehash,'checks':checks})
 families={stage:[v for v in archive.values() if v['record']['stage']==stage] for stage in pc.CAT}
 assert {k:len(v) for k,v in families.items()}=={'A':2,'B1':3,'B2':2,'C13':3}
 for sigma in [5000.,2000.]:
  rows=[];rejected=[];casechecks=[]
  for selected in itertools.product(*families.values()):
   cy=pc.cycle_for(selected);p=pc.best_rate(selected,cy,absorber_sigma=sigma);labels={e['record']['stage']:e['label'] for e in selected}
   ch=vi.result_checks(selected,cy,p);assert all(c['passed'] for c in ch),(labels,ch);casechecks.extend([dict(x,combination=labels) for x in ch])
   if not p.get('feasible'):rejected.append({'cases':labels,'reasons':p.get('reasons')});continue
   assert p['receiver_areal_mass_kg_m2']==sigma
   assert math.isclose(p['mass_breakdown_kg']['receiver'],4*math.pi*p['radius_m']**2*sigma,rel_tol=1e-12)
   summary={'cases':labels,'H_per_catalyst':{e['record']['stage']:e['record']['mixture_inputs']['hydrogen_per_catalyst'] for e in selected},'mW_kg':p['net_W_kg']*1000,'net_W':p['net_W'],'mass_kg':p['mass_kg'],'helium_per_s':p['helium_per_s'],'occupancy':p['sequential_occupancy_mass999'],'radiator_fraction':p['radiator_area_fraction'],'prepared_MeV_He':cy['U_prepared_MeV_per_He'],'fresh_C_H':cy['fresh_C_H_mass_ratio'],'fresh_C_kg_s':p['fresh_C_kg_s'],'fresh_H_kg_s':p['fresh_H_kg_s'],'H_return_kg_s':p['recovered_H_kg_s'],'CNO_return_kg_s':p['recovered_CNO_kg_s'],'gross_fuel_kg_s':p['gross_fuel_kg_s'],'pulse_interval_s':p['pulse_interval_s'],'H_processing_J_kg_for_10pct_export':p['H_processing_J_kg_for_10pct_export'],'surface_proxy_rise_K':p['max_material_plus_surface_prompt_rise_K']}
   rows.append({'summary':summary,'plant':p,'cycle':cy})
  rows.sort(key=lambda x:x['summary']['mW_kg'],reverse=True)
  best=rows[0];controlled={}
  for s in ['B1','B2']:
   comparisons=[]
   for row in rows:
    same=all(row['summary']['cases'][k]==best['summary']['cases'][k] for k in pc.CAT if k!=s)
    if same:comparisons.append(dict(row['summary'],relative_specific_power_to_best=row['summary']['mW_kg']/best['summary']['mW_kg']-1))
   controlled[s]=comparisons
  best_given={s:{str(h):next(row['summary'] for row in rows if row['summary']['H_per_catalyst'][s]==h) for h in sorted({row['summary']['H_per_catalyst'][s] for row in rows})} for s in ['B1','B2']}
  allchecks=checks+casechecks
  for check in allchecks:check['passed']=bool(check['passed'])
  for path,h in codehash.items():assert sha(path)==h,('implementation changed',path)
  output={'schema':1,'family':'median-carbon-128-prescribed-hot-flow','scope':__doc__,'comparison':{'receiver_areal_mass_kg_m2':sigma,'radius_m':1.1e9,'zones':128,'radiator_K':200,'heat_conversion':.6,'driver_efficiency':.6,'eventual_H_and_CNO_return':1,'energy_capture_fraction':.99,'maximum_main_window_occupancy':.8,'maximum_radiator_fraction':.9,'nuclear_rate_mode':'updated median carbon capture pairs; other190 original rates','source_preparation':'prescribed hot profiles; no source-mapped variant'},'ranked':[v['summary'] for v in rows],'top_detailed':rows[:10],'controlled_single_stage_comparisons':controlled,'best_conditional_on_mixture':best_given,'rejected':rejected,'checks':{'passed':int(sum(c['passed'] for c in allchecks)),'total':len(allchecks),'details':allchecks},'source_manifest':source,'code_sha256':codehash,'qualification':'Conditional integrated accounting; grey-opacity spherical source, collector and hardware assumptions retain their separate limitations. Exact eventual material return is a closed-reservoir accounting case, not perfect first-pass extraction. Omitted .1% tails remain retained and routed; main-window occupancy is not complete material clearance.'}
  write(a.output_dir/f'current-rate128-sigma{int(sigma)}.json',output)
  print(json.dumps({'sigma':sigma,'feasible':len(rows),'rejected':len(rejected),'best':best['summary'],'checks':output['checks']['passed']},indent=2),flush=True)
if __name__=='__main__':run()
