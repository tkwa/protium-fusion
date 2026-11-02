"""Archive and integrate the completed unchanged-source128 mesh sensitivity.
Keeps the source64 packet frozen and compares only coherent same-recipe families.
"""
from pathlib import Path
import json,sys,shutil,math
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'mixture-design'),str(ROOT/'component-design'),str(ROOT)]
import variant_integrate as vi
import pulse_diagnostics as pd
from fast_eos import FastEOS
import material_cycle as mc
sha=vi.sha
RAW=Path('/home/tkwa/src/protium-fusion/output/source128')
TAGS={'n14':('A','3'),'n15':('B1','2'),'c12':('B2','2'),'c13':('C13','4')}
def main():
 maps={'source128':{},'nominal128':{}};audits=[];eos=FastEOS()
 for cat,(tag,h) in TAGS.items():
  stem=RAW/f'{tag}-H{h}-source128-carbon-median';dest=P/'records'/stem.name;dest.parent.mkdir(parents=True,exist_ok=True)
  for ext in ['.json','.initial.json','.inputs.json']:shutil.copyfile(stem.with_suffix(ext),dest.with_suffix(ext))
  for family,rec in [('source128',dest.with_suffix('.json')),('nominal128',ROOT/f'mixture-design/rate-results/carbon128/{tag}-H{h}-fixed_energy-z128-k0p001-t1-rates-updated-median.json')]:
   d=vi.read(rec);inp=rec.with_suffix('.inputs.json');init=rec.with_suffix('.initial.json');i=vi.read(inp);tr=Path(d['trace_path']);duration=i.get('duration_s',i.get('physical_duration_s'))
   assert d['success'] and d['termination_reason']=='final_time' and d['zones']==128 and d['viscosity_mode']=='compression_cap' and abs(d['elapsed_s']/duration-1)<1e-12
   maps[family][cat]={'record':str(rec.relative_to(ROOT)),'initial':str(init.relative_to(ROOT)),'inputs':str(inp.relative_to(ROOT)),'trace':str(tr),'duration_s':duration,'initial_Y':i['initial_Y'],'eos':'FastEOS','result_sha256':sha(rec),'initial_profile_sha256':sha(init),'input_record_sha256':sha(inp),'raw_trace_sha256':sha(tr)}
   if family=='source128':
    errors=i['source_refinement']['remap_errors'];assert max(abs(x) for x in errors.values())<1e-10
    s,_=pd.reduce({'label':tag+' unchanged source128','family':'source128-carbon-median','result':str(rec),'trace':str(tr),'reaction_flows':False},eos,False);assert s['checks']['all_passed'];audits.append(s)
 for family,mp in maps.items():vi.write(P/(family+'-stage-map.json'),{'family':family+'-A3-B1H2-B2H2-C13H4-carbon-median','stage_map':mp,'scope':__doc__})
 vi.write(P/'trace-audit.json',{'cases':audits,'checks':{'count':sum(x['checks']['count'] for x in audits),'passed':sum(x['checks']['passed'] for x in audits),'all_passed':True},'reducer_sha256':sha(ROOT/'mixture-design/pulse_diagnostics.py')})
 receiver=vi.module('source128_receiver',ROOT/'mixture-design/receiver.py');plant=vi.module('source128_plant',ROOT/'mixture-design/plant_compare.py');receiver.self_check()
 codes=[Path(__file__),P/'refine_source.py',ROOT/'mixture-design/variant_integrate.py',ROOT/'mixture-design/plant_compare.py',ROOT/'mixture-design/receiver.py',ROOT/'component-design/material_cycle.py',ROOT/'cold-sink-redesign/station.py',ROOT/'open-loop-redesign/operating-cases.json',ROOT/'cold-sink-redesign/cases.json'];hashes={str(x):sha(x) for x in codes};checks=[];sources={};families={}
 for family in maps:
  selected,archive,cc,ss=vi.extract_explicit(P/(family+'-stage-map.json'),ROOT,RAW/'integration'/family,receiver,mc);sources.update(ss);checks+=cc;cy=plant.cycle_for(selected);pp={}
  for sig in [5000.,2000.]:
   out=plant.best_rate(selected,cy,radius=1.1e9,absorber_sigma=sig);checks+=vi.result_checks(selected,cy,out);pp[str(int(sig))]=out
  families[family]={'stage_inputs':archive,'cycle':cy,'capacity_plants':pp}
 oldpath=P.parent/'source-family-carbon-median/integration.json';old=vi.read(oldpath);sources[str(oldpath)]=sha(oldpath)
 def ratios(a,b,keys):return {k:a[k]/b[k]-1 for k in keys}
 compares={};stages={}
 for name,left,right in [('source128_vs_nominal128',families['source128'],families['nominal128']),('source128_vs_source64',families['source128'],old['families']['source64']),('nominal128_vs_nominal64',families['nominal128'],old['families']['nominal64'])]:
  compares[name]={'cycle_relative_changes':ratios(left['cycle'],right['cycle'],['U_prepared_MeV_per_He','H_nonnu_MeV_per_He','fresh_C_H_mass_ratio','gross_fuel_kg_per_He']),'capacity_relative_changes':{sig:ratios(left['capacity_plants'][sig],right['capacity_plants'][sig],['net_W','mass_kg','net_W_kg','helium_per_s','sequential_occupancy_mass999']) for sig in ['5000','2000']}}
 for cat in TAGS:
  n=families['nominal128']['stage_inputs'][cat]['selected']['record'];s=families['source128']['stage_inputs'][cat]['selected']['record'];o=old['families']['source64']['stage_inputs'][cat]['selected']['record']
  row={'source64_gain':o['deposited_nuclear_over_initial_energy'],'source128_gain':s['deposited_nuclear_over_initial_energy'],'nominal128_gain':n['deposited_nuclear_over_initial_energy'],'source128_target_yield':s['desired_product_per_initial_catalyst'],'source_mesh_relative_gain':s['deposited_nuclear_over_initial_energy']/o['deposited_nuclear_over_initial_energy']-1,'source_vs_nominal128_relative_gain':s['deposited_nuclear_over_initial_energy']/n['deposited_nuclear_over_initial_energy']-1,'physical_duration_s':s['elapsed_s'],'prepared_energy_J':s['initial_energy_J'],'remap':s['mixture_inputs']['source_refinement']}
  for refname,ref in [('nominal128',n),('source64',o)]:
   for key in ['mass_kg','initial_energy_J','elapsed_s']:
    err=s[key]/ref[key]-1;checks.append({'name':cat+' '+key+' match '+refname,'value':err,'passed':abs(err)<1e-10})
  stages[cat]=row
 for f,h in hashes.items():assert sha(f)==h
 for c in checks:c['passed']=bool(c['passed'])
 assert all(c['passed'] for c in checks)
 out={'schema':1,'scope':__doc__,'families':families,'comparisons':compares,'stage_comparison':stages,'checks':checks,'source_sha256':sources,'code_sha256':hashes,'parameters':old['parameters'],'qualification':'One spatial refinement with the original source entropy choices held fixed; this does not establish continuum convergence or3D cold assembly. Both receiver masses retain conditional99percent capture and hot-profile source timing. No old64 family is replaced.'}
 vi.write(P/'integration.json',out)
 print(json.dumps({'checks':len(checks),'stages':{k:{q:v[q] for q in ['source64_gain','source128_gain','nominal128_gain','source_mesh_relative_gain','source_vs_nominal128_relative_gain']} for k,v in stages.items()},'plants':{k:{sig:{q:v[q] for q in ['feasible','net_W_kg','net_W','mass_kg','helium_per_s','sequential_occupancy_mass999','fresh_C_H_mass_ratio']} for sig,v in f['capacity_plants'].items()} for k,f in families.items()},'comparisons':compares},indent=2))
if __name__=='__main__':main()
