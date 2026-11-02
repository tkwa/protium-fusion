"""Matched B1 source/rate envelope within the coherent source128 other stages.
Only B1 is varied; the rest retain their exact median-carbon source128 endpoints.
NACREII low/high are evaluation limits, not assigned probability quantiles.
"""
from pathlib import Path
import sys,json,copy,shutil
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'mixture-design'),str(ROOT/'component-design'),str(ROOT)]
import variant_integrate as vi
import material_cycle as mc
import pulse_diagnostics as pd
from fast_eos import FastEOS
sha=vi.sha
RAW=Path('/home/tkwa/src/protium-fusion/output/b1-source-robustness')
def main():
 basefile=P.parent/'source-family-carbon-median128/source128-reduced-inputs.json';plantfile=P.parent/'source-family-carbon-median128/integration.json';base=vi.read(basefile);baseplant=vi.read(plantfile);cases=vi.read(P/'source-variants.json')['cases'];receiver=vi.module('b1_envelope_receiver',ROOT/'mixture-design/receiver.py');pc=vi.module('b1_envelope_plant',ROOT/'mixture-design/plant_compare.py');codes=[Path(__file__),ROOT/'mixture-design/receiver.py',ROOT/'mixture-design/plant_compare.py',ROOT/'component-design/material_cycle.py',ROOT/'cold-sink-redesign/station.py'];hashes={str(p):sha(p) for p in codes};checks=[];results={};reduced={};audits=[]
 for h,stem in [('H2','B1-H2-source128-carbon-median-rates-carbon-n15-high-updated-median'),('H1p5','B1-H1p5-source128-carbon-rates-carbon-n15-high-updated-median')]:
  raw=RAW/'high'/stem;dest=P/'records'/stem
  for ext in ['.json','.initial.json','.inputs.json']:shutil.copyfile(raw.with_suffix(ext),dest.with_suffix(ext))
  d=vi.read(dest.with_suffix('.json'));assert d['success'] and d['termination_reason']=='final_time';i=vi.read(dest.with_suffix('.inputs.json'));cases[h+'_high']={'record':str(dest.with_suffix('.json').relative_to(ROOT)),'trace':str(raw.with_suffix('.npz')),'duration_s':i['duration_s'],'initial_Y':i['initial_Y'],'rate_update':d['rate_update']}
  s,_=pd.reduce({'label':h+'_source_high','family':'B1-source-rate-envelope','result':str(dest.with_suffix('.json')),'trace':str(raw.with_suffix('.npz')),'reaction_flows':False},FastEOS(),False);assert s['checks']['all_passed'];audits.append(s)
  other=ROOT/cases[h+'_carbon']['initial'];checks.append({'name':h+' high/control bit-identical source profile','passed':sha(dest.with_suffix('.initial.json'))==sha(other)})
 for label,spec in cases.items():
  rec=ROOT/spec['record'];d=vi.read(rec);assert d['zones']==128 and d['success'] and d['termination_reason']=='final_time';r=receiver.reduce_case(rec,ROOT.parents[1],radius=1e9,thermal_points=65536);assert r['usable_ballistic'] and r['completed_requested_duration'];st=mc.load_stage(rec,trace_path=Path(spec['trace']));entry={'name':rec.stem,'path':str(rec),'record':d,'receiver':r,'stage_data':st,'label':label};reduced[label]=entry
  selected=[copy.deepcopy(e) for e in base['cases'].values() if e['record']['stage']!='B1']+[entry];cy=pc.cycle_for(selected);plants={}
  for sig in [5000.,2000.]:
   pp=pc.best_rate(selected,cy,absorber_sigma=sig);checks+=vi.result_checks(selected,cy,pp);plants[str(int(sig))]=pp
  results[label]={'cycle':cy,'capacity_plants':plants,'source_gain':d['deposited_nuclear_over_initial_energy'],'target_yield':d['desired_product_per_initial_catalyst'],'source_prepared_energy_J':d['initial_energy_J'],'rate_update':d['rate_update']}
 for sig in ['5000','2000']:
  expected=baseplant['families']['source128']['capacity_plants'][sig];actual=results['H2_carbon']['capacity_plants'][sig]
  for key in ['net_W','mass_kg','net_W_kg','helium_per_s']:checks.append({'name':'exact source128 baseline replay '+sig+' '+key,'value':actual[key]/expected[key]-1,'passed':abs(actual[key]/expected[key]-1)<1e-11})
 for p,h in hashes.items():assert sha(p)==h
 for c in checks:c['passed']=bool(c['passed'])
 assert all(c['passed'] for c in checks)
 out={'scope':__doc__,'fixed_other_stages':{k:e for k,e in base['cases'].items() if e['record']['stage']!='B1'},'B1_cases':reduced,'results':results,'checks':checks,'high_curve_trace_audits':audits,'source_sha256':{str(basefile):sha(basefile),str(plantfile):sha(plantfile),str(P/'source-variants.json'):sha(P/'source-variants.json')},'code_sha256':hashes,'qualification':'Stage-local N15-rate sensitivity with exact other-stage source128 states. Three-dimensional cold contact remains unqualified. OriginalREACLIB N15 is the control, not NACREII adopted; no probability level assigned to high/low. Separate5000/2000receiver cases retain99percent energycapture only as a common accounting assumption.'}
 vi.write(P/'source128-envelope.json',out);print(json.dumps({'checks':len(checks),'results':{label:{'Q':r['source_gain'],'plants':{sig:{q:pp[q] for q in ['feasible','net_W_kg','sequential_occupancy_mass999','fresh_C_H_mass_ratio']} for sig,pp in r['capacity_plants'].items()}} for label,r in results.items()}},indent=2))
if __name__=='__main__':main()
