import json,sys,copy,hashlib,argparse,tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent));import material_cycle as g
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--output',type=Path,default=Path(__file__).with_name('material-cycle-validation.json'));ap.add_argument('--verify-raw-traces',action='store_true');args=ap.parse_args();R=args.root.resolve()
raw=json.loads((R/'membrane-receiver/pulse-loading.json').read_text());mapping=json.loads((R/'open-loop-redesign/corrected-stage-map.json').read_text());baseline=json.loads((R/'open-loop-redesign/station-inputs.json').read_text())['cycle']
stages={};checks=[];composition_map=json.loads((R/'component-design/material-cycle-baseline-map.json').read_text())
def check(n,c,e=None):
 checks.append({'name':n,'passed':bool(c),'error':e})
 if not c:raise AssertionError(n)
for k,v in mapping.items():
 d=json.loads((R/v).read_text());tag=d['stage'];trace=raw['stages'][tag]['sources']['raw_trace']['path']
 stages[k]=g.load_stage(R/v,trace_path=trace) if args.verify_raw_traces else g.load_stage(R/v,initial_Y=composition_map[k]['initial_Y'])
 check(k+' actual input normalization',abs(sum(stages[k]['input_counts'].get(n,0) for n in g.DEFAULT_STATES)-1)<1e-12)
 check(k+' recorded H reference',abs(stages[k]['protons_per_catalyst']-(2 if k=='c12' else 1))<1e-12)
 check(k+' hot baryons',abs(stages[k]['hot_baryon_residual'])<1e-10)
r=g.solve_supply_target(stages,eta_heat=.4)
for k in ['fresh_C_H_mass_ratio','fresh_carbon_kg_per_He','fresh_hydrogen_kg_per_He','nuclear_H_consumed_atoms_per_He','gross_H_input_atoms_per_He','unreturned_H_atoms_per_He','gross_fuel_kg_per_He','returned_fuel_kg_per_He','nuclear_Q_MeV_per_He','neutrino_MeV_per_He','H_nonnu_MeV_per_He','U_prepared_MeV_per_He','escaped_prompt_MeV_per_He','bulk_kinetic_MeV_per_He','hot_internal_MeV_per_He','delayed_decay_heat_MeV_per_He','preaux_net_MeV_per_He','minimum_CNO_return_fraction']:
 e=abs(r[k]/baseline[k]-1);check('baseline '+k,e<3e-8,e)
# A change of stage basis produces genuinely nonidentity F while preserving
# the same physically counted recipe and every extensive ledger quantity.
keys=list(stages);mix=np.eye(4)*.99
for j in range(4):mix[(j+1)%4,j]=.01
mixed={}
for j,k in enumerate(keys):
 s=copy.deepcopy(stages[k])
 for field in ['input_counts','hot_output','cold_output']:
  s[field]={}
  for i,kk in enumerate(keys):g.add(s[field],stages[kk][field],mix[i,j])
 for field in ['Q_hot_MeV','Q_cold_MeV','Q_total_MeV','nu_hot_MeV','nu_cold_MeV','nu_MeV','input_MeV','radiated_MeV','baryons_per_attempt']:
  s[field]=sum(mix[i,j]*stages[kk][field] for i,kk in enumerate(keys))
 hot=[stages[kk]['input_MeV']+stages[kk]['Q_hot_MeV']-stages[kk]['nu_hot_MeV']-stages[kk]['radiated_MeV'] for kk in keys]
 s['final_kinetic_fraction']=sum(mix[i,j]*hot[i]*stages[kk]['final_kinetic_fraction'] for i,kk in enumerate(keys))/sum(mix[i,j]*hot[i] for i in range(4))
 s['protons_per_catalyst']=s['input_counts'].get('p',0.);s['attempts_per_pulse']=s['mass_kg']/(s['baryons_per_attempt']*g.MU)
 mixed[k]=s
b=g.solve_cycle(stages,carbon_return=.9997,hydrogen_return=.99)
m=g.solve_cycle(mixed,carbon_return=.9997,hydrogen_return=.99)
for k in ['fresh_C_H_mass_ratio','nuclear_Q_MeV_per_He','H_nonnu_MeV_per_He','U_prepared_MeV_per_He','bulk_kinetic_MeV_per_He','hot_internal_MeV_per_He','preaux_net_MeV_per_He','gross_fuel_kg_per_He']:
 e=abs(m[k]/b[k]-1);check('nonidentity recipe '+k,e<1e-8,e)
# Add inertly carried excess H: return1 cancels it, return<1 charges actual loss.
extra=copy.deepcopy(stages)
for s in extra.values():
 s['input_counts']['p']+=2;s['cold_output']['p']+=2;s['hot_output']['p']+=2;s['baryons_per_attempt']+=2;s['protons_per_catalyst']+=2
 s['attempts_per_pulse']=s['mass_kg']/(s['baryons_per_attempt']*g.MU)
b=g.solve_cycle(stages,carbon_return=.9997,hydrogen_return=1)
e=g.solve_cycle(extra,carbon_return=.9997,hydrogen_return=1)
check('excess H fully returned does not dilute freshH',abs(e['fresh_hydrogen_kg_per_He']/b['fresh_hydrogen_kg_per_He']-1)<1e-9)
b=g.solve_cycle(stages,carbon_return=.9997,hydrogen_return=.99)
e=g.solve_cycle(extra,carbon_return=.9997,hydrogen_return=.99)
expected=.01*2*sum(b['stage_attempts_per_He'].values());actual=e['external_input_atoms_per_He']['p']-b['external_input_atoms_per_He']['p']
check('excess H imperfect return incurs correct loss',abs(actual-expected)<1e-9,actual-expected)
# No initial-composition inference by stage name.
try:g.load_stage(R/mapping['c12'])
except ValueError as e:check('missing actual composition rejected','Actual initial composition required' in str(e))
else:check('missing actual composition rejected',False)
# Multi-product explicit decay chain has correct baryons and Q.
o,q,nu=g.cold_decay({'be8':1.},decay_map={'be8':{'he4':2}},neutrino_MeV={'be8':0})
check('explicit multiple-daughter decay',o=={'he4':2.} and q>0 and nu==0)
# Physical-duration metadata supersedes sound-time labels; contradictory input
# composition and genuinely unfinished endpoints are rejected.
with tempfile.TemporaryDirectory(prefix='protium-material-check-') as td:
 p=Path(td)/'duration-check.json';q=json.loads((R/mapping['c12']).read_text());q['initial_Y']=stages['c12']['initial_Y'];q['mixture_inputs']={'duration_s':q['elapsed_s']};q['elapsed_sound_times']=3.9;p.write_text(json.dumps(q))
 case=g.load_stage(p);check('physical duration accepts fewer than four sound times',case['requested_physical_duration_s']==q['elapsed_s'])
 try:g.load_stage(p,initial_counts={'p':1,'c12':1})
 except ValueError:check('explicit composition disagreement rejected',True)
 else:check('explicit composition disagreement rejected',False)
 q['elapsed_s']*=.99;p.write_text(json.dumps(q))
 try:g.load_stage(p)
 except ValueError:check('unfinished physical duration rejected',True)
 else:check('unfinished physical duration rejected',False)
working=g.solve_supply_target(stages,eta_heat=.6)
out={'initial_composition_validation':'Direct raw trace' if args.verify_raw_traces else 'Archived trace-extracted composition; pass --verify-raw-traces to repeat extraction',
 'scope':'Independent reproduction, nonidentity recipe change-of-basis, excess-H return and explicit multi-daughter decay checks; synthetic mixtures are algebra tests, not pulse predictions.',
 'baseline40':r,'working60':working,'loaded_stages':stages,'checks':checks,'passed':sum(v['passed'] for v in checks),'total':len(checks),'module_sha256':hashlib.sha256(Path(g.__file__).read_bytes()).hexdigest()}
args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'passed':out['passed'],'total':len(checks),'working60':{k:working[k] for k in ['minimum_CNO_return_fraction','fresh_C_H_mass_ratio','H_nonnu_MeV_per_He','U_prepared_MeV_per_He','preaux_net_MeV_per_He']}},indent=2))
