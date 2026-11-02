"""Freeze one coherent four-stage source family and independently reduce its traces.
No discovery or replacement of previous maps. All sources must be completed64,
FastEOS/compression-cap, kappa.001, same median-carbon table and runtime kernel.
"""
from pathlib import Path
import sys,json,hashlib,shutil
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent.parent;sys.path[:0]=[str(ROOT/'mixture-design'),str(ROOT)]
import pulse_diagnostics as pd
from fast_eos import FastEOS
import pulse_model as pm
raw={
'n14':Path('/home/tkwa/src/protium-fusion/output/source-updated-four/AH3-N2000-local-entropy-fixed-energy-reference-floor-rates-updated-median'),
'n15':Path('/home/tkwa/src/protium-fusion/output/b1-source/updated/B1-H2-constant-contact-T0p9-rates-updated-median'),
'c12':Path('/home/tkwa/src/protium-fusion/output/b2-source/B2-H2-constant-contact-carbon-median'),
'c13':Path('/home/tkwa/src/protium-fusion/output/source-updated-four/C13H4-N2000-local-entropy-fixed-energy-rates-updated-median')}
failed=Path('/home/tkwa/src/protium-fusion/output/b2-source/B2-H2-N2000-carbon-median')
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
checks=[];stage_map={};stages={};rows=[];rate_identity=None;eos=FastEOS()
for cat,stem in raw.items():
 rec=stem.with_suffix('.json');inp=stem.with_suffix('.inputs.json');init=stem.with_suffix('.initial.json');trace=stem.with_suffix('.npz');d=json.loads(rec.read_text());i=json.loads(inp.read_text());p=json.loads(init.read_text());tag={'n14':'A','n15':'B1','c12':'B2','c13':'C13'}[cat];expectedh={'n14':3,'n15':2,'c12':2,'c13':4}[cat];duration=i.get('duration_s',i.get('physical_duration_s'))
 assert d['success'] and d['termination_reason']=='final_time' and abs(d['elapsed_s']/duration-1)<1e-12
 assert d['stage']==tag and d['zones']==64 and d['opacity_m2_kg']==.001 and d['viscosity_mode']=='compression_cap'
 assert abs(i['initial_Y']['p']/i['initial_Y'][cat]-expectedh)<1e-12
 rate=d['rate_update'];ident={k:rate[k] for k in ['revision','mode','quantile','table_sha256','kernel_sha256','interpolation','reverse_policy']};assert ident['mode']=='updated' and ident['quantile']=='median'
 if rate_identity is None:rate_identity=ident
 else:assert rate_identity==ident
 dest=P/'records'/stem.name;dest.parent.mkdir(parents=True,exist_ok=True)
 for ext in ['.json','.initial.json','.inputs.json']:shutil.copyfile(stem.with_suffix(ext),dest.with_suffix(ext))
 spec={'record':str(dest.with_suffix('.json').relative_to(ROOT)),'initial':str(dest.with_suffix('.initial.json').relative_to(ROOT)),'inputs':str(dest.with_suffix('.inputs.json').relative_to(ROOT)),'trace':str(trace),'duration_s':duration,'initial_Y':i['initial_Y'],'eos':'FastEOS','result_sha256':sha(rec),'initial_profile_sha256':sha(init),'input_record_sha256':sha(inp),'raw_trace_sha256':sha(trace)};stage_map[cat]=spec;stages[cat]={'stage':tag,'hydrogen_per_catalyst':expectedh,**spec,'initial_profile_exact':p,'rate_update':rate}
 s,h=pd.reduce(dict(label=f'{tag} source median64',family='source-median64',result=str(rec),trace=str(trace),reaction_flows=False),eos,False);assert s['checks']['all_passed']
 baseline=json.loads(Path(i['baseline_source']).read_text());delta=d['initial_energy_J']/baseline['initial_energy_J']-1;assert abs(delta)<1e-12
 checks.extend([dict(name=tag+' completed/profile/recipe/rate identity',passed=True),dict(name=tag+' prepared-energy match',value=delta,passed=True)])
 rows.append(s)
# Preserve the completed low-gain geometry control separately, never as a selected stage.
failed_dest=P.parent/'b2-source'/failed.name
for ext in ['.json','.initial.json','.inputs.json']:shutil.copyfile(failed.with_suffix(ext),failed_dest.with_suffix(ext))
s,_=pd.reduce(dict(label='B2 N2000 completed low-gain control',family='excluded-source-control',result=str(failed.with_suffix('.json')),trace=str(failed.with_suffix('.npz')),reaction_flows=False),eos,False);assert s['checks']['all_passed'];rows.append(s)
definition={'family':'source64-A3-B1H2-B2H2-C13H4-carbon-median','scope':'One explicit four-stage set: controlled A preheat, constant-contact-entropy B1/B2, localN2000 C13. All full27network64, same median carbon tables/kernel, FastEOS, opacity.001, original reference physicaldurations, matched prepared energies. Prescribedhotflow source tests, not3Dcoldassembly.','stage_map':stage_map,'stages':stages,'rate_identity':rate_identity,'excluded_control':{'record':str(failed_dest.with_suffix('.json').relative_to(ROOT)),'trace':str(failed.with_suffix('.npz')),'status':'Completed low-gain N2000 B2 entropy profile; retained as a physical profile sensitivity, not numericalfailure.'},'builder_sha256':sha(__file__)}
(P/'stage-map.json').write_text(json.dumps(definition,indent=2)+'\n');(P/'stage-paths.json').write_text(json.dumps({k:v['record'] for k,v in stage_map.items()},indent=2)+'\n');audit={'scope':pd.__doc__,'stage_map_sha256':sha(P/'stage-map.json'),'diagnostics_code_sha256':sha(ROOT/'mixture-design/pulse_diagnostics.py'),'code_sha256':sha(__file__),'cases':rows,'additional_checks':checks,'checks':{'count':sum(x['checks']['count'] for x in rows)+len(checks),'passed':sum(x['checks']['passed'] for x in rows)+sum(x['passed'] for x in checks),'all_passed':True}};(P/'trace-audit.json').write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n')
print(json.dumps({'family':definition['family'],'checks':audit['checks'],'cases':[{k:s[k] for k in ['label','nuclear_gain','solver_saved_step_peak_temperature_K','maximum_sampled_energy_residual_fraction_initial']} for s in rows]},indent=2))
