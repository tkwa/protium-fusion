"""Preserve the four B1 source/recipe/rate controls and independently audit traces."""
from pathlib import Path
import sys,json,hashlib,shutil
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'mixture-design'),str(ROOT)]
import pulse_diagnostics as pd
from fast_eos import FastEOS
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
base=Path('/home/tkwa/src/protium-fusion/output')
paths={'H2_carbon':base/'source128/B1-H2-source128-carbon-median','H2_low':base/'b1-source-robustness/h2low/B1-H2-source128-carbon-median-rates-carbon-n15-low-updated-median','H1p5_carbon':base/'b1-source-robustness/B1-H1p5-source128-carbon','H1p5_low':base/'b1-source-robustness/B1-H1p5-source128-low'}
eos=FastEOS();rows={};audit=[];checks=[]
for name,stem in paths.items():
 d=json.loads(stem.with_suffix('.json').read_text());inp=json.loads(stem.with_suffix('.inputs.json').read_text());profile=json.loads(stem.with_suffix('.initial.json').read_text());assert d['success'] and d['termination_reason']=='final_time' and d['zones']==128
 dest=P/'records'/stem.name
 for ext in ['.json','.initial.json','.inputs.json']:shutil.copyfile(stem.with_suffix(ext),dest.with_suffix(ext))
 yy=inp['initial_Y'];h=yy['p']/yy['n15'];duration=inp['duration_s'];assert abs(d['elapsed_s']/duration-1)<1e-12
 summary,_=pd.reduce({'label':name,'family':'B1-source-recipe-rate-robustness','result':str(stem.with_suffix('.json')),'trace':str(stem.with_suffix('.npz')),'reaction_flows':False},eos,False);assert summary['checks']['all_passed'];audit.append(summary)
 rows[name]={'record':str(dest.with_suffix('.json').relative_to(ROOT)),'initial':str(dest.with_suffix('.initial.json').relative_to(ROOT)),'inputs':str(dest.with_suffix('.inputs.json').relative_to(ROOT)),'trace':str(stem.with_suffix('.npz')),'duration_s':duration,'initial_Y':yy,'eos':'FastEOS','result_sha256':sha(stem.with_suffix('.json')),'initial_profile_sha256':sha(stem.with_suffix('.initial.json')),'input_record_sha256':sha(stem.with_suffix('.inputs.json')),'raw_trace_sha256':sha(stem.with_suffix('.npz')),'h_per_catalyst':h,'gain':d['deposited_nuclear_over_initial_energy'],'target_yield':d['desired_product_per_initial_catalyst'],'prepared_energy_J':d['initial_energy_J'],'rate_update':d['rate_update'],'peak_temperature_K':d['peak_temperature_K'],'energy_residual_fraction_initial':d['energy_residual_fraction_initial']}
for h in ['H2','H1p5']:
 a=rows[h+'_carbon'];b=rows[h+'_low'];checks.extend([{'name':h+' bit-identical initial profiles across rate curves','passed':a['initial_profile_sha256']==b['initial_profile_sha256']},{'name':h+' equal prepared energies across rate curves','passed':abs(a['prepared_energy_J']/b['prepared_energy_J']-1)<1e-12}])
assert all(x['passed'] for x in checks)
out={'scope':__doc__,'cases':rows,'trace_audits':audit,'additional_checks':checks,'checks':{'count':sum(x['checks']['count'] for x in audit)+len(checks),'passed':sum(x['checks']['passed'] for x in audit)+sum(x['passed'] for x in checks),'all_passed':True},'qualification':'Source-shaped spherical128 profiles only; no3D cold assembly. Low is the NACREII evaluation low curve, not an assigned credible quantile. The inherited contact-construction scope strings describe64 original-rate preparation; actual mesh and nuclear revision are the top-level result metadata and this packet. No source retuning between rates.','code_sha256':sha(__file__),'diagnostics_sha256':sha(ROOT/'mixture-design/pulse_diagnostics.py')}
(P/'source-variants.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'checks':out['checks'],'cases':{k:{q:v[q] for q in ['gain','target_yield','prepared_energy_J','peak_temperature_K']} for k,v in rows.items()}},indent=2))
