"""Archive and validate the single C13H2 polymer source128 discriminator."""
from pathlib import Path
import json,sys,shutil,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'mixture-design'),str(ROOT)]
import pulse_diagnostics as pd
from fast_eos import FastEOS
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def write(p,j):p.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
def main():
 stem=Path('/home/tkwa/src/protium-fusion/output/c13-polymer-source/C13-H2-source128-carbon-median')
 d=read(stem.with_suffix('.json'));i=read(stem.with_suffix('.inputs.json'));init=read(stem.with_suffix('.initial.json'))
 nominal=Path(i['nominal_record_expected']);n=read(nominal);duration=i['physical_duration_s']
 checks=[]
 def check(name,ok,value=None):
  checks.append(dict(name=name,passed=bool(ok),value=value))
  if not ok:raise AssertionError(name)
 check('completed requested physical duration',d['success'] and d['termination_reason']=='final_time' and abs(d['elapsed_s']/duration-1)<1e-12)
 check('same128 mesh and corrected viscosity',d['zones']==n['zones']==128 and d['viscosity_mode']==n['viscosity_mode']=='compression_cap')
 check('exact initial composition',d['initial_Y']==n['initial_Y'])
 check('same rates and opacity',d['rate_update']==n['rate_update'] and d['opacity_m2_kg']==n['opacity_m2_kg']==.001)
 for k in ['initial_energy_J','mass_kg','elapsed_s']:
  v=d[k]/n[k]-1;check('matched '+k,abs(v)<1e-10,v)
 check('no profile clipping',i['target_reference_T_factor']==1 and abs(i['prepared_energy_relative_error'])<1e-11)
 dest=P/'records'/stem.name;dest.parent.mkdir(exist_ok=True)
 for ext in ['.json','.initial.json','.inputs.json']:shutil.copyfile(stem.with_suffix(ext),dest.with_suffix(ext))
 trace=Path(d['trace_path']);s,_=pd.reduce({'label':'C13H2 polymer constant-contact source128','family':'polymer-source128-carbon-median','result':str(dest.with_suffix('.json')),'trace':str(trace),'reaction_flows':False},FastEOS(),False)
 check('trace conservation diagnostics',s['checks']['all_passed'])
 mp={'record':str(dest.with_suffix('.json').relative_to(ROOT)),'initial':str(dest.with_suffix('.initial.json').relative_to(ROOT)),'inputs':str(dest.with_suffix('.inputs.json').relative_to(ROOT)),'trace':str(trace),'duration_s':duration,'initial_Y':i['initial_Y'],'eos':'FastEOS','result_sha256':sha(dest.with_suffix('.json')),'initial_profile_sha256':sha(dest.with_suffix('.initial.json')),'input_record_sha256':sha(dest.with_suffix('.inputs.json')),'raw_trace_sha256':sha(trace)}
 write(P/'stage-map.json',{'family':'C13H2-polymer-source128-carbon-median-single-stage','stage_map':{'c13':mp},'scope':'Explicit replacement stage, other stages supplied by separate named family. No material/plant results inferred here.'})
 write(P/'audit.json',{'checks':checks,'trace':s,'source':mp,'nominal_record':str(nominal),'nominal_sha256':sha(nominal),'comparison':{'source_gain':d['deposited_nuclear_over_initial_energy'],'nominal_gain':n['deposited_nuclear_over_initial_energy'],'gain_relative_difference':d['deposited_nuclear_over_initial_energy']/n['deposited_nuclear_over_initial_energy']-1,'source_target_yield':d['desired_product_per_initial_catalyst'],'nominal_target_yield':n['desired_product_per_initial_catalyst'],'source_peak_temperature_K':d['peak_temperature_K'],'source_peak_temperature_scope':'solver accepted-step peak; sampled pressure/radiation maxima are separate'},'code_sha256':{str(p):sha(p) for p in [Path(__file__),P/'run_source.py',ROOT/'mixture-design/pulse_diagnostics.py']}})
 print(json.dumps({'checks':len(checks),'trace_checks':s['checks'],'comparison':read(P/'audit.json')['comparison'],'record':str(dest.with_suffix('.json'))},indent=2))
if __name__=='__main__':main()
