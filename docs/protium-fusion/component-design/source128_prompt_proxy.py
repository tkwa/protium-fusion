"""Replay existing finite-depth prompt proxy on exact source128 traces.
Only missing H/catalyst display metadata is adapted; raw physics arrays unchanged.
"""
from pathlib import Path
import json,hashlib,sys,math
P=Path(__file__).parent;sys.path.insert(0,str(P));import mixture_radiation as model
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
src=P/'source-family-carbon-median128/integration.json';f=json.loads(src.read_text())['families']['source128'];scratch=Path('/home/tkwa/src/protium-fusion/output/source128-prompt-proxy');scratch.mkdir(parents=True,exist_ok=True)
rows={};checks=[];inputhash={str(src):sha(src),str(Path(model.__file__)):sha(model.__file__)}
for cat,item in f['stage_inputs'].items():
 selected=item['selected'];stage=selected['stage_data'];record=Path(item['record_path']);trace=Path(selected['receiver']['source_paths']['trace']);d=json.loads(record.read_text());inputhash[str(record)]=sha(record);inputhash[str(trace)]=sha(trace)
 d.setdefault('mixture_inputs',{})['hydrogen_per_catalyst']=stage['protons_per_catalyst']
 stem=scratch/cat;stem.with_suffix('.json').write_text(json.dumps(d,indent=2)+'\n');link=stem.with_suffix('.npz')
 if link.is_symlink():assert link.resolve()==trace.resolve()
 elif link.exists():assert sha(link)==sha(trace)
 else:link.symlink_to(trace)
 r=model.calculate(stem);r['actual_source_record']=str(record);r['actual_source_trace']=str(trace);r['source_family']='source128 exact completed hot profile';r['metadata_adapter']='Only missing mixture_inputs.hydrogen_per_catalyst display value added from actual stage input counts; no array or physical scalar changed.'
 for name,c in r['cases'].items():
  m=c['materials']['carbon'];c['bulk_carbon_k50_cp1800_reference_rise_K']=m['resolved_surface_no_conduction_sensible_rise_K']/1.8+m['below1keV_all_surface_conduction_peak_K']*math.sqrt(100*1000/(50*1800))
  c['shallow_lowk4_cp1800_reference_rise_K']=m['resolved_surface_no_conduction_sensible_rise_K']/1.8+m['below1keV_all_surface_conduction_peak_K']*math.sqrt(100*1000/(4*1800))
  c['rescaling_scope']='Resolved no-conduction depth term scales as1/cp; only unresolved soft surface term scales as1/sqrt(k cp). This remains a homogeneous first-interaction proxy, not a solved layered damaged skin.'
 assert r['checks']['all_passed'];rows[stage['label']]=r
out={'scope':__doc__,'radius_m':model.R,'source_sha256':inputhash,'code_sha256':sha(__file__),'model_dependency_sha256':{str(p):sha(p)for p in[model.ROOT/'component-design/receiver_radiation.py',model.ROOT/'membrane-receiver/spectral_proxy.py']},'stages':rows,'qualification':['Actual source128 grey energy/time/outer-state histories; frequency shapes remain explicit proxy alternatives, not new spectral radiation transport.','Prompt and slow material peaks are separated in time; do not add their isolated maxima as simultaneous.','All source records/traces and the one metadata-only adapter are recorded.','k50/cp1800 and shallow k4/cp1800 are required-property sensitivities, not measured irradiated assemblies.','The first-interaction finite-depth result does not certify electron transport, secondary scattering, fluorescence or layered thermal contact.']}
(P/'source128-prompt-proxy.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({n:{'fluence_J_m2':r['prompt_fluence_J_m2'],'checks':r['checks'],'proxy_rises':{k:{'bulk50':v['bulk_carbon_k50_cp1800_reference_rise_K'],'skin4':v['shallow_lowk4_cp1800_reference_rise_K']}for k,v in r['cases'].items()}}for n,r in rows.items()},indent=2))
