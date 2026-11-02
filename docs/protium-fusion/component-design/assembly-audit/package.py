"""Archive compact assembly records and recompute independent comparison checks."""
from pathlib import Path
import json,hashlib,shutil,math,sys
import numpy as np
P=Path(__file__).resolve().parent;RAW=Path('/home/tkwa/src/protium-fusion/output/assembly-audit')
from initial_screen import quant,load,sha
NAMES=['B1-backward-feos128','B1-backward-feos128-refined','B2-backward-feos128','B2-backward-feos128-refined','B2-backward-feos256']
checks=[];records={};manifest={};destination=P/'records';destination.mkdir(exist_ok=True)
for name in NAMES:
 source=RAW/(name+'.json');d=json.loads(source.read_text());trace=RAW/(name+'.npz');checks.extend([d['numerics']['success'],math.isclose(d['numerics']['elapsed_s'],d['numerics']['requested_duration_s'],rel_tol=1e-12),d['energy']['max_total_energy_residual_fraction']<1e-8,d['minimum_shell_width_m']>0,d['contact_crossing_max_relative_T_error']<1e-6,sha(trace)==d['trace_sha256']]);shutil.copyfile(source,destination/source.name);manifest[str(source)]={'sha256':sha(source),'archive':'records/'+source.name};manifest[str(trace)]={'sha256':sha(trace),'scope':'Large trace retained in persistent research output.'};records[name]=d
comparisons={}
for stage in ['B1','B2']:
 a=records[stage+'-backward-feos128'];b=records[stage+'-backward-feos128-refined'];comparisons[stage+'_time_EOS_refinement']={key:[x/y-1 for x,y in zip(b[key],a[key])] for key in ['contact_crossing_time_mass_quantiles_s','contact_crossing_radius_mass_quantiles_m']}
a=records['B2-backward-feos128-refined'];b=records['B2-backward-feos256'];comparisons['B2_128_to256']={key:[x/y-1 for x,y in zip(b[key],a[key])] for key in ['contact_crossing_time_mass_quantiles_s','contact_crossing_radius_mass_quantiles_m']}
for label in ['k0','k0p001','k0p01']:
 name='B1-forward-'+label
 for ext in ['.json','.initial.json','.inputs.json']:
  source=RAW/(name+ext);shutil.copyfile(source,destination/source.name);manifest[str(source)]={'sha256':sha(source),'archive':'records/'+source.name}
 trace=RAW/(name+'.npz');manifest[str(trace)]={'sha256':sha(trace),'scope':'Large trace retained in persistent research output.'};d=json.loads((RAW/(name+'.json')).read_text());records[name]=d
 checks.extend([d['success'],d['termination_reason']=='final_time',not d['nuclear_enabled'],d['artificial_viscosity_linear']==0,d['artificial_viscosity_quadratic']==0,abs(d['energy_residual_fraction_initial'])<1e-8])
control=records['B1-forward-k0'];w=np.array(control['profile']['mass_fraction_of_total']);target=load('B1')[0];fcomp={}
for label in ['k0','k0p001','k0p01']:
 d=records['B1-forward-'+label];entry={'prepared_J':d['initial_energy_J'],'radiation_escaped_J':d['radiation_escaped_J'],'radiation_escaped_fraction_prepared':d['radiation_escaped_J']/d['initial_energy_J'],'energy_residual_fraction_initial':d['energy_residual_fraction_initial'],'state_relative_to_no_radiation':{}}
 for field in ['T_K','rho_kg_m3','velocity_m_s','outer_radius_m']:
  actual=np.array(d['profile'][field]);ref=np.array(control['profile'][field]);rel=actual/ref-1;entry['state_relative_to_no_radiation'][field]={'mass_quantiles_relative':quant(rel,w),'max_absolute_relative':float(max(abs(rel))),'mass_weighted_absolute_change':float(w@abs(actual-ref))}
 entry['temperature_relative_to_original_source128_mass_quantiles']=quant(np.array(d['profile']['T_K'])/np.array(target['T_K'])-1,w);fcomp[label]=entry
 checks.append(math.isclose(d['initial_energy_J'],control['initial_energy_J'],rel_tol=1e-14))
summary={'scope':'Conditional assembly audit: kinematic and EOS checks, reversible FEOS pressure-work continuation, actual controlled FastEOS30us preparation radiation and separate A/C13 omission screens. These do not supply a 3D cold-source trajectory.','code_sha256':sha(__file__),'backward_cases':{name:{k:d[k] for k in ['stage','numerics','energy','contact_active_mass_fraction','contact_crossed_mass_fraction','contact_not_crossed_active_mass_fraction','contact_crossing_time_mass_quantiles_s','contact_crossing_radius_mass_quantiles_m','reverse_compression_max_mass_fraction','return_to_initial_max_dimensionless_error']} for name,d in records.items() if 'backward' in name},'refinement_comparisons':comparisons,'forward_B1_comparison':fcomp,'validation':{'passed':sum(checks),'count':len(checks),'all_passed':all(checks)}}
assert all(checks);(P/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
for name in ['initial-screen.json','omission-screen.json','other-stages-screen.json','B1-precursor-EOS-bridge.json','characteristic-refinement.json']:
 manifest[str(P/name)]={'sha256':sha(P/name),'scope':'Current assembly result.'}
for stage in ['A','C13']:
 path=RAW/(stage+'-ballistic-isentropic-screen.npz');manifest[str(path)]={'sha256':sha(path),'scope':'Large diagnostic trace retained in persistent research output.'}
(P/'manifest.json').write_text(json.dumps({'files':manifest,'summary_sha256':sha(P/'summary.json')},indent=2)+'\n');print(json.dumps(summary['validation']))
