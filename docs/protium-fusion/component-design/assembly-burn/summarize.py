"""Validate and compare exact delivered-profile hot continuations without rerunning."""
from pathlib import Path
import hashlib,json,math
import numpy as np
P=Path(__file__).resolve().parent
ROOT=P.parent.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 checks=[]
 def check(name,ok,detail=None):
  checks.append({'name':name,'passed':bool(ok),'detail':detail})
  if not ok:raise AssertionError(name)
 prep=json.loads((P/'prepared-cases.json').read_text())['cases']
 bridge=json.loads((P/'bridge.json').read_text())
 basepath=ROOT/'component-design/source-family-carbon-median128/records/B1-H2-source128-carbon-median.json'
 baseline=json.loads(basepath.read_text())
 def scalars(d):
  return {k:d[k] for k in ['initial_energy_J','initial_internal_energy_J','initial_kinetic_energy_J','nuclear_energy_J','weak_energy_escaped_upper_J','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','final_internal_J','final_kinetic_J','radiation_escaped_J','energy_residual_fraction_initial','elapsed_s','wall_seconds']}
 rows={}
 for label in ['k0','k0p001','k0p01']:
  stem='B1-delivered-'+label;p=P/'records'/(stem+'.json');d=json.loads(p.read_text());inp=d['mixture_inputs'];profile=json.loads((P/'inputs'/(stem+'.initial.json')).read_text());tr=np.load(d['trace_path']);z=int(tr['zones']);nv=int(tr['variables_per_zone']);s=tr['state'].reshape(z,nv,-1);species=tr['species'].tolist();m=tr['cell_mass_kg'];y=np.array([inp['initial_Y'].get(n,0.) for n in species]);E=inp['preparation_initial_energy_J'];ledger=d['assembly_energy_ledger'];rawprep=json.loads(Path(inp['preparation_record']).read_text());mass=d['mass_kg']
  check(label+': source record hash',sha(inp['preparation_record'])==inp['preparation_record_sha256'])
  check(label+': source trace hash',sha(Path(inp['preparation_record']).with_suffix('.npz'))==inp['preparation_trace_sha256'])
  check(label+': hot trace hash',sha(d['trace_path'])==d['trace_sha256'])
  check(label+': exact delivered profile',all(profile[k]==rawprep['profile'][k] for k in ['outer_radius_m','velocity_m_s','T_K','mass_fraction_of_total']))
  check(label+': exact mass grid',np.max(np.abs(m/mass-np.array(profile['mass_fraction_of_total'])))<1e-16)
  check(label+': hot initial T',np.max(np.abs(s[:,27,0]*d['T0_K']/profile['T_K']-1))<1e-13)
  check(label+': hot initial radius',np.max(np.abs(s[:,28,0]*d['R0_m']/profile['outer_radius_m']-1))<1e-13)
  check(label+': hot initial velocity',np.max(np.abs(s[:,29,0]*float(tr['velocity_unit_cm_s'])/100-np.array(profile['velocity_m_s'])))<1e-8)
  check(label+': actual H2 composition',np.max(np.abs(s[:,:27,0]-y))<1e-15 and math.isclose(y[species.index('p')]/y[species.index('n15')],2.))
  check(label+': physical duration complete',d['success'] and d['termination_reason']=='final_time' and abs(d['elapsed_s']/inp['duration_s']-1)<1e-12 and abs(tr['time_s'][-1]/inp['duration_s']-1)<1e-12)
  check(label+': identical hot settings',d['zones']==128 and d['radiation_diffusion'] and d['opacity_m2_kg']==.001 and d['viscosity_mode']=='compression_cap' and d['artificial_viscosity_linear']==.5 and d['artificial_viscosity_quadratic']==2.)
  check(label+': preparation energy closure',abs((d['initial_energy_J']+inp['preparation_escaped_radiation_J'])/E-1)<1e-9)
  check(label+': no burn energy renormalization',abs(d['initial_energy_J']/inp['burn_start_energy_J']-1)<1e-12)
  check(label+': combined energy closure',abs(ledger['combined_energy_residual_fraction_preparation'])<1e-5)
  check(label+': baryon conservation',d['max_cell_baryon_residual']<1e-10)
  check(label+': negligible negative abundance',d['minimum_mass_fraction']>-1e-12)
  check(label+': same rate metadata',d['rate_update']['quantile']=='median')
  row=scalars(d);row.update(assembly_energy_ledger=ledger,output_mass_fractions=d['output_mass_fractions'],record_path=str(p),record_sha256=sha(p),preparation_radiation_fraction=inp['preparation_escaped_radiation_J']/E,preparation_record_sha256=inp['preparation_record_sha256'],trace_path=d['trace_path'],trace_sha256=d['trace_sha256'],rate_update=d['rate_update'])
  rows[label]=row
 d0=rows['k0'];base_dep=baseline['nuclear_energy_J']-baseline['weak_energy_escaped_upper_J'];d0dep=d0['nuclear_energy_J']-d0['weak_energy_escaped_upper_J']
 for label,row in rows.items():
  dep=row['nuclear_energy_J']-row['weak_energy_escaped_upper_J'];row['deposited_nuclear_relative_to_delivered_no_prep_radiation']=dep/d0dep-1;row['deposited_nuclear_relative_to_original_hot_reference']=dep/base_dep-1;row['conversion_relative_to_delivered_no_prep_radiation']=row['desired_product_per_initial_catalyst']/d0['desired_product_per_initial_catalyst']-1;row['conversion_relative_to_original_hot_reference']=row['desired_product_per_initial_catalyst']/baseline['desired_product_per_initial_catalyst']-1
 result={'scope':'Fixed profiles and settings; these calculations isolate hot burn response to whole-fluid preparation endpoints. Preparation is an auxiliary spherical model, and its FEOS/FastEOS bridge is explicit. No cold-launch, merger, 3D stability, or plant qualification is inferred.','code_sha256':sha(__file__),'bridge_path':str(P/'bridge.json'),'bridge_sha256':sha(P/'bridge.json'),'bridge_relative_FastEOS_minus_FEOS_total':bridge['precursor_FastEOS_prepared_J']/bridge['precursor_FEOS_direct_total_J']-1,'original_hot_reference':{'path':str(basepath),'sha256':sha(basepath),**scalars(baseline)},'cases':rows,'checks':checks,'checks_passed':sum(q['passed'] for q in checks),'checks_total':len(checks)}
 (P/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'checks':len(checks),'cases':{k:{q:v[q] for q in ['deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','deposited_nuclear_relative_to_delivered_no_prep_radiation','deposited_nuclear_relative_to_original_hot_reference','peak_temperature_K']} for k,v in rows.items()}},indent=2))
if __name__=='__main__':main()
