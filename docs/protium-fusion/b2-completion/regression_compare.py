"""Compare the four saved small stages against corrected-viscosity reruns."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--scratch',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',default='compression_cap');a=p.parse_args()
old_dir=a.root/'reassessment-2049/small-stages';rows=[];checks=[];hashes={}
def hashed(path):
 hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
def check(stage,name,value,limit,passed):
 checks.append(dict(stage=stage,check=name,value=value,limit=limit,passed=bool(passed)))
for stage in ['A','B1','C13','He']:
 radius=3 if stage=='He' else 30
 old_path=old_dir/f'scaled-hot-{stage}-64-R{radius}.json';new_path=a.scratch/f'regression-{stage}-{a.mode}.json';trace_path=new_path.with_suffix('.npz');input_path=new_path.with_suffix('.inputs.json')
 old=json.loads(old_path.read_text());new=json.loads(new_path.read_text());inputs=json.loads(input_path.read_text())
 for path in [old_path,new_path,trace_path,input_path]:hashed(path)
 fields=['initial_energy_J','nuclear_energy_J','integrated_nuclear_energy_J','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','final_internal_J','final_kinetic_J','weak_energy_escaped_upper_J','radiation_escaped_J','peak_temperature_K','energy_residual_fraction_initial','elapsed_s']
 values={}
 for key in fields:
  x=old.get(key);y=new.get(key)
  values[key]=dict(legacy=x,corrected=y,relative_change=None if x is None or x==0 else y/x-1)
 products={key:dict(legacy=x,corrected=new['output_mass_fractions'][key],difference=new['output_mass_fractions'][key]-x) for key,x in old['output_mass_fractions'].items()}
 nuclear_agreement=(new['integrated_nuclear_energy_J']-new['nuclear_energy_J'])/new['initial_energy_J']
 with np.load(trace_path) as z:
  nn=len(z['species']);nz=int(z['zones']);nv=int(z['variables_per_zone']);state=z['state'].reshape(nz,nv,-1)
  rr=np.concatenate([np.zeros((1,state.shape[-1])),state[:,nn+1,:]])
  vv=np.concatenate([np.zeros((1,state.shape[-1])),state[:,nn+2,:]])
  dv=np.diff(rr*rr*vv,axis=0)
  widths=np.diff(rr,axis=0)
  bulk=np.maximum(-3*dv*np.diff(rr,axis=0)/np.diff(rr**3,axis=0),0)
  radial=np.maximum(-np.diff(vv,axis=0),0)
  compress=np.minimum(radial,bulk) if a.mode=='compression_cap' else bulk/3
  wrong=int(np.count_nonzero((compress>0)&(dv>0)))
  mass_numbers=np.array([float(''.join(c for c in k if c.isdigit())) if any(c.isdigit() for c in k) else 1 for k in z['species']])
  trace_fractions=state[:,:nn,:]*mass_numbers[None,:,None]
  trace_baryon=np.max(np.abs(np.sum(trace_fractions,axis=1)-1))
  trace_minimum=float(np.min(trace_fractions))
 check(stage,'full requested duration',new['elapsed_sound_times'],4,new['success'] and new['termination_reason']=='final_time' and abs(new['elapsed_sound_times']-4)<1e-12)
 check(stage,'energy residual relative to initial',new['energy_residual_fraction_initial'],5e-5,abs(new['energy_residual_fraction_initial'])<5e-5)
 check(stage,'independent integrated nuclear energy agrees with isotope mass difference',nuclear_agreement,5e-6,abs(nuclear_agreement)<5e-6)
 check(stage,'final cell baryon residual',new['max_cell_baryon_residual'],1e-10,new['max_cell_baryon_residual']<1e-10)
 check(stage,'stored trace baryon residual',float(trace_baryon),1e-10,trace_baryon<1e-10)
 check(stage,'minimum final mass fraction',new['minimum_mass_fraction'],-1e-12,new['minimum_mass_fraction']>=-1e-12)
 check(stage,'minimum stored trace mass fraction',trace_minimum,-1e-12,trace_minimum>=-1e-12)
 check(stage,'initial energy unchanged',new['initial_energy_J']/old['initial_energy_J']-1,1e-12,abs(new['initial_energy_J']/old['initial_energy_J']-1)<1e-12)
 check(stage,'mass unchanged',new['mass_kg']/old['mass_kg']-1,1e-12,abs(new['mass_kg']/old['mass_kg']-1)<1e-12)
 check(stage,'positive physical shell widths',float(np.min(widths)),0,np.min(widths)>0)
 check(stage,'no expanding cell receives scalar artificial pressure in stored trace',wrong,0,wrong==0)
 rows.append(dict(stage=stage,meaning='He denotes the helium-4 triple-alpha factory at3m reference radius, not helium-3 fuel' if stage=='He' else '',success=new['success'],termination_reason=new['termination_reason'],values=values,products=products,inputs=inputs,corrected_metadata_path=str(new_path),corrected_trace_path=str(trace_path)))
for path in [Path(__file__),a.scratch/'regression_runner.py',a.scratch/'regression_sign_audit.py',a.scratch/'model/radial_completion.py',a.scratch/'regression-legacy-sign.json']:hashed(path)
out=dict(scope='Narrow regression of A,B1,C13 and the3m helium factory only. Frozen hot inflows,64zones,full network,radiation and4sound times retained. Corrected continuous compression-cap scalar viscosity; no physical driver or multidimensional validation.',legacy_tolerances=dict(rtol=1e-6,atol=1e-11),corrected_tolerances=dict(rtol=1e-7,atol=1e-12),limitations=['The continuous cap limits the old radial-gradient viscosity when local volume compression is weaker; it preserves homologous compression and the planar-shock limit. This changes the physical approximation in addition to eliminating negative expansion heating.','Old and new integrations have different tolerances and the corrected run carries two passive diagnostic state variables. Changes are not a strict same-tolerance one-line subtraction.','Stored trace checks sample the history; nonnegative viscous heat also follows algebraically from the corrected compression measure.','The helium factory remains a small positive carbon yield; a successful numerical integration does not establish favorable complete-cycle economics.'],legacy_sign_audit=json.loads((a.scratch/'regression-legacy-sign.json').read_text()),stages=rows,checks=checks,checks_passed=sum(c['passed'] for c in checks),checks_total=len(checks),hashes=hashes)
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(checks_passed=out['checks_passed'],checks_total=len(checks),stages=[dict(stage=r['stage'],gain=r['values']['deposited_nuclear_over_initial_energy'],yield_=r['values']['desired_product_per_initial_catalyst'],energy_residual=r['values']['energy_residual_fraction_initial']) for r in rows]),indent=2))
