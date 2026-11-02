"""Independent saved-trace energy audit and matched nonideal/control comparison.
No hydrodynamic integration. Optionally archive numerically unchanged result records
with the inherited solver's generic no-Coulomb limitation corrected in metadata.
"""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(HERE),str(ROOT)]
from nonideal_eos import NonidealEOS,validate
import pulse_model as pm
CASES=[('A',1),('A',3),('B1',1),('B1',2),('B2',2),('B2',2.5),('C13',1),('C13',4)]
BASE=Path('/home/tkwa/src/protium-fusion/output')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tag(h):return format(h,'.8g').replace('.','p')
def audit(path,d):
 f=np.load(path.with_suffix('.npz'));nn=len(f['species']);nz=int(f['zones']);nv=int(f['variables_per_zone']);z=f['state'].T.reshape(-1,nz,nv);m=f['cell_mass_kg'];r=z[:,:,nn+1]*d['R0_m']*100;rho=(m*1000)/(4*np.pi/3*np.diff(np.c_[np.zeros(len(z)),r**3],axis=1));T=z[:,:,nn]*d['T0_K'];y=np.maximum(z[:,:,:nn],0);norm=y@pm.AA;ye=(y@pm.ZZ)/norm;yi=y.sum(axis=2)/norm
 c=d['mixture_inputs']['EOS_stress'];eos=NonidealEOS(c['frozen_charge_weighted_Z_two_thirds'],strength=c['strength']);state=eos.state(rho,T,ye,yi)
 assert np.min(state.dp_drho)>0,'Nonpositive sampled isothermal compressibility'
 U=(state.e*m).sum(axis=1)*1e-4;v=z[:,:,nn+2]*float(f['velocity_unit_cm_s'])*.01;nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];K=.5*(v*v*nm).sum(axis=1)
 unit=float(f['energy_unit_erg_g'])*1e-4;weak=(z[:,:,nn+3]*m).sum(axis=1)*unit;rad=(z[:,:,nn+4]*m).sum(axis=1)*unit;Q=(z[:,:,nn+5]*m).sum(axis=1)*unit
 E0=K[0]+U[0];res=(K+U+weak+rad-Q)/E0-1
 dynmom=(y@(pm.ZZ**(5/3)))/(y@pm.ZZ);uc,pc=eos.correction(rho,ye);deltaU=(uc*(dynmom/c['frozen_charge_weighted_Z_two_thirds']-1)*m).sum(axis=1)*1e-4
 checks={'initial_internal_relative':U[0]/d['initial_internal_energy_J']-1,'initial_kinetic_relative':K[0]/d['initial_kinetic_energy_J']-1,'final_internal_relative':U[-1]/d['final_internal_J']-1,'final_kinetic_relative':K[-1]/d['final_kinetic_J']-1,'final_nuclear_relative':Q[-1]/d['integrated_nuclear_energy_J']-1,'final_radiation_relative':rad[-1]/d['radiation_escaped_J']-1,'final_record_energy_residual_difference':res[-1]-d['energy_residual_fraction_initial'],'physical_duration_relative':float(f['time_s'][-1])/d['mixture_inputs']['duration_s']-1}
 assert max(abs(v) for v in checks.values())<1e-8,checks
 return {'checks':{k:float(v) for k,v in checks.items()},'saved_trace_energy_residual_max_abs_initial':float(max(abs(res))),'minimum_pressure_Pa':float(np.min(state.p))*.1,'minimum_isothermal_dp_drho_m2_s2':float(np.min(state.dp_drho))*1e-4,'minimum_adiabatic_sound_speed_squared_m2_s2':float(np.min(state.gamma1*state.p/rho))*1e-4,'minimum_cv_J_kg_K':float(np.min(state.c_v))*1e-4,'source_diagnostic_minima':d['EOS_stress_diagnostics'],'frozen_moment_omission_energy_scale':{'largest_abs_if_replaced_by_instantaneous_true_moment_over_prepared':float(max(abs(deltaU))/E0),'final_if_replaced_by_instantaneous_true_moment_over_prepared':float(deltaU[-1]/E0),'dynamic_moment_over_frozen_min_max':[float(np.min(dynmom/c['frozen_charge_weighted_Z_two_thirds'])),float(np.max(dynmom/c['frozen_charge_weighted_Z_two_thirds']))],'scope':'Postprocessed strong-coupling energy difference only; not a corrected trajectory or a complete reaction-screening sensitivity.'},'source_trace_sha256':sha(path.with_suffix('.npz'))}
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=BASE/'nonideal-runs/nonideal_comparison.json');p.add_argument('--allow-partial',action='store_true');p.add_argument('--archive-dir',type=Path);a=p.parse_args();rows=[];missing=[]
 for stage,h in CASES:
  name=f'{stage}-H{tag(h)}-fixed_energy-z64-k0p001-t1';path=BASE/'nonideal-runs'/(name+'-coulomb1.json')
  if not path.exists():missing.append(name);continue
  matches=list((BASE/'mixture-runs').rglob(name+'.json'));assert len(matches)==1,matches;control=matches[0];b=json.loads(control.read_text());d=json.loads(path.read_text())
  assert d['success'] and d['termination_reason']=='final_time',name
  metrics={k:{'control':b[k],'stress':d[k],'relative_change':d[k]/b[k]-1} for k in ['initial_energy_J','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','radiation_escaped_J','final_kinetic_J','final_internal_J']}
  assert abs(metrics['initial_energy_J']['relative_change'])<1e-10
  row={'stage':stage,'H_per_catalyst':h,'control_path':str(control),'stress_path':str(path),'control_sha256':sha(control),'stress_sha256':sha(path),'comparison':metrics,'audit':audit(path,d)}
  if a.archive_dir:
   a.archive_dir.mkdir(parents=True,exist_ok=True);out=a.archive_dir/path.name
   d['archived_raw_result_path']=str(path);d['archived_raw_result_sha256']=sha(path);d['metadata_archive_script_sha256']=sha(__file__)
   d['limitations']=[('No reaction-screening correction; EOS uses a frozen-effective-charge ion-sphere hydrodynamic stress rather than a precision evolving-composition Coulomb model' if x=='no screening or Coulomb EOS correction' else x) for x in d['limitations']]
   d['metadata_scope']='Numerical values unchanged from raw result; inherited generic no-Coulomb limitation corrected to describe the process-local EOS override.'
   out.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');row['archived_record']=str(out);row['archived_record_sha256']=sha(out)
  rows.append(row)
 if missing and not a.allow_partial:raise RuntimeError(f'Missing expected cases: {missing}')
 out={'scope':__doc__,'code_sha256':sha(__file__),'EOS_module_sha256':sha(HERE/'nonideal_eos.py'),'algebra_validation':validate(),'cases':rows,'missing':missing,'all_complete':not missing,'interpretation':'A hydrodynamic free-energy stress at matched total initial energy and physical duration. Frozen charge moment and unchanged nuclear screening prevent interpreting this as a complete physically self-consistent Coulomb plasma calculation.'}
 a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'complete':not missing,'cases':[{'stage':r['stage'],'H':r['H_per_catalyst'],'Q':r['comparison']['deposited_nuclear_over_initial_energy'],'trace_energy_error':r['audit']['saved_trace_energy_residual_max_abs_initial']} for r in rows]},indent=2))
if __name__=='__main__':main()
