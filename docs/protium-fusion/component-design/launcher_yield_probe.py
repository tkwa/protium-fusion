"""Fixed-total-energy hot-profile entropy sensitivity; no delivery-model claim."""
from pathlib import Path
import argparse,json,hashlib,sys,time
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT)]
from radial_completion import simulate,initial
from fast_eos import FastEOS
ap=argparse.ArgumentParser();ap.add_argument('--stage',default='C13',choices=['C13','B2','A']);ap.add_argument('--factor',type=float,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
src={'C13':ROOT/'reassessment-2049/small-stages/scaled-hot-C13-64-R30-initial.json','B2':ROOT/'b2-completion/B2-64-compression-cap.initial.json','A':ROOT/'reassessment-2049/small-stages/scaled-hot-A-64-R30-initial.json'}[a.stage]
p=json.loads(src.read_text());m=np.array(p['mass_fraction_of_total']);r=np.r_[0.,p['outer_radius_m']];rho=1.1309733552923255e13*m/(4*np.pi/3*np.diff(r**3));v=np.array(p['velocity_m_s']);T=np.array(p['T_K']);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
y,_,_=initial(a.stage)
import pulse_model as pm
e=FastEOS();ye=pm.ZZ@y;yi=sum(y);U=float(m@e.state(rho/1000,T,ye,yi).e)*1e-4;K=float(.5*nm@(v*v));Unew=float(m@e.state(rho/1000,T*a.factor,ye,yi).e)*1e-4
if Unew>=K+U:raise ValueError('No kinetic energy left at requested thermal state')
scale=np.sqrt((K+U-Unew)/K);p['T_K']=(T*a.factor).tolist();p['velocity_m_s']=(v*scale).tolist();p['description']='Uniform hot-profile temperature perturbation; velocity scale restorestotalpreparedenergy. No cold-to-hot preparation claim.'
a.out.mkdir(parents=True,exist_ok=True);tag=format(a.factor,'.4g').replace('.','p');stem=a.out/f'{a.stage}-Tfactor{tag}';stem.with_suffix('.initial.json').write_text(json.dumps(p,indent=2)+'\n')
inputs={'stage':a.stage,'temperature_factor':a.factor,'kinetic_velocity_scale':float(scale),'old_specific_internal_J_kg':U,'old_specific_kinetic_J_kg':K,'new_specific_internal_J_kg':Unew,'transferred_kinetic_to_internal_J_kg':Unew-U,'profile_source':str(src),'profile_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'Uniform initialhotentropyprofile sensitivity atfixedtotalpreparedenergy; not a modelof whereorwhenmergerheat isdeposited.'}
stem.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
result=simulate(a.stage,3e8,1e8,30,zones=64,nuclear=True,reduced_a=a.stage=='A',radiation=True,initial_profile=p,final_time_units=4,rtol=1e-7,atol=1e-12,viscosity_mode='compression_cap',save_trace=stem.with_suffix('.npz'),trace_samples=80)
result['launcher_sensitivity']=inputs;stem.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','wall_seconds']}),flush=True)
