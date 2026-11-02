"""Full-network64 C13 response to nonuniform FEOS merger-entropy profile.
Hot-density-profile sensitivity at fixed prepared energy; not a 3D driver simulation.
"""
from pathlib import Path
import argparse,json,sys,hashlib
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent;sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT)]
from radial_completion import simulate,initial
from fast_eos import FastEOS
import pulse_model as pm
ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=2000);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=True)
src=ROOT/'inputs/pulse-profiles/scaled-hot-C13-64-R30-initial.json';entropy=P/'launcher-feos.json';p=json.loads(src.read_text());j=json.loads(entropy.read_text());row=next(r for r in j['stages']['C13']['local_speed_entropy']['fixed_N_profiles'] if r['effective_parallel_directions']==a.N)
m=np.array(p['mass_fraction_of_total']);r=np.r_[0.,p['outer_radius_m']];M=4*np.pi/3*1e8*30**3;rho=M*m/(4*np.pi/3*np.diff(r**3));v=np.array(p['velocity_m_s']);T=np.array(p['T_K']);mapped=T.copy();mapped[row['active_cell_indices']]=row['mapped_T_K'];nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]]
y,target,cats=initial('C13');e=FastEOS();ye=pm.ZZ@y;yi=sum(y);U=float(m@e.state(rho/1000,T,ye,yi).e)*1e-4;K=float(.5*nm@(v*v));Un=float(m@e.state(rho/1000,mapped,ye,yi).e)*1e-4;scale=float(np.sqrt((K+U-Un)/K));assert np.isfinite(scale)
p['T_K']=mapped.tolist();p['velocity_m_s']=(v*scale).tolist();p['description']='Local-speed FEOS contact-entropy mapping, C13N%d, then FastEOS totalpreparedenergy restored by uniformvelocityrescale. Inner0.1%%seed retainsreferenceT. No delivered3Dassembly.'%a.N
stem=a.out/f'C13-N{a.N}-local-entropy-fixed-energy';stem.with_suffix('.initial.json').write_text(json.dumps(p,indent=2)+'\n')
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
inputs={'stage':'C13','effective_angular_regions':a.N,'mass_kg':M,'initial_Y':dict(zip(pm.SPECIES,y.tolist())),'nuclear_network_species':pm.SPECIES,'nuclear_mass_convention':'pynucastro atomic mass-excess table in pulse_model; numberabundances perbaryon. FEOS CH usesC13=13,H=1u inertialmasses. FastEOS ye,yi fromsameintegerA composition.','profile_source':str(src),'profile_sha256':sha(src),'entropy_mapping_sha256':sha(entropy),'calculation_sha256':sha(__file__),'radial_module_sha256':sha(ROOT/'b2-completion/radial_completion.py'),'old_specific_internal_J_kg':U,'old_specific_kinetic_J_kg':K,'new_specific_internal_J_kg':Un,'velocity_scale':scale,'prepared_energy_relative_error':(Un+K*scale**2)/(U+K)-1,'opacity_m2_kg':.001,'reference_time_unit_definition':'R30m dividedbyinitialreferencecs atT3e8K,rho1e8kgm3,originalC13+H composition; unchangedfrombaseline.','scope':'Fixed hotdensity profile; FEOSentropy createsTshape, FastEOSe definesequalinputenergy. Not a selfconsistentcold-tohottrajectory or FEOSburnEOS.'}
stem.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
print(json.dumps({'started':str(stem),'velocity_scale':scale,'oldU':U,'newU':Un,'same_energy_error':inputs['prepared_energy_relative_error']}),flush=True)
res=simulate('C13',3e8,1e8,30,zones=64,nuclear=True,reduced_a=False,radiation=True,kappa=.001,initial_profile=p,final_time_units=4,rtol=1e-7,atol=1e-12,viscosity_mode='compression_cap',save_trace=stem.with_suffix('.npz'),trace_samples=120)
res['launcher_mapped_inputs']=inputs;stem.with_suffix('.json').write_text(json.dumps(res,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:res[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','wall_seconds']}),flush=True)
