"""Nested 64-to-128 source profile sensitivity; no contact-energy retuning.
Preserves each piecewise-constant density/temperature field by enclosed mass,
exact total pulse mass and prepared energy, composition and physical duration.
"""
from pathlib import Path
import argparse,hashlib,json,sys,math
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT/'mixture-design'),str(ROOT)]
import radial_completion as radial
import pulse_model as pm
from fast_eos import FastEOS
import updated_network
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--cat',required=True,choices=['n14','n15','c12','c13']);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--prepare-only',action='store_true');a=ap.parse_args()
 mapping=ROOT/'component-design/source-family-carbon-median/stage-map.json';spec=json.loads(mapping.read_text())['stage_map'][a.cat]
 src=ROOT/spec['record'];init=ROOT/spec['initial'];infile=ROOT/spec['inputs'];b=json.loads(src.read_text());p=json.loads(init.read_text());inp=json.loads(infile.read_text())
 for f,key in [(src,'result_sha256'),(init,'initial_profile_sha256'),(infile,'input_record_sha256')]:assert sha(f)==spec[key]
 assert b['zones']==64 and b['success'] and b['termination_reason']=='final_time' and b['viscosity_mode']=='compression_cap'
 m0=np.asarray(p['mass_fraction_of_total']);q0=np.r_[0.,np.cumsum(m0)];q0[-1]=1.;r0=np.r_[0.,p['outer_radius_m']];v0=np.r_[0.,p['velocity_m_s']];T0=np.asarray(p['T_K']);mass=b['mass_kg']
 # Split each existing cell at the midpoint in its original mass-radius label.
 # This is exactly the linspace^3 grid for the saved64 family and remains nested.
 mid=(.5*(q0[:-1]**(1/3)+q0[1:]**(1/3)))**3
 q=np.sort(np.r_[q0,mid]);m=np.diff(q);r=np.interp(q,q0,r0**3)**(1/3);r[0]=0.;parent=np.repeat(np.arange(64),2)
 T=T0[parent];v=np.interp(r,r0,v0)[1:];rho=mass*m/(4*np.pi/3*np.diff(r**3));rho0=mass*m0/(4*np.pi/3*np.diff(r0**3))
 y=np.array([inp['initial_Y'].get(n,0.) for n in pm.SPECIES]);assert abs(pm.AA@y-1)<1e-12
 eos=FastEOS();ye=float(pm.ZZ@y);yi=float(sum(y));u=eos.state(rho/1000,T,ye,yi).e*1e-4;u0=eos.state(rho0/1000,T0,ye,yi).e*1e-4
 U=mass*float(m@u);U0=mass*float(m0@u0);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];Kraw=.5*mass*float(nm@(v*v));target=b['initial_energy_J'];scale=math.sqrt((target-U)/Kraw);v*=scale;K=.5*mass*float(nm@(v*v))
 errors={'density_max_relative':float(max(abs(rho/rho0[parent]-1))),'temperature_max_relative':float(max(abs(T/T0[parent]-1))),'internal_energy_relative_to64':float(U/U0-1),'total_energy_relative_to64':float((U+K)/target-1),'mass_relative_to64':float(sum(m)-1),'kinetic_energy_relative_to64':float(K/b['initial_kinetic_energy_J']-1)}
 assert errors['density_max_relative']<1e-10 and max(abs(x) for x in errors.values())<1e-10
 profile={'outer_radius_m':r[1:].tolist(),'velocity_m_s':v.tolist(),'T_K':T.tolist(),'mass_fraction_of_total':m.tolist(),'description':'Nested128 remap of exact completed64 source profile. Piecewise rho/T fixed by mass; globally scaled interpolated radial velocity restores original total energy. No contact entropy or recipe retuning.'}
 duration=spec['duration_s'];_,desired,_=radial.initial(b['stage']);catY=float(y[pm.SPECIES.index(a.cat)]);radial.initial=lambda stage:(y.copy(),desired,catY)
 rate=updated_network.install(radial,quantile='median');assert rate['kernel_sha256']==b['rate_update']['kernel_sha256'];assert rate['table_sha256']==b['rate_update']['table_sha256']
 newinp=dict(inp);newinp.update(zones=128,duration_s=duration,physical_duration_s=duration,rate_update=rate,initial_Y=inp['initial_Y'],source_refinement={'source_record':str(src),'source_record_sha256':sha(src),'source_initial':str(init),'source_initial_sha256':sha(init),'source_inputs':str(infile),'source_inputs_sha256':sha(infile),'stage_map_sha256':sha(mapping),'runner_sha256':sha(__file__),'parent_cell_for_each_new_cell':parent.tolist(),'velocity_scale_after_interpolation':scale,'initial_energy_J':U+K,'initial_internal_J':U,'initial_kinetic_J':K,'remap_errors':errors,'no_source_retuning':True,'scope':__doc__})
 a.out.parent.mkdir(parents=True,exist_ok=True);write(a.out.with_suffix('.initial.json'),profile);write(a.out.with_suffix('.inputs.json'),newinp)
 print(json.dumps({'prepared':str(a.out),'velocity_scale':scale,'remap_errors':errors}),flush=True)
 if a.prepare_only:return
 center=eos.state(b['rho0_kg_m3']/1000,b['T0_K'],ye,yi);tu=b['R0_m']*100/float(np.sqrt(center.gamma1*center.p/(b['rho0_kg_m3']/1000)))
 result=radial.simulate(b['stage'],b['T0_K'],b['rho0_kg_m3'],b['R0_m'],zones=128,nuclear=True,reduced_a=False,radiation=b['radiation_diffusion'],kappa=b['opacity_m2_kg'],initial_profile=profile,final_time_units=duration/tu,rtol=b['rtol'],atol=b['atol'],viscosity_mode='compression_cap',save_trace=a.out.with_suffix('.npz'),trace_samples=400)
 result.update(initial_Y=inp['initial_Y'],mixture_inputs=newinp,rate_update=rate,trace_path=str(a.out.with_suffix('.npz')));write(a.out.with_suffix('.json'),result)
 print(json.dumps({'case':str(a.out),**{k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','peak_temperature_K','wall_seconds']}}),flush=True)
if __name__=='__main__':main()
