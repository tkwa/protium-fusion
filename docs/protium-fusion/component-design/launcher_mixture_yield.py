"""Full-network64 response to selected-recipe shell-dependent merger entropy.
Each process overrides only initialcomposition; physicalduration andenergy matchitsbaseline.
"""
from pathlib import Path
import argparse,json,sys,hashlib
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent;sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT)]
import radial_completion as radial
from fast_eos import FastEOS
import pulse_model as pm
ap=argparse.ArgumentParser();ap.add_argument('--case',choices=['C13H4','AH3'],required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=True)
mapfile=P/'launcher-mixture-entropy.json';d=json.loads(mapfile.read_text())['cases'][a.case];stage=d['stage'];src=Path(d['profile_source']);ip=Path(d['inputs_source']);baseinputs=json.loads(ip.read_text());p=json.loads(src.read_text());basepath=src.with_name(src.name.replace('.initial.json','.json'));baseline=json.loads(basepath.read_text());y=np.array([d['initial_Y'].get(s,0.) for s in pm.SPECIES]);_,target,_=radial.initial(stage);cat='c13' if stage=='C13' else'n14';catidx=pm.SPECIES.index(cat);radial.initial=lambda stage:(y.copy(),target,float(y[catidx]));assert abs(pm.AA@y-1)<1e-14
m=np.array(p['mass_fraction_of_total']);r=np.r_[0,p['outer_radius_m']];rho=baseline['mass_kg']*m/(4*np.pi/3*np.diff(r**3));v=np.array(p['velocity_m_s']);T=np.array(p['T_K']);mapped=T.copy()
for cell in d['cells']:mapped[cell['index']]=cell['mapped_T_K']
nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];e=FastEOS();ye=float(pm.ZZ@y);yi=float(sum(y));U=float(m@e.state(rho/1000,T,ye,yi).e)*1e-4;K=float(.5*nm@(v*v));Un=float(m@e.state(rho/1000,mapped,ye,yi).e)*1e-4;scale=float(np.sqrt((K+U-Un)/K));assert np.isfinite(scale)
p['T_K']=mapped.tolist();p['velocity_m_s']=(v*scale).tolist();p['description']=f'{a.case} N2000localFEOSentropyprofile, FastEOSequalpreparedenergy, originaldensity, inner0.1%seedT unchanged.'
center=e.state(baseline['rho0_kg_m3']/1000,baseline['T0_K'],ye,yi);cs=float(np.sqrt(center.gamma1*center.p/(baseline['rho0_kg_m3']/1000)));tunit=baseline['R0_m']*100/cs;duration=d['physical_duration_s'];stem=a.out/f'{a.case}-N2000-local-entropy-fixed-energy'
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
inputs={'stage':stage,'case':a.case,'initial_Y':d['initial_Y'],'entropy_mapping_source':str(mapfile),'entropy_mapping_sha256':sha(mapfile),'profile_source':str(src),'profile_sha256':sha(src),'baseline_source':str(basepath),'baseline_sha256':sha(basepath),'code_sha256':sha(__file__),'radial_module_sha256':sha(ROOT/'b2-completion/radial_completion.py'),'old_specific_internal_J_kg':U,'old_specific_kinetic_J_kg':K,'new_specific_internal_J_kg':Un,'velocity_scale':scale,'prepared_energy_relative_error':(Un+K*scale**2)/(U+K)-1,'physical_duration_s':duration,'time_unit_s':tunit,'opacity_m2_kg':.001,'FEOS_material_reference':d['cold_reference'],'scope':'Full27network64prescribedhotdensityresponse; no3Dassembly. FEOSentropy mapsT, FastEOSsetsenergy/burnEOS; C13usesintegerisotopemasses13+4H, A14+3H.'}
stem.with_suffix('.initial.json').write_text(json.dumps(p,indent=2)+'\n');stem.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n');print(json.dumps({'started':a.case,'velocity_scale':scale,'Uold':U,'Unew':Un,'duration_s':duration}),flush=True)
res=radial.simulate(stage,baseline['T0_K'],baseline['rho0_kg_m3'],baseline['R0_m'],zones=64,nuclear=True,reduced_a=False,radiation=True,kappa=.001,initial_profile=p,final_time_units=duration/tunit,rtol=1e-7,atol=1e-12,viscosity_mode='compression_cap',save_trace=stem.with_suffix('.npz'),trace_samples=120)
res['initial_Y']=d['initial_Y'];res['launcher_mixture_inputs']=inputs;stem.with_suffix('.json').write_text(json.dumps(res,indent=2,allow_nan=False)+'\n');print(json.dumps({k:res[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','wall_seconds']}),flush=True)
