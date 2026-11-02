"""B1 constructive constant transverse-collision energy, followed by mapped hot flow.
One shock entropy is imposed by a common normal collision speed, then each parcel
isentropically reaches its prescribed radial density. Outer parcels can rarefy.
This tests an entropy-shaped hot initial condition, not a delivered 3D assembly.
"""
from pathlib import Path
import argparse,ctypes as c,os,json,math,hashlib,sys
import numpy as np
from scipy.optimize import root,brentq
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT)]
import radial_completion as radial
from fast_eos import FastEOS
import pulse_model as pm
ap=argparse.ArgumentParser();ap.add_argument('--map',type=Path,default=P/'entropy-map.json');ap.add_argument('--case',default='H2-rho1062p5');ap.add_argument('--target-factor',type=float,default=.9);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--library-cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--no-burn',action='store_true');a=ap.parse_args();a.map=a.map.resolve();a.out=a.out.resolve();a.out.parent.mkdir(parents=True,exist_ok=True)
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
d=json.loads(a.map.read_text())['cases'][a.case];recipe=d['recipe'];j=recipe['entity'];lib=a.library_cache.resolve()/'FEOS/Code/libfeos.so';os.chdir(a.cache.resolve())
l=c.CDLL(str(lib));D=c.c_double;I=c.c_int;DP=c.POINTER(D);IP=c.POINTER(I)
l.FEOS_Initialize.argtypes=[I,I];l.FEOS_Init_Mat.argtypes=[I,I,I,I,DP,I,IP,IP,DP];l.FEOS_Get_EOS.argtypes=[I,I,I,D,D,IP,DP,DP,DP,DP,DP,DP];l.FEOS_Initialize(4,0)
for jj in range(1,5):
 ne=I();ni=I();tt=D();l.FEOS_Init_Mat(jj,9200+jj,0,0,None,0,c.byref(ne),c.byref(ni),c.byref(tt))
AR=7.5657e-16
def st(rho,T):
 b=I();p=D();e=D();s=D();f=D();q=(D*2)();qt=D();l.FEOS_Get_EOS(j,0,0,rho/1000,T/11604.51812,c.byref(b),c.byref(p),c.byref(e),c.byref(s),c.byref(f),q,c.byref(qt));return dict(p=p.value*.1+AR*T**4/3,e=e.value*1e-4+AR*T**4/rho,s=s.value*1e-4/11604.51812+4*AR*T**3/(3*rho),q=list(q))
def shock(q):
 rho0=recipe['rho'];s0=st(rho0,recipe['T'])
 def f(x):
  rho,T=np.exp(x);s=st(rho,T);return [(s['e']-s0['e']-q)/q,(.5*(s['p']+s0['p'])*(1/rho0-1/rho)-q)/q]
 rr=root(f,np.log([rho0*4,max(1e3,q/1e4)]));rho,T=np.exp(rr.x);err=float(max(abs(np.array(f(rr.x)))))
 if err>1e-5 or rho<=rho0:raise ValueError('contact shock root failed')
 return dict(rho_kg_m3=float(rho),T_K=float(T),**st(rho,T),q_J_kg=float(q),error=err)
def quant(x,w):
 x=np.asarray(x);o=np.argsort(x);return np.interp([.1,.5,.9],np.cumsum(w[o])-.5*w[o],x[o]).tolist()
src=Path(d['reference_profile']);ref=Path(d['reference_metadata']);baseline=json.loads(ref.read_text());profile=json.loads(src.read_text());m=np.array(profile['mass_fraction_of_total']);r=np.r_[0,profile['outer_radius_m']];rho=baseline['mass_kg']*m/(4*np.pi/3*np.diff(r**3));T=np.array(profile['T_K']);v=np.array(profile['velocity_m_s']);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];active=np.array([x['index'] for x in d['required_cells']]);w=m[active]/sum(m[active]);qref=np.array([x['required_shock_q_J_kg_by_target_Tfactor'][str(a.target_factor)] for x in d['required_cells']]);q=float(w@qref);contact=shock(q)
mapped=T.copy()
for idx in active:mapped[idx]=math.exp(brentq(lambda x:st(rho[idx],math.exp(x))['s']-contact['s'],math.log(1),math.log(3e9)))
y=np.array([d['initial_Y'].get(s,0.) for s in pm.SPECIES]);_,target,_=radial.initial('B1');ci=pm.SPECIES.index('n15');radial.initial=lambda stage:(y.copy(),target,float(y[ci]));assert abs(pm.AA@y-1)<1e-14
eos=FastEOS();ye=float(pm.ZZ@y);yi=float(sum(y));oldU=eos.state(rho/1000,T,ye,yi).e*1e-4;newU=eos.state(rho/1000,mapped,ye,yi).e*1e-4;U=float(m@oldU);Un=float(m@newU);K=float(.5*nm@(v*v));scale=float(np.sqrt((K+U-Un)/K));newv=v*scale;localK=.25*(newv*newv+np.r_[0.,newv[:-1]**2]);energy=localK+newU
# Local cold energy allocation is a preparation diagnostic, not an advected invariant:
# pressure work can move energy between parcels in a physical assembly trajectory.
fractions=q/energy[active];assert np.max(fractions)<1
angles=np.degrees(np.arcsin(np.sqrt(fractions)));cold_speed=np.sqrt(2*energy[active]);radial_speed=np.sqrt(2*(energy[active]-q))
profile['T_K']=mapped.tolist();profile['velocity_m_s']=newv.tolist();profile['description']='B1 constant FEOS contact-shock entropy; nonseed parcels rarefy/compress isentropically to original radial density; FastEOS matched total prepared energy; no external heater.'
center=eos.state(baseline['rho0_kg_m3']/1000,baseline['T0_K'],ye,yi);cs=float(np.sqrt(center.gamma1*center.p/(baseline['rho0_kg_m3']/1000)));tunit=baseline['R0_m']*100/cs;duration=d['physical_duration_s']
inputs=dict(stage='B1',case=a.case,initial_Y=d['initial_Y'],physical_duration_s=duration,duration_s=duration,entropy_mapping_source=str(a.map),entropy_mapping_sha256=sha(a.map),profile_source=str(src),profile_sha256=sha(src),baseline_source=str(ref),baseline_sha256=sha(ref),code_sha256=sha(__file__),radial_module_sha256=sha(ROOT/'b2-completion/radial_completion.py'),network_source_sha256=sha(ROOT/'pulse_model.py'),FEOS_library_sha256=sha(lib),FEOS_material_reference=recipe,contact=contact,target_factor_used_to_choose_one_contact_q=a.target_factor,temperature_ratio_mass_q10_q50_q90=quant(mapped[active]/T[active],w),temperature_ratio_minmax=[float(min(mapped[active]/T[active])),float(max(mapped[active]/T[active]))],old_specific_internal_J_kg=U,old_specific_kinetic_J_kg=K,new_specific_internal_J_kg=Un,velocity_scale=scale,prepared_energy_relative_error=(Un+K*scale**2)/(U+K)-1,contact_heat_J=q*baseline['mass_kg']*sum(m[active]),contact_heat_fraction_prepared=q*sum(m[active])/(K+U),transverse_collision_speed_each_stream_m_s=math.sqrt(2*q),relative_transverse_pair_speed_m_s=2*math.sqrt(2*q),parcel_energy_allocation_scope='Cold scalar allocation E_i=final hot specific internal+radial kinetic, zero cold binding reference; summed energy matches prepared input. Angle is a diagnostic for symmetric equal-energy colliding streams. Real pressure-work/geometry can redistribute local energy; no direct3Dlaunch solution claimed.',pair_half_angle_degrees_mass_q10_q50_q90=quant(angles,w),pair_half_angle_degrees_minmax=[float(min(angles)),float(max(angles))],cold_total_speed_m_s_mass_q10_q50_q90=quant(cold_speed,w),cold_post_transverse_radial_speed_m_s_mass_q10_q50_q90=quant(radial_speed,w),postshock_rarefaction_mass_fraction=float(sum(m[active][rho[active]<contact['rho_kg_m3']])),cells=[dict(index=int(idx),initial_T_K=float(mapped[idx]),cold_energy_allocation_J_kg=float(energy[idx]),half_angle_degrees=float(angle),contact_q_J_kg=q) for idx,angle in zip(active,angles)],time_unit_s=tunit,opacity_m2_kg=.001,scope='Full27network64 originalrates; actualB1density; fixedpreparedenergyandphysicalduration. The contact heating comes from transverse packet KE, not an added externalpreheater. Constant shock entropy is a controlledgeometrycandidate; coldgrain3Dmerger and sourcewaveform remain unqualified.')
a.out.with_suffix('.initial.json').write_text(json.dumps(profile,indent=2)+'\n');a.out.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n');print(json.dumps({k:inputs[k] for k in ['contact','temperature_ratio_mass_q10_q50_q90','velocity_scale','contact_heat_fraction_prepared','transverse_collision_speed_each_stream_m_s','pair_half_angle_degrees_mass_q10_q50_q90','postshock_rarefaction_mass_fraction']}),flush=True)
if not a.no_burn:
 result=radial.simulate('B1',baseline['T0_K'],baseline['rho0_kg_m3'],baseline['R0_m'],zones=64,nuclear=True,reduced_a=False,radiation=True,kappa=.001,initial_profile=profile,final_time_units=duration/tunit,rtol=1e-7,atol=1e-12,viscosity_mode='compression_cap',save_trace=a.out.with_suffix('.npz'),trace_samples=400)
 result['initial_Y']=d['initial_Y'];result['launcher_mixture_inputs']=inputs;a.out.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','wall_seconds']}),flush=True)
