"""B2 H2 uniform contact-shock entropy from controlled equal normal collision speed.
Uses the archived CH2 FEOS material and median carbon-rate update. The generated
hot profile retains radial density and total energy; physical3Ddelivery remains open.
"""
from pathlib import Path
import argparse,ctypes as c,os,json,math,hashlib,sys
import numpy as np
from scipy.optimize import root,brentq
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT/'mixture-design'),str(ROOT)]
import radial_completion as radial
from fast_eos import FastEOS
import pulse_model as pm
import updated_network
ap=argparse.ArgumentParser();ap.add_argument('--library-cache',type=Path,required=True);ap.add_argument('--target-factor',type=float,default=1.);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out=a.out.resolve();a.out.parent.mkdir(parents=True,exist_ok=True)
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
lib=a.library_cache.resolve()/'FEOS/Code/libfeos.so';db=a.library_cache.resolve()/'FEOS/EOS-Data/FEOS_Material-DB.dat';mapping=json.loads((P.parent/'launcher-feos.json').read_text());assert sha(db)==mapping['material_database_sha256'];os.chdir(db.parent)
l=c.CDLL(str(lib));D=c.c_double;I=c.c_int;DP=c.POINTER(D);IP=c.POINTER(I);l.FEOS_Initialize.argtypes=[I,I];l.FEOS_Init_Mat.argtypes=[I,I,I,I,DP,I,IP,IP,DP];l.FEOS_Get_EOS.argtypes=[I,I,I,D,D,IP,DP,DP,DP,DP,DP,DP];l.FEOS_Initialize(2,0)
for j,n in [(1,9001),(2,9002)]:
 ne=I();ni=I();tt=D();l.FEOS_Init_Mat(j,n,0,0,None,0,c.byref(ne),c.byref(ni),c.byref(tt))
AR=7.5657e-16

def st(rho,T):
 b=I();p=D();e=D();s=D();f=D();q=(D*2)();qt=D();l.FEOS_Get_EOS(1,0,0,rho/1000,T/11604.51812,c.byref(b),c.byref(p),c.byref(e),c.byref(s),c.byref(f),q,c.byref(qt));return dict(p=p.value*.1+AR*T**4/3,e=e.value*1e-4+AR*T**4/rho,s=s.value*1e-4/11604.51812+4*AR*T**3/(3*rho),q=list(q))
def shock(q):
 rho0=944.2;s0=st(rho0,300.)
 def f(x):
  rho,T=np.exp(x);s=st(rho,T);return [(s['e']-s0['e']-q)/q,(.5*(s['p']+s0['p'])*(1/rho0-1/rho)-q)/q]
 r=root(f,np.log([rho0*4,max(1e4,q/8000)]));rho,T=np.exp(r.x);error=float(max(abs(np.array(f(r.x)))));assert error<1e-5 and rho>rho0
 return dict(rho_kg_m3=float(rho),T_K=float(T),**st(rho,T),q_J_kg=float(q),error=error)
def quant(x,w):
 x=np.asarray(x);i=np.argsort(x);return np.interp([.1,.5,.9],np.cumsum(w[i])-.5*w[i],x[i]).tolist()
ref=Path('/home/tkwa/src/protium-fusion/output/mixture-runs/grid/B2-H2-fixed_energy-z64-k0p001-t1');baseline=json.loads(ref.with_suffix('.json').read_text());bi=json.loads(ref.with_suffix('.inputs.json').read_text());profile=json.loads(ref.with_suffix('.initial.json').read_text());m=np.array(profile['mass_fraction_of_total']);r=np.r_[0,profile['outer_radius_m']];rho=baseline['mass_kg']*m/(4*np.pi/3*np.diff(r**3));T=np.array(profile['T_K']);v=np.array(profile['velocity_m_s']);active=np.where(np.cumsum(m)-.5*m>=.001)[0];w=m[active]/sum(m[active]);targetS=float(w@np.array([st(rho[i],T[i]*a.target_factor)['s'] for i in active]));q=math.exp(brentq(lambda x:shock(math.exp(x))['s']-targetS,math.log(1e7),math.log(1e11)));contact=shock(q);mapped=T.copy()
for i in active:mapped[i]=math.exp(brentq(lambda x:st(rho[i],math.exp(x))['s']-contact['s'],math.log(10),math.log(3e9)))
y=np.array([bi['initial_Y'].get(s,0.) for s in pm.SPECIES]);_,target,_=radial.initial('B2');ci=pm.SPECIES.index('c12');radial.initial=lambda stage:(y.copy(),target,float(y[ci]));rate_update=updated_network.install(radial,quantile='median')
eos=FastEOS();ye=float(pm.ZZ@y);yi=float(sum(y));oldU=eos.state(rho/1000,T,ye,yi).e*1e-4;newU=eos.state(rho/1000,mapped,ye,yi).e*1e-4;U=float(m@oldU);Un=float(m@newU);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];K=float(.5*nm@(v*v));scale=float(np.sqrt((K+U-Un)/K));nv=v*scale;localK=.25*(nv*nv+np.r_[0.,nv[:-1]**2]);energy=newU+localK;assert max(q/energy[active])<1;angles=np.degrees(np.arcsin(np.sqrt(q/energy[active])));profile['T_K']=mapped.tolist();profile['velocity_m_s']=nv.tolist();profile['description']='B2H2 constant contact-shock entropy, from common normal collision speed; same radial densities and total prepared energy; median carbon rates.'
center=eos.state(baseline['rho0_kg_m3']/1000,baseline['T0_K'],ye,yi);cs=float(np.sqrt(center.gamma1*center.p/(baseline['rho0_kg_m3']/1000)));tu=baseline['R0_m']*100/cs;duration=bi['duration_s']
inputs=dict(stage='B2',hydrogen_per_catalyst=2,initial_Y=bi['initial_Y'],duration_s=duration,physical_duration_s=duration,contact=contact,target_reference_T_factor=a.target_factor,constant_contact_entropy_J_kg_K=targetS,temperature_ratio_mass_q10_q50_q90=quant(mapped[active]/T[active],w),temperature_ratio_minmax=[float(min(mapped[active]/T[active])),float(max(mapped[active]/T[active]))],old_specific_internal_J_kg=U,new_specific_internal_J_kg=Un,old_specific_kinetic_J_kg=K,velocity_scale=scale,prepared_energy_relative_error=(Un+K*scale**2)/(U+K)-1,contact_heat_J=q*baseline['mass_kg']*sum(m[active]),contact_heat_fraction_prepared=q*sum(m[active])/(U+K),transverse_collision_speed_each_stream_m_s=math.sqrt(2*q),pair_half_angle_degrees_mass_q10_q50_q90=quant(angles,w),pair_half_angle_degrees_minmax=[float(min(angles)),float(max(angles))],postshock_rarefaction_mass_fraction=float(sum(m[active][rho[active]<contact['rho_kg_m3']])),cold_reference=dict(composition='C12H2',rho_kg_m3=944.2,T_K=300,bulk_modulus_Pa=2e9),rate_update=rate_update,FEOS_library_sha256=sha(lib),FEOS_material_database_sha256=sha(db),profile_source=str(ref.with_suffix('.initial.json')),profile_sha256=sha(ref.with_suffix('.initial.json')),baseline_source=str(ref.with_suffix('.json')),baseline_sha256=sha(ref.with_suffix('.json')),code_sha256=sha(__file__),network_source_sha256=sha(ROOT/'pulse_model.py'),radial_module_sha256=sha(ROOT/'b2-completion/radial_completion.py'),scope='Controlled collision-energy geometry, not temperature clipping. FEOSdefinesonecontactshockentropy andeachparcelmapsisentropicallytooriginalhotdensity. Anglesuseanexplicitlocalcoldenergyallocationequaltohotinternal+radialKE; actualpressurework/3Dgeometryremainsunresolved. Noexternalultrafastheateroradditionalinputenergy.')
a.out.with_suffix('.initial.json').write_text(json.dumps(profile,indent=2)+'\n');a.out.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n');print(json.dumps({k:inputs[k] for k in ['contact','temperature_ratio_mass_q10_q50_q90','velocity_scale','contact_heat_fraction_prepared','transverse_collision_speed_each_stream_m_s','pair_half_angle_degrees_mass_q10_q50_q90']}),flush=True)
res=radial.simulate('B2',baseline['T0_K'],baseline['rho0_kg_m3'],baseline['R0_m'],zones=64,nuclear=True,radiation=True,kappa=.001,initial_profile=profile,final_time_units=duration/tu,rtol=1e-7,atol=1e-12,viscosity_mode='compression_cap',save_trace=a.out.with_suffix('.npz'),trace_samples=400)
res['initial_Y']=bi['initial_Y'];res['mixture_inputs']=inputs;res['rate_update']=rate_update;res['trace_path']=str(a.out.with_suffix('.npz'));a.out.with_suffix('.json').write_text(json.dumps(res,indent=2,allow_nan=False)+'\n');print(json.dumps({k:res[k] for k in ['success','termination_reason','elapsed_s','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','peak_temperature_K','energy_residual_fraction_initial','wall_seconds']}),flush=True)
