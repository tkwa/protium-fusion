"""B1 H1.5 controlled contact source at128; carbon or low-N15 rate control.
Cold physical payload is separate hydrazine and nitrogen, with75%/25% of N
atoms respectively; the original homogenized FEOS contact surrogate is retained.
"""
from pathlib import Path
import argparse,hashlib,json,sys,math
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT/'b2-completion'),str(ROOT/'mixture-design'),str(ROOT)]
import radial_completion as radial
from fast_eos import FastEOS
import pulse_model as pm
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def wr(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--curve',choices=['carbon','low'],required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 stem=Path('/home/tkwa/src/protium-fusion/output/b1-source-robustness/B1-H1p5-contact64');inp=json.loads(stem.with_suffix('.inputs.json').read_text());p=json.loads(stem.with_suffix('.initial.json').read_text());base=Path(inp['baseline_source']);b=json.loads(base.read_text());y=np.array([inp['initial_Y'].get(n,0.) for n in pm.SPECIES]);ye=float(pm.ZZ@y);yi=float(sum(y));eos=FastEOS()
 m0=np.asarray(p['mass_fraction_of_total']);q0=np.r_[0.,np.cumsum(m0)];q0[-1]=1.;r0=np.r_[0.,p['outer_radius_m']];v0=np.r_[0.,p['velocity_m_s']];T0=np.asarray(p['T_K']);mid=(.5*(q0[:-1]**(1/3)+q0[1:]**(1/3)))**3;q=np.sort(np.r_[q0,mid]);m=np.diff(q);r=np.interp(q,q0,r0**3)**(1/3);r[0]=0.;parent=np.repeat(np.arange(64),2);T=T0[parent];v=np.interp(r,r0,v0)[1:];rho=b['mass_kg']*m/(4*np.pi/3*np.diff(r**3));rho0=b['mass_kg']*m0/(4*np.pi/3*np.diff(r0**3));u=eos.state(rho/1000,T,ye,yi).e*1e-4;U=b['mass_kg']*float(m@u);nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];K0=.5*b['mass_kg']*float(nm@(v*v));scale=math.sqrt((b['initial_energy_J']-U)/K0);v*=scale;K=.5*b['mass_kg']*float(nm@(v*v));errors={'rho_max_relative':float(max(abs(rho/rho0[parent]-1))),'T_max_relative':float(max(abs(T/T0[parent]-1))),'prepared_energy_relative':float((U+K)/b['initial_energy_J']-1)};assert max(abs(x) for x in errors.values())<1e-10
 profile={'outer_radius_m':r[1:].tolist(),'velocity_m_s':v.tolist(),'T_K':T.tolist(),'mass_fraction_of_total':m.tolist(),'description':'B1H1.5 constant-contact source, nested128 remap with matched total prepared energy; no source retuning between rate controls.'}
 _,target,_=radial.initial('B1');radial.initial=lambda stage:(y.copy(),target,float(y[pm.SPECIES.index('n15')]))
 if a.curve=='carbon':
  import updated_network as update
  rate=update.install(radial,quantile='median')
 else:
  import updated_network_n15 as update
  rate=update.install(radial,quantile='median',n15_curve='low')
 inputs=dict(inp);inputs.update(hydrogen_per_catalyst=1.5,zones=128,rate_update=rate,source_recipe={'N_atoms_from_hydrazine':.75,'N_atoms_from_N2':.25,'hydrazine_payload_mass_fraction':17/22,'N2_payload_mass_fraction':5/22,'isotope_convention':'N15 and protium; molecular mass34 forN2H4 and30 forN2. Atom-normalized mixture isH/N1.5.','phase_scope':'Separate condensed hydrazine and sealed cryogenicN2 packet populations; not a room-temperature pureN2H3 compound. FEOSrho1000,T250K is a homogenized contact surrogate, not the physical common starting phase. Actual capsule/support/IR and homogenization require design.'},refinement={'parent_profile':str(stem.with_suffix('.initial.json')),'parent_profile_sha256':sha(stem.with_suffix('.initial.json')),'parent_inputs_sha256':sha(stem.with_suffix('.inputs.json')),'runner_sha256':sha(__file__),'velocity_scale_after_interpolation':scale,'errors':errors,'source_contact_not_retuned':True},duration_s=inp['physical_duration_s'])
 a.out.parent.mkdir(parents=True,exist_ok=True);wr(a.out.with_suffix('.initial.json'),profile);wr(a.out.with_suffix('.inputs.json'),inputs);print(json.dumps({'prepared':str(a.out),'curve':a.curve,'remap':errors,'contact_q_J_kg':inp['contact']['q_J_kg']}),flush=True)
 center=eos.state(b['rho0_kg_m3']/1000,b['T0_K'],ye,yi);tu=b['R0_m']*100/float(np.sqrt(center.gamma1*center.p/(b['rho0_kg_m3']/1000)))
 result=radial.simulate('B1',b['T0_K'],b['rho0_kg_m3'],b['R0_m'],zones=128,nuclear=True,reduced_a=False,radiation=True,kappa=.001,initial_profile=profile,final_time_units=inputs['duration_s']/tu,rtol=1e-7,atol=1e-12,viscosity_mode='compression_cap',save_trace=a.out.with_suffix('.npz'),trace_samples=400)
 result.update(initial_Y=inp['initial_Y'],mixture_inputs=inputs,rate_update=rate,trace_path=str(a.out.with_suffix('.npz')));wr(a.out.with_suffix('.json'),result)
 print(json.dumps({'curve':a.curve,**{k:result[k] for k in ['success','termination_reason','deposited_nuclear_over_initial_energy','desired_product_per_initial_catalyst','energy_residual_fraction_initial','wall_seconds']}}),flush=True)
if __name__=='__main__':main()
