"""Thirty-microsecond B1 whole-fluid preparation sensitivity from one precursor.
Compare FastEOS controls with and without grey radiation; nuclear and artificial
viscosity disabled. This tests nonlinear transport, not a cold-grain contact path.
"""
from pathlib import Path
import argparse,json,sys,math,hashlib,time
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent.parent
sys.path[:0]=[str(ROOT),str(ROOT/'b2-completion')]
from initial_screen import load,sha
from fast_eos import FastEOS
import pulse_model as pm
import radial_completion as radial

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--opacity',type=float,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
 source=Path('/home/tkwa/src/protium-fusion/output/assembly-audit/B1-backward-feos128.npz');d=np.load(source);duration=30e-6;t=d['time_s'];j=int(np.argmin(abs(t-duration)));assert abs(t[j]-duration)<1e-14
 original,inp,base,*_=load('B1');m=d['mass_kg'];profile={'outer_radius_m':d['radius_m'][j,1:].tolist(),'velocity_m_s':(-d['reverse_velocity_m_s'][j,1:]).tolist(),'T_K':d['temperature_K'][j].tolist(),'mass_fraction_of_total':(m/sum(m)).tolist(),'description':'Same whole-fluid precursor30microseconds before target source hot state from reversible FEOS continuation; forward FastEOS radiation sensitivity. Contact boundaries not resolved.'}
 y=np.array([inp['initial_Y'].get(n,0.) for n in pm.SPECIES]);_,desired,_=radial.initial('B1');cat=float(y[pm.SPECIES.index('n15')]);radial.initial=lambda stage:(y.copy(),desired,cat)
 eos=FastEOS();ye=float(pm.ZZ@y);yi=float(sum(y));center=eos.state(base['rho0_kg_m3']/1000,base['T0_K'],ye,yi);tu=base['R0_m']*100/math.sqrt(center.gamma1*center.p/(base['rho0_kg_m3']/1000))
 inputs={'scope':__doc__,'source_trace':str(source),'source_trace_sha256':sha(source),'precursor_index':j,'duration_s':duration,'opacity_m2_kg':a.opacity,'radiation':a.opacity>0,'viscosity_linear':0.,'viscosity_quadratic':0.,'nuclear':False,'initial_Y':inp['initial_Y'],'code_sha256':sha(__file__),'radial_code_sha256':sha(ROOT/'b2-completion/radial_completion.py'),'EOS_code_sha256':sha(ROOT/'fast_eos.py'),'EOS_table_sha256':sha(ROOT/'electron-eos.npz'),'EOS_scope':'The reversible predecessor uses FEOS; this forward pair both uses FastEOS. The no-radiation control quantifies that model switch before radiation differences are interpreted.'}
 a.out.with_suffix('.initial.json').write_text(json.dumps(profile,indent=2)+'\n');a.out.with_suffix('.inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
 result=radial.simulate('B1',base['T0_K'],base['rho0_kg_m3'],base['R0_m'],zones=len(m),nuclear=False,radiation=a.opacity>0,kappa=a.opacity if a.opacity else .001,final_time_units=duration/tu,rtol=1e-8,atol=1e-12,art_linear=0.,art_quad=0.,initial_profile=profile,viscosity_mode='compression_cap',save_trace=a.out.with_suffix('.npz'),trace_samples=400)
 result['assembly_forward_inputs']=inputs;a.out.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['success','termination_reason','elapsed_s','energy_residual_fraction_initial','peak_temperature_K','radiation_escaped_J','initial_energy_J','wall_seconds']},indent=2),flush=True)
if __name__=='__main__':main()
