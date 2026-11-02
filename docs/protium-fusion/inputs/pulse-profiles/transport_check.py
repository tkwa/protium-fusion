"""Independent frozen-profile opacity/thermalization checks for small stages."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
sys.path[:0]=[str(a.root),str(a.root/'transport/pulse')]
from microphysics import *
import pulse_model as pm
rows=[]
for p in sorted(a.output.glob('scaled-hot-*-64-R*.json')):
 if p.name.endswith(('-initial.json','-inputs.json','-failure.json')):continue
 d=json.loads(p.read_text());trace=p.with_suffix('.npz');raw=np.load(trace)
 nz=d['zones'];nv=int(raw['variables_per_zone']);nn=nv-5;z=raw['state'].reshape(nz,nv,-1);nt=z.shape[-1];t=raw['time_s']
 assert raw['species'].tolist()==pm.SPECIES
 m=d['mass_kg']*np.diff(np.linspace(0,1,nz+1)**3);r=np.vstack([np.zeros(nt),z[:,nn+1]*d['R0_m']]);dr=np.diff(r,axis=0)
 rho=m[:,None]/(4*np.pi/3*np.diff(r**3,axis=0));T=z[:,nn]*d['T0_K'];Y=z[:,:nn].transpose(0,2,1);Y=np.maximum(Y,0);Y/=np.sum(Y*pm.AA,axis=-1)[...,None]
 inc=np.maximum(-m[:,None]*np.sum(np.diff(Y,axis=1)*pm.MASS,axis=-1)*pm.MEV/pm.MU,0);w=np.zeros_like(T);w[:,:-1]+=.5*inc;w[:,1:]+=.5*inc
 ne=rho/MU*np.sum(Y*pm.ZZ,axis=-1);z2=rho/MU*np.sum(Y*pm.ZZ**2,axis=-1);eta,nm,np_=fermions(T,ne);ks,valid=poutanen(T,eta,nm+np_,rho);af,em,aplanck=freefree(T,ne,z2)
 tau=np.sqrt(3*af*(af+ks*rho))*dr
 var=np.minimum(d['elapsed_s'],np.minimum(T/np.maximum(abs(np.gradient(T,t,axis=1)),1e-99),rho/np.maximum(abs(np.gradient(rho,t,axis=1)),1e-99)))
 q=lambda x:weighted_quantile(x,w)
 assert np.all(dr>0) and np.all(rho>0) and np.all(np.diff(t)>0)
 rows.append({'stage':d['stage'],'metadata':p.name,'metadata_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'trace_sha256':hashlib.sha256(trace.read_bytes()).hexdigest(),'samples':nt,'minimum_sampled_outer_radius_m':float(r[-1].min()),'peak_sampled_density_kg_m3':float(rho.max()),'Poutanen_outside_domain_burn_weight_fraction':float(w[~valid].sum()/w.sum()),'burn_weighted_q01_q50_q99':{'rho_kg_m3':q(rho),'T_K':q(T),'Compton_opacity_m2_kg':q(ks),'effective_thermalization_depth_per_cell_x5':q(tau),'effective_depth_both_coefficients_divided100':q(.01*tau),'absorption_time_over_local_evolution_both_coefficients_divided100':q(100/(C*af*var))}})
 print(d['stage'],rows[-1]['burn_weighted_q01_q50_q99'],flush=True)
out={'scope':'Independent frozen saved-profile diagnostic. No coupled transport yield prediction. Free-free coefficient suppression is a stress test, not a certified bound; invalid scattering-fit cells are separately quantified.','microphysics_sha256':hashlib.sha256((a.root/'transport/pulse/microphysics.py').read_bytes()).hexdigest(),'rows':rows}
(a.output/'transport-check.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
