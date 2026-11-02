"""Positive nuclear-photon source and outward first-collision bypass screen.
Reconstructs full network flows on saved selected-mixture traces; no hydro rerun.
"""
from pathlib import Path
import hashlib,json,sys,math
import numpy as np
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion');sys.path.insert(0,str(ROOT))
import pulse_model as pm
from nuclear_screen import KB
KAPPAS=[.0001,.0003,.001,.003]
def run(path):
 path=Path(path);d=json.loads(path.read_text());f=np.load(path.with_suffix('.npz'));nn=len(f['species']);nz=int(f['zones']);nv=int(f['variables_per_zone']);st=f['state'].T.reshape(-1,nz,nv);mass=f['cell_mass_kg'];t=f['time_s'];rad=st[:,:,nn+1]*d['R0_m'];r=np.c_[np.zeros(len(t)),rad];dr=np.diff(r);rho=mass/(4*np.pi/3*np.diff(r**3));T=st[:,:,nn]*d['T0_K'];y=np.maximum(st[:,:,:nn],0)
 nreact=np.array([len(r.reactants) for r in pm.RATES]);isweak=np.array([bool(r.weak_type) for r in pm.RATES]);capture=np.array([len(r.products)==1 and len(r.reactants)>=2 and not r.weak_type and pm.Q[j]>0 for j,r in enumerate(pm.RATES)])
 positive=(pm.Q>0)&(~isweak)&(nreact>=2);ncap=[];nall=[];net=[];neg_mean_power=0.;mean_energy_max=0.
 for ti,ri,yi in zip(T,rho/1000,y):
  t9=ti/1e9;tf=np.stack([np.ones(nz),1/t9,t9**(-1/3),t9**(1/3),t9,t9**(5/3),np.log(t9)],axis=1);dtf=np.stack([np.zeros(nz),-1/t9,-t9**(-1/3)/3,t9**(1/3)/3,t9,5*t9**(5/3)/3,np.ones(nz)],axis=1)
  terms=np.exp(tf@pm.COEFF.T);dterms=terms*(dtf@pm.COEFF.T);rates=np.stack([np.bincount(pm.RATE_INDEX,weights=a,minlength=pm.NR) for a in terms]);deriv=np.stack([np.bincount(pm.RATE_INDEX,weights=a,minlength=pm.NR) for a in dterms]);lnslope=np.divide(deriv,rates,out=np.zeros_like(deriv),where=rates>0)
  yp=np.c_[yi,np.ones(nz)];flow=rates*pm.PREF*ri[:,None]**pm.DEXP*np.prod(yp[:,pm.REACT],axis=2)
  ecm=(KB*ti[:,None]/pm.MEV)*(1.5*(nreact-1)+lnslope);gamma=pm.Q+np.maximum(ecm,0)
  power=flow*gamma*pm.MEV/pm.MU
  neg_mean_power+=float(np.sum(abs(flow[:,positive]*np.minimum(ecm[:,positive],0))))
  live=(flow[:,positive]>np.max(flow[:,positive])*1e-12)
  if np.any(live):mean_energy_max=max(mean_energy_max,float(np.max(gamma[:,positive][live])))
  ncap.append(np.sum(power[:,capture],axis=1)*mass);nall.append(np.sum(power[:,positive],axis=1)*mass);net.append((flow@pm.Q)*pm.MEV/pm.MU@mass)
 ncap=np.array(ncap);nall=np.array(nall);net=np.array(net);columns=rho*dr;outer=np.cumsum(columns[:,::-1],axis=1)[:,::-1]-columns
 results={}
 for label,source in [('radiative_captures',ncap),('all_positive_nonweak_reactions_as_photons',nall)]:
  gross=float(np.trapz(source.sum(axis=1),t));rows=[]
  for k in KAPPAS:
   edge=np.exp(-k*outer);x=k*columns;uniform=np.divide(-np.expm1(-x),x,out=np.ones_like(x),where=x>0)
   escaped=float(np.trapz(np.sum(source*edge,axis=1),t));hom=float(np.trapz(np.sum(source*edge*uniform,axis=1),t));coarse=float(np.trapz(np.sum((source*edge)[::2],axis=1),t[::2]))
   rows.append({'kappa_first_interaction_m2_kg':k,'outer_cell_edge_all_outward_first_collision_bypass_J':escaped,'fraction_of_initial_plus_nuclear':escaped/(d['initial_energy_J']+d['nuclear_energy_J']),'uniform_radial_column_source_sensitivity_J':hom,'quadrature_every_other_sample_relative':coarse/escaped-1})
  results[label]={'gross_photon_source_J':gross,'gross_per_net_nuclear_J':gross/d['nuclear_energy_J'],'outermost_cell_source_J':float(np.trapz(source[:,-1],t)),'rows':rows}
 return {'source':str(path),'sha256':{s:hashlib.sha256(path.with_suffix(s).read_bytes()).hexdigest() for s in ['.json','.npz']},'saved_time_points':len(t),'max_relevant_Q_plus_reacting_pair_energy_MeV':mean_energy_max,'negative_mean_energy_diagnostic_unweighted_sum':neg_mean_power,'net_rate_quadrature_relative_to_saved_integrated_Q':float(np.trapz(net,t)/d['integrated_nuclear_energy_J']-1),'results':results,'minimum_total_radial_column_kg_m2':float(np.min(np.sum(columns,axis=1))),'hot_weak_whole_Q_J_upper':d['weak_energy_escaped_upper_J']}
def main():
 spectrum=json.loads(Path(__file__).with_name('spectrum.json').read_text());rows={s:run(v['source']) for s,v in spectrum['stages'].items()}
 aggregates=[{'kappa_first_interaction_m2_kg':k,'all_positive_as_photons_direct_bypass_MeV_per_He':sum(v['results']['all_positive_nonweak_reactions_as_photons']['rows'][i]['outer_cell_edge_all_outward_first_collision_bypass_J']*spectrum['stages'][s]['attempts_per_He']/spectrum['stages'][s]['attempts_per_pulse']/pm.MEV for s,v in rows.items())} for i,k in enumerate(KAPPAS)]
 for a in aggregates:a['all_positive_as_photons_D1p10_lab_bypass_MeV_per_He']=1.10*a['all_positive_as_photons_direct_bypass_MeV_per_He']
 return {'aggregate_cycle':aggregates,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'nuclear_screen_sha256':hashlib.sha256((ROOT/'nuclear_screen.py').read_bytes()).hexdigest(),'network_sha256':hashlib.sha256((ROOT/'pulse_model.py').read_bytes()).hexdigest(),'stages':rows,'scope':__doc__,'thermal_energy':'Maxwellian reacting-pair mean E=kT[3(nreact−1)/2+dln(rate)/dlnT], with analytical Reaclib derivatives. Negative fit means clipped only in photon-energy addendum and recorded; no rates changed.','limits':['Cross sections are explicit photon first-interaction scenarios, not solved dense/hot relativistic scattering coefficients.','Outer-edge placement and all-outward direction maximize direct bypass within each piecewise resolved source. Unresolved subcell production or fast source evolution remains a model limitation.','Direct first-collision attenuation does not equal full thermalization; a scattered gamma or pair cascade may still carry energy outward.','Uniform-in-column emissivity is a separate sensitivity, not a bound.','All positive nonweak reaction energy assigned to photons is an intentionally loose source envelope; includes captures balanced by reverse reactions.'], 'source_spectrum_sha256':hashlib.sha256(Path(__file__).with_name('spectrum.json').read_bytes()).hexdigest()}
if __name__=='__main__':
 o=main();Path(__file__).with_name('hot-bypass.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({s:{'net_quad':v['net_rate_quadrature_relative_to_saved_integrated_Q'],'capture':v['results']['radiative_captures'],'all':v['results']['all_positive_nonweak_reactions_as_photons']} for s,v in o['stages'].items()},indent=2))
