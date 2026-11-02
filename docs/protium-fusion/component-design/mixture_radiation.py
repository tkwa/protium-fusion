#!/usr/bin/env python3
"""Finite-depth photon energy-transfer proxy on newly saved mixture traces.
No hydrodynamic rerun. Spectrum remains conditional on a grey escape history.
"""
import argparse,hashlib,importlib.util,json,math
from pathlib import Path
import numpy as np
from scipy.signal import fftconvolve
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion')
def load(p,n):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
sp=load(ROOT/'membrane-receiver/spectral_proxy.py','sp');rr=load(ROOT/'component-design/receiver_radiation.py','rr')
R=1.1e9;AREA=4*math.pi*R**2;SIG=5.670374419e-8

def calculate(stem):
 stem=Path(stem);d=json.loads(stem.with_suffix('.json').read_text());f=np.load(stem.with_suffix('.npz'));nn=len(f['species']);nz=int(f['zones']);nv=int(f['variables_per_zone']);z=f['state'].T.reshape(-1,nz,nv);m=f['cell_mass_kg'];t=f['time_s'];cum=(z[:,:,nn+4]@m)*float(f['energy_unit_erg_g'])*1e-4;df=np.maximum(np.diff(cum),0)/AREA
 T=z[:,-1,nn]*d['T0_K'];r=z[:,-1,nn+1]*d['R0_m'];kt=.5*(T[:-1]+T[1:])*sp.KB_KEV;rm=.5*(r[:-1]+r[1:]);L=np.diff(cum)/np.diff(t);teff=np.maximum(L,1e-200)/(4*math.pi*rm**2*SIG);teff=teff**.25*sp.KB_KEV
 grid=np.unique(np.r_[np.geomspace(1,1000,6000),rr.TABLES['Al'][:,0],rr.TABLES['Al'][:,0]*(1-1e-10)])
 al=rr.coeff('Al',grid);oxide=al*(2*.0269815385/.101960077)+sp.coeff('O',grid)*(3*.015999/.101960077);carbon=sp.coeff('C',grid)
 props={'Al2O3':(3900.,124.771/.101960077,10.,oxide),'carbon':(2000.,1000.,100.,carbon)}
 cases={};checks=[abs(cum[-1]/d['radiation_escaped_J']-1)<1e-9]
 for name,kts in [('outer_T_quarter',kt*.25),('outer_T_half',kt*.5),('outer_T',kt),('outer_T_double',kt*2),('surface_flux_equivalent_Teff',teff)]:
  spectrum,low,high=sp.planck_parts(kts,df,grid);norm=(np.trapz(spectrum,grid)+low+high)/sum(df);checks.append(abs(norm-1)<1e-5)
  lowgrid=np.geomspace(1e-7,1,1800);lowbins=np.trapz(sp.spectrum(lowgrid,kts),lowgrid,axis=1)*df
  arrival=t-(r-r[0])/299792458.;ft=np.linspace(arrival[0],arrival[-1],12001);dt=ft[1]-ft[0];lowq=np.diff(np.interp(ft,arrival,np.r_[0,np.cumsum(lowbins)]))/dt;fullq=np.diff(np.interp(ft,arrival,np.r_[0,np.cumsum(df)]))/dt;weights=np.sqrt(np.arange(len(lowq))+1)-np.sqrt(np.arange(len(lowq)))
  mats={}
  for mat,(rho,cp,k,coef) in props.items():
   mu,men=coef;depth=np.r_[0.,np.geomspace(1e-8,1.,160)];adiabatic=np.trapz(spectrum[None,:]*men[None,:]*np.exp(-rho*depth[:,None]*mu[None,:]),grid,axis=1)/cp
   lowrise=float(max(fftconvolve(lowq,weights)[:len(lowq)])*2*math.sqrt(dt)/(math.sqrt(math.pi*k*rho*cp)))
   fullrise=float(max(fftconvolve(fullq,weights)[:len(fullq)])*2*math.sqrt(dt)/(math.sqrt(math.pi*k*rho*cp)))
   thick=50e-6;sig=rho*thick;transfer=float(np.trapz(spectrum*men/mu*(-np.expm1(-sig*mu)),grid));allinter=float(np.trapz(spectrum*(-np.expm1(-sig*mu)),grid))
   mats[mat]={'rho_kg_m3':rho,'cp_J_kg_K':cp,'k_W_m_K':k,'depth_m':depth.tolist(),'resolved_no_conduction_sensible_rise_K':adiabatic.tolist(),'resolved_surface_no_conduction_sensible_rise_K':float(adiabatic[0]),'below1keV_all_surface_conduction_peak_K':lowrise,'resolved_plus_soft_upper_within_homogeneous_first_interaction_model_K':float(adiabatic[0])+lowrise,'all_photons_surface_conduction_peak_K':fullrise,'tighter_homogeneous_model_upper_K':min(float(adiabatic[0])+lowrise,fullrise),'first50um_energy_transfer_J_m2':transfer,'first50um_all_interaction_energy_J_m2':allinter,'first50um_fraction_of_incident_energy_transfer':transfer/sum(df),'resolved_all_first_interaction_surface_sensible_K':float(np.trapz(spectrum*mu,grid)/cp),'scope':'mu_en first-interaction local retention; depth profile has no thermal conduction and bounds homogeneous heat diffusion for this assigned profile. Multiple photon scattering, fluorescence, electron transport, changing opacity and layered thermal contact are not solved. Soft tail uses actual time bins and surface conduction; adding independent maxima is conservative within these assumptions.'}
   checks.extend([0<=transfer<=allinter<=sum(df)*1.00001,np.all(np.diff(adiabatic)<=1e-10)])
  # Explicit oxide-on-carbon depth deposition: uncollided transmission through 50um oxide.
  spec_carbon=spectrum*np.exp(-oxide[0]*3900*50e-6)
  carbon_front=float(np.trapz(spec_carbon*carbon[1],grid)/1000)
  cases[name]={'weighted_kT_keV':float(df@kts/sum(df)),'below1keV_fraction':float(low/sum(df)),'above1MeV_fraction':float(high/sum(df)),'normalization_ratio':float(norm),'materials':mats,'50um_oxide_on_carbon_first_interaction_max_no_conduction_K':max(mats['Al2O3']['resolved_surface_no_conduction_sensible_rise_K'],carbon_front),'carbon_immediately_behind50um_oxide_K':carbon_front}
 return {'source_stem':str(stem),'source_sha256':{s:hashlib.sha256(stem.with_suffix(s).read_bytes()).hexdigest() for s in ['.json','.npz']},'grey_opacity_m2_kg':d['opacity_m2_kg'],'H_per_catalyst':d['mixture_inputs']['hydrogen_per_catalyst'],'prompt_fluence_J_m2':float(sum(df)),'prompt_total_J':float(cum[-1]),'saved_times':len(t),'duration_s':float(t[-1]-t[0]),'weighted_Teff_over_Touter':float(df@(teff/kt)/sum(df)),'cases':cases,'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':bool(all(checks))}}
def main(paths):
 out={Path(p).name:calculate(p) for p in paths};return {'radius_m':R,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'import_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'component-design/receiver_radiation.py',ROOT/'membrane-receiver/spectral_proxy.py']},'traces':out,'scope':__doc__,'color_temperature_scope':'Outer T spectral shape is conditional. Quarter, half and double are sensitivity choices, not confidence limits. Flux-equivalent Teff assumes a blackbody photosphere exactly at the outer radius; it is another conditional spectrum, not independent transport.','limitations':['A grey luminosity and total energy alone cannot give a rigorous narrow bound on escaping soft X rays.','mu_en local-retention calculation is a first-interaction deposition model, not full photon/electron Monte Carlo.','Homogeneous depth-profile maximum is a bound only for the chosen heating profile and constant material properties; actual layered conduction/fluorescence can redistribute it.','Pulse angular smearing and Doppler effects omitted; saved radial-ray arrival bins drive soft-tail conduction.','Material k,rho,cp are explicit reference values. Composition and temperature dependent opacity and multiphase response remain to be modeled.']}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stems',nargs='+');p.add_argument('--output',default=str(Path(__file__).with_name('mixture-radiation.json')));a=p.parse_args();o=main(a.stems);Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({s:{'F':v['prompt_fluence_J_m2'],'checks':v['checks'],'cases':{k:{'kT':q['weighted_kT_keV'],'soft_fraction':q['below1keV_fraction'],'oxide_rise':q['materials']['Al2O3']['resolved_plus_soft_upper_within_homogeneous_first_interaction_model_K'],'carbon_rise':q['materials']['carbon']['resolved_plus_soft_upper_within_homogeneous_first_interaction_model_K']} for k,q in v['cases'].items()}} for s,v in o['traces'].items()},indent=2))
