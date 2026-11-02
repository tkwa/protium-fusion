#!/usr/bin/env python3
"""Trace-derived opacity scales; no hydro rerun and no calibrated replacement EOS.
Poutanen2017 Compton mean fits; nondegenerate Kramers free-free diagnostic g=1.
"""
import argparse,hashlib,json,math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from pynucastro.eos.electron_eos import ElectronEOS
import pynucastro as pn
SIGT=6.6524587321e-29;MU=1.66053906660e-27;KB=1.380649e-23;HBAR=1.054571817e-34;ME=9.1093837139e-31;C=299792458.;HC=6.62607015e-27;KBC=1.380649e-16;NA=6.02214076e23
EE=ElectronEOS()
def lamb(kt,eta):
 # Table1: eq36 degeneracy fit2–300keV; eta<-10 uses eq41 nondegenerate fit.
 if kt<2 or kt>300 or eta>7:return None
 if eta<-10:return 1+(kt/41.5)**.90
 xi=math.exp(.682*eta-.0454*eta*eta)
 f1=1+.24*xi+.0043*xi*xi;f2=1+.05*xi-.0067*xi*xi;f3=1-.037*xi+.0031*xi*xi
 return f1*(1+(kt/(43.3*f2))**(.885*f3))
def state(rhosi,T,ye,z2,yi,cell_width_m):
 rho=rhosi/1000; comp=SimpleNamespace(zbar=ye,abar=1.,ye=ye)
 es,ps=EE.pe_state(rho,T,comp,compute_derivs=False);kt=T*8.617333262e-8;theta=kt/510.99895
 em=float(es.eta);ep=-em-2/theta;lm=lamb(kt,em);lp=lamb(kt,ep)
 kappat=SIGT*ye/MU
 kc=None if lm is None or lp is None else SIGT*(es.n/lm+ps.n/lp)*1e6/rhosi
 ne=es.n*1e6;pf=HBAR*(3*math.pi**2*ne)**(1/3);ef=(math.sqrt((pf*C)**2+(ME*C*C)**2)-ME*C*C)/KB
 # Classical full-ionization monochromatic absorption; density in cgs.
 x=np.geomspace(1e-4,100,1600);nu=x*KBC*T/HC
 af=3.692e8*T**(-.5)*(rho*NA*ye)*(rho*NA*z2)*nu**(-3)*(-np.expm1(-x))
 kf=af/rho*.1;weight=x**4*np.exp(-x)/(-np.expm1(-x))**2;norm=np.trapz(weight,x)
 krff=norm/np.trapz(weight/kf,x)
 krmix=None if kc is None else norm/np.trapz(weight/(kf+kc),x)
 alpha4=float(np.interp(4,x,af))*100 # m^-1
 scalpha=0 if kc is None else kc*rhosi
 tauff4=alpha4*cell_width_m;taueff=math.sqrt(tauff4*(tauff4+scalpha*cell_width_m))
 return {'rho_kg_m3':rhosi,'T_K':T,'kT_keV':kt,'Ye':ye,'sum_Z2Y':z2,'ion_Y':yi,'electron_degeneracy_eta':em,'electron_Fermi_temperature_K':ef,'T_over_TF':T/ef,'Z7_ion_coupling_indicator':7**(5/3)*2.307077552e-28/(3/(4*math.pi*ne))**(1/3)/(KB*T),'plasma_frequency_energy_over_kT':HBAR*math.sqrt(ne*1.602176634e-19**2/(8.8541878128e-12*ME))/(KB*T),'pair_number_per_net_electron':float(ps.n/(rho*NA*ye)),'Thomson_net_electron_opacity_m2_kg':kappat,'Compton_Poutanen_fit_m2_kg':kc,'classical_ff_Rosseland_g1_m2_kg':float(krff),'grey_Compton_plus_spectral_classical_ff_Rosseland_m2_kg':None if krmix is None else float(krmix),'classical_ff_qualitative_validity':bool(em<-.5 and kt<50),'cell_width_m':cell_width_m,'classical_ff_absorption_tau_at4kT':tauff4,'effective_thermalization_tau_at4kT':taueff,'scope':'Kramers ff numbers are unsuppressed classical reference; no dense-plasma/degeneracy/relativistic correction. Effective tau is a cell diagnostic, not a spectral photosphere solution.'}
def trace(path):
 path=Path(path);d=json.loads(path.with_suffix('.json').read_text());f=np.load(path.with_suffix('.npz'));sp=[pn.Nucleus(s) for s in f['species']];A=np.array([n.A for n in sp]);Z=np.array([n.Z for n in sp]);nn=len(sp);nz=int(f['zones']);nv=int(f['variables_per_zone']);st=f['state'].T.reshape(-1,nz,nv);mass=f['cell_mass_kg'];r=st[:,:,nn+1]*d['R0_m'];dr=np.diff(np.c_[np.zeros(len(st)),r],axis=1);rho=mass/(4*math.pi/3*np.diff(np.c_[np.zeros(len(st)),r**3],axis=1));T=st[:,:,nn]*d['T0_K'];y=np.maximum(st[:,:,:nn],0);norm=y@A;ye=(y@Z)/norm;yi=y.sum(axis=2)/norm;z2=(y@(Z*Z))/norm
 rad=np.sum(st[:,:,nn+4]*mass,axis=1);qidx=[int(np.searchsorted(rad,frac*rad[-1])) for frac in [.05,.5,.95]]
 hot=np.unravel_index(np.argmax(T),T.shape);dense=np.unravel_index(np.argmax(rho),rho.shape);pts={(0,0),(0,nz-1),(len(st)-1,0),(len(st)-1,nz-1),hot,dense}
 for i in qidx:pts.update({(min(i,len(st)-1),0),(min(i,len(st)-1),nz-1)})
 rows=[]
 for i,j in sorted(pts):
  row=state(float(rho[i,j]),float(T[i,j]),float(ye[i,j]),float(z2[i,j]),float(yi[i,j]),float(dr[i,j]));row.update(time_s=float(f['time_s'][i]),cell=int(j),trace_index=int(i));rows.append(row)
 return {'source_stem':str(path),'input_sha256':{p.suffix:hashlib.sha256(p.read_bytes()).hexdigest() for p in [path.with_suffix('.npz'),path.with_suffix('.json')]},'all_saved_range':{'T_K':[float(T.min()),float(T.max())],'rho_kg_m3':[float(rho.min()),float(rho.max())],'Ye':[float(ye.min()),float(ye.max())],'cell_old_kappa_tau':[float((.001*rho*dr).min()),float((.001*rho*dr).max())]},'selected_states':rows}
def main():
 stems={'A':'/var/tmp/tkwa-chatgpt/protium-b2-regression/regression-A-compression_cap','B1':'/var/tmp/tkwa-chatgpt/protium-b2-regression/regression-B1-compression_cap','C13':'/var/tmp/tkwa-chatgpt/protium-b2-regression/regression-C13-compression_cap','B2':'/var/tmp/tkwa-chatgpt/protium-b2-completion/B2-64-compression-cap'}
 traces={k:trace(v) for k,v in stems.items()};mix={}
 for label,A,Z,x0 in [('B2',12,6,2),('C13',13,6,1),('A',14,7,1),('B1',15,7,1)]:
  mix[label]=[{'H_per_catalyst':x,'Ye':(Z+x)/(A+x),'capture_density_factor':x/(A+x)**2,'ff_density_factor':(Z+x)*(Z*Z+x)/(A+x)**2,'ideal_thermal_particles_per_baryon':(Z+1+2*x)/(A+x),'Thomson_m2_kg':SIGT/MU*(Z+x)/(A+x)} for x in sorted(set([x0,2,4,8,A,Z**1.5,32]))]
 checks=[abs(lamb(10,-100)-(1+(10/41.5)**.9))<1e-12, lamb(1,0) is None, lamb(50,8) is None]
 for d in traces.values():
  for r in d['selected_states']:
   checks.extend([r['Thomson_net_electron_opacity_m2_kg']>0,r['classical_ff_Rosseland_g1_m2_kg']>0, r['pair_number_per_net_electron']>=0])
 return {'checks':{'count':len(checks),'passed':int(sum(checks)),'all_passed':all(checks)},'scope':__doc__,'sources':{'Compton_fit':'https://doi.org/10.3847/1538-4357/835/2/119','primary_ff_application':'https://doi.org/10.1051/0004-6361/202451835','relativistic_Gaunt_reference':'https://arxiv.org/abs/1502.06788','dense_plasma_opacity':'https://doi.org/10.1103/PhysRevE.102.043209'},'traces':traces,'fixed_rho_T_mixture_factors':mix,'Compton_fit_domain':'2≤kT≤300keV, eta≤7; eta<-10 uses nondegenerate limit, within table uses full fit. Both electron and positron contributions evaluated with actual EOS populations and physical eta_plus=-eta_minus-2/theta. Pair contribution is a componentwise extension of electron fit, flagged diagnostic.','limitations':['Grey Rosseland sensitivity does not establish LTE, spectral emission or multigroup transport.','Classical ff opacity has no degeneracy/ion-correlation/relativistic/plasma-frequency correction; only comparison scale.','The source EOS includes equilibrium electron-positron pairs but fully ionized ideal ions and no Coulomb/ionization correction.','Sparse representative states establish scale, not extrema of fitted opacity throughout trajectories.'], 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--output',default=str(Path(__file__).with_name('opacity-interface.json')));p=a.parse_args();out=main();Path(p.output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({s:{'range':d['all_saved_range'],'selected_kappa':[r['Compton_Poutanen_fit_m2_kg'] for r in d['selected_states']],'eta':[r['electron_degeneracy_eta'] for r in d['selected_states']],'pair':[r['pair_number_per_net_electron'] for r in d['selected_states']]} for s,d in out['traces'].items()},indent=2))
