"""Independent units, limiting-formula and diagnostic-integral checks."""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--sources',type=Path,required=True);ap.add_argument('--results',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();sys.path.insert(0,str(a.research_root/'phase-two/pulse-transport'))
from microphysics import conductivity_envelope
from pynucastro.neutrino_cooling.sneut5_mod import sneut5
checks=[]
def ck(n,v,limit):
 checks.append({'name':n,'absolute_error':float(abs(v)),'limit':limit,'passed':bool(abs(v)<=limit)})
 if abs(v)>limit:raise AssertionError(n)
abar=1/(1/3+1/12+1/168);ye=1/3+2/12+26/168
expected=773682819.5195054;q=sneut5(1e8,2e9,abar=abar,zbar=abar*ye)
ck('published pynucastro three-species cgs example',q/expected-1,1e-12)
KB=1.380649e-23;ME=9.1093837139e-31;HB=1.054571817e-34;C=299792458.;EC=2.307077552e-28
for ne in [1e25,1e30,1e35]:
 T=1e7;zeff=5.;pf=HB*(3*np.pi**2*ne)**(1/3);mf=np.sqrt(ME**2+(pf/C)**2);vf=pf/mf
 nu_full=4*np.pi*EC**2*zeff*ne/(pf*pf*vf);nu_compact=4*zeff*EC**2*mf/(3*np.pi*HB**3)
 ck('degenerate collision cancellation '+str(ne),nu_full/nu_compact-1,1e-12)
 _,kc,kd=conductivity_envelope(np.array([T]),np.array([ne]),np.array([zeff*ne]),lnlambda=1.)
 expected_k=np.pi**2*KB**2*T*ne/(3*mf*nu_compact)
 ck('relativistic-degenerate conductivity '+str(ne),kd[0]/expected_k-1,2e-8)
_,kc,_=conductivity_envelope(np.array([1e7]),np.array([1e25]),np.array([1e25]),lnlambda=10.)
coefficient=kc[0]*(.419*.225/.4)*10/(1e7)**2.5
ck('classical hydrogen conductivity SI normalization',coefficient/1.84e-10-1,.02)
source=json.loads(a.sources.read_text());summary={}
for key in source['cases']:
 p=a.results/(key+'.json');d=json.loads(p.read_text());trace=a.results/(key+'.diagnostic.npz');z=np.load(trace);t=z['time_s'];m=z['mass_kg'];dt=np.diff(t);q=z['neutrino_components_W_kg'];power=z['conductive_face_power_W'];T=z['T_K'];r=z['radius_bounds_m'];rho=z['rho_kg_m3']
 ck(key+' diagnostic hash',0 if hashlib.sha256(trace.read_bytes()).hexdigest()==d['diagnostic_trace_sha256'] else 1,0)
 lum=np.einsum('i,itj->tj',m,q);E=np.sum((lum[:-1]+lum[1:])/2*dt[:,None],axis=0)
 ck(key+' component energy integral',max(abs(E[i]/v-1) for i,v in enumerate(d['thermal_neutrinos']['components_J'].values()) if v),1e-12)
 face=np.sum((abs(power[:,:-1])+abs(power[:,1:]))/2*dt[None,:])
 ck(key+' conductive interface integral',face/d['conduction']['absolute_heat_crossing_all_interfaces_J']-1,1e-12)
 cv=z['cv_total_J_m3_K'];vol=4*np.pi/3*np.diff(r**3,axis=0);cell=np.vstack([np.zeros(len(t)),power])-np.vstack([power,np.zeros(len(t))]);local=abs(cell)/(cv*vol*T);integral=np.sum((local[:,:-1]+local[:,1:])/2*dt[None,:],axis=1)
 ck(key+' local temperature perturbation',max(integral)/d['conduction']['max_integrated_local_abs_deltaT_over_T_total_Cv']-1,1e-12)
 ck(key+' density geometry',float(np.max(abs(rho*vol/m[:,None]-1))),1e-12)
 ck(key+' energy transfer closes',float(np.max(abs(cell.sum(axis=0)))/max(1.,np.max(np.sum(abs(cell),axis=0)))),1e-12)
 ck(key+' entropy sign',min(0.,float(np.min(power*(1/T[1:]-1/T[:-1])))),1e-5)
 rc=.75*np.diff(r**4,axis=0)/np.diff(r**3,axis=0);grad=abs(np.diff(T,axis=0)/np.diff(rc,axis=0));width=np.diff(r,axis=0);L=np.minimum(T[:-1],T[1:])/np.maximum(grad,1e-300)
 summary[key]={'minimum_zone_width_m':float(width.min()),'median_zone_width_m_unweighted':float(np.median(width)),'minimum_neighbor_temperature_scale_m':float(L.min())}
out={'scope':__doc__,'checks':checks,'passed':sum(c['passed'] for c in checks),'total':len(checks),'geometric_scales':summary,'hydrogen_conductivity_coefficient_SI':float(coefficient),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'passed':out['passed'],'total':len(checks),'hydrogen_coefficient_SI':coefficient},indent=2))
