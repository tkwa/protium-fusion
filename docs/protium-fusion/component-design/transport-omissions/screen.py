"""Frozen-trajectory thermal-neutrino and electron-conduction significance screen.

Neutrinos: published Itoh/Timmes sneut5 via pynucastro 3.1.0, cgs input/output.
Conduction: classical Lorentz and relativistic-degenerate electron-ion limits,
with no e-e resistance, compared as a generous sensitivity envelope. This is
not a calibrated interpolation or a coupled radiation/conduction hydro run.
"""
import argparse,json,sys,hashlib,time,math,warnings
from pathlib import Path
import numpy as np
from scipy.special import kv
ap=argparse.ArgumentParser();ap.add_argument('--research-root',type=Path,required=True);ap.add_argument('--sources',type=Path,required=True);ap.add_argument('--out-dir',type=Path,required=True);ap.add_argument('--case',default='all');a=ap.parse_args();r=a.research_root.resolve();sys.path.insert(0,str(r));sys.path.insert(0,str(r/'transport/pulse'))
import pulse_model as pm
from fast_eos import FastEOS
from microphysics import conductivity_envelope,weighted_quantile
from pynucastro.neutrino_cooling import sneut5_mod as nu
import pynucastro
KB=1.380649e-23;ME=9.1093837139e-31;MU=1.66053906660e-27;C=299792458.;HB=1.054571817e-34;EC=2.307077552e-28;ARAD=7.565733250033928e-16
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
a.out_dir.mkdir(parents=True,exist_ok=True);sources=json.loads(a.sources.read_text());eos=FastEOS();start=time.monotonic()

def components(rho,T,ye,yi,z2,*,moments=False,floor=False):
 out=np.zeros((*T.shape,5));low=T<1e7
 for index in np.ndindex(T.shape):
  temp=float(T[index]);density=float(rho[index]/1000)
  if temp<1e7 and not floor:continue
  temp=max(temp,1e7)
  zz=float(z2[index]/ye[index]) if moments else float(ye[index]/yi[index]);aa=zz/float(ye[index])
  with np.errstate(over='ignore',under='ignore',invalid='raise',divide='raise'):
   total,c=nu.sneut5(density,temp,abar=aa,zbar=zz,full_output=True)
  vals=np.array([c.spair,c.splas,c.sphot,c.sbrem,c.sreco],float)*1e-4
  if np.any(vals<0) or not np.isfinite(vals).all():raise ValueError(('invalid neutrino component',index,vals))
  if abs(vals.sum()-total*1e-4)>1e-10*max(1.,total*1e-4):raise ValueError('component sum mismatch')
  out[index]=vals
 return out

def case(key,entry):
 record=Path(entry['record']);trace=Path(entry['trace']);d=json.loads(record.read_text());z=np.load(trace);times=z['time_s'];nz=int(z['zones']);nv=int(z['variables_per_zone']);state=z['state'].reshape(nz,nv,len(times));names=d['state_variables'];nn=names.index('T/T0');m=z['cell_mass_kg'];T=state[:,nn]*d['T0_K'];rb=np.vstack([np.zeros(len(times)),state[:,names.index('r/R0')]*d['R0_m']]);rho=m[:,None]/(4*np.pi/3*np.diff(rb**3,axis=0));Y=np.maximum(np.moveaxis(state[:,:nn],1,-1),0.);Y/=np.sum(Y*pm.AA,axis=-1,keepdims=True);ye=Y@pm.ZZ;yi=Y.sum(axis=-1);z2=Y@(pm.ZZ**2);dt=np.diff(times);duration=times[-1]-times[0];E0=d['initial_energy_J'];Enuc=d['nuclear_energy_J']
 q=components(rho,T,ye,yi,z2);lum=np.sum(m[:,None,None]*q,axis=0);energy=np.trapezoid(lum,times,axis=0)
 # Midpoint states test time interpolation without changing the accepted trace.
 mid=lambda v:(v[:,:-1]+v[:,1:])/2
 rhom=np.sqrt(rho[:,:-1]*rho[:,1:]);Tm=np.sqrt(T[:,:-1]*T[:,1:]);qm=components(rhom,Tm,mid(ye),mid(yi),mid(z2));lm=np.sum(m[:,None,None]*qm,axis=0);E_mid=np.sum(lm*dt[:,None],axis=0)
 # Match charge density and sum Z² n_i for the effective-ion brem sensitivity.
 qz=components(rho,T,ye,yi,z2,moments=True);Ez=np.trapezoid(np.sum(m[:,None,None]*qz,axis=0),times,axis=0)
 coldmask=T<1e7;floor_cold=np.zeros_like(T)
 if coldmask.any():floor_cold[coldmask]=components(rho[coldmask],T[coldmask],ye[coldmask],yi[coldmask],z2[coldmask],floor=True).sum(axis=-1)
 Ecold=float(np.trapezoid(np.sum(m[:,None]*floor_cold,axis=0),times))
 # Deliberately loose independent envelope: peak recorded T at every density,
 # maximal rate across observed Ye/mean-charge corners, all mass for whole time.
 peakT=max(float(T.max()),d.get('peak_temperature_K',0));maxq=0.;maxstate=None
 for rr in np.geomspace(float(rho.min()),float(rho.max()),64):
  for yy in [float(ye.min()),float(ye.max())]:
   for zz in [float((ye/yi).min()),float((ye/yi).max()),float((z2/ye).max())]:
    v=float(nu.sneut5(rr/1000,peakT,abar=zz/yy,zbar=zz))*1e-4
    if v>maxq:maxq=v;maxstate=[rr,peakT,yy,zz]
 Ebox=maxq*m.sum()*duration
 # Electron-ion conduction limits at lnLambda=1. Additional resistance omitted.
 ne=rho/MU*ye;z2ni=rho/MU*z2;kenv,kclass,kdeg=conductivity_envelope(T,ne,z2ni,lnlambda=1.)
 pf=HB*(3*np.pi**2*ne)**(1/3);mf=np.sqrt(ME**2+(pf/C)**2);vf=pf/mf;EF=(mf-ME)*C*C;theta_deg=KB*T/EF;theta_rel=KB*T/(ME*C*C)
 # Maxwell-Boltzmann LTE pair inventory is a diagnostic upper estimate when
 # net electrons are degenerate; pair scattering is deliberately omitted.
 nzero=(ME*C/HB)**3/np.pi**2*theta_rel*kv(2,1/theta_rel)
 minus=.5*(ne+np.sqrt(ne*ne+4*nzero*nzero));plus=np.divide(nzero*nzero,minus,out=np.zeros_like(minus),where=minus>0);boost=1+2*plus/ne
 kenv*=boost
 es=eos.state(rho/1000,T,ye,yi);cvtotal=rho*es.c_v*1e-4
 max_ratio=0.
 for temperature in np.geomspace(1e7,peakT,16):
  for density in np.geomspace(float(rho.min()),float(rho.max()),32):
   for yy in [float(ye.min()),float(ye.max())]:
    for zz in [float((ye/yi).min()),float((ye/yi).max())]:
     loss=float(nu.sneut5(density/1000,temperature,abar=zz/yy,zbar=zz))*1e-4
     energy_state=float(eos.state(density/1000,temperature,yy,yy/zz).e)*1e-4
     max_ratio=max(max_ratio,loss/energy_state)
 Eint_max=float(np.max(np.sum(m[:,None]*es.e*1e-4,axis=0)));E_energy_bound=max_ratio*Eint_max*duration
 qtotal=q.sum(axis=-1);cooltime=np.divide(cvtotal*T,rho*qtotal,out=np.full_like(T,np.inf),where=qtotal>0)
 outside=(rho/1000*ye)<1.;Eout=float(np.trapezoid(np.sum(m[:,None]*qtotal*outside,axis=0),times))
 ld=np.log(rho/1000*ye);lt=np.log(T);ve=np.exp(eos.ee.ev(ld,lt))*0.1 # erg/cm³ -> J/m³
 cve=ve*eos.ee.ev(ld,lt,dy=1)/T
 if np.any(cve<=0):raise ValueError('Nonpositive electron heat capacity')
 chi=kenv/cvtotal;chie=kenv/cve;vthermal=C*np.sqrt(1-(1+1.5*theta_rel)**-2);vchar=np.maximum(vthermal,vf);mfp=3*kenv/(cve*vchar)
 rc=.75*np.diff(rb**4,axis=0)/np.diff(rb**3,axis=0);dr=np.diff(rc,axis=0);grad=np.diff(T,axis=0)/dr
 kface=np.maximum(kenv[:-1],kenv[1:]);flux=-kface*grad;power=4*np.pi*rb[1:-1]**2*flux
 absolute_face_energy=np.trapezoid(abs(power),times,axis=1);total_cross=float(absolute_face_energy.sum())
 cellpower=np.vstack([np.zeros(len(times)),power])-np.vstack([power,np.zeros(len(times))])
 entropy=power*(1/T[1:]-1/T[:-1]);entropy_min=float(entropy.min());closure=float(np.max(abs(cellpower.sum(axis=0)))/max(1.,np.max(np.sum(abs(cellpower),axis=0))))
 def integral(v):return np.trapezoid(v,times,axis=1)
 dTfrac=integral(abs(cellpower)/(m[:,None]*es.c_v*1e-4*T));dTefrac=integral(abs(cellpower)/(cve*4*np.pi/3*np.diff(rb**3,axis=0)*T))
 localgrad=np.maximum(np.vstack([np.zeros(len(times)),abs(grad)]),np.vstack([abs(grad),np.zeros(len(times))]));LT=np.divide(T,localgrad,out=np.full_like(T,np.inf),where=localgrad>0);width=np.diff(rb,axis=0)
 ldif=np.sqrt(2*integral(chi));ldife=np.sqrt(2*integral(chie));gammac=EC*(4*np.pi*ne/3)**(1/3)/(KB*T);dom=np.argmax(Y*pm.ZZ,axis=-1);Gamma=gammac*pm.ZZ[dom]**(5/3)
 # Physical impact-parameter screen; correlation physics is NOT certified by it.
 pchar=mf*vchar;zmean=z2/ye;screen=np.sqrt(KB*T/(4*np.pi*EC*(ne+z2ni)));ai=(3/(4*np.pi*rho/MU*yi))**(1/3);bmax=np.maximum(screen,ai);bmin=np.maximum(zmean*EC/(pchar*vchar),HB/(2*pchar));clog=.5*np.log1p((bmax/bmin)**2)
 weights=m[:,None]*np.r_[dt[0]/2,(dt[:-1]+dt[1:])/2,dt[-1]/2][None,:]
 def quant(v):return weighted_quantile(v,weights,(.01,.5,.99))
 idx=np.unravel_index(np.argmax(mfp/LT),T.shape);jmax=int(np.argmax(dTfrac));cidx=np.unravel_index(np.argmax(abs(power)),power.shape)
 out={'case':key,'stage':d['stage'],'duration_s':duration,'saved_samples':len(times),'zones':nz,'mass_kg':float(m.sum()),'prepared_J':E0,'nuclear_J':Enuc,
  'state_ranges':{'rho_kg_m3':[float(rho.min()),float(rho.max())],'T_K':[float(T.min()),float(T.max())],'record_peak_T_K':peakT,'T_over_TF_mass_time_quantiles':quant(theta_deg),'T_over_TF_range':[float(theta_deg.min()),float(theta_deg.max())],'Gamma_dominant_range':[float(Gamma.min()),float(Gamma.max())],'pair_carrier_multiplier_max':float(boost.max()),'approximate_impact_log_range':[float(clog.min()),float(clog.max())]},
  'thermal_neutrinos':{'components_J':dict(zip(['pair','plasma','photo','brems','recombination'],energy.tolist())),'total_J':float(energy.sum()),'fraction_prepared':float(energy.sum()/E0),'fraction_nuclear':float(energy.sum()/Enuc),'midpoint_interpolated_J':float(E_mid.sum()),'midpoint_vs_trapezoid_relative':float(E_mid.sum()/energy.sum()-1),'charge_square_matched_effective_ion_J':float(Ez.sum()),'effective_ion_relative_change':float(Ez.sum()/energy.sum()-1),'T_below1e7_floor_diagnostic_J':Ecold,'all_mass_peak_T_envelope_J':float(Ebox),'peak_T_envelope_fraction_prepared':float(Ebox/E0),'envelope_maximizing_state_rho_T_Ye_Z':maxstate,'minimum_instantaneous_thermal_cooling_time_s':float(cooltime.min()),'energy_limited_envelope_J':E_energy_bound,'energy_limited_envelope_fraction_prepared':E_energy_bound/E0,'maximum_grid_loss_per_internal_energy_s_inverse':max_ratio,'maximum_trace_internal_J':Eint_max,'rhoYe_below1_g_cm3_fraction_loss':Eout/energy.sum()},
  'conduction':{'model':'max(classical Lorentz,relativistic-degenerate e-ion) at lnLambda1, multiplied by diagnostic pair-carrier enhancement; no e-e resistance. Max adjacent conductivity at each interface. Sensitivity envelope, not calibrated interpolation.',
   'absolute_heat_crossing_all_interfaces_J':total_cross,'fraction_prepared':total_cross/E0,'fraction_nuclear':total_cross/Enuc,'100fold_conductivity_fraction_prepared':100*total_cross/E0,
   'max_integrated_local_abs_deltaT_over_T_total_Cv':float(dTfrac.max()),'mass_fraction_local_deltaT_over_T_above1pct':float(m[dTfrac>.01].sum()/m.sum()),'max_integrated_local_abs_deltaTe_over_T_electron_Cv':float(dTefrac.max()),
   'total_Cv_diffusion_length_m_mass_quantiles':weighted_quantile(ldif,m,(.01,.5,.99)),'electron_Cv_diffusion_length_m_mass_quantiles':weighted_quantile(ldife,m,(.01,.5,.99)),
   'maximum_diffusion_length_over_minimum_zone_width':float(np.max(ldif/width.min(axis=1))),'max_effective_heat_mfp_m':float(mfp.max()),'max_mfp_over_L_T':float(np.max(mfp/LT)),'max_mfp_over_zone_width':float(np.max(mfp/width)),'max_mfp_over_L_T_state':{'cell':int(idx[0]),'time_s':float(times[idx[1]]),'rho_kg_m3':float(rho[idx]),'T_K':float(T[idx]),'T_over_TF':float(theta_deg[idx]),'Gamma':float(Gamma[idx]),'mean_free_path_m':float(mfp[idx]),'L_T_m':float(LT[idx])},
   'largest_local_temperature_integral_cell':jmax,'largest_face_power_state':{'face':int(cidx[0]),'time_s':float(times[cidx[1]]),'power_W':float(power[cidx]),'adjacent_T_K':[float(T[cidx[0],cidx[1]]),float(T[cidx[0]+1,cidx[1]])]},'discrete_energy_closure_relative':closure,'minimum_entropy_production_W_K':entropy_min},
  'checks':{'positive_finite_neutrino_components':True,'positive_EOS_heat_capacity':True,'conduction_energy_conservation':closure<1e-12,'conduction_nonnegative_entropy':entropy_min>=-1e-6},
  'source_sha256':{str(p):sha(p) for p in [record,trace,Path(nu.__file__),r/'transport/pulse/microphysics.py',r/'fast_eos.py',r/'electron-eos.npz']},'code_sha256':sha(__file__),'pynucastro_version':pynucastro.__version__}
 (a.out_dir/(key+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 # Compact moments preserve every microphysical input consumed by this screen.
 np.savez_compressed(a.out_dir/(key+'.diagnostic.npz'),time_s=times,mass_kg=m,T_K=T,rho_kg_m3=rho,Ye=ye,Yi=yi,Z2Y=z2,dominant_Z=pm.ZZ[dom],radius_bounds_m=rb,neutrino_components_W_kg=q,neutrino_mid_components_W_kg=qm,conductive_face_power_W=power,conductivity_envelope_W_m_K=kenv,cv_total_J_m3_K=cvtotal,cv_electron_J_m3_K=cve)
 out['diagnostic_trace_sha256']=sha(a.out_dir/(key+'.diagnostic.npz'))
 (a.out_dir/(key+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({'case':key,'nu_over_prepared':out['thermal_neutrinos']['fraction_prepared'],'nu_peak_envelope':Ebox/E0,'conduction_over_prepared':total_cross/E0,'max_local_dT_T':float(dTfrac.max()),'mfp_L':float(np.max(mfp/LT)),'elapsed_s':time.monotonic()-start}),flush=True)
 return out

if __name__=='__main__':
 for key,entry in sources['cases'].items():
  if a.case in ['all',key]:case(key,entry)
