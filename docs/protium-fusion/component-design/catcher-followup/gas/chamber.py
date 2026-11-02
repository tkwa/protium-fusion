"""Actual-dispatch H/He gas loading and conditional diffuse-cavity response.
No hydro or pumping hardware is solved. Source routing and capture remain inputs.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
from scipy.signal import fftconvolve
from scipy.integrate import quad,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
P=Path(__file__).resolve().parent;KB=1.380649e-23

def quantile_time(s,q,species,R):
 m=np.asarray(s['HHe_cell_mass_kg'][species]);v=np.array(s['velocity_bounds_m_s']);r=np.array(s['radius_bounds_m']);cs=np.r_[0,np.cumsum(m)];target=np.asarray(q)*cs[-1];i=np.minimum(len(m)-1,np.searchsorted(cs,target,side='right')-1);f=(target-cs[i])/m[i];vel=v[i]+f*np.diff(v)[i];pos=r[i]+f*np.diff(r)[i]
 return s['endpoint_time_s']+(R-pos)/vel

def emission(d,sp,dt,points):
 end=370*86400+dt;grid=np.arange(0,end+dt,dt);out=np.zeros(len(grid));R=d['R_m'];curves={};q=(np.arange(points)+.5)/points
 for cat,s in d['stages'].items():
  tt=quantile_time(s,q,sp,R);amount=sum(s['HHe_cell_mass_kg'][sp])/points;curves[cat]={'arrival_quantile_midpoints_s':tt.tolist(),'atomic_kg_per_quadrature_point':amount}
  for e in d['events']:
   if e['cat']!=cat:continue
   x=(e['focus_s']+tt)/dt;i=np.ceil(x).astype(int);ok=(i>=0)&(i<len(out));np.add.at(out,i[ok],amount)
 return grid,out,curves

def uniform_response(t,source,Q,lam,desorb):
 if desorb==0:
  kernel=np.exp(-lam*t);initial=Q/lam*kernel
 else:
  mu=1/desorb;assert abs(mu-lam)>1e-10*max(mu,lam)
  kernel=mu/(mu-lam)*(np.exp(-lam*t)-np.exp(-mu*t));b=Q/(lam-mu);initial=(Q/lam-b)*np.exp(-lam*t)+b*np.exp(-mu*t)
 return initial+fftconvolve(source,kernel)[:len(t)]

def firstflight_I(x):
 if x<=0:return 1.
 if x>=2:return 0.
 if x<1e-4:return 1+4*x/3+3*x*x/4
 return (x*x-2*x-2*(1+x*x)*math.log(abs(1-x)))/(2*x*x)

def firstflight_K(z):
 if z==0:return 1.
 if z<.1:
  return quad(lambda w:2*w**3*math.exp(-w*w)*firstflight_I(z*w),0,min(12,2/z),epsabs=2e-9)[0]
 return quad(lambda x:2*x**3/z**4*math.exp(-(x/z)**2)*firstflight_I(x),0,2,points=[1.],epsabs=2e-9,epsrel=2e-8)[0]

def firstflight_S(z):
 z=np.array(z,dtype=float);out=np.ones_like(z);nz=z>0;b=np.zeros_like(z);b[nz]=4/z[nz]**2;regular=nz&(b>.01);small=nz&~regular
 out[regular]=1-2/b[regular]+(1+2/b[regular])*np.exp(-b[regular]);u=b[small];out[small]=u*u/6-u**3/12+u**4/40-u**5/180
 return out

def firstflight_pressure_kernel(z):
 # 2/3 of the surviving kinetic energy in kT per emitted gas particle.
 z=np.asarray(z,dtype=float);out=np.full_like(z,4/3);nz=z>0;b=np.zeros_like(z);b[nz]=4/z[nz]**2;reg=nz&(b>.03);small=nz&~reg
 out[reg]=2/3*(2-6/b[reg]+(b[reg]+4+6/b[reg])*np.exp(-b[reg]));u=b[small];out[small]=2/3*(u**3/12-u**4/20+u**5/60-u**6/252+u**7/1344)
 return out

def fields_stats(t,columns,pressures,d):
 total=columns.sum(axis=0);pressure=pressures.sum(axis=0);queries={};R=d['R_m'];cum=np.r_[0.,np.cumsum(total[:-1]*np.diff(t))];qv=[.001,.01,.1,.5,.9,.99,.999]
 for cat,s in d['stages'].items():
  mm=np.array(s['source_mass_weights']);vv=np.array(s['source_velocity_bounds_m_s']);cs=np.r_[0,np.cumsum(mm)];catrows=[]
  for q in qv:
   i=min(len(mm)-1,np.searchsorted(cs,q,side='right')-1);v=vv[i]+(q-cs[i])/mm[i]*(vv[i+1]-vv[i]);L=R/v;shots=[e for e in d['events'] if e['cat']==cat and e['focus_s']>=L and e['focus_s']<=t[-1]];end=np.array([e['focus_s'] for e in shots]);start=end-L
   col=(np.interp(end,t,cum)-np.interp(start,t,cum))/L
   catrows.append({'source_mass_quantile':q,'velocity_m_s':v,'flight_s':L,'minimum_column_kg_m2':float(min(col)),'median_column_kg_m2':float(np.median(col)),'maximum_column_kg_m2':float(max(col)),'worst_focus_s':float(end[np.argmax(col)]),'full_retention_specific_heat_1cm_rho1000_J_kg':float(3*max(col)*v*v/(8*1000*.01))})
  queries[cat]=catrows
 return {'mean_radial_column_kg_m2':float(np.trapz(total,t)/(t[-1]-t[0])),'peak_radial_column_kg_m2':float(max(total)),'minimum_radial_column_kg_m2':float(min(total)),'peak_column_time_s':float(t[np.argmax(total)]),'mean_pressure_Pa':float(np.trapz(pressure,t)/(t[-1]-t[0])),'peak_pressure_Pa':float(max(pressure)),'source_profile_pellet_encounters':queries}

def main(dt=300,points=64,emission_override=None,Q_override=None,suffix=""):
 d=json.loads((P/'inputs.json').read_text());R=d['R_m'];A=4*math.pi*R*R;V=A*R/3;Q0=np.array([d['plant']['recovered_H_kg_s'],d['plant']['helium_per_s']*d['atomic_mass_kg']['he4']]);Q=Q0 if Q_override is None else np.asarray(Q_override);emit=emission if emission_override is None else emission_override;tg,h,curveh=emit(d,'p',dt,points);_,he,curvehe=emit(d,'he4',dt,points);source=[h,he];rows=[];massatom=np.array([d['atomic_mass_kg']['p'],d['atomic_mass_kg']['he4']]);curve={'H':curveh,'He':curvehe};stationary=[]
 for T,Hatoms in [(1000,2),(300,2),(2000,2),(1000,1)]:
  masses=massatom*np.array([Hatoms,1]);vbar=np.sqrt(8*KB*T/(math.pi*masses));tau=4*R/(3*vbar);M=Q*tau;col=M*3/A;p=M/masses*KB*T/V
  stationary.append({'T_K':T,'H_atoms_per_gas_particle':Hatoms,'species':['H2' if Hatoms==2 else 'H','He'],'Q_atomic_kg_s':Q.tolist(),'gas_particle_mass_kg':masses.tolist(),'vbar_m_s':vbar.tolist(),'unit_capture_residence_s':tau.tolist(),'unit_capture_inventory_kg':M.tolist(),'unit_capture_columns_kg_m2':col.tolist(),'unit_capture_pressures_Pa':p.tolist()})
  for g in ([1,.1,.01] if T==1000 and Hatoms==2 else [1]):
   for td in ([0,3600,86400,30*86400] if T==1000 and Hatoms==2 else [0]):
    M=np.array([uniform_response(tg,source[j],Q[j],g/tau[j],td) for j in range(2)]);columns=3*M/A;pressures=M/masses[:,None]*KB*T/V;stats=fields_stats(tg,columns,pressures,d);rows.append({'model':'well_mixed_volume_balance','T_K':T,'H_atoms_per_gas_particle':Hatoms,'effective_wall_capture_g':g,'material_release_time_s':td,'cavity_injected_fraction_by_H_He':(Q/Q0).tolist(),'unit_capture_residence_s':tau.tolist(),'initial_inventory_kg':(Q*tau/g).tolist(),'initial_material_release_reservoir_kg':(Q*td).tolist(),**stats})
 # Physical collisionless first-flight response: all-wall diffuse cosine emission,
 # flux-weighted Maxwellian speeds, actual spherical chords, unit removal at return.
 zgrid=np.r_[0,np.geomspace(1e-4,1e3,600)];kg=np.array([firstflight_K(z) for z in zgrid]);interp=PchipInterpolator(zgrid,kg);T=1000;masses=massatom*np.array([2,1]);vth=np.sqrt(2*KB*T/masses);cols=[];pressures=[];checks=[];kernels={}
 for j in range(2):
  z=tg*vth[j]/R;K=interp(z);S=firstflight_S(z);Pk=firstflight_pressure_kernel(z);intPk=np.r_[0,cumulative_trapezoid(Pk,tg)];intK=np.r_[0,cumulative_trapezoid(K,tg)];intS=np.r_[0,cumulative_trapezoid(S,tg)];vbar=2*vth[j]/math.sqrt(math.pi);IK=4*R/vbar;IS=IK/3
  ckint=quad(firstflight_K,0,2,epsrel=2e-7)[0]+quad(firstflight_K,2,100,epsrel=2e-7)[0];checks.append({'name':f'species{j} analytic steady radial column kernel','relative_error':ckint/(2*math.sqrt(math.pi))-1});assert abs(checks[-1]['relative_error'])<1e-5
  C=Q[j]/A*np.maximum(0,IK-intK)+fftconvolve(source[j],K)[:len(tg)]/A;M=Q[j]*np.maximum(0,IS-intS)+fftconvolve(source[j],S)[:len(tg)];cols.append(C);pressure=(Q[j]*np.maximum(0,IS-intPk)+fftconvolve(source[j],Pk)[:len(tg)])/masses[j]*KB*T/V;pressures.append(pressure);kernels[str(j)]={'dimensionless_age':zgrid.tolist(),'column_kernel':kg.tolist(),'inventory_survival':firstflight_S(zgrid).tolist(),'kinetic_pressure_kernel':firstflight_pressure_kernel(zgrid).tolist()}
 rows.append({'model':'collisionless_spherical_firstflight_g1','T_K':1000,'H_atoms_per_gas_particle':2,'effective_wall_capture_g':1.,'material_release_time_s':0.,'cavity_injected_fraction_by_H_He':(Q/Q0).tolist(),**fields_stats(tg,np.array(cols),np.array(pressures),d),'interpretation':'Pressure is the volume-averaged trace of the kinetic stress tensor, with a separate energy-weighted survival kernel; it is not a local wall/plume value. Gas molecules are ballistic and absorbed at their first return. Background at t0 represents a steady constant source before t0; then actual chronological emissions replace it. The radial column uses exact spherical chord geometry. Pellet exposure uses quasi-static radial columns averaged over each flight; thermal gas speeds are small compared with selected pellet speeds.'})
 out={'schema':1,'family':d['family'],'cavity_source_mass_rate_atomic_kg_s':Q.tolist(),'cavity_source_fraction_H_He':(Q/Q0).tolist(),'radius_m':R,'dt_s':dt,'arrival_quantile_points':points,'source_time_discretization':'Each release quadrature point is assigned to the next grid boundary (delay<dt), and pellet-flight integrals use left-continuous intervals. This prevents any outgoing gas from a pulse leaking numerically into its own earlier incoming flight.','stationary':stationary,'cases':rows,'source_arrival_quadrature':curve,'checks':checks,'source_sha256':{'inputs.json':hashlib.sha256((P/'inputs.json').read_bytes()).hexdigest()},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':['Every case uses actual final350 stage outputs and focus times. Full release and f=1 are response normalizations, not measured desorption fractions. Multiply H and He responses separately by their own cavity injection fractions.','Capture g includes aperture area, molecular transmission and irreversible routing into the closed return/product line. A hot wall bounce is not removal.','Well-mixed balance is a conditional diffuse model, with finite initial steady reservoirs. It is not a local plume bound. The separate first-flight kernel tests transient spherical geometry at g1 only.','Gas collision cross sections, charged-plume fields, side-channel backstream, coating ablation and dust are separate screens.','All mass is retained somewhere in the station: a source delay, pumping, heating or ionization does not create permanent loss.']};(P/f'chamber-dt{dt}-q{points}{suffix}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');np.savez_compressed(P/f'chamber-dt{dt}-q{points}{suffix}.npz',time_s=tg,firstflight_column_H_He=np.array(cols),firstflight_volume_pressure_H_He=np.array(pressures),source_H_He_kg=np.array(source));print(json.dumps([{k:r[k] for k in ['model','T_K','H_atoms_per_gas_particle','effective_wall_capture_g','material_release_time_s','mean_radial_column_kg_m2','peak_radial_column_kg_m2','peak_pressure_Pa']} for r in rows],indent=2))
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--dt',type=float,default=300);ap.add_argument('--points',type=int,default=64);args=ap.parse_args();main(args.dt,args.points)
