"""Bounded collisional column of old outward tails seen by new main-pulse grains.

Old cell mass is initially uniform in spatial volume and velocity is linear
in radius. Both flows subsequently ballistic; no fields or new plasma response.
Energy transfer assumes all intercepted relative kinetic energy thermalizes in
one grain. Incoming central0.1% belongs to a separately assembled seed.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--schedule',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--quadrature-order',type=int,default=12);a=ap.parse_args();R=a.root.resolve();schedule=json.loads(a.schedule.read_text());loading=json.loads((R/'membrane-receiver/pulse-loading.json').read_text());cache=json.loads((R/'open-loop-redesign/collector-inputs.json').read_text());radius=loading['radius_m'];incoming={v['record']['stage']:v for v in cache['stages'].values()};gx,gw=leggauss(a.quadrature_order);rows={};checks=[]
for key,pair in schedule['worst_adjacent_tail_pairs'].items():
 old,new=key.split('→');d=incoming[new]['initial'];rb=np.r_[0.,d['outer_radius_m']];vb=np.r_[0.,-np.asarray(d['velocity_m_s'])];weights=np.asarray(d['mass_fraction_of_total']);cum=np.r_[0,np.cumsum(weights)]
 out=loading['stages'][old]['cells'];m=np.asarray(out['mass_kg']);ro=np.asarray(out['radius_bounds_m']);vo=np.asarray(out['velocity_bounds_m_s']);sl=np.diff(vo)/np.diff(ro);intercept=vo[:-1]-sl*ro[:-1];F=pair['focus_separation_s'];samples=[]
 for q in np.r_[np.geomspace(.001,.1,30),np.linspace(.1,1,40)]:
  i=min(len(weights)-1,np.searchsorted(cum,q,side='right')-1);f=(q-cum[i])/weights[i];rend=(rb[i]**3+f*(rb[i+1]**3-rb[i]**3))**(1/3);u=vb[i]+(vb[i+1]-vb[i])*(rend-rb[i])/(rb[i+1]-rb[i]);L=F-(radius-rend)/u
  if L<=0:raise ValueError('Main incoming packet predates previous shot')
  cross=(rend+u*F-ro)/(u+vo);knots=np.unique(np.r_[L,cross[(cross>L)&(cross<F)],F]);column=heat=impulse=0.
  for lo,hi in zip(knots[:-1],knots[1:]):
   t=(hi+lo)/2+(hi-lo)/2*gx;r=rend+u*(F-t)
   ii=np.searchsorted(ro+vo*((lo+hi)/2),rend+u*(F-(lo+hi)/2),side='right')-1
   if ii<0 or ii>=len(m):continue
   orig=(r-intercept[ii]*t)/(1+sl[ii]*t);vel=intercept[ii]+sl[ii]*orig
   rho=(3*m[ii]*orig**2/(ro[ii+1]**3-ro[ii]**3))/(1+sl[ii]*t)/(4*np.pi*r**2)
   rel=u+vel;fac=(hi-lo)/2;column+=float(gw@(rho*rel))*fac;heat+=float(gw@(.5*rho*rel**3))*fac;impulse+=float(gw@(rho*rel**2))*fac
  samples.append({'initial_enclosed_mass_fraction':float(q),'incoming_velocity_m_s':float(u),'launch_after_previous_s':float(L),'relative_intercept_column_kg_m2':column,'intercept_heat_J_m2':heat,'intercept_momentum_kg_m_s_m2':impulse})
 worst=max(samples,key=lambda x:x['intercept_heat_J_m2']);rho_s=1000.
 rows[key]={'focus_separation_s':F,'maximum_heat_sample':worst,'max_column_kg_m2':max(x['relative_intercept_column_kg_m2'] for x in samples),
 'grain_heating_J_kg_max_at_density1000':{str(rad):3/(4*rho_s*rad)*worst['intercept_heat_J_m2'] for rad in [1e-6,1e-5,1e-4,.001]},
 'minimum_grain_radius_for100kJ_kg_heat_m':3/(4*rho_s)*worst['intercept_heat_J_m2']/1e5,'samples':samples}
 checks.append({'name':key+' positive outward density/energy','passed':all(s['relative_intercept_column_kg_m2']>=0 and s['intercept_heat_J_m2']>=0 for s in samples)})
out={'scope':__doc__,'quadrature_order':a.quadrature_order,'rows':rows,'criterion':'100kJ/kg is an illustrative incoming-grain heating allowance; actual volatile composition and pulse-preheat design set the applicable limit.',
 'central_tail_scope':'Central0.1%incomingmass excluded by explicit fast-seed architecture; remaining outward mass is retained, possibly re-entrained, never discarded.',
 'checks':checks,'sources_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.schedule,R/'membrane-receiver/pulse-loading.json',R/'open-loop-redesign/collector-inputs.json']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:{'column':v['max_column_kg_m2'],'heat_at_1mm':v['grain_heating_J_kg_max_at_density1000']['0.001'],'minradius':v['minimum_grain_radius_for100kJ_kg_heat_m']} for k,v in rows.items()},indent=2))
