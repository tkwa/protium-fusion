"""Exact current-family endpoint and finite-flight H/He surface spectra.
Mass means hydro baryon mass A*m_u. No range or retention assumption is made.
"""
from pathlib import Path
import json,math,hashlib,importlib.util,sys
import numpy as np
P=Path(__file__).resolve().parent;R=P.parents[2];CHECKOUT=R.parents[1]
MU=1.66053906660e-27;YEAR=365.25*86400
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def loadmod(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def setup():
 rs=loadmod('catcher_receiver_stopping',R/'mixture-design/receiver_stopping.py');pm,st,decay,weak,hashes=rs.setup(CHECKOUT)
 return rs,pm,st,decay,weak,hashes

def packets(record,trace,pm,decay,radius=1.1e9,order=64):
 d=json.loads(record.read_text());z=np.load(trace);nz=int(z['zones']);nv=int(z['variables_per_zone']);q=z['state'].reshape(nz,nv,-1)[:,:,-1];names=d['state_variables'];assert names[:pm.NN]==pm.SPECIES
 m=z['cell_mass_kg'];y=np.maximum(q[:,:pm.NN],0);v=np.r_[0.,q[:,names.index('v/cs0')]*float(z['velocity_unit_cm_s'])/100];r=np.r_[0.,q[:,names.index('r/R0')]*d['R0_m']];assert np.all(np.diff(v)>0) and np.all(np.diff(r)>0)
 g,w=np.polynomial.legendre.leggauss(order);dv=np.diff(v);dr=np.diff(r);a=dr/dv;b=radius-r[:-1]+a*v[:-1]
 speed=.5*(v[1:,None]+v[:-1,None])+dv[:,None]*g/2;rr=.5*(r[1:,None]+r[:-1,None])+dr[:,None]*g/2;mass=m[:,None]*np.broadcast_to(w/2,speed.shape);time=(radius-rr)/speed;proper=time*np.sqrt(1-(speed/pm.C)**2)
 yy=np.stack([decay(y[i],proper[i]) for i in range(nz)]);area=4*np.pi*radius**2;n=(mass[:,:,None]/pm.MU*yy/area).reshape(-1,pm.NN);endpoint=(m[:,None]/pm.MU*y/area).sum(axis=0);cold=(m[:,None]/pm.MU*np.stack([decay(row,1e8) for row in y])/area).sum(axis=0)
 eu=(.5*pm.MU*speed**2/pm.MEV).ravel();dvweight=(dv[:,None]*np.broadcast_to(w/2,speed.shape)).ravel();jac=(speed**2/b[:,None]).ravel();smooth=float(np.sum(m/dv*np.diff(v**3)/6))/area
 checks={'mass_closure_relative':float(np.sum(n*pm.AA)*pm.MU/(d['mass_kg']/area)-1),'smooth_kinetic_closure_relative':float(np.sum(n*pm.AA*eu[:,None])*pm.MEV/smooth-1),'endpoint_mass_fractions_L1':float(np.sum(abs(endpoint*pm.AA*pm.MU*area/d['mass_kg']-np.array([d['output_mass_fractions'][s] for s in pm.SPECIES]))))}
 assert max(abs(x) for x in checks.values())<1e-9
 return {'record':d,'n':n,'eu':eu,'speed':speed.ravel(),'arrival':time.ravel()+d['elapsed_s'],'arrival_after_endpoint':time.ravel(),'mass_fluence':mass.ravel()/area,'dvweight':dvweight,'arrival_jacobian':jac,'cell_index':np.repeat(np.arange(nz),order),'endpoint':endpoint,'cold':cold,'area':area,'checks':checks,'smooth_K_J_m2':smooth,'vmax':v[-1]}

def quantile(x,w):
 k=np.argsort(x);c=np.cumsum(w[k]);return {str(p):float(x[k[min(np.searchsorted(c,p*c[-1]),len(k)-1)]]) for p in [.001,.01,.1,.5,.9,.99,.999]}
def main():
 rs,pm,st,decay,weak,hashes=setup();mp=R/'component-design/final-family/delivered-B1-C13H2-stage-map.json';sp=R/'component-design/final-family-dispatch/sigma350/input.json';mapping=json.loads(mp.read_text());schedule=json.loads(sp.read_text());stages={};checks=[];refs={'A':96.1,'B1':1.92,'B2':20.7,'C13':30.3}
 for cat,decl in mapping['stage_map'].items():
  record=R/decl['record'];trace=Path(decl['trace']);assert sha(record)==decl['result_sha256'] and sha(trace)==decl['raw_trace_sha256'];x=packets(record,trace,pm,decay,order=32);fine=packets(record,trace,pm,decay,order=64);d=x['record'];stage=d['stage'];wall=x['n'].sum(axis=0);period=schedule['plant']['pulse_interval_s'][cat];light={}
  for s in ['p','he4']:
   j=pm.SPECIES.index(s);count=x['n'][:,j];E=pm.AA[j]*x['eu'];nm=float(count.sum());flux=count/x['dvweight']*x['arrival_jacobian'];im=int(np.argmax(flux));light[s]={'nuclei_m2':nm,'mass_kg_m2_baryon':nm*pm.AA[j]*pm.MU,'mass_ug_m2_baryon':nm*pm.AA[j]*pm.MU*1e9,'incident_kinetic_J_m2':float(count@E*pm.MEV),'endpoint_nuclei_m2':float(x['endpoint'][j]),'fully_decayed_nuclei_m2':float(x['cold'][j]),'energy_MeV_quantiles_number':quantile(E,count),'arrival_s_quantiles_number':quantile(x['arrival'],count),'peak_number_flux_m2_s':float(flux[im]),'peak_number_flux_time_s':float(x['arrival'][im]),'packet_nuclei_m2':count.tolist()}
  hhe=sum(z['mass_ug_m2_baryon'] for z in light.values());total=d['mass_kg']/x['area']*1e9
  row={'stage':stage,'catalyst':cat,'record_path':str(record),'trace_path':str(trace),'record_sha256':sha(record),'trace_sha256':sha(trace),'zones':d['zones'],'quadrature_order':32,'period_s':period,'pulses_per_year':YEAR/period,'total_mass_fluence_kg_m2_baryon':d['mass_kg']/x['area'],'total_mass_ug_m2_baryon':total,'H_plus_He_ug_m2_baryon':hhe,'H_plus_He_endpoint_ug_m2_baryon':sum(x['endpoint'][pm.SPECIES.index(s)]*pm.AA[pm.SPECIES.index(s)]*pm.MU*1e9 for s in ['p','he4']),'user_rounded_H_plus_He_ug_m2':refs[stage],'species':light,'whole_pulse_species_wall_nuclei_m2':{s:float(wall[j]) for j,s in enumerate(pm.SPECIES)},'whole_pulse_species_endpoint_nuclei_m2':{s:float(x['endpoint'][j]) for j,s in enumerate(pm.SPECIES)},'whole_pulse_species_fully_decayed_nuclei_m2':{s:float(x['cold'][j]) for j,s in enumerate(pm.SPECIES)},'smooth_kinetic_fluence_J_m2':x['smooth_K_J_m2'],'nodal_kinetic_fluence_J_m2':d['final_kinetic_J']/x['area'],'packet_energy_MeV_u':x['eu'].tolist(),'packet_velocity_m_s':x['speed'].tolist(),'packet_arrival_after_ignition_s':x['arrival'].tolist(),'packet_cell_index':x['cell_index'].tolist(),'packet_mass_kg_m2_baryon':x['mass_fluence'].tolist(),'packet_delta_v_weight_m_s':x['dvweight'].tolist(),'packet_arrival_jacobian_m_s2':x['arrival_jacobian'].tolist(),'checks':dict(x['checks'],quadrature32_vs64_species_L1=float(np.sum(abs(wall-fine['n'].sum(axis=0)))/wall.sum()))}
  stages[stage]=row;checks.append({'stage':stage,'passed':max(abs(v) for v in row['checks'].values())<1e-6});print(json.dumps({'stage':stage,'total_ug_m2':total,'HHe_ug_m2':hhe,'light':{s:{k:z[k] for k in ['nuclei_m2','mass_ug_m2_baryon','incident_kinetic_J_m2','energy_MeV_quantiles_number']} for s,z in light.items()}}),flush=True)
 result={'scope':'Exact coherent delivered-B1/C13H2 source128 family; radial normal incidence on intact graphite, R1.1e9m,350kg/m2. Bulk-speed collisionless projection from actual endpoint cells with finite-flight weak decays; daughter inherits parent velocity. Hydro baryon-mass convention A*m_u; p=H1,he4=He4. Charge-state, beta/recoil transport and post-endpoint pressure acceleration are separate. No retention or release probability inferred.','radius_m':1.1e9,'area_m2':4*np.pi*(1.1e9)**2,'graphite_density_kg_m3':2000.,'slab_sigma_kg_m2':350.,'slab_thickness_m':.175,'nuclear_species':pm.SPECIES,'mass_numbers':pm.AA.tolist(),'charges':pm.ZZ.tolist(),'weak_channels':weak,'source_map':str(mp),'source_map_sha256':sha(mp),'cadence_source':str(sp),'cadence_source_sha256':sha(sp),'code_sha256':sha(__file__),'helper_sha256':hashes,'stages':stages,'checks':checks}
 assert all(c['passed'] for c in checks);write(P/'surface-spectra.json',result)
if __name__=='__main__':main()
