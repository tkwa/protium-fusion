"""Cached receiver/source reducer for arbitrary saved uniform-mixture runs.
No hydro integration. --runs recursively discovers final record JSONs.
Actual trace composition drives the source EOS. Named
radiation accumulator is reduced without rescaling to the record. Outflow uses
mass-conserving uniform-dm/dv cells; nodal and smooth energies remain separate.
"""
from pathlib import Path
import argparse,hashlib,json,math,sys
import numpy as np
from scipy.signal import fftconvolve
RREF=1e9;C=299792458.
PROPS={'k_W_m_K':100.,'rho_kg_m3':2000.,'cp_J_kg_K':1000.}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def temp_rise(q,dt):
 n=np.arange(len(q),dtype=float);kernel=2*math.sqrt(dt)*(np.sqrt(n+1)-np.sqrt(n))
 return fftconvolve(q,kernel)[:len(q)]/math.sqrt(math.pi*np.prod(list(PROPS.values())))

def reps(r,v):
 rc=.75*np.diff(r**4)/np.diff(r**3)
 vc=v[:-1]+np.diff(v)*(rc-r[:-1])/np.diff(r)
 return rc,vc

def quantile(v,m,f,weight='mass'):
 slope=m/np.diff(v);bins=m if weight=='mass' else slope*np.diff(v**3)/6
 cum=np.r_[0.,np.cumsum(bins)];target=f*cum[-1];i=int(np.clip(np.searchsorted(cum,target,side='right')-1,0,len(m)-1))
 speed=v[i]+(target-cum[i])/slope[i] if weight=='mass' else (v[i]**3+6*(target-cum[i])/slope[i])**(1/3)
 return float(speed),i

def timing(ri,vi,ro,vo,m,radius):
 result={}
 for name,f,weight in [('mass999',.001,'mass'),('energy99',.01,'energy')]:
  a,i=quantile(vi,m,f,weight);b,j=quantile(vo,m,f,weight)
  rin=np.interp(a,vi,ri);rout=np.interp(b,vo,ro)
  result[name]={'incoming_s':radius/a,'outgoing_s':radius/b,'occupied_s':radius/a+radius/b,'finite_radius_incoming_s':float((radius-rin)/a),'finite_radius_outgoing_s':float((radius-rout)/b),'omitted_slow_mass_on_each_leg':.001 if name=='mass999' else None}
 _,a=reps(ri,vi);_,b=reps(ro,vo);result['full_resolved']={'incoming_s':float(radius/min(a)),'outgoing_s':float(radius/min(b)),'occupied_s':float(radius/min(a)+radius/min(b))}
 return result

def reduce_case(record_path,root,radius=RREF,thermal_points=65536):
 record_path=record_path.resolve();root=root.resolve()
 base=root/'docs/protium-fusion';sys.path.insert(0,str(base))
 from fast_eos import FastEOS
 import pulse_model as pm
 d=json.loads(record_path.read_text());initial_path=record_path.with_suffix('.initial.json');inputs_path=record_path.with_suffix('.inputs.json');trace_path=Path(d.get('trace_path',record_path.with_suffix('.npz')))
 initial=json.loads(initial_path.read_text());requested_duration=float(json.loads(inputs_path.read_text())['duration_s']);tr=np.load(trace_path);nz=int(tr['zones']);nv=int(tr['variables_per_zone']);tt=tr['time_s'];state=tr['state'].reshape(nz,nv,len(tt));names=d['state_variables'];nn=names.index('T/T0');assert names[:nn]==list(tr['species'])==pm.SPECIES
 m=np.asarray(tr['cell_mass_kg']);area=4*math.pi*radius**2;mi=d['mass_kg']*np.asarray(initial['mass_fraction_of_total']);assert len(mi)==nz and np.max(np.abs(mi/m-1))<1e-10
 ri=np.r_[0.,initial['outer_radius_m']];vi=np.r_[0.,np.abs(initial['velocity_m_s'])]
 ro=np.r_[0.,state[:,names.index('r/R0'),-1]*d['R0_m']]
 vo=np.r_[0.,state[:,names.index('v/cs0'),-1]*float(tr['velocity_unit_cm_s'])/100]
 hashes={str(p):sha(p) for p in [record_path,initial_path,inputs_path,trace_path,base/'fast_eos.py',base/'electron-eos.npz',base/'pulse_model.py',base/'nuclear_screen.py']}
 out={'case':record_path.stem,'stage':d['stage'],'source_paths':{'record':str(record_path),'initial':str(initial_path),'inputs':str(inputs_path),'trace':str(trace_path)},'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'radius_m':radius,'thermal_points':thermal_points,'source_success':d['success'],'source_termination':d['termination_reason'],'source_elapsed_s':d['elapsed_s'],'requested_physical_duration_s':requested_duration,'completed_requested_duration':d['success'] and d['termination_reason']=='final_time' and math.isclose(d['elapsed_s'],requested_duration,rel_tol=1e-8,abs_tol=1e-14),'source_hydro_energy_residual_fraction_initial':d['energy_residual_fraction_initial'],'source_zones':nz,'initial_Y':d.get('initial_Y',{}),'hydrogen_per_catalyst':d.get('mixture_inputs',{}).get('hydrogen_per_catalyst'),'warnings':[]}
 Y0=state[:,:nn,0];b0=Y0@pm.AA;ye0=Y0@pm.ZZ;yi0=Y0.sum(axis=1)
 if d.get('initial_Y'):
  declared=np.array([d['initial_Y'].get(k,0.) for k in pm.SPECIES]);assert np.max(np.abs(Y0-declared))<1e-12
 assert np.max(abs(b0-1))<1e-10
 eos=FastEOS();rhoi=m/(4*np.pi/3*np.diff(ri**3));Ui=m*eos.state(rhoi/1000,np.array(initial['T_K']),ye0,yi0).e*1e-4
 initialK=.25*np.sum(m*(vi[:-1]**2+vi[1:]**2));initialU=float(Ui.sum())
 checks={'physical_duration_relative_to_requested':d['elapsed_s']/requested_duration-1,'trace_endpoint_relative_to_record':float(tt[-1]/d['elapsed_s']-1),'initial_mass_relative':float(sum(m)/d['mass_kg']-1),'initial_internal_relative':initialU/d['initial_internal_energy_J']-1,'initial_kinetic_relative':initialK/d['initial_kinetic_energy_J']-1}
 assert max(abs(checks[x]) for x in ['trace_endpoint_relative_to_record','initial_mass_relative','initial_internal_relative','initial_kinetic_relative'])<1e-8,checks
 if np.any(np.diff(vi)<=0) or np.any(np.diff(ri)<=0):out['usable_ballistic']=False;out['warnings'].append('Nonmonotone incoming profile.');out['checks']=checks;return out
 rc,vc=reps(ri,vi);launch=-(radius-rc)/vc;edge=np.r_[launch[0]-.5*(launch[1]-launch[0]),.5*(launch[1:]+launch[:-1]),launch[-1]+.5*(launch[-1]-launch[-2])];assert np.all(np.diff(edge)>0)
 centroidK=.5*m*vc*vc;centroidK*=initialK/centroidK.sum();localK=.25*m*(vi[:-1]**2+vi[1:]**2)
 out['source_rating']={'peak_kinetic_electric_W':float(np.max(centroidK/np.diff(edge))/.6),'peak_total_prepared_electric_W':float(np.max((centroidK+Ui)/np.diff(edge))/.6),'peak_total_exact_local_halfcell_allocation_W':float(np.max((localK+Ui)/np.diff(edge))/.6),'electric_storage_J':float((initialK+initialU)/.6),'initial_electrons_per_baryon_min_max':[float(min(ye0)),float(max(ye0))],'initial_ions_per_baryon_min_max':[float(min(yi0)),float(max(yi0))],'scope':'60% preparation efficiency. Actual initial composition and EOS; launch-time bins retain historical centroid allocation for comparison, with exact half-cell kinetic allocation separately reported. No switching/focusing overhead is inferred.'}
 ridx=names.index('radiated_energy/uunit');cum=m@state[:,ridx,:]*float(tr['energy_unit_erg_g'])*1e-4
 checks['prompt_accumulator_relative']=float(cum[-1]/d['radiation_escaped_J']-1) if d['radiation_escaped_J'] else float(cum[-1]);checks['minimum_prompt_increment_fraction']=float(min(np.diff(cum))/max(1.,cum[-1]))
 assert abs(checks['prompt_accumulator_relative'])<1e-8 and checks['minimum_prompt_increment_fraction']>-1e-7
 def prompt_screen(time):
  time=time-time[0];grid=np.linspace(0,time[-1],8193);dt=grid[1]-grid[0];q=np.diff(np.interp(grid,time,cum))/(dt*area);rise=temp_rise(q,dt)
  return {'peak_surface_rise_K':float(max(rise)),'peak_flux_W_m2':float(max(q)),'fluence_J_m2':float((cum[-1]-cum[0])/area),'duration_s':float(time[-1])}
 out['prompt']={'energy_J':float(cum[-1]),'emission_time_surface_screen':prompt_screen(tt),'radial_ray_arrival_surface_screen':prompt_screen(tt-state[:,names.index('r/R0'),:][-1]*d['R0_m']/C),'scope':'Named mass-weighted radiation accumulator, no renormalization. Surface-deposition thermal upper reference; spectrum-dependent depth, angular smearing and material phase changes are separate.'}
 if np.any(np.diff(vo)<=0) or np.any(vo[1:]<=0) or np.any(np.diff(ro)<=0):out['usable_ballistic']=False;out['warnings'].append('Nonmonotone/nonoutgoing endpoint; ballistic projection is invalid.');out['checks']=checks;return out
 out['usable_ballistic']=True
 Y=np.maximum(state[:,:nn,-1],0.);bary=Y@pm.AA;rho=m/(4*np.pi/3*np.diff(ro**3));u=eos.state(rho/1000,state[:,names.index('T/T0'),-1]*d['T0_K'],(Y@pm.ZZ)/bary,Y.sum(axis=1)/bary).e*1e-4
 Knode=.25*np.sum(m*(vo[:-1]**2+vo[1:]**2));Pnode=.5*np.sum(m*(vo[:-1]+vo[1:]));slope=m/np.diff(vo);Ksm=float(np.sum(slope*np.diff(vo**3)/6));Usm=float(m@u)
 checks['final_kinetic_nodal_relative']=float(Knode/d['final_kinetic_J']-1);checks['final_internal_relative']=float(Usm/d['final_internal_J']-1);assert abs(checks['final_kinetic_nodal_relative'])<1e-6 and abs(checks['final_internal_relative'])<1e-6
 drdv=np.diff(ro)/np.diff(vo);b=radius-ro[:-1]+drdv*vo[:-1];offset=tt[-1]-drdv
 peakK=.5*slope*vo[1:]**4/b/area;peakP=slope*vo[1:]**3/b/area;peakM=slope*vo[1:]**2/b/area
 out['timing_at_radius']=timing(ri,vi,ro,vo,m,radius)
 tq=[]
 for f in [.99,.5,.01]:
  vel,j=quantile(vo,m,f,'energy');tq.append(float(offset[j]+b[j]/vel))
 # Per-cell internal energy is an advected upper scenario, not rescaled kinetic mass.
 tmax=tq[-1]*8;dt=tmax/thermal_points;grid=(np.arange(thermal_points)+.5)*dt;q=np.zeros(thermal_points);qu=q.copy();pressure=q.copy();massflux=q.copy()
 for j in range(nz):
  speed=b[j]/(grid-offset[j]);ok=(speed>=vo[j])&(speed<vo[j+1])&(grid>tt[-1]);dm=slope[j]*speed[ok]**2/b[j]/area
  massflux[ok]+=dm;q[ok]+=.5*dm*speed[ok]**2;qu[ok]+=dm*u[j];pressure[ok]+=dm*speed[ok]
 rise=temp_rise(q,dt);riseU=temp_rise(q+qu,dt);integrated=float(sum(q)*dt*area);checks['thermal_grid_kinetic_relative_vs_smooth']=integrated/Ksm-1
 out['ballistic']={'mass_kg':float(sum(m)),'mass_fluence_kg_m2':float(sum(m)/area),'kinetic_J_exact_nodal':float(Knode),'kinetic_J_smooth_mass_conserving':Ksm,'smooth_vs_nodal_kinetic_relative':Ksm/Knode-1,'kinetic_fluence_J_m2_exact_nodal':float(Knode/area),'kinetic_fluence_J_m2_smooth':Ksm/area,'momentum_fluence_Pa_s':float(Pnode/area),'peak_kinetic_flux_W_m2':float(max(peakK)),'peak_momentum_flux_Pa':float(max(peakP)),'peak_mass_flux_kg_m2_s':float(max(peakM)),'kinetic_peak_arrival_s':float(tt[-1]+(radius-ro[1:][np.argmax(peakK)])/vo[1:][np.argmax(peakK)]),'arrival_1_50_99_percent_energy_s':tq,'maximum_speed_m_s':float(max(vo)),'maximum_MeV_u':float(.5*pm.MU*max(vo)**2/pm.MEV),'endpoint_internal_J':Usm,'internal_to_kinetic_ratio':Usm/Knode,'minimum_nonzero_nodal_speed_m_s':float(vo[1]),'literal_full_clearance':'Infinite for uniform central-cell dm/dv extending to zero. Full_resolved is a finite shell-representative diagnostic; the retained slow core needs its own route.'}
 out['thermal']={'properties':PROPS,'kinetic_only_peak_rise_K':float(max(rise)),'kinetic_plus_advected_endpoint_internal_peak_rise_K':float(max(riseU)),'material_peak_time_s':float(grid[np.argmax(riseU)]),'grid_duration_s':float(tmax),'grid_points':thermal_points,'scope':'Semi-infinite constant-property carbon, no surface radiation/cooling. Actual mass-conserving velocity histogram; per-cell residual internal energy advected as an upper scenario. Subsequent pressure acceleration or radiation changes its partition. Delayed nuclear decay is not included in this material or prompt screen.'}
 out['energy_channels_J']={'prepared':float(initialK+initialU),'outgoing_kinetic_nodal':float(Knode),'endpoint_internal':Usm,'escaped_prompt_radiation':float(cum[-1]),'future_decay':'Separate nuclear ledger required; no fixed-stage beta inventory assumed.'}
 out['checks']=checks
 if not out['completed_requested_duration']:out['warnings'].append('Short/guarded/failed source trajectory; do not rank as a completed-stage design.')
 if abs(checks['thermal_grid_kinetic_relative_vs_smooth'])>.003:out['warnings'].append('Thermal grid energy error exceeds 0.3%; refine before using temperature extrema.')
 out['scope']='Cached generalized receiver requirements, conditional on ballistic continuation of this actual mixture. Complete mass retained; no mass normalization to nodal kinetic energy. Not a validated apparatus or permanent matter-loss calculation.'
 return out

def self_check():
 # Exact constant-flux half-space solution checks FFT normalization and time indexing.
 q=7.;dt=.03;n=1000
 expected=2*q*np.sqrt(dt*np.arange(1,n+1))/math.sqrt(math.pi*np.prod(list(PROPS.values())))
 err=float(np.max(np.abs(temp_rise(np.full(n,q),dt)/expected-1)))
 assert err<1e-12
 # Independently integrate a finite-radius ballistic shell with Gaussian quadrature.
 # This checks dm/dt, kinetic-power and pressure Jacobians simultaneously.
 v0,v1,r0,r1,R,m=2.,5.,3.,7.,100.,11.
 drdv=(r1-r0)/(v1-v0);b=R-r0+drdv*v0;off=-drdv;slope=m/(v1-v0)
 lo,hi=(R-r1)/v1,(R-r0)/v0
 x,w=np.polynomial.legendre.leggauss(128);t=(hi+lo)/2+(hi-lo)*x/2;v=b/(t-off);flux=slope*v*v/b
 values=[float(np.dot(w,flux*v**power)*(hi-lo)/2) for power in [0,1,2]]
 targets=[m,m*(v0+v1)/2,m*(v0*v0+v0*v1+v1*v1)/3]
 errors=[a/b-1 for a,b in zip(values,targets)];assert max(abs(x) for x in errors)<1e-12
 return {'constant_flux_halfspace_max_relative':err,'ballistic_mass_momentum_twice_kinetic_integral_relative':errors,'passed':True}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--runs',type=Path);ap.add_argument('--record',type=Path,action='append');ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--radius',type=float,default=RREF);ap.add_argument('--thermal-points',type=int,default=65536);ap.add_argument('--force',action='store_true');a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 paths=a.record or [p for p in sorted(a.runs.rglob('*.json')) if not any(p.name.endswith(x) for x in ['.initial.json','.inputs.json','.pretrace.json'])]
 validation=self_check();summary=[]
 for path in paths:
  try:
   d=json.loads(path.read_text())
   if 'trace_path' not in d:continue
   target=a.output_dir/(path.stem+'.receiver.json');old=json.loads(target.read_text()) if target.exists() else None
   valid=not a.force and old and old.get('code_sha256')==sha(Path(__file__)) and old.get('radius_m')==a.radius and old.get('thermal_points')==a.thermal_points and all(Path(p).exists() and sha(Path(p))==h for p,h in old.get('input_sha256',{}).items())
   if valid:result=old;status='cached'
   else:result=reduce_case(path,a.root,a.radius,a.thermal_points);target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');status='reduced'
   summary.append({'case':path.stem,'status':status,'completed':result['completed_requested_duration'],'usable_ballistic':result['usable_ballistic'],'warnings':result['warnings']})
  except Exception as e:
   summary.append({'case':path.stem,'status':'error','error':repr(e)})
 (a.output_dir/'index.json').write_text(json.dumps({'scope':'Only final source-record files found at this invocation; the running grid may add cases. Successful completion is reported separately.','code_sha256':sha(Path(__file__)),'validation':validation,'cases':summary},indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
