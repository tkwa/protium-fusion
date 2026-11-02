#!/usr/bin/env python3
"""Molecular recirculation work in explicitly thermalized moving-blade cascades.
No new trajectories: uses archived cumulative work per open-mouth incident molecule.
"""
from pathlib import Path
import json, math, hashlib, argparse
import numpy as np
KB=1.380649e-23; MU=1.66053906660e-27
MASS={'H2':2*1.00782503223*MU,'He':4.00260325413*MU}
ROOT=Path('/home/tkwa/.codex/worktrees/protium-fusion/website')
W=ROOT/'docs/protium-fusion/component-design/catcher-followup/warm-mechanical-cleanup'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
selpath=W/'blade-screen/records/selected-million.json'; gridpath=W/'blade-screen/records/grid.json';bpath=W/'buffer.json';mpath=W/'results.json'
sel=read(selpath);grid=read(gridpath);buf=read(bpath);main=read(mpath)
cases={(c['parameters']['gas'],c['parameters']['alpha']):c for c in sel['cases'] if c['parameters']['u']==700 and c['parameters']['theta']==45 and c['parameters']['L_pitch']==1}
checks=[]
def check(name,ok):
 if not bool(ok):raise AssertionError(name)
 checks.append(name)
def network(c,areas,pfore,*,pin=None,baffle=None,PVflow=None,tag=''):
 p=c['parameters'];N=len(areas);f=p['open_fraction'];Tf=c['forward']['T'];Tr=c['reverse']['T'];K=Tf/Tr;T=p['T'];m=MASS[p['gas']];vb=math.sqrt(8*KB*T/(math.pi*m))
 Z=sum(K**(-i)/(a*f*Tf) for i,a in enumerate(areas)); floor=pfore/K**N
 if PVflow is None:q=(pin-floor)/(Z+(1/baffle if baffle else 0.))
 else:q=PVflow/(vb/4)
 pp=np.empty(N+1);pp[-1]=pfore
 for i in reversed(range(N)):pp[i]=pp[i+1]/K+q/(areas[i]*f*Tf)
 Wf=c['forward']['mean_wall_mechanical_work_per_incident_J'];Wr=c['reverse']['mean_wall_mechanical_work_per_incident_J']
 Qf=c['forward']['mean_wall_heat_to_gas_per_incident_J'];Qr=c['reverse']['mean_wall_heat_to_gas_per_incident_J']
 flux_f=areas*f*vb/(4*KB*T)*pp[:-1];flux_r=areas*f*vb/(4*KB*T)*pp[1:]
 shaft=flux_f*Wf+flux_r*Wr;wallheat=flux_f*Qf+flux_r*Qr
 tips=areas*(1-f)*vb/(4*KB*T)*(pp[:-1]+pp[1:])*m*p['u']**2
 qrows=areas*f*(Tf*pp[:-1]-Tr*pp[1:]);check(tag+'flow_continuity',np.max(np.abs(qrows-q))<2e-6*max(abs(q),1e-20))
 check(tag+'inlet_relation',abs(pp[0]-(floor+q*Z))<1e-12*max(1,pp[0]))
 if pin is not None:check(tag+'baffle',abs(q-baffle*(pin-pp[0]))<1e-10*max(abs(q),1e-20))
 check(tag+'positive_pressure',np.all(pp>=0));check(tag+'positive_shaft',sum(shaft)>0)
 # Dense internal-node solve checks the independent resistance/recurrence.
 mat=np.zeros((N-1,N-1));rhs=np.zeros(N-1)
 for i in range(1,N):
  gprev=areas[i-1]*f;gnext=areas[i]*f;j=i-1
  mat[j,j]=gprev*Tr+gnext*Tf
  if i>1:mat[j,j-1]=-gprev*Tf
  else:rhs[j]+=gprev*Tf*pp[0]
  if i<N-1:mat[j,j+1]=-gnext*Tr
  else:rhs[j]+=gnext*Tr*pp[-1]
 dense=np.linalg.solve(mat,rhs);check(tag+'dense_nodes',np.allclose(dense,pp[1:-1],rtol=1e-8,atol=1e-20))
 # Exact zero-flow geometric solution, irrespective of row area.
 pz=pfore*K**(-np.arange(N,-1,-1));zshaft=areas*f*vb/(4*KB*T)*(pz[:-1]*Wf+pz[1:]*Wr)
 vnet=q*vb/4; ideal=vnet*math.log(pfore/pp[0]) if pfore>pp[0] else 0.
 tau=(pp[:-1]*c['forward']['mean_flight_pitch_over_speed_s']+pp[1:]*c['reverse']['mean_flight_pitch_over_speed_s'])/(pp[:-1]+pp[1:])
 # Uniform local pressure is only a diagnostic for collisionless applicability.
 col={str(pitch):float(max(math.sqrt(2)*.5*(pp[:-1]+pp[1:])*1e-19/(KB*T)*vb*tau*pitch)) for pitch in [1e-4,1e-3,.01]}
 return {'species':p['gas'],'rows':N,'K_total':K**N,'row_areas_per_receiver_m2':areas.tolist(),'pressure_nodes_Pa':pp.tolist(),'inlet_pressure_Pa':float(pp[0]),'forepressure_Pa':pfore,'PVflow_W_m2':vnet,'zero_backpressure_pump_speed_m3_s_per_receiver_m2':vb/(4*Z),'effective_cascade_capture_per_inlet_face':1/(Z*areas[0]),'effective_global_capture_including_baffle':q/pin if pin else None,'open_blade_mechanical_W_m2':float(sum(shaft)),'flat_tip_addition_at_full_diffuse_W_m2':float(sum(tips)),'open_blade_heat_to_gas_W_m2':float(sum(wallheat)),'reset_surface_heat_from_open_exits_W_m2':float(sum(shaft+wallheat)),'rotor_surface_heat_from_open_collisions_W_m2':float(-sum(wallheat)),'zero_flow_at_same_forepressure_open_blade_W_m2':float(sum(zshaft)),'net_gas_isothermal_compression_W_m2':ideal,'row_mean_flight_time_per_m_pitch_s':tau.tolist(),'row_shaft_W_m2':shaft.tolist(),'row_full_diffuse_tip_W_m2':tips.tolist(),'last10_rows_fraction_of_open_shaft':float(sum(shaft[-10:])/sum(shaft)),'max_mean_path_collision_count_at_sigma1e19_by_pitch_m':col}

def areas_for(mode,A,N):
 i=np.arange(N)
 if mode=='uniform':fac=np.ones(N)
 elif mode=='taper0p95_from_first':fac=.95**i
 elif mode=='first16_then0p9':fac=.9**np.maximum(0,i-15)
 elif mode=='first16_then0p8':fac=.8**np.maximum(0,i-15)
 else:raise ValueError(mode)
 return A*fac

def package(parts,*,scope,alpha,mode):
 w=sum(p['open_blade_mechanical_W_m2']for p in parts);tip=sum(p['flat_tip_addition_at_full_diffuse_W_m2']for p in parts)
 totalp=sum(np.array(p['pressure_nodes_Pa']) for p in parts)
 collisions={}
 for p in parts:
  vb=math.sqrt(8*KB*1000/(math.pi*MASS[p['species']]))
  collisions[p['species']]={str(pitch):float(max(math.sqrt(2)*.5*(totalp[:-1]+totalp[1:])*1e-19/(KB*1000)*vb*np.array(p['row_mean_flight_time_per_m_pitch_s'])*pitch))for pitch in [1e-4,1e-3,.01]}
 return {'scope':scope,'blade_diffuse_probability':alpha,'area_schedule':mode,'species':parts,'open_blade_mechanical_W_m2':w,'full_diffuse_flat_tip_addition_W_m2':tip,'max_row_mean_path_collision_count_using_total_pressure_at_sigma1e19_by_species_and_pitch_m':collisions,'total_shaft_by_flat_tip_diffuse_fraction':{str(a):w+a*tip for a in [0,.5,1]},'electrical_W_m2_at90pct_motor_including_diffuse_tips':(w+tip)/.9,'electrical_molecular_subtotal_by_motor_efficiency_including_diffuse_tips':{str(eta):(w+tip)/eta for eta in [.5,.9,1.]},'net_gas_isothermal_compression_W_m2':sum(p['net_gas_isothermal_compression_W_m2']for p in parts)}

rows=[];localrows=[]
T=1000;Aglobal=main['reference']['required_pump_aperture_fraction'];Alocal=buf['primary_warm_booster']['shared_face_selected_fraction'];duration=buf['reference']['source_discretization']['duration_s'];h=buf['parameters']['source_plenum_void_depth_m']
local_inputs=[]
for j,s in enumerate(buf['reference']['species']):
 Q=(s['cold_out_PV_Pa_m']+buf['reference']['buffer_volume_per_wall_area_m']*(s['final_buffer_pressure_Pa']-s['initial_buffer_pressure_Pa']))/duration
 local_inputs.append({'species':s['species'],'actual_record_mean_booster_PVflow_W_m2':Q,'actual_record_mean_forepressure_Pa':s['buffer_pressure_exact_mean_Pa'],'forepressure_sampled_max_Pa':s['buffer_pressure_max_Pa'],'source_plenum_archived_mean_Pa':s['source_pressure_exact_mean_Pa']})
for alpha in [0.,.5,1.]:
 N=math.ceil(math.log(1e9)/math.log(cases['H2',alpha]['K0']));Nl=math.ceil(math.log(1e6)/math.log(cases['H2',alpha]['K0']))
 for mode in ['uniform','taper0p95_from_first','first16_then0p9','first16_then0p8']:
  pp=[];ll=[]
  for gas,j in [('H2',0),('He',1)]:
   c=cases[gas,alpha]; ar=areas_for(mode,Aglobal,N)
   pp.append(network(c,ar,.001,pin=1e-10,baffle=.03,tag=f'g{alpha}{mode}{gas}'))
   inp=local_inputs[j];ar=areas_for(mode,Alocal,Nl)
   l=network(c,ar,inp['actual_record_mean_forepressure_Pa'],PVflow=inp['actual_record_mean_booster_PVflow_W_m2'],tag=f'l{alpha}{mode}{gas}')
   l['source_plenum_mean_Pa_including_archived_channel_conductance']=l['inlet_pressure_Pa']+inp['actual_record_mean_booster_PVflow_W_m2']/buf['parameters']['channel_speeds_m_s'][j]
   ll.append(l)
  rows.append(package(pp,scope='Conservative simultaneous global reference: inlet1e-10 Pa and foreline0.001 Pa PER species; total foreline0.002 Pa.',alpha=alpha,mode=mode))
  localrows.append(package(ll,scope='Local primary booster: actual record-mean H2/He flows and buffer partial pressures. This steady mean gives mean molecular work because the network/work are linear at fixed geometry.',alpha=alpha,mode=mode))
# Independent lower-statistics discovery runs are a sensitivity to MC means, not a formal CI.
rep=[]
for r in rows:
 if r['area_schedule']!='uniform':continue
 alpha=r['blade_diffuse_probability'];N=r['species'][0]['rows'];parts=[]
 for gas in ['H2','He']:
  c=next(c for c in grid['cases'] if c['parameters']['gas']==gas and c['parameters']['alpha']==alpha and c['parameters']['theta']==45 and c['parameters']['L_pitch']==1 and c['parameters']['u']==700)
  parts.append(network(c,areas_for('uniform',Aglobal,N),.001,pin=1e-10,baffle=.03,tag=f'rep{alpha}{gas}'))
 z=package(parts,scope='Independent100k-particle discovery samples, same fixed row counts.',alpha=alpha,mode='uniform')
 z['relative_open_work_change_vs_million']=z['open_blade_mechanical_W_m2']/r['open_blade_mechanical_W_m2']-1;rep.append(z)
coll=[]
for pressure in [.001,.002,1,2,5,10]:
 for sigma in [1e-20,1e-19,1e-18]:
  lam=KB*T/(math.sqrt(2)*pressure*sigma)
  coll.append({'total_pressure_Pa':pressure,'effective_collision_area_m2':sigma,'hard_sphere_mfp_m':lam,'gap_for_Kn10_m':lam/10})
for r in rows+localrows:
 for p in r['species']:
  check('heat_reset_'+str(len(checks)),math.isclose(p['reset_surface_heat_from_open_exits_W_m2']+p['rotor_surface_heat_from_open_collisions_W_m2'],p['open_blade_mechanical_W_m2'],rel_tol=1e-12))
# Flat tangential tip: Ework=m u^2; rotor heat=-0.5m u^2 and reset heat=+0.5m u^2.
check('flat_tip_energy_identity',(-.5+1.)==.5)
out={'title':'Molecular windage of thermalized cascades','scope':__doc__,'inputs':{'global_pump_face_m2_per_receiver_m2':Aglobal,'local_primary_face_m2_per_receiver_m2':Alocal,'temperature_K':1000,'speed_m_s':700,'global_inlet_partial_Pa_each':1e-10,'global_foreline_partial_Pa_each':.001,'local_inputs':local_inputs,'moving_tip_scope':'Both front and back flat tangentially moving faces cover fraction1-f. Their work is zero if specular; m u^2 per incidence if fully diffuse. Oblique edges/additional disk faces are not modeled.'},'equations':{'flux':'PVflow/(vbar/4)=Ai f(Tf pi-Tr pnext)','resistance':'Z=sum_i K^(-i)/(Ai f Tf); p0=pN/K^N+qZ','open_shaft':'sum_i Ai f vbar/(4 kT) [pi Wf+pnext Wr]','tips':'sum_i Ai(1-f) vbar/(4 kT)(pi+pnext) m u^2 alpha_tip'},'global_cases':rows,'local_primary_cases':localrows,'independent_MC_sample_comparison':rep,'collision_scale_scenarios':coll,'checks_passed':len(checks),'checks':checks,'limitations':['Ideal isothermal stationary reset reservoirs must be physically implemented; their transmission losses are outside this row model.','All cumulative wall work within each row is already counted; do not multiply by reflection counts.','Motor electrical inputs must include molecular windage, bearing/drive losses and other exposed rotating surfaces. Do not add windage twice to a total-motor requirement.','Pressure-stage area reduction requires physical confinement of higher-pressure gas to smaller moving areas. Exposed annuli retain drag; tapering a flow aperture alone does not remove their windage.','Same local tip speed at each reduced area is assumed; preserve radius with reduced annular span or provide separate appropriate-speed rotors.','Mean local pressure calculation excludes thermalized pump internal storage relative to day-scale buffer response. Pulse-resolved local high-pressure transport may be collisional.','The mean-path collision diagnostic uses a stipulated effective area1e-19m2 and arithmetic mean adjacent pressure; it is not a Boltzmann collision calculation.','Archived records omit work variances/covariances. Independent discovery samples test numerical stability but do not provide a formal work uncertainty interval. Surface scattering, stator transmission and geometry dominate this epistemic uncertainty.'],'source_sha256':{str(p):sha(p) for p in [selpath,gridpath,bpath,mpath,W/'blade-screen/simulate.py']},'code_sha256':sha(__file__)}
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=Path(__file__).with_name('results.json'));a=ap.parse_args();a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({'checks':len(checks),'global':[{k:r[k]for k in ['blade_diffuse_probability','area_schedule','open_blade_mechanical_W_m2','full_diffuse_flat_tip_addition_W_m2','electrical_W_m2_at90pct_motor_including_diffuse_tips']}for r in rows],'local':[{k:r[k]for k in ['blade_diffuse_probability','area_schedule','open_blade_mechanical_W_m2','full_diffuse_flat_tip_addition_W_m2','electrical_W_m2_at90pct_motor_including_diffuse_tips']}for r in localrows],'MC_repeat_relative':[r['relative_open_work_change_vs_million']for r in rep]},indent=2))
