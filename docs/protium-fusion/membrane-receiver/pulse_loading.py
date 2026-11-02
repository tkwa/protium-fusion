"""Species-resolved wall input from saved corrected64 trajectories; no hydro.

Exact nodal half-cell accounting and a separate continuous uniform-dm/dv
reconstruction are both retained. The latter conserves mass and species;
its kinetic quadrature difference is reported, never removed by scaling mass.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True);ap.add_argument('--radius',type=float,default=1.1e9)
a=ap.parse_args();sys.path.insert(0,str(a.root))
import pulse_model as pm
from fast_eos import FastEOS,KB_CGS,NA,ARAD_CGS
MU=pm.MU;EV=pm.MEV/1e6;ME=9.1093837139e-31;C=pm.C
AREA=4*np.pi*a.radius**2;NN=pm.NN
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
cache_path=a.root/'open-loop-redesign/collector-inputs.json'
cache=json.loads(cache_path.read_text())
W=np.zeros((NN,NN));channels=[]
for j,r in enumerate(pm.RATES):
 if r.weak_type:
  assert len(r.reactants)==1 and r.dens_exp==0
  parent=pm.NUCS.index(r.reactants[0]);rate=r.eval(1e8)
  W[:,parent]+=pm.STOICH[:,j]*rate
  channels.append({'name':r.fname,'parent':pm.SPECIES[parent],
   'products':[str(n).lower() for n in r.products],'rate_s_inverse':rate,
   'half_life_s':float(np.log(2)/rate),'weak_type':r.weak_type})
eig,V=np.linalg.eig(W);invV=np.linalg.inv(V)
assert np.max(np.abs(eig.imag))<1e-14
eig=eig.real;V=V.real;invV=invV.real
def decay(y,t):
 t=np.asarray(t)
 z=(np.exp(t[...,None]*eig)*(invV@y))@V.T
 return np.maximum(z,0.)
for t in [1.,100.,1000.]:
 y=np.arange(1,NN+1,dtype=float);y/=pm.AA@y
 assert np.max(np.abs(decay(y,t)-expm(W*t)@y))<1e-13

def packets(mass,Y,v,r,t0,order):
 gx,gw=np.polynomial.legendre.leggauss(order)
 vel=.5*(v[1:,None]+v[:-1,None])+np.diff(v)[:,None]*gx/2
 rr=.5*(r[1:,None]+r[:-1,None])+np.diff(r)[:,None]*gx/2
 weights=mass[:,None]*np.broadcast_to(gw/2,vel.shape)
 flight=(a.radius-rr)/vel;gamma=1/np.sqrt(1-(vel/C)**2)
 ya=np.stack([decay(Y[i],flight[i]/gamma[i]) for i in range(len(mass))])
 return vel,weights,flight+t0,ya

def reduce_stage(cat,item):
 d=item['record'];paths=item['source_paths'];record_path=Path(paths['record']);trace_path=Path(paths['raw_trace'])
 assert sha(record_path)==item['source_sha256']['record']
 tr=np.load(trace_path);time=tr['time_s'];nz=int(tr['zones']);nv=int(tr['variables_per_zone'])
 assert list(tr['species'])==pm.SPECIES
 state=tr['state'].reshape(nz,nv,len(time));last=state[:,:,-1]
 rawY=last[:,:NN];Y=np.maximum(rawY,0.)
 m=np.asarray(tr['cell_mass_kg']);v=np.r_[0,last[:,NN+2]*float(tr['velocity_unit_cm_s'])/100]
 r=np.r_[0,last[:,NN+1]*d['R0_m']]
 assert np.all(np.diff(v)>0) and np.all(np.diff(r)>0)
 T=last[:,NN]*d['T0_K'];rho=m/(4*np.pi/3*np.diff(r**3))
 ye=(Y@pm.ZZ)/(Y@pm.AA);yi=Y.sum(axis=1)/(Y@pm.AA)
 eos=FastEOS().state(rho/1000,T,ye,yi)
 u=eos.e*1e-4;ui=1.5*KB_CGS*NA*T*yi*1e-4;ur=ARAD_CGS*T**4/(rho/1000)*1e-4
 # Exact staggered kinetic quadrature, retaining the zero-speed half-cell.
 nm=np.zeros(nz+1);ny=np.zeros((nz+1,NN))
 for i in range(nz):
  for j in [i,i+1]:nm[j]+=m[i]/2;ny[j]+=m[i]/2*Y[i]
 ny/=nm[:,None]
 K=.5*np.sum(nm*v*v);P=np.sum(nm*v)
 node_times=np.full(nz+1,np.inf);node_times[1:]=time[-1]+(a.radius-r[1:])/v[1:]
 node_arrival=np.stack([decay(ny[i],(node_times[i]-time[-1])/np.sqrt(1/(1-(v[i]/C)**2))) if i else ny[i] for i in range(nz+1)])
 # Uniform dm/dv within each original cell: exact analytic kinetic moments.
 slope=m/np.diff(v);ke=slope*np.diff(v**3)/6;momentum=slope*np.diff(v**2)/2
 drdv=np.diff(r)/np.diff(v);b=a.radius-r[:-1]+drdv*v[:-1];offset=time[-1]-drdv
 peakq=.5*slope*v[1:]**4/b/AREA;peakp=slope*v[1:]**3/b/AREA
 c=np.r_[0,np.cumsum(ke)];quant=[]
 for f in [.99,.5,.01]:
  k=min(nz-1,np.searchsorted(c,f*c[-1],side='right')-1)
  vel=(v[k]**3+6*(f*c[-1]-c[k])/slope[k])**(1/3)
  quant.append(float(offset[k]+b[k]/vel))
 # Exact integrals with two bounded quadrature orders for in-flight decay.
 p8=packets(m,Y,v,r,time[-1],8);p16=packets(m,Y,v,r,time[-1],16)
 vel,weights,arrival,ya=p8
 pre_counts=np.sum(m[:,None]*Y,axis=0)/(MU*AREA)
 after_counts=np.sum(weights[:,:,None]*ya,axis=(0,1))/(MU*AREA)
 after16=np.sum(p16[1][:,:,None]*p16[3],axis=(0,1))/(MU*AREA)
 number_pre=pre_counts.sum();number_after=after_counts.sum()
 cell_beta_Q=np.sum(weights[:,:,None]*(Y[:,None,:]-ya)*pm.MASS,axis=(0,1,2))/MU*pm.MEV/AREA
 # Flux samples retain both sides of each cell boundary. No hidden kernel smoothing.
 tmin=float(time[-1]+(a.radius-r[-1])/v[-1]);tmax=quant[-1]*5
 boundaries=time[-1]+(a.radius-r[1:])/v[1:]
 grid=np.unique(np.r_[np.geomspace(tmin*(1+1e-10),tmax,480),boundaries*(1-1e-9),boundaries*(1+1e-9)])
 grid=grid[(grid>=tmin)&(grid<=tmax)]
 mf=np.zeros(len(grid));power=mf.copy();pressure=mf.copy();nf=mf.copy()
 for i in range(nz):
  vi=b[i]/(grid-offset[i]);ok=(vi>=v[i])&(vi<v[i+1])
  dm_dt=slope[i]*vi[ok]**2/b[i]/AREA
  mf[ok]+=dm_dt;power[ok]+=.5*vi[ok]**2*dm_dt;pressure[ok]+=vi[ok]*dm_dt
  gamma=1/np.sqrt(1-(vi[ok]/C)**2)
  yy=decay(Y[i],(grid[ok]-time[-1])/gamma)
  nf[ok]+=dm_dt/MU*np.sum(yy,axis=1)
 # Prompt radiation is the exact mass-weighted saved accumulator.
 ri=d['state_variables'].index('radiated_energy/uunit')
 cum=m@state[:,ri,:]*float(tr['energy_unit_erg_g'])*1e-4
 assert abs(cum[-1]/d['radiation_escaped_J']-1)<1e-10
 prompt_arrival=time+(a.radius-state[:,NN+1,:][-1]*d['R0_m'])/C
 electrons_pre=float(pre_counts@pm.ZZ);electrons_after=float(after_counts@pm.ZZ)
 electron_KE=float(np.sum(.5*ME*vel[:,:,None]**2*weights[:,:,None]/MU*ya*pm.ZZ)/AREA)
 raw_mass_error=float(np.max(np.abs(rawY@pm.AA-1)))
 negative_correction=float(np.sum(m[:,None]*(Y-rawY)*pm.AA)/sum(m))
 # A broad atomic binding upper bound: each of Z electrons bound by <14 Z² eV.
 binding_upper=float(after_counts@(14.*pm.ZZ**3)*EV)
 out={'stage':d['stage'],'radius_m':a.radius,'area_m2':AREA,
  'sources':{k:{'path':v,'sha256':sha(v)} for k,v in paths.items()},
  'endpoint':{'time_s':float(time[-1]),'zones':nz,'requested_sound_times':d['elapsed_sound_times'],
   'success':d['success'],'termination_reason':d['termination_reason'],
   'internal_to_kinetic_ratio':d['final_internal_J']/K,'outward_strictly_monotone':True,
   'pressure_conversion_uniform_speed_factor':float(np.sqrt(1+d['final_internal_J']/K)),
   'pressure_scope':'Global residual internal energy envelope only; individual shell acceleration is unresolved. B2 has a separately saved eight-time sensitivity.'},
  'cells':{'mass_kg':m.tolist(),'radius_bounds_m':r.tolist(),'velocity_bounds_m_s':v.tolist(),
   'raw_Y_endpoint':rawY.tolist(),'transport_Y_endpoint':Y.tolist(),
   'specific_internal_J_kg':u.tolist(),'specific_ion_internal_J_kg':ui.tolist(),
   'specific_radiation_internal_J_kg':ur.tolist(),'specific_electron_internal_J_kg':(u-ui-ur).tolist()},
  'exact_nodal':{'mass_kg':nm.tolist(),'velocity_m_s':v.tolist(),'Y_endpoint':ny.tolist(),
   'arrival_time_s':[None if not np.isfinite(t) else float(t) for t in node_times],
   'Y_at_wall':node_arrival.tolist(),'stationary_origin_mass_fraction':float(nm[0]/sum(m)),
   'kinetic_J':K,'momentum_kg_m_s':float(P),'kinetic_reconstruction_relative_error':K/d['final_kinetic_J']-1,
   'scope':'Each cell supplies half its mass and isotope inventory to each bounding node; exact saved nodal kinetic energy. Node0 has no finite ballistic arrival. Discrete impulses have no unique instantaneous power.'},
  'totals':{'mass_fluence_kg_m2':float(sum(m)/AREA),'nuclei_fluence_m2':float(number_after),
   'nuclei_before_flight_m2':float(number_pre),'species_before_flight_m2':pre_counts.tolist(),
   'species_at_wall_m2':after_counts.tolist(),'kinetic_fluence_exact_nodal_J_m2':K/AREA,
   'kinetic_fluence_smooth_J_m2':float(sum(ke)/AREA),'momentum_fluence_exact_nodal_Pa_s':float(P/AREA),
   'momentum_fluence_smooth_Pa_s':float(sum(momentum)/AREA),
   'prompt_fluence_J_m2':float(cum[-1]/AREA),'residual_internal_fluence_J_m2':d['final_internal_J']/AREA,
   'inflight_decay_total_Q_upper_J_m2':float(cell_beta_Q)},
  'smooth':{'model':'uniform dm/dv within each cell; original cell isotope composition retained; finite initial radius included',
   'mass_per_velocity_kg_s_m':slope.tolist(),'kinetic_discrepancy_vs_nodal':float(sum(ke)/K-1),
   'peak_kinetic_power_W_m2':float(max(peakq)),'peak_momentum_flux_Pa':float(max(peakp)),
   'peak_power_arrival_s':float(time[-1]+(a.radius-r[1:][np.argmax(peakq)])/v[1:][np.argmax(peakq)]),
   'arrival_1_50_99_percent_energy_s':quant,
   'old_collector_energy_only_normalization_factor':float(K/sum(ke)),
   'waveform':{'time_s':grid.tolist(),'mass_flux_kg_m2_s':mf.tolist(),'kinetic_power_W_m2':power.tolist(),
     'momentum_flux_Pa':pressure.tolist(),'nuclei_flux_m2_s':nf.tolist()}},
  'spectral_packets':{'quadrature':'8-point Gauss-Legendre in velocity per cell, mass-conserving',
   'cell_index':np.repeat(np.arange(nz),8).tolist(),'velocity_m_s':vel.ravel().tolist(),
   'mass_fluence_kg_m2':(weights/AREA).ravel().tolist(),'arrival_time_s':arrival.ravel().tolist(),
   'kinetic_MeV_per_u':(.5*MU*vel**2/pm.MEV).ravel().tolist(),
   'Y_at_wall':ya.reshape(-1,NN).tolist()},
  'prompt':{'emission_time_s':time.tolist(),'radial_ray_arrival_time_s':prompt_arrival.tolist(),
   'cumulative_fluence_J_m2':(cum/AREA).tolist(),
   'outer_temperature_K':(state[-1,NN,:]*d['T0_K']).tolist(),
   'outer_radius_m':(state[-1,NN+1,:]*d['R0_m']).tolist(),
   'scope':'Grey radiation integral, no photon spectrum. Radial-ray delay; angular source-size smearing is not resolved.'},
  'electron_bookkeeping':{'fully_stripped_electrons_before_flight_m2':electrons_pre,
   'fully_stripped_neutralizing_electrons_at_wall_m2':electrons_after,
   'comoving_electron_kinetic_energy_scale_J_m2':electron_KE,
   'complete_atomic_recombination_upper_J_m2':binding_upper,
   'surface_potential_energy_scale_J_m2_per_V':electrons_after*EV,
   'scope':'Electron kinetic energy is a component-scale diagnostic, not an extra addition to hydrodynamic K. Actual ion charge and sheath potential require plasma transport. Recombination radiation or carried ionization energy must not be counted twice.'},
  'checks':{'raw_baryon_max_error':raw_mass_error,'negative_Y_clipped_mass_fraction_for_transport':negative_correction,
   'transport_composition_is_not_renormalized':True,'8_vs_16_decay_species_L1_relative_nuclei':float(np.sum(np.abs(after_counts-after16))/number_after),
   'nodal_mass_relative_error':float(sum(nm)/sum(m)-1),
   'nodal_isotope_relative_L1':float(np.sum(np.abs(nm@ny-m@Y))/np.sum(m@Y)),
   'smooth_packet_mass_relative_error':float(weights.sum()/sum(m)-1),
   'smooth_packet_kinetic_relative_error':float(np.sum(.5*weights*vel**2)/sum(ke)-1),
   'internal_reconstruction_relative_error':float(m@u/d['final_internal_J']-1)}}
 return out

out={'schema_version':1,'radius_m':a.radius,'species':pm.SPECIES,'A':pm.AA.tolist(),'Z':pm.ZZ.tolist(),
 'units':{'Y':'nuclei per baryon; number fluence = mass_fluence / atomic_mass_unit * Y',
  'kinetic_MeV_per_u':'nonrelativistic energy per baryonic atomic-mass unit; multiply by A for energy per projectile'},
 'atomic_mass_unit_kg':MU,'MeV_J':pm.MEV,'weak_decay_channels':channels,
 'decay_scope':'Same one-body beta rates as the saved network, evolved separately in every velocity packet using proper flight time. Nuclear daughter keeps parent bulk velocity; beta recoil broadening and charged decay-particle transport are separate. In-flight decay Q is not added to the ion kinetic waveform.',
 'source_hashes':{str(p):sha(p) for p in [Path(__file__),cache_path,a.root/'pulse_model.py',a.root/'fast_eos.py',a.root/'electron-eos.npz']},
 'stages':{}}
for cat in ['c13','n14','n15','c12']:
 out['stages'][cache['stages'][cat]['record']['stage']]=reduce_stage(cat,cache['stages'][cat])
out['sources']=[{'url':'https://reaclib.jinaweb.org/','use':'The frozen local ReacLib beta-rate coefficients used for finite-flight daughter inventories.'},
 {'url':'https://ncsp.llnl.gov/sites/ncsp/files/2021-05/BROWN1.pdf','use':'BNL nuclear-chart reference cross-check for O14/O15/N13 lifetime scale.'},
 {'url':'https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html','use':'Atomic ionization energies; shared packet uses a deliberately broader14 Z^3 eV all-electron binding upper bound.'}]
a.output.write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n')
print(json.dumps({k:{'totals':v['totals'],'smooth':{x:v['smooth'][x] for x in ['kinetic_discrepancy_vs_nodal','peak_kinetic_power_W_m2','peak_momentum_flux_Pa','arrival_1_50_99_percent_energy_s']},'checks':v['checks']} for k,v in out['stages'].items()},indent=2))
