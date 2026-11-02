"""Read-only kinematics/geometry from frozen corrected64 input profiles."""
import argparse,json,hashlib,math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1],help='docs/protium-fusion root')
parser.add_argument('--output',type=Path,default=Path(__file__).with_name('incoming-timing.json'))
args=parser.parse_args()
ROOT=args.root.resolve()
checks=[]
def check(name,condition):
 checks.append({'name':name,'passed':bool(condition)})
 if not condition:raise AssertionError(name)
R=1.1e9;NREF=1000
paths=[ROOT/'open-loop-redesign'/n for n in ['collector-inputs.json','collector.json','operating-cases.json']]
cache,collector,ops=[json.loads(p.read_text()) for p in paths]
stages={}
for cat,item in cache['stages'].items():
 d=item['record'];p=item['initial'];M=d['mass_kg'];w=np.asarray(p['mass_fraction_of_total']);m=M*w
 r=np.r_[0.,p['outer_radius_m']];u=np.r_[0.,-np.asarray(p['velocity_m_s'])]
 check(d['stage']+' ordered initial radii and speeds',np.all(np.diff(r)>0) and np.all(np.diff(u)>0))
 nm=np.r_[.5*(m[:-1]+m[1:]),.5*m[-1]];K=.5*nm@u[1:]**2;vrms=math.sqrt(2*K/M)
 cum=np.r_[0,np.cumsum(w)];quant={}
 for f in [.01,.05,.5,.95,.99,.999]:
  mass_inside=1-f;i=min(len(m)-1,np.searchsorted(cum,mass_inside,side='right')-1)
  frac=(mass_inside-cum[i])/w[i]
  rr=(r[i]**3+frac*(r[i+1]**3-r[i]**3))**(1/3)
  uu=u[i]+(u[i+1]-u[i])*(rr-r[i])/(r[i+1]-r[i])
  quant[str(f)]=(R-rr)/uu
 rc=.75*np.diff(r**4)/np.diff(r**3);uc=u[:-1]+np.diff(u)*(rc-r[:-1])/np.diff(r)
 contact=[]
 for density in [1000.,2000.]:
  def rho_at(t):return m/(4*np.pi/3*np.diff((r+u*t)**3))
  t=brentq(lambda t:rho_at(t).max()-density,0,1,xtol=1e-16)
  rr=r+u*t
  check(d['stage']+' contact density '+str(density),abs(rho_at(t).max()/density-1)<1e-9)
  contact.append({'max_density_kg_m3':density,'lookback_s':t,'outer_radius_m':rr[-1],
   'mean_density_kg_m3':M/(4*np.pi/3*rr[-1]**3)})
 rating=collector['stages'][cat]['source_rating']['peak_total_prepared_electric_W']/1.1
 rows=[]
 for n in [100,1000,10000,100000]:
  cap=math.acos(1-2/n);neighbor=math.sqrt(4*np.pi/n)
  # Cone-average transverse energy in a parallel beam made locally radial.
  esector=vrms**2/n*(1-2/(3*n))
  epair=.5*vrms**2*math.sin(neighbor/2)**2
  check(d['stage']+' cap solid angle '+str(n),abs(2*np.pi*(1-math.cos(cap))/(4*np.pi/n)-1)<1e-11)
  rows.append({'launchers':n,'equal_area_cap_halfangle_deg':math.degrees(cap),
   'typical_neighbor_angle_deg':math.degrees(neighbor),'source_sector_radius_m':R*math.sin(cap),
   'footprint_radius_at_2000_density_backmap_m':contact[-1]['outer_radius_m']*math.sin(cap),
   'mass_per_launcher_kg':M/n,'peak_prepared_electric_W_per_launcher':rating/n,
   'cone_average_transverse_energy_J_kg':esector,
   'cone_energy_fraction_initial_kinetic':esector*M/d['initial_kinetic_energy_J'],
   'cone_energy_fraction_initial_internal':esector*M/d['initial_internal_energy_J'],'neighbor_pair_merger_energy_J_kg':epair})
 stages[d['stage']]={'source_cat':cat,'mass_kg':M,'reference_radius_m':d['R0_m'],
  'initial_internal_energy_J':d['initial_internal_energy_J'],'initial_kinetic_energy_J':d['initial_kinetic_energy_J'],
  'hot_initial_outer_radius_m':r[-1],'v_rms_m_s':vrms,'fastest_speed_m_s':u[-1],
  'mass_flight_quantiles_s':quant,'middle90_launch_window_s':quant['0.95']-quant['0.05'],
  'central_cell_centroid_flight_s':float((R-rc[0])/uc[0]),
  'ballistic_focus_times_from_hot_snapshot_s':{'min':float(np.min(r[1:]/u[1:])),
   'max':float(np.max(r[1:]/u[1:]))},
  'packing_density_backmaps':contact,'count_comparisons':rows,
  'timing_from_frozen_collector_at_R':{k:{a:b*1.1 for a,b in row.items()} for k,row in collector['stages'][cat]['timing_at_R1e9'].items()},
  'peak_prepared_electric_W':rating,
  'snapshot_error_budget_examples':[{'allowed_position_error_m':x,'angle_rad':x/R,
    'uncorrected_timing_s':x/vrms,'uncorrected_fractional_velocity_error':x/R} for x in [1.,10.,30.]]}

c13=stages['C13'];M=c13['mass_kg'];window=c13['middle90_launch_window_s'];v=c13['v_rms_m_s']
grains=[]
for radius in [.001,.01,.1]:
 mg=4*np.pi/3*1000*radius**3;total=M/mg
 grains.append({'radius_m':radius,'density_kg_m3':1000,'grain_mass_kg':mg,
  'physical_grains_per_shot':total,'physical_grains_per_launcher_at1000':total/1000,
  'middle90_mean_rate_per_launcher_s':.9*total/(1000*window)})
bundle=[]
for k in [1,64,1000,1e6]:
 P=1000*k;mass=M/P
 bundle.append({'packets_per_launcher':k,'packets_per_shot':P,'mass_per_packet_kg':mass,
  'solid_equivalent_radius_at1000kg_m3_m':(3*mass/(4*np.pi*1000))**(1/3)})
scenario=ops['cases']['high_source_rating_surface'];intervals=scenario['pulse_intervals_s']
u=v;ww=1.1e9/94.33837343481528
overlap=[{'shot_separation_s':delta,'crossing_radius_m':u*ww*delta/(u+ww),
 'crossing_time_after_first_shot_s':u*delta/(u+ww),
 'crossing_within_receiver':u*ww*delta/(u+ww)<R} for delta in [60.,600.,1800.,intervals['n15']]]
out={'scope':'Kinematic and conditional geometric screen, no discrete packets or angular hydrodynamics simulated.',
 'radius_m':R,'reference_launcher_count':NREF,
 'sources':{str(p.relative_to(ROOT)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths},
 'stages':stages,'C13_illustrative_grains':grains,'C13_illustrative_packet_bundles':bundle,
 'current_operating_case':{'name':'high_source_rating_surface','pulse_intervals_s':intervals,
  'all_stage_mean_interval_s':1/sum(1/x for x in intervals.values()),'sequential_occupancy':scenario['sequential_occupancy']},
 'example_opposing_shot_crossings':{'incoming_speed_m_s':u,'outgoing_speed_m_s':ww,
  'minimum_center_to_center_separation_s_to_avoid_this_pair':R/u+R/ww,'cases':overlap}}
for row in overlap:
 t=row['crossing_time_after_first_shot_s'];delta=row['shot_separation_s']
 check('opposing crossing '+str(delta),abs(ww*t-u*(delta-t))/max(1.,ww*t)<1e-12)
# Analytic homogeneous-flow limit: r(q)=r_out*(1-q)^(1/3).
rtest=np.array([0.,1.,2.]);utest=3*rtest;wtest=np.diff(rtest**3)/8.;cumtest=np.r_[0,np.cumsum(wtest)]
for f in [.05,.5,.95,.999]:
 mi=1-f;i=min(1,np.searchsorted(cumtest,mi,side='right')-1);frac=(mi-cumtest[i])/wtest[i]
 rr=(rtest[i]**3+frac*(rtest[i+1]**3-rtest[i]**3))**(1/3)
 check('uniform-density quantile '+str(f),abs(rr/(2*(1-f)**(1/3))-1)<1e-12)
out['validation']={'passed':sum(x['passed'] for x in checks),'total':len(checks),'checks':checks}
out['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'output':str(args.output),'checks':out['validation']['passed'],'stages':list(stages)},indent=2))
