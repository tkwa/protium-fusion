"""Broad cleanup magnet interaction, cycle-energy and slow-tail requirements."""
from pathlib import Path
import json,csv,math,hashlib
import numpy as np
from scipy.optimize import brentq
P=Path(__file__).resolve().parent;C=P.parents[1]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
src=C/'final-family-dispatch/sigma350/input.json';events=src.with_name('events.csv');sched=src.with_name('schedule.json');spec=P.parent/'stopping/surface-spectra.json'
d=json.loads(src.read_text());ss=json.loads(spec.read_text());ev=list(csv.DictReader(events.open()));horizon=json.loads(sched.read_text())['count_horizon_s'];coil=json.loads((P/'coil-array.json').read_text());rate=len(ev)/horizon
def duty(full=False,ramp=0.):
 spans=[]
 for row in ev:
  t=float(row['focus_s']);end=t+d['timing'][row['cat']]['full_resolved']['outgoing_s']*1.1 if full else float(row['main_capture_end_s'])
  spans.append((max(0,t-1-ramp),min(horizon,end+10+ramp)))
 out=[]
 for lo,hi in sorted(spans):
  if out and lo<=out[-1][1]:out[-1][1]=max(hi,out[-1][1])
  else:out.append([lo,hi])
 return 1-sum(b-a for a,b in out)/horizon
rows=[];e=1.602176634e-19;mu=1.66053906660e-27;mu0=4*math.pi*1e-7
def ramp_case(U,t):
 live=duty(False,t);me=9.1093837139e-31;eps=8.8541878188e-12;v=2299.93667;r=.001;L=.01;dk=2*r*v/3;entrance=.5*v/4;ve=math.sqrt(200*e/me);k=3.5e-21*ve
 def S(nu):
  si=math.sqrt(dk*nu)*math.tanh(L*math.sqrt(nu/dk));return si*entrance/(si+entrance)
 nu=brentq(lambda n:4*.5*live*S(n)/v-.01,1e-12,1e8);ne=nu/k;w=math.sqrt(80*e/me)/r;Breq=max(.2,math.sqrt(4*me*ne/eps),me/e*(w+ne*e*e/(2*eps*me*w)))
 return {'each_ramp_s':t,'minimum_mean_ramp_converter_W_m2':U/t,'converter_kg_m2_at1kW_kg':U/t/1000,'main_mass999_live_fraction_including_ramps':live,'g_with_original_cloud_unchanged':.01*live/duty(False),'required_ne_m3_to_restore_g001':ne,'required_Bmin_T_to_restore_g001_in_same_rotation_screen':Breq,'scope':'Ramp energy/rating is evaluated at the archived base field. Restoringg=.01 at lowerduty requires re-solving cable/current/shield sizing at the listednewB, not silently retainingtheold hardware.'}
for r in coil['rows']:
 if r['status']!='finite_fixed_point' or r['Je_A_mm2']!=30 or r['pitch_m']not in[10.,50.,100.] or r['shield_radial_column_kg_m2']!=2000:continue
 U=r['return_field_energy_J_m2'];p=r['pitch_m'];dep=r['coil_center_depth_below_hot_face_m'];I=r['current_A_per_rail'];Bint=mu0*I/(2*math.pi)*math.log(1/math.tanh(math.pi*dep/p))
 rows.append({'pitch_m':p,'field_energy_J_m2':U,'unsigned_charge_and_discharge_throughput_W_m2':2*U*rate,'cycle_loss_W_m2':{str(f):f*U*rate for f in[.0001,.001,.01,.1,1.]},'cold_fraction_of_cycle_loss':{'scope':'Loss location matters: dissipation at65K also incurs1/COP refrigerator work. Room-temperature recovery electronics are separately cooled.','refrigerator_multiplier_for_cold_loss':1/coil['parameters']['COP65to200']},'ramp_ratings':[ramp_case(U,t)for t in[10.,100.,1000.]],'independent_energy_store_kg_m2_at100kJ_kg':U/1e5,'integrated_tangential_B_Tm_on_ray_above_rail':Bint,'1MeV_proton_rigidity_Tm':math.sqrt(2*mu*1e6*e)/e,'off_field_fraction_required_for0p01rad_1MeV_proton_reference':.01*math.sqrt(2*mu*1e6*e)/e/Bint,'local_module_energy_J_at_length100p':U*p*(100*p),'ampere_turns_per_rail':I,'series_turns_reference_at20kA_supply':math.ceil(I/20000)})
tails=[]
for stage,x in ss['stages'].items():
 cat=x['catalyst'];cut=d['timing'][cat]['mass999']['outgoing_s']*1.1;tm=np.array(x['packet_arrival_after_ignition_s']);late=tm>cut;m=np.array(x['packet_mass_kg_m2_baryon']);v=np.array(x['packet_velocity_m_s']);dv=np.array(x['packet_delta_v_weight_m_s']);jac=np.array(x['packet_arrival_jacobian_m_s2']);pressure=m/dv*jac*v
 tails.append({'stage':stage,'cut_after_focus_s':cut,'full_resolved_time_s':d['timing'][cat]['full_resolved']['outgoing_s']*1.1,'smooth_quadrature_mass_fraction_later_than_mass999_cut':float(m[late].sum()/m.sum()),'peak_ballistic_ram_pressure_Pa_sampled':float(pressure.max()),'species_number_fraction_in_late_packets':{name:float(np.array(vv['packet_nuclei_m2'])[late].sum()/sum(vv['packet_nuclei_m2']))for name,vv in x['species'].items() if name in['p','he4']}})
out={'scope':__doc__,'code_sha256':sha(__file__),'source_sha256':{str(p):sha(p)for p in[src,events,sched,spec,P/'coil-array.json']},'event_rate_s':rate,'events_per_year':rate*365.25*86400.,'live_mass999_with10s_recovery':duty(False),'live_full_resolved_with10s_recovery':duty(True),'field_Bmin_T':coil['parameters']['Bmin_front_T'],'magnetic_pressure_at_Bmin_Pa':coil['parameters']['Bmin_front_T']**2/(2*mu0),'rows':rows,'tails':tails,'conclusions':['Static broad fields are not a drop-in addition to the ballistic receiver. Outgoing charged ions are strongly turned; a locally quasineutral plasma does not remove collective magnetic stresses. Neutral trajectories require their actual charge-state distribution.','The cycle rows are conditional switching requirements, not a completed converter or AC-loss design. Off during main outgoing mass leaves charged slow tails that need a separate field-guided collection/routing calculation. Waiting for every resolved tail has almost no global cleanup live time because successiveB1 tails overlap.','Turning off only the electron population does not remove the magnet field. An alternative is a compact closed-flux magnet that keeps its integrated transverse field below the actual particle-transport acceptance, with its cold/hot pole placement and radiation/thermal mass accounted.','The99.9%main-mass timing is not a loss allocation. Diverted tails remain inside the closed station and may be recovered; their destination, gas release and collector heating must be explicitly modeled.','One-percent per-cycle loss means0.01 timesstoredenergy for the combined discharge/charge cycle. This is not independently applying1% to each leg. Recovery, switch/conductorAC loss and cryogenic dissipation must be separated.','A putative300kV collection cathode addressesion destination afterionization; it doesnot solve broadfield interference, cloudconfinement or molecular admission.']}
assert abs(out['live_mass999_with10s_recovery']-.7970058470499236)<1e-12
(P/'switching.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
