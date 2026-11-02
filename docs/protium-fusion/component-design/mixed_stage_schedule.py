"""Chronological dispatch, finite stocks and conditional central-seed timing.

Main delivery/capture windows retain the frozen mass99.9% criterion. No material
is discarded: omitted incoming core is assigned a faster assembled seed, and
outgoing tails stay in the material inventory pending collection. The seed's
hydrodynamic equivalence and tail interaction remain separate qualifications.
"""
import argparse,json,hashlib,csv,math
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--events-csv',type=Path,required=True);ap.add_argument('--days',type=float,default=370);a=ap.parse_args();R=a.root.resolve()
timing=json.loads((R/'launcher-reference/incoming-timing.json').read_text());cache=json.loads((R/'open-loop-redesign/collector-inputs.json').read_text());ops=json.loads((R/'open-loop-redesign/operating-cases.json').read_text())['cases']['high_source_rating_surface'];cycle=json.loads((R/'open-loop-redesign/station-inputs.json').read_text())['cycle'];material=json.loads((R/'component-design/material-cycle-validation.json').read_text());stages=material['loaded_stages'];H=ops['helium_per_s'];Radius=ops['radius_m'];DAY=86400.;horizon=a.days*DAY
phases={'c12':.37,'c13':.05,'n14':.12,'n15':.5};raw=[];CAT={s['stage']:k for k,s in cycle['stage_pulse_data'].items()}
for cat,P in ops['pulse_intervals_s'].items():
 tag=stages[cat]['label'];w=timing['stages'][tag]['timing_from_frozen_collector_at_R']['mass999']
 for n in range(int(horizon/P)+2):
  nominal=(n+phases[cat])*P
  if nominal>horizon:continue
  raw.append({'cat':cat,'stage':tag,'nominal_focus_s':nominal,'incoming_s':w['incoming_s'],'outgoing_s':w['outgoing_s']})
raw.sort(key=lambda x:x['nominal_focus_s']);end=0.;events=[]
for j in raw:
 start=max(j['nominal_focus_s']-j['incoming_s'],end);focus=start+j['incoming_s'];end=focus+j['outgoing_s']
 events.append({**j,'main_launch_start_s':start,'focus_s':focus,'main_capture_end_s':end,'dispatch_delay_s':focus-j['nominal_focus_s']})
checks=[]
def check(name,value):
 checks.append({'name':name,'passed':bool(value)})
 if not value:raise AssertionError(name)
check('every main window disjoint',all(b['main_launch_start_s']>=p['main_capture_end_s']-1e-7 for p,b in zip(events,events[1:])))
check('none dispatched before requested nominal focus',all(j['dispatch_delay_s']>=-1e-7 for j in events))
stats={}
for cat in stages:
 es=[j for j in events if j['cat']==cat];tt=np.array([j['focus_s'] for j in es]);stats[cat]={'stage':stages[cat]['label'],'pulse_count':len(es),'target_count':horizon/ops['pulse_intervals_s'][cat],
  'min_focus_interval_s':float(np.diff(tt).min()) if len(tt)>1 else None,'max_focus_interval_s':float(np.diff(tt).max()) if len(tt)>1 else None,
  'max_dispatch_delay_s':max([j['dispatch_delay_s'] for j in es] or [0])}
# Exact input isotope counts and cold output amounts, finite atomic stock.
MU=1.6605390666e-27;U=931.49410242
# Obtain atomic masses indirectly from the already validated reference values;
# nuclear-mass corrections here are below0.1%, so explicit baryonic kg used.
def massA(name):return 1 if name in ['p','n'] else int(''.join(x for x in name if x.isdigit()))
fresh={n:material['working60']['external_input_atoms_per_He'].get(n,0.)*H*massA(n)*MU for n in ['c12','c13','n14','n15','p']}
stockcases=[]
for hold_days in [1,30,365]:
 changes=[]
 for e in events:
  s=stages[e['cat']];attempt=s['attempts_per_pulse']
  draw={k:-s['input_counts'].get(k,0.)*attempt*massA(k)*MU for k in fresh}
  ret={k:s['cold_output'].get(k,0.)*attempt*massA(k)*MU*(cycle['minimum_CNO_return_fraction'] if k!='p' else 1-cycle['hydrogen_loss_per_cold_stage_output']) for k in fresh}
  changes.append((e['main_launch_start_s'],draw));changes.append((e['main_capture_end_s']+hold_days*DAY,ret))
 changes.sort(key=lambda x:x[0]);running={k:0. for k in fresh};minimum=running.copy();maximum=running.copy();last=0.
 for t,ds in changes:
  if t>horizon:break
  for k in fresh:
   running[k]+=fresh[k]*(t-last);maximum[k]=max(maximum[k],running[k]);running[k]+=ds.get(k,0.);minimum[k]=min(minimum[k],running[k]);maximum[k]=max(maximum[k],running[k])
  last=t
 init={k:max(0.,-v) for k,v in minimum.items()}
 stockcases.append({'hold_days_after_main_capture':hold_days,'minimum_starting_atomic_stock_baryonic_kg':init,'total_minimum_starting_stock_kg':sum(init.values()),
 'maximum_ready_stock_kg_by_species':{k:maximum[k]+init[k] for k in fresh},'horizon_end_stock_kg_by_species':{k:running[k]+init[k] for k in fresh},
 'steady_gross_immobilized_fuel_kg':ops['all_stage_flow_kg_s']*hold_days*DAY,
 'scope':'Historical maximum-loss boundary case; finite370day history, exact integer dispatch and complete delayed-return batches; annual rotating-skin operation also needs a prefilled age distribution to reach steady state.'})
# Incoming central mass with long ballistic launch history: exact cumulative
# mass/velocity integrals for uniform rho per cell and linear v(r).
seed={}
for cat,x in cache['stages'].items():
 d=x['record'];p=x['initial'];m=np.asarray(p['mass_fraction_of_total'])*d['mass_kg'];rr=np.r_[0.,p['outer_radius_m']];vv=np.r_[0.,-np.asarray(p['velocity_m_s'])]
 target=.001*d['mass_kg'];remain=target;ks=0.;seeds=[]
 for i,mm in enumerate(m):
  take=min(mm,remain)
  if take<=0:break
  rout=(rr[i]**3+take/mm*(rr[i+1]**3-rr[i]**3))**(1/3);s=(vv[i+1]-vv[i])/(rr[i+1]-rr[i]);b=vv[i]-s*rr[i]
  F=lambda z:s*s*z**5/5+s*b*z**4/2+b*b*z**3/3
  ki=.5*3*mm/(rr[i+1]**3-rr[i]**3)*(F(rout)-F(rr[i]));ks+=ki;seeds.append({'cell':i,'mass_kg':take,'outer_radius_m':rout,'max_inward_m_s':s*rout+b,'kinetic_J':ki});remain-=take
 seed[d['stage']]={'core_mass_fraction':.001,'core_mass_kg':target,'core_initial_kinetic_J':ks,'core_K_fraction_total_prepared':ks/d['initial_energy_J'],
 'core_outer_radius_m':seeds[-1]['outer_radius_m'],'core_outer_velocity_m_s':seeds[-1]['max_inward_m_s'],'cells':seeds,
 'faster_seed_extra_energy_bound':{str(speed):{'gross_carrier_K_J':.5*target*speed*speed,'fraction_existing_prepared':.5*target*speed*speed/d['initial_energy_J']} for speed in [1e6,2e6,3e6]},
 'scope':'These gross K scales must first be compared with existing core internal-energy allowance; no additive cost is asserted when collision heating replaces that prepared energy.'}
# How much old outward material remains inside R at the next main window/focus?
loading=json.loads((R/'membrane-receiver/pulse-loading.json').read_text());pairs={}
def remaining(tag,t):
 c=loading['stages'][tag]['cells'];m=np.asarray(c['mass_kg']);rb=np.asarray(c['radius_bounds_m']);vb=np.asarray(c['velocity_bounds_m_s']);result=0.
 for mm,r0,r1,v0,v1 in zip(m,rb[:-1],rb[1:],vb[:-1],vb[1:]):
  s=(v1-v0)/(r1-r0);orig=(Radius-(v0-s*r0)*t)/(1+s*t)
  f=np.clip((orig**3-r0**3)/(r1**3-r0**3),0,1);result+=mm*f
 return float(result)
for p,e in zip(events,events[1:]):
 key=p['stage']+'→'+e['stage'];gap=e['focus_s']-p['focus_s'];oldmain=remaining(p['stage'],e['main_launch_start_s']-p['focus_s']);oldfocus=remaining(p['stage'],gap)
 row={'focus_separation_s':gap,'outgoing_tail_at_next_launch_kg':oldmain,'outgoing_tail_at_next_focus_kg':oldfocus,'tail_mass_over_next_pulse_at_focus':oldfocus/stages[e['cat']]['mass_kg']}
 if key not in pairs or row['tail_mass_over_next_pulse_at_focus']>pairs[key]['tail_mass_over_next_pulse_at_focus']:pairs[key]=row
with a.events_csv.open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(events[0]));writer.writeheader();writer.writerows(events)
firstC=next(j['focus_s'] for j in events if j['stage']=='C13');excerpt=[j for j in events if abs(j['focus_s']-firstC)<2*DAY]
result={'scope':__doc__,'days':a.days,'radius_m':Radius,'phase_choices':phases,'dispatch_rule':'Nominal periodic jobs sorted by requested focus; delay each main launch until preceding main capture ends. Positive20%mean slack absorbsrare jobs; nominal stage rates unchanged in long run.',
 'events_csv':a.events_csv.name,'event_count':len(events),'target_mean_all_stage_interval_s':1/sum(1/v for v in ops['pulse_intervals_s'].values()),
 'occupied_fraction':sum(j['incoming_s']+j['outgoing_s'] for j in events)/horizon,'max_dispatch_delay_s':max(j['dispatch_delay_s'] for j in events),'stage_stats':stats,
 'excerpt_near_first_C13':excerpt,'stock_cases':stockcases,'central_seed_mass_and_energy':seed,'worst_adjacent_tail_pairs':pairs,
 'tail_return_requirement':'All residual outgoing material stays in material inventory; harmless passage, re-entrainment into later pulse, or recovery must conserve it. No discard credit is taken. Separate collision-column calculation required.',
 'checks':checks,'source_sha256':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/'launcher-reference/incoming-timing.json',R/'open-loop-redesign/collector-inputs.json',R/'open-loop-redesign/operating-cases.json',R/'open-loop-redesign/station-inputs.json',R/'membrane-receiver/pulse-loading.json']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['event_count','occupied_fraction','max_dispatch_delay_s','stage_stats','stock_cases','worst_adjacent_tail_pairs']},indent=2))
