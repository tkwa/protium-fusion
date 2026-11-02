"""Chronological dispatch and finite stock screen for an archived mixture plant.

99.9% main windows are serialized. For each species, the full omitted 0.1% pulse mass is conservatively withheld
through the finite study horizon as retained slow-tail inventory.
Long solid holds and the local hydrogen gas hold are distinct. Central seed
assembly and the eventual wall-to-source geometry remain engineering inputs.
"""
import argparse,json,csv,hashlib,math
from pathlib import Path
import material_cycle as mc
ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--events-csv',type=Path,required=True);ap.add_argument('--days',type=float,default=370);a=ap.parse_args()
d=json.loads(a.input.read_text());cy=d['cycle'];p=d['plant'];st=d['loaded_stages'];H=p['helium_per_s'];DAY=86400.;horizon=a.days*DAY
phases={'c12':.37,'c13':.05,'n14':.12,'n15':.5};raw=[];windows={}
for cat,s in st.items():
 w=d['timing'][cat];scale=p['radius_m']/w['radius_m'];windows[cat]={k:w['mass999'][k]*scale for k in ['incoming_s','outgoing_s','occupied_s']}
 period=p['pulse_interval_s'][cat]
 for n in range(int(horizon/period)+2):
  nominal=(n+phases[cat])*period
  if nominal<=horizon:raw.append({'cat':cat,'stage':s['label'],'nominal_focus_s':nominal})
raw.sort(key=lambda x:x['nominal_focus_s']);end=0.;events=[]
for j in raw:
 w=windows[j['cat']];start=max(j['nominal_focus_s']-w['incoming_s'],end);focus=start+w['incoming_s'];end=focus+w['outgoing_s']
 events.append({**j,'incoming_s':w['incoming_s'],'outgoing_s':w['outgoing_s'],'main_launch_start_s':start,'focus_s':focus,'main_capture_end_s':end,'dispatch_delay_s':focus-j['nominal_focus_s']})
endh=max(horizon,events[-1]['main_capture_end_s']);checks=[]
def ck(n,b):
 checks.append({'name':n,'passed':bool(b)})
 if not b:raise AssertionError(n)
ck('disjoint main windows',all(b['main_launch_start_s']>=c['main_capture_end_s']-1e-7 for c,b in zip(events,events[1:])))
ck('nonnegative dispatch delay',min(e['dispatch_delay_s'] for e in events)>=-1e-7)
stats={}
for cat,s in st.items():
 es=[e for e in events if e['cat']==cat];times=[e['focus_s'] for e in es];dt=[y-x for x,y in zip(times,times[1:])]
 stats[cat]={'stage':s['label'],'count':len(es),'target_count':horizon/p['pulse_interval_s'][cat], 'nominal_period_s':p['pulse_interval_s'][cat],'min_interval_s':min(dt) if dt else None,'max_interval_s':max(dt) if dt else None,'max_dispatch_delay_s':max(e['dispatch_delay_s'] for e in es)}
 ck(cat+' integer count within one expected shot',abs(len(es)-horizon/p['pulse_interval_s'][cat])<1)
names=[*mc.DEFAULT_STATES,'p'];fresh={n:cy['external_input_atoms_per_He'].get(n,0)*H*mc.atomic_kg(n) for n in names};stock=[]
for solid_days in [1,30,365.25]:
 changes=[];tail={n:0. for n in names}
 for e in events:
  s=st[e['cat']];aa=s['attempts_per_pulse'];draw={n:-s['input_counts'].get(n,0)*aa*mc.atomic_kg(n) for n in names}
  out={n:s['cold_output'].get(n,0)*aa*mc.atomic_kg(n) for n in names}
  changes.append((e['main_launch_start_s'],draw))
  withheld={n:min(out[n],.001*1.01*s['mass_kg']) for n in names}
  changes.append((e['main_capture_end_s']+solid_days*DAY,{n:out[n]-withheld[n] for n in names if n!='p'}))
  changes.append((e['main_capture_end_s']+3600,{'p':out['p']-withheld['p']}))
  for n in names:tail[n]+=withheld[n]
 changes.sort(key=lambda x:x[0]);run={n:0. for n in names};minimum=run.copy();maximum=run.copy();last=0.
 for t,change in changes:
  if t>endh:break
  for n in names:
   run[n]+=fresh[n]*(t-last);maximum[n]=max(maximum[n],run[n]);run[n]+=change.get(n,0);minimum[n]=min(minimum[n],run[n]);maximum[n]=max(maximum[n],run[n])
  last=t
 for n in names:run[n]+=fresh[n]*(endh-last)
 ready={n:max(0,-minimum[n]) for n in names}
 stock.append({'solid_hold_days':solid_days,'H_local_hold_s':3600,'minimum_starting_ready_atomic_kg':ready,'total_starting_ready_kg':sum(ready.values()),'max_ready_atomic_kg':{n:maximum[n]+ready[n] for n in names},'endpoint_ready_atomic_kg':{n:run[n]+ready[n] for n in names},'uncredited_slow_tail_atomic_kg':tail,'scope':'Finite chronological upper stock requirement: For every species, up to the full omitted 0.1% of pulse mass (plus a 1% atomic/baryonic mass margin) remains uncredited throughout the horizon; these independent upper bounds overallocate retained tails, with no fresh-loss debit. A steady rotating-skin plant also needs its resident age distribution.'})
firstC=next(e['focus_s'] for e in events if e['stage']=='C13')
with a.events_csv.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(events[0]));w.writeheader();w.writerows(events)
out={'scope':__doc__,'selected_cases':p['selected_cases'],'days_requested':a.days,'count_horizon_s':horizon,'dispatch_end_s':endh,'radius_m':p['radius_m'],'phase_choices':phases,'event_count':len(events),'occupied_fraction_requested_horizon':sum(e['incoming_s']+e['outgoing_s'] for e in events)/horizon,'continuous_occupancy':p['sequential_occupancy_mass999'],'stage_stats':stats,'max_dispatch_delay_s':max(e['dispatch_delay_s'] for e in events),'stock_cases':stock,'excerpt_near_first_C13':[e for e in events if abs(e['focus_s']-firstC)<2*DAY],'events_csv':a.events_csv.name,'checks':checks,'sources_sha256':{a.input.name:mc.hashfile(a.input)},'code_sha256':mc.hashfile(__file__),'module_sha256':mc.hashfile(mc.__file__)}
a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps({k:out[k] for k in ['event_count','occupied_fraction_requested_horizon','continuous_occupancy','stage_stats','max_dispatch_delay_s','stock_cases']},indent=2))
