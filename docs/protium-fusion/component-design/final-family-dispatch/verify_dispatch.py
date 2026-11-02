"""Independent stdlib reconstruction of archived dispatch and source alignment."""
from pathlib import Path
import csv,json,math,hashlib
P=Path(__file__).resolve().parent;checks=[]
def ck(n,b):
 checks.append({'name':n,'passed':bool(b)})
 if not b:raise AssertionError(n)
def eq(x,y):return math.isclose(x,y,rel_tol=1e-12,abs_tol=1e-7)
for sigma in [250,350]:
 S=P/f'sigma{sigma}';d=json.loads((S/'input.json').read_text());q=json.loads((S/'schedule.json').read_text());p=d['plant'];cy=d['cycle'];ev=list(csv.DictReader((S/'events.csv').open()));H=p['helium_per_s'];h=q['count_horizon_s'];prev=0;occupation=0;names={}
 for e in ev:
  cat=e['cat'];w=d['timing'][cat];scale=p['radius_m']/w['radius_m'];inc=w['mass999']['incoming_s']*scale;out=w['mass999']['outgoing_s']*scale;nom=float(e['nominal_focus_s']);start=max(nom-inc,prev);focus=start+inc;end=focus+out
  if not all(eq(x,float(e[k])) for k,x in [('main_launch_start_s',start),('focus_s',focus),('main_capture_end_s',end),('incoming_s',inc),('outgoing_s',out),('dispatch_delay_s',focus-nom)]):raise AssertionError((sigma,e))
  assert start>=prev-1e-7;prev=end;occupation+=inc+out;names[cat]=names.get(cat,0)+1
 ck(f'{sigma} all greedy chronological events',True);ck(f'{sigma} complete event count',len(ev)==q['event_count']);ck(f'{sigma} finite occupancy',eq(occupation/h,q['occupied_fraction_requested_horizon']));ck(f'{sigma} completion horizon',eq(max(h,prev),q['dispatch_end_s']));stationary=0.
 for cat,s in d['loaded_stages'].items():
  period=1/(H*cy['stage_pulse_data'][cat]['pulses_per_He']);ck(f'{sigma} {cat} cadence',eq(period,p['pulse_interval_s'][cat]));times=[float(e['focus_s']) for e in ev if e['cat']==cat];stats=q['stage_stats'][cat];ck(f'{sigma} {cat} count',len(times)==stats['count']);ck(f'{sigma} {cat} rounding',abs(len(times)-h/period)<1);ck(f'{sigma} {cat} period',eq(period,stats['nominal_period_s']));ck(f'{sigma} {cat} mixture',s['input_counts']==cy['stage_pulse_data'][cat]['input_counts']);w=d['timing'][cat];stationary+=(w['mass999']['incoming_s']+w['mass999']['outgoing_s'])*p['radius_m']/w['radius_m']/period
 ck(f'{sigma} stationary occupancy',eq(stationary,p['sequential_occupancy_mass999']));ck(f'{sigma} eventual CNO return',cy['carbon_return']==1);ck(f'{sigma} eventual H return',cy['hydrogen_return']==1);ck(f'{sigma} source no extra miss',d['selected_row']['additional_missed_hot_fraction_HplusUprecursor']==0);ck(f'{sigma} column',d['selected_row']['receiver_sigma_kg_m2']==sigma);ck(f'{sigma} cold routing',d['selected_row']['missed_hot_routing']=='cold_equipment')
print(json.dumps({'checks_passed':len(checks),'checks':checks,'scope':'Chronological reconstruction, actual-composition/cadence and source-policy checks. Stored stock cases are the unchanged existing material scheduler; this audit does not validate physical tail trajectories or isotope extraction hardware.'},indent=2))
