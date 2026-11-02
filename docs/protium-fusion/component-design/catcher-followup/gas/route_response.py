"""Apply actual packet-specific ideal mobile routing to the chronological gas model."""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import chamber as c
P=Path(__file__).resolve().parent;r=json.loads((P/'route-inputs.json').read_text());d=json.loads((P/'inputs.json').read_text());tag=sys.argv[1] if len(sys.argv)>1 else 'cap0.1_eta100000.0';dt=float(sys.argv[2]) if len(sys.argv)>2 else 150.;stages=r['cases'][tag]
def emitted(d,sp,dt,points):
 end=370*86400+dt;t=np.arange(0,end+dt,dt);out=np.zeros(len(t));curves={}
 for cat,s in stages.items():
  tt=np.array(s['arrival_after_ignition_s']);raw=np.array(s['species'][sp]['leaked_atomic_kg_per_packet']);amount=np.maximum(raw,0.);curves[cat]={'routing_case':tag,'route_inputs_sha256':hashlib.sha256((P/'route-inputs.json').read_bytes()).hexdigest(),'quadrature_packet_count':len(tt),'total_front_atomic_kg':float(sum(amount)),'raw_fraction_of_full_stage_species':s['species'][sp]['total_front_fraction'],'roundoff_negative_weight_correction_kg':float(sum(amount)-sum(raw)),'raw_minimum_exit_probability':float(min(s['species'][sp]['packet_front_fraction']))}
  for e in d['events']:
   if e['cat']!=cat:continue
   i=np.ceil((e['focus_s']+tt)/dt).astype(int);ok=(i>=0)&(i<len(out));np.add.at(out,i[ok],amount[ok])
 return t,out,curves
Q=[sum(sum(max(x,0.) for x in s['species'][sp]['leaked_atomic_kg_per_packet'])/s['species'][sp]['period_s'] for s in stages.values()) for sp in ['p','he4']]
c.main(dt,128,emission_override=emitted,Q_override=Q,suffix='-'+tag)
outpath=P/f'chamber-dt{dt}-q128-{tag}.json';out=json.loads(outpath.read_text());out['routing_code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();out['routing_input_sha256']=hashlib.sha256((P/'route-inputs.json').read_bytes()).hexdigest();out['arrival_quadrature_scope']='All4096 archived packet nodes per stage with species-specific routed weights; arrival_quantile_points is a baseline API parameter unused by this explicit override. Tiny negative finite-cap probability roundoff is projected to zero only in the auxiliary emitted-gas weights; per-stage absolute corrections are recorded. Original routing data and all hydro abundances are unchanged.';outpath.write_text(json.dumps(out,indent=2)+'\n')
