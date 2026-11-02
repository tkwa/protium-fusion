"""Independent totals and time-grid check for the explicit routing override."""
from pathlib import Path
import json,math,hashlib
P=Path(__file__).resolve().parent
d=json.loads((P/'inputs.json').read_text());routes=json.loads((P/'route-inputs.json').read_text());checks=[]
def check(name,condition,value=None):
 checks.append({'name':name,'passed':bool(condition),'value':value})
 assert condition,(name,value)
for tag in ['cap0.1_eta100000.0','cap0.01_eta100000.0']:
 r=json.loads((P/f'chamber-dt150.0-q128-{tag}.json').read_text())
 for j,sp in enumerate(['p','he4']):
  total=0
  for cat,s in routes['cases'][tag].items():
   z=s['species'][sp];raw=sum(z['leaked_atomic_kg_per_packet']);mass=sum(max(x,0.) for x in z['leaked_atomic_kg_per_packet']);full=d['stages'][cat]['cold_H_atomic_kg' if sp=='p' else 'cold_He_atomic_kg']
   check(f'{tag}/{cat}/{sp} raw routed mass',abs(raw/full/z['total_front_fraction']-1)<1e-6,raw/full/z['total_front_fraction']-1)
   check(f'{tag}/{cat}/{sp} bounded probability roundoff',min(z['packet_front_fraction'])>=-1e-10 and max(z['packet_front_fraction'])<=1 and (mass-raw)/full<1e-10,{'minimum_probability':min(z['packet_front_fraction']),'emission_correction_kg':mass-raw,'correction_per_full_species_mass':(mass-raw)/full})
   check(f'{tag}/{cat}/{sp} positive and causal',mass>=0 and min(s['arrival_after_ignition_s'])>0)
   total+=mass/z['period_s']
  check(f'{tag}/{sp} steady source',math.isclose(total,r['cavity_source_mass_rate_atomic_kg_s'][j],rel_tol=1e-12))
 check(tag+' exact chamber source hash',r['code_sha256']==hashlib.sha256((P/'chamber.py').read_bytes()).hexdigest())
 check(tag+' routing source hash',r['routing_code_sha256']==hashlib.sha256((P/'route_response.py').read_bytes()).hexdigest())
a=json.loads((P/'chamber-dt150.0-q128-cap0.1_eta100000.0.json').read_text());b=json.loads((P/'chamber-dt300.0-q128-cap0.1_eta100000.0.json').read_text());errors=[]
for x,y in zip(a['cases'],b['cases']):
 error=max(abs(i['maximum_column_kg_m2']/j['maximum_column_kg_m2']-1) for cat in x['source_profile_pellet_encounters'] for i,j in zip(x['source_profile_pellet_encounters'][cat],y['source_profile_pellet_encounters'][cat]))
 peak=abs(x['peak_radial_column_kg_m2']/y['peak_radial_column_kg_m2']-1)
 errors.append({'model':x['model'],'g':x['effective_wall_capture_g'],'delay_s':x['material_release_time_s'],'max_encounter_relative_change':error,'peak_relative_change':peak})
 check('routed refinement '+str(len(errors)),max(error,peak)<.001,max(error,peak))
out={'checks_passed':len(checks),'checks':checks,'refinement':errors,'source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['inputs.json','route-inputs.json','chamber-dt150.0-q128-cap0.1_eta100000.0.json','chamber-dt300.0-q128-cap0.1_eta100000.0.json']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'routing-validation.json').write_text(json.dumps(out,indent=2)+'\n');print('routing checks',len(checks),'max grid error',max(max(x['max_encounter_relative_change'],x['peak_relative_change']) for x in errors))
