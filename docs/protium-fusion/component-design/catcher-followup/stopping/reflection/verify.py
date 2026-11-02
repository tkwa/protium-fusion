"""Independent packet identities and provenance checks; does not rerun transport."""
from pathlib import Path
import json,hashlib,math
P=Path(__file__).resolve().parent;checks=[]
def ck(name,value):checks.append({'check':name,'passed':bool(value)})
def close(a,b,tol=3e-12):return abs(a-b)<=tol*max(abs(a),abs(b),1e-100)
single=json.loads((P/'single-collision.json').read_text());j=json.loads((P/'comparison.json').read_text());src=json.loads((P/'sources/primary-tables.json').read_text());surf=json.loads((P.parent/'surface-spectra.json').read_text())
for name,d in src['datasets'].items():
 ck(name+' grid',len(d['rows'])==84 and len({(r['E_eV'],r['angle_deg'])for r in d['rows']})==84)
 for r in d['rows']:
  ck(name+f" E{r['E_eV']} theta{r['angle_deg']} units",0<=r['RE']<=r['RN']<=1 and all(0<=v<=1 for v in r['energy_fraction_quantiles_10_30_50_70_90']))
for name,val in j['inputs_sha256'].items():ck('inputhash '+name,hashlib.sha256(Path(name).read_bytes()).hexdigest()==val)
for name,val in single['input_sha256'].items():ck('kernelinputhash '+name,hashlib.sha256(Path(name).read_bytes()).hexdigest()==val)
ck('kernelcodehash',single['code_sha256']==hashlib.sha256((P/'single_collision.py').read_bytes()).hexdigest())
ck('quadraturechecks',single['checks']['number']==14==single['checks']['passed'])
for row in j['stages']:
 name=row['stage']+' '+row['species'];x=surf['stages'][row['stage']];N=sum(x['species'][row['species']]['packet_nuclei_m2']);ck(name+' count',close(N,row['nuclei_m2']));ck(name+' cadence',close(row['mean_incident_flux_m2_s'],N/x['period_s']))
 for key,s in row['scenarios'].items():
  ck(name+key+' branches',close(s['R'],s['low_primary_R']+s['He_5to20keV_primary_R']+s['high_model_R']))
  ck(name+key+' energy',close(s['reflected_KE_J_m2'],N*s['R']*s['mean_first_reflected_E_eV']*1.602176634e-19))
  for cut in [1,100]:
   lo,hi=s[f'exit_below_{cut}eV_fraction_bracket'];ck(name+key+f' cut{cut}',0<=lo<=hi<=s['all_reflected_eventually_warm_fraction_ceiling'])
for species,m in j['cadence_weighted'].items():
 rows=[r for r in j['stages']if r['species']==species];flux=sum(r['mean_incident_flux_m2_s']for r in rows);ck(species+' meanflux',close(flux,m['mean_incident_flux_m2_s']))
 for key,v in m['scenarios'].items():
  expected=sum(r['mean_incident_flux_m2_s']*r['scenarios'][key]['R']for r in rows)/flux;ck(species+key+' weightedR',close(expected,v['R']))
  ck(species+key+' meanenergy',close(v['mean_reflected_power_W_m2'],flux*v['R']*v['mean_first_reflected_E_eV']*1.602176634e-19))
 ck(species+' nuclearthreshold',close(m['effective_backward_cross_section_barn_to_add_0p001_if_all_column_events_escape']*m['column_event_probability_per_barn'],.001))
ck('no bundled report PDF or page image',not any(p.suffix.lower()in['.pdf','.png','.gif']for p in P.rglob('*')if p.is_file()))
result={'scope':__doc__,'number':len(checks),'passed':sum(c['passed']for c in checks),'checks':checks};(P/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result['passed'], '/',result['number']);assert all(c['passed']for c in checks)
files=[p for p in sorted(P.rglob('*'))if p.is_file()and p.name!='manifest.json'and '__pycache__'not in p.parts]
(P/'manifest.json').write_text(json.dumps({'scope':'New reflection subpacket only; frozen parent archives unchanged.','files':{str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files}},indent=2)+'\n')
