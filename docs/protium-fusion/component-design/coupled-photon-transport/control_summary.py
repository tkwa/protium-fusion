from pathlib import Path
import json,math,hashlib,gzip,statistics
P=Path(__file__).resolve().parent;O=P/'packet';s=json.loads((O/'source128-spectrum.json').read_text());H=s['thermal_boundary_MeV_per_He'];line={x['energy_MeV']:x['MeV_per_He'] for x in s['primary_photon_lines']};w={2.5438523:1.1*line[2.312593]/H,5.657344:1.1*(line[.511]+s['main_beta_kinetic_MeV_per_He'])/H};base=json.loads((P/'final-S500-isotropic.json').read_text());b={x['incident_energy_MeV']:x for x in base['records']};out=[]
for name in ['final-S500-cut1','final-S500-R1e4','final-S500-surface']:
 p=P/(name+'.json');d=json.loads(p.read_text());delta=0;var=0;rs=[]
 for r in d['records']:
  E=r['incident_energy_MeV'];diff=r['leakage_fraction']-b[E]['leakage_fraction'];se=math.hypot(r['sampling_standard_error'],b[E]['sampling_standard_error']);delta+=w[E]*diff;var+=(w[E]*se)**2;rs.append({k:v for k,v in r.items() if k!='events'})
 row={'name':name,'source_sha256':d['code_sha256'],'raw_path':str(p),'raw_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'parameters':d['parameters'],'per_energy':rs,'source128_weighted_leakage_difference_fraction_HplusU':delta,'difference_sampling_standard_error':math.sqrt(var),'difference_sampling_80pct_interval':[delta-1.2815515655446004*math.sqrt(var),delta+1.2815515655446004*math.sqrt(var)],'scope':'Independent-seed comparison against production S500; combines only the two dominant broad delayed-energy channels. Prompt source unchanged. This statistical comparison does not establish a mathematical discretization-error bound.'};out.append(row)
 ev={'parameters':d['parameters'],'source_sha256':d['code_sha256'],'records':[{'energy_MeV':r['incident_energy_MeV'],'values':[[e['deposited_MeV'],e['escaped_calorimetric_MeV'],e['vacuum_deposit_MeV'],e['residual_inside_kinetic_MeV']] for e in r['events']]} for r in d['records']]}
 with (O/'events'/f'{name}.json.gz').open('wb') as f:
  with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as g:g.write(json.dumps(ev,separators=(',',':')).encode())
(O/'controls.json').write_text(json.dumps({'source128_spectrum_sha256':hashlib.sha256((O/'source128-spectrum.json').read_bytes()).hexdigest(),'rows':out},indent=2)+'\n')
print(json.dumps([{k:v for k,v in x.items() if k not in ['per_energy','parameters']} for x in out],indent=2))
