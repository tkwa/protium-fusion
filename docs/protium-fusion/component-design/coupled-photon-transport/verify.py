"""Rebuild sufficient statistics from archived independent event-energy arrays."""
from pathlib import Path
import json,gzip,math,hashlib
from reweight import evaluate
P=Path(__file__).resolve().parent

def main():
 b=json.loads((P/'basis.json').read_text());checks=[];Nall=0;emax=0.;cache={}
 def check(name,v):
  checks.append({'name':name,'passed':bool(v)})
  if not v:raise AssertionError(name)
 for S,node in b['columns'].items():
  for r in node['responses']:
   p=P/r['event_archive'];key=str(p)
   if key not in cache:
    cache[key]=json.loads(gzip.decompress(p.read_bytes()))
    rawname=r['raw_name'];check(rawname+':archivehash',hashlib.sha256(p.read_bytes()).hexdigest()==b['raw_inputs'][rawname]['event_archive_sha256'])
   z=cache[key]['records'][r['event_archive_record']];a=z['values'];E=z['calorimetric_energy_MeV'];N=len(a);Nall+=N;leak=[v[1]/E for v in a];mean=math.fsum(leak)/N;se=math.sqrt(math.fsum((v-mean)**2 for v in leak)/(N*(N-1)));er=max(abs(math.fsum(v)/E-1) for v in a);emax=max(emax,er);tag=f'{S}:{r["particle"]}:{r["energy_MeV"]}'
   check(tag+':N',N==r['histories']);check(tag+':mean',abs(mean-r['leakage_fraction'])<1e-13);check(tag+':SE',abs(se-r['sampling_standard_error'])<1e-13);check(tag+':energy',er<1e-7);check(tag+':no_unfinished_energy',all(v[3]==0 for v in a))
 results=json.loads((P/'weighted-results.json').read_text())
 for fam in ['current128','source128']:
  s=json.loads((P/(fam+'-spectrum.json')).read_text());H=s['thermal_boundary_MeV_per_He'];prompt=s['cycle']['escaped_prompt_MeV_per_He']
  for row in results[fam]['rows']:
   S=row['sigma_kg_m2']
   for scenario in ['nominal_lines_hardbeta','broad_D1p10']:
    x=evaluate(b,s,S,scenario,'EM20',hot_bypass_ppm=6.4);y=evaluate(b,s,S,scenario,'all_to_cold',hot_bypass_ppm=6.4);old=row['scenarios'][scenario+':EM20']
    check(f'{fam}:{S}:{scenario}:saved',abs(x['leakage_fraction']-old['leakage_fraction'])<1e-14)
    for pr in ['EM20','all_to_cold']:
     for c in [0.,.5,1.]:
      q=evaluate(b,s,S,scenario,pr,c,6.4);check(f'{fam}:{S}:{scenario}:{pr}:partition{c}',abs(q['hot_capture_fraction']+q['cold_equipment_fraction']+q['external_escape_fraction']-1)<1e-14)
    check(f'{fam}:{S}:{scenario}:prompt_once',abs(y['leakage_fraction']-x['leakage_fraction']-(prompt-x['prompt_leakage_MeV_per_He'])/H)<1e-14)
    check(f'{fam}:{S}:{scenario}:uplift_available',x['source_uplift_drawn_from_existing_bulk_MeV_per_He']<s['cycle']['bulk_kinetic_MeV_per_He'])
    check(f'{fam}:{S}:{scenario}:sampling_target',x['sampling_80pct_half_width']<.001)
 controls=json.loads((P/'controls.json').read_text());Nc=0
 for r in controls['rows']:
  d=json.loads(gzip.decompress((P/'events'/(r['name']+'.json.gz')).read_bytes()))
  for z,rec in zip(d['records'],r['per_energy']):
   E=z['energy_MeV'];a=z['values'];N=len(a);Nc+=N;vals=[v[1]/E for v in a];mean=math.fsum(vals)/N;check(r['name']+str(E)+':mean',abs(mean-rec['leakage_fraction'])<1e-14);check(r['name']+str(E)+':energy',max(abs(math.fsum(v)/E-1) for v in a)<1e-7)
 out={'checks_passed':sum(x['passed'] for x in checks),'checks_total':len(checks),'production_histories':Nall,'control_histories':Nc,'maximum_production_relative_event_energy_error':emax,'checks':checks};(P/'validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
if __name__=='__main__':main()
