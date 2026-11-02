"""Verify archived response means, sampling errors, closure and basis preservation."""
from pathlib import Path
import json,gzip,math,hashlib
P=Path(__file__).resolve().parent
basis=json.loads((P/'combined-basis.json').read_text());base=json.loads((P/'../coupled-photon-transport/basis.json').resolve().read_text());checks=[]
def ck(name,ok):
 checks.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
cache={};histories=0;maxerr=0.
for col,c in basis['columns'].items():
 for r in c['responses']:
  p=(P/r['event_archive']).resolve()
  if p not in cache:
   with gzip.open(p,'rt') as f:cache[p]=json.load(f)
  e=cache[p]['records'][r['event_archive_record']];values=e['values'];E=e['calorimetric_energy_MeV'];n=len(values);a=[v[1]/E for v in values];mean=math.fsum(a)/n;se=(math.fsum((v-mean)**2 for v in a)/(n*(n-1)))**.5;err=max(abs(math.fsum(v)/E-1) for v in values);tag=f'{col}:{r["particle"]}:{E}'
  ck(tag+' count',n==r['histories']);ck(tag+' mean',abs(mean-r['leakage_fraction'])<1e-13);ck(tag+' sampling SE',abs(se-r['sampling_standard_error'])<1e-13);ck(tag+' energy',err<1e-7);ck(tag+' all tracks complete',all(v[3]==0 for v in values));histories+=n;maxerr=max(maxerr,err)
for col,c in base['columns'].items():
 new=json.loads(json.dumps(basis['columns'][col]))
 for r in new['responses']:r['event_archive']=r['event_archive'].removeprefix('../coupled-photon-transport/')
 ck('original column '+col+' unchanged',new==c)
print(json.dumps({'checks_passed':len(checks),'histories':histories,'max_event_relative_energy_error':maxerr,'checks':checks},indent=2))
