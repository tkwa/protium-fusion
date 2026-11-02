"""Bounded independent integration/API and provenance checks for the supplement."""
from pathlib import Path
import json,hashlib,math
from bank_cycle import bank_cycle
P=Path(__file__).resolve().parent
r=json.loads((P/'results.json').read_text());d=json.loads((P/'inputs.json').read_text());l=json.loads((P/'lineage.json').read_text());checks=[]
def ck(name,ok,value=None):
 checks.append({'name':name,'passed':bool(ok),'value':value});assert ok,(name,value)
for c in r['two_bank_cycles']:
 a=bank_cycle(standby_K=c['standby_temperature_K'],fraction_Carnot=c['fraction_of_Carnot'],electrical_warming=True)
 for ak,ckey in [('active_bank_radiation_W_m2','active_bank_radiation_electric_W_m2'),('inactive_sensible_cooldown_W_m2','inactive_sensible_cooldown_electric_average_W_m2'),('inactive_radiation_during_cooldown_W_m2','inactive_radiation_during_cooldown_electric_average_W_m2'),('inactive_standby_radiation_W_m2','inactive_standby_radiation_electric_average_W_m2'),('inactive_total_W_m2','total_inactive_electric_average_with_electrical_warming_W_m2'),('both_banks_peak_5K_cooling_W_m2','both_banks_5K_end_cooling_load_W_m2'),('both_banks_peak_electric_W_m2','both_banks_peak_electric_W_m2_local_module')]:
  ck(f'API {c["standby_temperature_K"]}/{c["fraction_of_Carnot"]}/{ak}',abs(a[ak]-c[ckey])<1e-12,abs(a[ak]-c[ckey]))
for path,h in d['source_sha256'].items():
 f=Path(path);ck(f'original frozen source {f.name}',hashlib.sha256(f.read_bytes()).hexdigest()==h)
for case in l['cases']:
 for row in case['species']:
  I=row['incident_I_kg'];L=row['direct_local_L_kg'];C=row['captured_from_cavity_integral_C_kg'];delta=row['delta_cavity_G_kg'];ck(f'finite lineage {case["effective_capture_H_He"]}/{row["species"]}',abs(L+C+delta-I)/I<2e-12)
for c in r['checks']+l['checks']:ck('saved '+c['name'],c['passed'])
for T in [200.,200.189,210.]:
 a=bank_cycle(sink_K=T);b=bank_cycle(sink_K=T,electrical_warming=False)
 expected=4.8666935800053074e-5*1000*(T-5)/3600
 ck(f'warming branch{T}',abs(a['inactive_total_W_m2']-b['inactive_total_W_m2']-expected)<1e-14)
 ck(f'source bypass{T}',a['source_hot_to_cold_bypass_increment_W_m2']==0.)
out={'checks':checks,'passed':sum(x['passed'] for x in checks),'total':len(checks),'scope':'API consistency, finite atom closure and unchanged source provenance. These checks do not qualify hardware or source populations.','source_sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['inputs.json','results.json','regional.json','lineage.json','bank_cycle.py']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':out['passed'],'total':out['total']}))
