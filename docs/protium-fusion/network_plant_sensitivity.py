"""Selected engineering sensitivities at the 140 T recycling reference."""
from pathlib import Path
import json,gzip
from scipy.optimize import brentq
from network_plant import plant
ROOT=Path(__file__).parent
energy=json.loads((ROOT/'network-ion-budget.json').read_text());family=json.loads(gzip.decompress((ROOT/'network-channel-family.json.gz').read_bytes()))
variants=[('reference',{}),('8 GPa hoop / 4 GPa compression',{'sigma_hoop':8e9,'sigma_comp':4e9}),
 ('90% ion heating',{'heater':.9}),('90% direct conversion',{'direct':.9}),
 ('99.8% coherent end transfer',{'end_in':.998,'end_out':.998}),
 ('1 kg/W conversion equipment',{'conversion_kg_W':1.}),
 ('0.03% annual mass replacement',{'annual_replacement_fraction':3e-4}),
 ('1 MeV additional electrical load per primary',{'extra_processing_MeV_per_primary':1.}),
 ('2 MeV additional electrical load per primary',{'extra_processing_MeV_per_primary':2.})]
rows=[]
for name,kw in variants:
 r=plant(energy,family,**kw);r['variant']=name;rows.append(r)
 print(name,round(r['net_W_kg']*1e3,6),flush=True)
limit=brentq(lambda x:plant(energy,family,extra_processing_MeV_per_primary=x)['net_W_kg']-.001,0.,3.,xtol=1e-10)
(ROOT/'network-plant-sensitivity.json').write_text(json.dumps({'scope':__doc__,'cases':rows,'additional_electrical_MeV_per_primary_at_1mW_kg':limit},indent=2,allow_nan=False)+'\n')
print('additional load at 1 mW/kg',limit)
