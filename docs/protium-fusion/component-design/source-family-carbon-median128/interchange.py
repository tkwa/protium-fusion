"""Schema-only source128 interchange for existing spectrum/receiver tools.
No new material, receiver, or plant calculation and no physical-field changes.
"""
from pathlib import Path
import json,copy,hashlib
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
src=P/'integration.json';d=json.loads(src.read_text());family=d['families']['source128'];cases={};manifest={};labels={};hmap={};adapters=[]
for cat,row in family['stage_inputs'].items():
 e=copy.deepcopy(row['selected']);tag=e['record']['stage'];Y=e['stage_data']['initial_Y'];h=Y['p']/Y[cat];label=f'source128-{tag}-H{h:g}';e['label']=label;labels[tag]=label;hmap[tag]=h
 if 'hydrogen_per_catalyst' not in e['record']['mixture_inputs']:
  e['record']['mixture_inputs']['hydrogen_per_catalyst']=h;adapters.append({'case':label,'field':'record.mixture_inputs.hydrogen_per_catalyst','value':h,'scope':'Derived exactly from archived initialY; metadata-only addition for interface compatibility.'})
 else:assert abs(e['record']['mixture_inputs']['hydrogen_per_catalyst']-h)<1e-12
 cases[label]=e;manifest[label]={'record':e['path'],'stage':tag,'hydrogen_per_catalyst':h,'source_sha256':row['input_sha256']}
common={'integration_source':str(src),'integration_sha256':sha(src),'code_sha256':sha(__file__),'scope':__doc__,'metadata_adapters':adapters,'source_preparation':'Unretuned source128 FEOS entropy-mapped hot states, no3Dassemblyqualification.'}
write(P/'source128-reduced-inputs.json',{'schema':1,**common,'cases':cases,'source_manifest':{'family':'source128-A3-B1H2-B2H2-C13H4-carbon-median','cases':manifest},'checks':d['checks']})
for sig,p in family['capacity_plants'].items():
 c=family['cycle'];summary={'cases':labels,'H_per_catalyst':hmap,'mW_kg':p['net_W_kg']*1000,'net_W':p['net_W'],'mass_kg':p['mass_kg'],'helium_per_s':p['helium_per_s'],'occupancy':p['sequential_occupancy_mass999'],'radiator_fraction':p['radiator_area_fraction'],'prepared_MeV_He':c['U_prepared_MeV_per_He'],'fresh_C_H':c['fresh_C_H_mass_ratio'],'pulse_interval_s':p['pulse_interval_s']}
 write(P/f'source128-sigma{sig}.json',{'schema':1,'family':'source128-A3-B1H2-B2H2-C13H4-carbon-median',**common,'comparison':dict(d['parameters'],zones=128,receiver_areal_mass_kg_m2=float(sig),energy_capture_fraction=.99),'top_detailed':[{'summary':summary,'plant':p,'cycle':c}],'ranked':[summary],'qualification':d['qualification']})
 print(json.dumps({'sigma':sig,'B1_interval_s':p['pulse_interval_s']['n15'],'mW_kg':p['net_W_kg']*1000}))
