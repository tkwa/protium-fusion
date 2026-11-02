"""Replay existing component kernels on final delivered-B1/C13-CH2 source family.
Only metadata schema, sigma grid and destinations change. No new hydrodynamics.
"""
from pathlib import Path
import json,hashlib
P=Path(__file__).parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
source=P/'final-family/delivered-B1-C13H2/integration.json';d=json.loads(source.read_text())
adapter=P/'final-receiver-thermal-input.json'
a={'families':{'source128':{'stage_inputs':d['stage_inputs'],'cycle':d['hot_cycle'],'capacity_plants':{k.replace('.0',''):v for k,v in d['capacity_plants_preparation_external_no_recovery'].items()}}},'source_family':d['family'],'source_file':str(source),'source_sha256':sha(source),'metadata_adapter':'Exact hot stage arrays/records/caches; capacity-plant key spelling changed for the existing thermal interface. The precursor preparation-radiation increment remains separately allocated to cold equipment by the plant ledger.'}
adapter.write_text(json.dumps(a,indent=2,allow_nan=False)+'\n')
operations=[
 ('source128_thin_thermal.py',{"P/'source-family-carbon-median128/integration.json'":"P/'final-receiver-thermal-input.json'",'for k,cp in[(50.,1800.),(25.,1800.),(100.,1800.)]:':'for k,cp in[(50.,1800.)]:','for sigma in[500.,750.,1000.]:':'for sigma in[200.,250.,350.,500.]:',"P/'source128-thin-thermal.json'":"P/'final-receiver-thermal.json'",'source128 carbon median, exact prescribed source profile remap':'delivered-B1/C13-CH2 source128 carbon median, exact archived hot source family'}),
 ('source128_prompt_proxy.py',{"P/'source-family-carbon-median128/integration.json'":"P/'final-receiver-thermal-input.json'",'/output/source128-prompt-proxy':'/output/final-receiver-prompt-proxy',"P/'source128-prompt-proxy.json'":"P/'final-receiver-prompt-proxy.json'"}),
 ('source128_prompt_radiance.py',{"P/'source-family-carbon-median128/integration.json'":"P/'final-receiver-thermal-input.json'","P/'source128-prompt-radiance.json'":"P/'final-receiver-prompt-radiance.json'"}),
 ('source128_thermal_combine.py',{"P/'source128-thin-thermal.json'":"P/'final-receiver-thermal.json'","P/'source128-prompt-proxy.json'":"P/'final-receiver-prompt-proxy.json'","P/'source128-prompt-radiance.json'":"P/'final-receiver-prompt-radiance.json'","P/'source128-thermal-combined.json'":"P/'final-receiver-combined.json'"}),
]
provenance=[]
for name,changes in operations:
    original=P/name;text=original.read_text()
    for before,after in changes.items():
        assert text.count(before)==1,(name,before,text.count(before))
        text=text.replace(before,after)
    exec(compile(text,str(original),'exec'),{'__file__':str(original),'__name__':'final_component_replay'})
    provenance.append({'parent_script':name,'parent_sha256':sha(original),'substitutions':changes})
(P/'final-receiver-thermal-replay.json').write_text(json.dumps({'code_sha256':sha(__file__),'family_source_sha256':sha(source),'family':d['family'],'operations':provenance},indent=2)+'\n')
