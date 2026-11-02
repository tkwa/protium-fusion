"""Independent finite-skin routine and saved-ledger checks; no hydro rerun."""
from pathlib import Path
import json,hashlib,importlib.util,math
import numpy as np
P=Path(__file__).resolve().parent;root=P.parents[2];base=root/'docs/protium-fusion'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
s=importlib.util.spec_from_file_location('independent_stopping',base/'mixture-design/receiver_stopping.py');rs=importlib.util.module_from_spec(s);s.loader.exec_module(rs);pm,st,decay,_,_=rs.setup(root)
checks=[]
for name in ['surface-damage.json','surface-damage-current128.json']:
 d=json.loads((P/name).read_text());assert d['code_sha256']==sha(P/'surface_damage.py')
 for path,expected in d['input_sha256'].items():assert sha(root/path if not Path(path).is_absolute() else path)==expected,path
 # Use the pre-existing membrane finite-thickness implementation, independently
 # of the new depth-bin integrator, for the most demanding C13 spectrum.
 stage=next(x for x in d['stages'] if x['stage']=='C13');record=next(Path(p) for p in d['input_sha256'] if Path(p).name==stage['case']+'.json')
 rec,n,eu,_,_,_,_=rs.spectrum(record,pm,decay,d['radius_m'],64)
 for factor in [.5,1.,2.]:
  r=next(x for x in stage['rows'] if x['electronic_stopping_multiplier']==factor)
  for thick in [50.,100.]:
   rr=st.run_case('carbon',thick*1000,n,eu[:,None]*pm.AA,pm.AA,pm.ZZ,d['area_m2'],factor)
   got=next(x for x in r['skins'] if x['thickness_um']==thick)
   err=rr['elastic_transfer_J_m2']/got['elastic_transfer_J_m2']-1;dp=rr['NRT_proxy_dpa_Ed40_eV']/got['NRT_proxy_dpa_per_pulse']['40.0']-1
   assert abs(err)<1e-12 and abs(dp)<1e-12
   checks.append({'family':name,'electronic_multiplier':factor,'skin_um':thick,'elastic_transfer_relative':err,'NRT40_relative':dp})
 for a in d['annual']:
  for skin in a['skins']:
   for ed,dpa in skin['NRT_proxy_dpa_year'].items():
    direct=sum(next(x for x in next(r for r in stage['rows'] if r['electronic_stopping_multiplier']==a['electronic_stopping_multiplier'])['skins'] if x['thickness_um']==skin['thickness_um'])['NRT_proxy_dpa_per_pulse'][ed]*365.25*86400/stage['period_s'] for stage in d['stages'])
    assert math.isclose(direct,dpa,rel_tol=1e-12)
 assert d['prompt_reference_matches_selected']==(name=='surface-damage.json')
 assert all(x['process_W_m2_at_MJ_kg']['2.0']>=0 for x in d['restoration'])
print(json.dumps({'finite_skin_comparisons':checks,'source_hashes_and_annual_dose_identity':'pass'},indent=2))
(P/'surface-damage-validation.json').write_text(json.dumps({'code_sha256':sha(__file__),'input_sha256':{name:sha(P/name)for name in ['surface-damage.json','surface-damage-current128.json']},'finite_skin_comparisons':checks,'source_hashes_and_annual_dose_identity':'pass'},indent=2)+'\n')
