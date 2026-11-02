"""Apply the unchanged delayed-spectrum reducer to the exact source128 family."""
from pathlib import Path
import hashlib,importlib.util,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
CODE=ROOT/'component-design/current-photon128/current_photon128.py'
spec=importlib.util.spec_from_file_location('current_photons',CODE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
FAMILY=ROOT/'component-design/source-family-carbon-median128'
FILES=[FAMILY/'source128-reduced-inputs.json',FAMILY/'source128-sigma5000.json',ROOT/'component-design/current-photon128/sources.json']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def main():
 sp=m.extract(ROOT,*[read(p) for p in FILES]);sp['family']='median-carbon128 A3/B1H2/B2H2/C13H4 source-mapped inflow'
 sp['input_file_sha256']={str(p):sha(p) for p in FILES+[CODE]};sp['source_adapter_sha256']=sha(__file__)
 (HERE/'spectrum.json').write_text(json.dumps(sp,indent=2,allow_nan=False)+'\n')
 rows=[{'sigma_kg_m2':s,'scenarios':{mode:m.evaluate(sp,s,mode=mode) for mode in ['nominal_lines','broad_rest','broad_actual_bulk','broad_D1p10']},'broad_plus_carried6p4ppm':m.evaluate(sp,s,hot_bypass_ppm=6.4)} for s in [100,200,300,400,500,600,750,1000,1250,1500,2000,3000,5000]]
 out={'family':sp['family'],'spectrum_sha256':sha(HERE/'spectrum.json'),'code_sha256':sha(__file__),'rows':rows,'carried_reserve_scope':'6.4ppm is an inherited hot direct-bypass sensitivity, not a recomputed source128 nuclear-flow bound. All existing localized-cascade, source-position and secondary-transport assumptions remain.'}
 (HERE/'containment.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
 print(json.dumps({'family':sp['family'],'HplusU_MeV_He':sp['thermal_boundary_MeV_per_He'],'max_endpoint_bulk_D':sp['max_endpoint_bulk_D'],'minimum_D1p10_thermal_sigma':min(v['D1p10_endpoint_thermal_tail_screen']['margin_sigma'] for v in sp['stages'].values()),'broad_pluscarry':[{'sigma':r['sigma_kg_m2'],'loss':r['broad_plus_carried6p4ppm']['leakage_upper_fraction_thermal']} for r in rows]},indent=2))
if __name__=='__main__':main()
