"""Join exact-profile rate comparisons across saved64/128 meshes; no reruns."""
from pathlib import Path
import json,hashlib,argparse
ap=argparse.ArgumentParser();ap.add_argument('--runs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();records={}
for family in ['nominal64','extension64','nominal128','root256']:
 for p in (a.runs/family).glob('*-rates-updated-median.json'):
  d=json.loads(p.read_text());bpath=Path(d['rate_update']['baseline_record']);b=json.loads(bpath.read_text());name=p.name.split('-fixed')[0]
  records.setdefault(name,{})[str(d['zones'])]={'updated_record':str(p.resolve()),'updated_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'baseline_record':str(bpath),'old_gain':b['deposited_nuclear_over_initial_energy'],'updated_gain':d['deposited_nuclear_over_initial_energy'],'old_target':b['desired_product_per_initial_catalyst'],'updated_target':d['desired_product_per_initial_catalyst'],'energy_residual_fraction_initial':d['energy_residual_fraction_initial'],'physical_duration_s':d['elapsed_s'],'requested_duration_s':d['mixture_inputs']['duration_s'],'same_initial_energy':d['initial_energy_J']==b['initial_energy_J'],'completion':d['termination_reason']}
for name,row in records.items():
 for lower,upper in [('64','128'),('128','256')]:
  if lower in row and upper in row:row['updated_'+upper+'_over_'+lower+'_gain_change']=row[upper]['updated_gain']/row[lower]['updated_gain']-1
out={'scope':'Exact-profile within-mesh rate comparisons; initial profiles are nested/refined baselines across meshes. Nonmonotone numerical changes are retained, not extrapolated away.','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':records}
a.out.write_text(json.dumps(out,indent=2)+'\n')
print('case old64 updated64 old128 updated128 old256 updated256')
for n,row in records.items():print(n,*[f"{row.get(z,{}).get(k,float('nan')):.9f}" for z in ['64','128','256'] for k in ['old_gain','updated_gain']])
