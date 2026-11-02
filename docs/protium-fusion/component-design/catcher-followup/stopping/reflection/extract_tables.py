"""Extract numerical coefficients from official TRIM downloads; raw reports stay in scratch."""
from pathlib import Path
import hashlib,json,sys
P=Path(__file__).resolve().parent
raw=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/home/tkwa/src/protium-fusion/output/catcher-followup/reflection')
names=['H_on_C.txt','H_on_C_1.00.txt','He_on_C.txt','He_on_C_0.00.txt']
out={'scope':'Primary BCA simulation tables, not hot irradiated laminate measurements. Numerical factual extraction only.','datasets':{}}
for name in names:
 p=raw/name; lines=p.read_text().splitlines(); rows=[]
 for i,line in enumerate(lines):
  try:a=[float(v)for v in line.split()]
  except ValueError:continue
  if len(a)<7 or a[0] not in[1,2] or a[2]!=6 or not 11<a[3]<13:continue
  nums=[]
  for line2 in lines[i+1:]:
   if not line2.strip():continue
   nums.append([float(v)for v in line2.split()])
   if len(nums)==6:break
  row={'E_eV':a[4],'angle_deg':a[5],'RN':a[6],'RE':a[7] if len(a)>7 else a[6]*sum(nums[0][:5])/5/a[4],'energy_fraction_quantiles_10_30_50_70_90':[v/a[4]for v in nums[0][:5]],'cosine_quantiles_conditional_on_energy': [v[:5]for v in nums[1:]]}
  if len(nums[0])==7:row['energy_fraction_bounds']=[v/a[4]for v in nums[0][5:]]
  rows.append(row)
 assert len(rows)==84
 out['datasets'][name]={'url':'https://www.eirene.de/old_eirene/'+name,'raw_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size,'rows':rows}
out['documentation']={'url':'https://www.eirene.de/old_eirene/trim.pdf','raw_sha256':hashlib.sha256((raw/'eirene-trim.pdf').read_bytes()).hexdigest(),'scope':'Raw energy quantiles and bounds are in eV; this extract divides by incident E to save fractions. Five equal-probability energy and conditional cosine quantiles; RN is particle reflection and RE is reflected incident-energy fraction. Old-format RE above is the documented five-quantile estimator.'}
out['checks']={'datasets':4,'rows_per_dataset':84,'unique_pairs':all(len({(r['E_eV'],r['angle_deg'])for r in d['rows']})==84 for d in out['datasets'].values()),'probabilities':all(0<=r['RE']<=r['RN']<=1 for d in out['datasets'].values()for r in d['rows'])}
(P/'sources/primary-tables.json').write_text(json.dumps(out,indent=2)+'\n')
print(out['checks'])
