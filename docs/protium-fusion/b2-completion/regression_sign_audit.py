"""Sign-only check of scalar artificial-viscosity heating on frozen traces."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--traces',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
rows=[]
for stage in ['A','B1','C13','He']:
 path=a.traces/f'scaled-hot-{stage}-64-R{3 if stage=="He" else 30}.npz'
 with np.load(path) as z:
  n=int(z['zones']);nv=int(z['variables_per_zone']);nn=len(z['species']);state=z['state'].reshape(n,nv,-1)
  r=np.concatenate([np.zeros((1,state.shape[-1])),state[:,nn+1,:]])
  v=np.concatenate([np.zeros((1,state.shape[-1])),state[:,nn+2,:]])
  # Positive scale factors cancel in both sign tests. Masses are uniform in
  # original radial coordinate cubed, rather than equal across cells.
  negative_heat=(np.diff(v,axis=0)<0)&(np.diff(r*r*v,axis=0)>0)
  fraction=np.diff(np.linspace(0,1,n+1)**3)@negative_heat
  rows.append(dict(stage=stage,trace=str(path),trace_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),stored_snapshots=state.shape[-1],snapshots_with_negative_viscous_heat=int(np.count_nonzero(fraction)),maximum_mass_fraction_with_negative_viscous_heat=float(np.max(fraction)),final_mass_fraction_with_negative_viscous_heat=float(fraction[-1]),scope='Sign-only diagnostic on saved legacy trace; no integration of erroneous heat.'))
a.output.write_text(json.dumps(rows,indent=2)+'\n')
