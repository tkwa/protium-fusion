"""Saved-sample hot-mass coverage of failed and constructive B2 source profiles."""
from pathlib import Path
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parent;raw=Path('/home/tkwa/src/protium-fusion/output/b2-source');rows=[]
for name in ['B2-H2-N2000-carbon-median','B2-H2-constant-contact-carbon-median']:
 stem=raw/name;d=json.loads(stem.with_suffix('.json').read_text());z=np.load(stem.with_suffix('.npz'));t=z['time_s'];a=z['state'].reshape(int(z['zones']),int(z['variables_per_zone']),len(t));nn=len(z['species']);T=a[:,nn,:]*d['T0_K'];w=z['cell_mass_kg']/sum(z['cell_mass_kg']);fractions={}
 for threshold in [1e8,3e8,5e8,1e9]:
  f=w@(T>threshold);k=int(np.argmax(f));fractions[str(threshold)]={'max_mass_fraction':float(f[k]),'time_at_max_s':float(t[k])}
 meanT=w@T;k=int(np.argmax(meanT));rows.append({'case':name,'sample_count':len(t),'peak_massmean_T_K':float(meanT[k]),'time_at_peak_mean_s':float(t[k]),'solver_step_peak_T_K':d['peak_temperature_K'],'mass_above_T_K':fractions,'metadata_sha256':hashlib.sha256(stem.with_suffix('.json').read_bytes()).hexdigest(),'trace_sha256':hashlib.sha256(stem.with_suffix('.npz').read_bytes()).hexdigest()})
out={'scope':'Maxima over saved samples, not all solver states. Temperature alone is not a reaction-rate or burn-yield calculation. This illustrates spatial coverage; it does not isolate a unique causal parcel.','cases':rows,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(P/'temperature-coverage.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
