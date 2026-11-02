"""Downsample saved raw traces for a compact, explicitly labeled archive."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();rows=[]
for p in sorted(a.source.glob('scaled-hot-*.npz')):
 if p.name.endswith('-compact.npz'):continue
 raw=np.load(p);z=int(raw['zones']);nv=int(raw['variables_per_zone']);nn=nv-5;nt=len(raw['time_s']);state=raw['state'].reshape(z,nv,nt)
 peak=int(np.unravel_index(np.argmax(state[:,nn,:]),(z,nt))[1]);ix=np.unique(np.r_[np.rint(np.linspace(0,nt-1,min(32,nt))).astype(int),np.arange(max(0,nt-8),nt),peak]);dest=a.output/(p.stem+'-compact.npz')
 np.savez_compressed(dest,time_s=raw['time_s'][ix],state=raw['state'][:,ix],species=raw['species'],variables_per_zone=raw['variables_per_zone'],zones=raw['zones'],original_sample_indices=ix,full_trace_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),sampling_note='32 evenly spaced stored-index samples plus last eight and maximum stored temperature. Not the full solver trace.')
 rows.append({'file':dest.name,'source':str(p),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'original_samples':nt,'archived_samples':len(ix),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
(a.output/'compact-traces.json').write_text(json.dumps({'scope':'Downsampled visualization/audit archive; full raw trace is needed to reproduce original quadrature diagnostics.','rows':rows},indent=2)+'\n')
print('Saved',len(rows),'compact traces')
