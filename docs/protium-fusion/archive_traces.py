"""Archive only the time-series quantities used by downstream screens.

Full local hydro traces can be regenerated from compression_probe.py and the
saved input parameters. This compact archive avoids storing every isotope in
every cell for every checkpoint in the website repository.
"""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).parent
archive=ROOT/'compression-traces.json'
out=json.loads(archive.read_text()) if archive.exists() else {
    'scope':'Decimated hydro checkpoints; extrema are sampled, not rigorous bounds.', 'traces':{}}
for path in sorted(ROOT.glob('compression-*-burning.json')):
    trace=path.with_suffix('.npz')
    if not trace.exists(): continue
    d=json.loads(path.read_text());a=np.load(trace)
    nz=int(a['zones']);nv=int(a['variables_per_zone']);nn=nv-5
    z=a['state'].reshape(nz,nv,-1)
    weights=np.array(d['profile']['mass_fraction_of_total'])
    rad=weights@z[:,nn+4,:]
    if rad[-1]>0: rad*=d['radiation_escaped_J']/rad[-1]
    r=np.vstack([np.zeros(z.shape[-1]),z[:,nn+1]*d['R0_m']])
    rho=d['mass_kg']*weights[:,None]/(4*np.pi/3*np.diff(r**3,axis=0))
    out['traces'][path.stem]={'source_npz_sha256':hashlib.sha256(trace.read_bytes()).hexdigest(),
      'zones':nz,'species':a['species'].tolist(),'time_s':a['time_s'].tolist(),
      'cumulative_radiation_escaped_J':rad.tolist(),
      'sampled_peak_density_kg_m3':float(rho.max())}
archive.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print('Archived',len(out['traces']),'compact traces')
