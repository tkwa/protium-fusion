"""Mass/solid-angle-weighted Gaussian-quantile timing errors on 3D patches.

This fixes the empirical error distribution while permuting its spatial pattern.
The finite 128-point distribution has finite tails; it is not a worst-case bound
over all distributions with a given standard deviation.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import ndtri
from packet_geometry import directions

p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=128)
p.add_argument('--seed',type=int,default=12345);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();assert a.n%2==0
z=(np.arange(a.n//2)+.5)/(a.n//2);phi=np.arange(a.n//2)*np.pi*(3-np.sqrt(5))
v=np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
sites=np.stack([v,-v],axis=1).reshape(-1,3)
sample=directions(524288);distance,index=cKDTree(sites).query(sample,workers=1)
weights=np.bincount(index,minlength=a.n)/len(sample)
errors=ndtri((np.arange(a.n)+.5)/a.n)
np.random.default_rng(a.seed).shuffle(errors)
errors-=weights@errors;errors/=np.sqrt(weights@(errors**2))
np.savetxt(a.out,errors,fmt='%.17g')
report=dict(n=a.n,seed=a.seed,angular_weights=weights.tolist(),errors=errors.tolist(),
            mean=float(weights@errors),rms=float(np.sqrt(weights@(errors**2))),
            minimum=float(errors.min()),maximum=float(errors.max()),
            direction_sine_rms=float(np.sqrt(np.mean(1-(1-distance**2/2)**2))))
a.out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print({k:report[k] for k in ['n','seed','mean','rms','minimum','maximum','direction_sine_rms']})
