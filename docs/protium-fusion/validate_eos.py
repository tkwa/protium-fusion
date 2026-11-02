"""Compare the saved EOS interpolant with its source at off-grid states."""
import json
from pathlib import Path
import numpy as np
import pynucastro as p
from fast_eos import FastEOS

slow=p.StellarEOS();fast=FastEOS()
comp=p.Composition(['p','n14']);comp['p']=1/15;comp['n14']=14/15
rng=np.random.default_rng(20261028)
points=list(zip(10**rng.uniform(-4,6.9,32),10**rng.uniform(5.2,9.45,32)))
points += [(1e5,3e8),(1e5,1.45e9),(1e3,4.5e8),(1e7,5e8),(1e7,1.6e9)]
rows=[]
for rho,temp in points:
    a=slow.pe_state(rho,temp,comp);b=fast.pe_state(rho,temp,comp)
    rows.append({'rho_g_cm3':rho,'T_K':temp,
                 **{k:float(abs(getattr(b,k)/getattr(a,k)-1))
                    for k in ['p','e','c_v','gamma1']}})
Path(__file__).with_name('eos-validation-extended.json').write_text(json.dumps(rows,indent=2)+'\n')
print({k:max(r[k] for r in rows) for k in ['p','e','c_v','gamma1']})
