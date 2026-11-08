"""Compare saved Noh results with exact finite-volume density averages."""
import json
from pathlib import Path
import numpy as np
from analyze import snapshot

base=Path('/var/tmp/tkwa-chatgpt/protium-stability')
out=[]
cases=[(n,f'noh-{n}x16-v3') for n in [128,256,512]]
if (base/'runs/noh-512x16-hlle/execution.json').exists():cases.append((512,'noh-512x16-hlle'))
for n,name in cases:
    case=base/'runs'/name
    t,r,theta,data=snapshot([sorted(case.glob('*.vtk'))[-1]])
    rho=data['dens'];center=.75*np.diff(r**4)/np.diff(r**3)
    left,right=r[:-1],r[1:]
    middle=np.clip(t/3,left,right)
    exact=(64*(middle**3-left**3)+(right+t)**3-(middle+t)**3)/(right**3-left**3)
    mask=(center>.06)&(center<.5);weight=np.diff(r**3)
    error=np.sum(weight[mask]*abs(rho[0,mask]-exact[mask]))/np.sum(weight[mask]*exact[mask])
    out.append(dict(name=name,zones=n,time=t,cell_average_relative_L1_density_error=float(error),
        postshock_mean_density=float(np.mean(rho[0,(center>.06)&(center<.15)])),
        max_angular_range=float(np.ptp(rho,axis=0).max())))
(base/'validation/noh-validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
